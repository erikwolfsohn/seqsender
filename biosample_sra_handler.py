#!/usr/bin/env python3

###########################	Description	##################################
# Functions for NCBI input/output handling
################################################################################

import pandas as pd
from typing import Set, Dict, Any, Tuple, List
import os
from lxml import etree
import sys
import re
import ncbi_handler
import file_handler
import upload_log

from settings import BIOSAMPLE_REGEX, SRA_REGEX

# Check raw reads files listed in metadata file
def check_raw_read_files(submission_name: str, submission_dir: str, metadata: pd.DataFrame) -> Set[str]:
	# Pop off the end directories of submission_dir 'submission_files/SRA'
	raw_reads_path_default = os.path.join(os.path.split(os.path.split(submission_dir)[0])[0], "raw_reads")
	# Separate samples stored in local and cloud
	local_df = metadata[metadata["sra-file_location"] == "local"]
	file_columns = [col for col in local_df.columns if re.match("sra-file_[1-9]\d*", col)]
	validated_files = set()
	invalid_raw_files = []
	for index, row in local_df.iterrows():
		# If multiple files check each one
		for column_name in file_columns:
			file = row[column_name]
			# At least one file required for SRA but extra file columns can have empty slots
			if (file is None or file.strip() == "") and column_name != "sra-file_1":
				continue
			file = file.strip()
			# If full path don't use rawreads folderfil
			if os.path.isabs(file):
				file_path = file
			else:
				# If just file name use raw reads default path
				file_path = os.path.join(raw_reads_path_default, file)
			if os.path.isfile(file_path) == False:
				invalid_raw_files.append(f"Error: Raw read files for {row['sra-sample_name']} could not be found. Field '{column_name}' must either be the full file path, or if just the file name it must be stored at '<submission_dir>/<submission_name>/raw_reads/<sra-file>'.")
			else:
				validated_files.add(file_path)
	if invalid_raw_files:
		for error_msg in invalid_raw_files:
			print(error_msg, file=sys.stderr)
		sys.exit(1)
	return validated_files

# Create files for optional manual submission to repositories biosample and sra
def create_manual_submission_files(database: str, submission_dir: str, metadata: pd.DataFrame, config_dict: Dict[str, Any]) -> None:
	if "SRA" in database:
		metadata_regex = "^sra-|^organism$|^collection_date$"
		rename_columns = {"sra-library_name":"sra-library_ID"}
		drop_columns = ["sra-file_location", "sra-loader", "sra-file"]
		column_ordered = ["sample_name","library_ID"]
		prefix = "sra-"
		# Create SRA specific fields
		if 'sra-title' not in metadata.columns or metadata['sra-title'].isnull().any() or (metadata['sra-title'] == '').any():
			metadata["sra-title"] = config_dict["Description"]["Title"]
		filename_cols = [col for col in metadata.columns.tolist() if re.match("sra-file_[1-9]\d*", col)]
		# Correct index for filename column
		for col in filename_cols:
			# Remove 0 index
			if col == "sra-file_1":
				rename_columns[col] = "sra-filename"
			else:
				rename_columns[col] = col.replace("sra-file_", "sra-filename")
	elif "BIOSAMPLE" in database:
		metadata_regex = "^bs-|^organism$|^collection_date$"
		rename_columns = {"bs-description":"sample_title","bioproject":"bioproject_accession"}
		drop_columns = ["bs-package"]
		column_ordered = ["sample_name"]
		prefix = "bs-"
	else:
		print("Error: create_manual_submission_files function only for databases SRA/BioSample. Not '{database}'.", file=sys.stderr)
		sys.exit(1)
	# Filter only needed columns
	database_df = metadata.filter(regex=metadata_regex).copy()
	database_df = database_df.drop_duplicates()
	# Rename columns for manual submission
	database_df = database_df.rename(columns={key: value for key, value in rename_columns.items() if key in database_df})
	# Remove columns not required for manual submission
	database_df = database_df.drop(columns=[col for col in drop_columns if col in database_df.columns])
	# Remove database prefix
	database_df.columns = database_df.columns.str.replace(prefix,"")
	# Order columns
	columns_no_order = list(filter(lambda x: (x not in column_ordered)==True, database_df.columns))
	column_ordered = column_ordered + columns_no_order
	database_df = database_df.reindex(columns=column_ordered)
	file_handler.save_csv(df=database_df, file_path=submission_dir, file_name="metadata.tsv", sep="\t")

# Create submission XML
def create_submission_xml(organism: str, database: str, submission_name: str, config_dict: Dict[str, Any], metadata: pd.DataFrame, failed_seqs_auto_removed: bool = True) -> bytes:
	# Submission XML header
	root = etree.Element("Submission")
	description = etree.SubElement(root, "Description")
	title = etree.SubElement(description, "Title")
	title.text = config_dict["Description"]["Title"]
	comment = etree.SubElement(description, "Comment")
	comment.text = config_dict["Description"]["Comment"]
	# Description info including organization and contact info
	organization = etree.SubElement(description, "Organization", type=config_dict["Description"]["Organization"]["Type"], role=config_dict["Description"]["Organization"]["Role"])
	org_name = etree.SubElement(organization, "Name")
	org_name.text = config_dict["Description"]["Organization"]["Name"]
	contact = etree.SubElement(organization, "Contact", email=config_dict["Description"]["Organization"]["Submitter"]["Email"])
	name = etree.SubElement(contact, "Name")
	first_name = etree.SubElement(name, "First")
	first_name.text = config_dict["Description"]["Organization"]["Submitter"]["Name"]["First"]
	last_name = etree.SubElement(name, "Last")
	last_name.text = config_dict["Description"]["Organization"]["Submitter"]["Name"]["Last"]
	if config_dict["Specified_Release_Date"]:
		release_date = etree.SubElement(description, "Hold", release_date=config_dict["Specified_Release_Date"])
	# XML actions
	if "BIOSAMPLE" in database:
		database_df = metadata.filter(regex=BIOSAMPLE_REGEX).copy()
		database_df = database_df.drop_duplicates()
		for index, row in database_df.iterrows():
			action = etree.SubElement(root, "Action")
			add_data = etree.SubElement(action, "AddData", target_db="BioSample")
			data = etree.SubElement(add_data, "Data", content_type="xml")
			xmlcontent = etree.SubElement(data, "XmlContent")
			biosample = etree.SubElement(xmlcontent, "BioSample", schema_version="2.0")
			sampleid = etree.SubElement(biosample, "SampleId")
			spuid = etree.SubElement(sampleid, "SPUID", spuid_namespace=config_dict["Spuid_Namespace"])
			spuid.text = row["bs-sample_name"]
			descriptor = etree.SubElement(biosample, "Descriptor")
			title = etree.SubElement(descriptor, "Title")
			title.text = row["bs-description"]
			organismxml = etree.SubElement(biosample, "Organism")
			organismname = etree.SubElement(organismxml, "OrganismName")
			organismname.text = row["organism"]
			if pd.notnull(row["bioproject"]) and row["bioproject"].strip() != "":
				bioproject = etree.SubElement(biosample, "BioProject")
				primaryid = etree.SubElement(bioproject, "PrimaryId", db="BioProject")
				primaryid.text = row["bioproject"]
			package = etree.SubElement(biosample, "Package")
			package.text = config_dict["BioSample_Package"]
			# Attributes
			attributes = etree.SubElement(biosample, "Attributes")
			# Remove columns with bs-prefix that are not attributes
			biosample_cols = [col for col in database_df.columns.tolist() if (col.startswith('bs-')) and (col not in ["bs-sample_name", "bs-package", "bs-description"])]
			for col in biosample_cols:
				attribute_value = row[col]
				if pd.notnull(attribute_value) and attribute_value.strip() != "":
					attribute = etree.SubElement(attributes, "Attribute", attribute_name=col.replace("bs-",""))
					attribute.text = row[col]
			# Add collection date to Attributes
			attribute = etree.SubElement(attributes, "Attribute", attribute_name="collection_date")
			attribute.text = row["collection_date"]
			identifier = etree.SubElement(add_data, "Identifier")
			spuid = etree.SubElement(identifier, "SPUID", spuid_namespace=config_dict["Spuid_Namespace"])
			spuid.text = row["bs-sample_name"]
	if "SRA" in database:
		database_df = metadata.filter(regex=(SRA_REGEX)).copy()
		database_df = database_df.drop_duplicates()
		file_columns = [col for col in database_df.columns if re.match("sra-file_[1-9]\d*", col)]
		for index, row in database_df.iterrows():
			action = etree.SubElement(root, "Action")
			addfiles = etree.SubElement(action, "AddFiles", target_db="SRA")
			for column_name in file_columns:
				# extra file columns can be empty but not first file column
				if column_name == "sra-file_1" and (row[column_name] is None or row[column_name].strip() == ""):
					print(f"Error: metadata must contain a file for {row['sra-sample_name']} in column sra-file_1", file=sys.stderr)
					sys.exit(1)
				elif column_name != "sra-file_1" and (row[column_name] is None or row[column_name].strip() == ""):
					continue
				if row["sra-file_location"].strip().lower() == "cloud":
					file = etree.SubElement(addfiles, "File", cloud_url = row[column_name].strip())
				elif row["sra-file_location"].strip().lower() == "local":
					file = etree.SubElement(addfiles, "File", file_path = os.path.basename(row[column_name].strip()))
				else:
					print("Error: Metadata field file_location must be either cloud or local. Field currently contains: " + row["sra-file_location"].strip().lower(), file=sys.stderr)
					sys.exit(1)
				datatype = etree.SubElement(file, "DataType")
				datatype.text = "generic-data"
			# Remove columns with sra- prefix that are not attributes
			sra_cols = [col for col in database_df.columns.tolist() if col.startswith('sra-') and not re.match("(sra-sample_name|sra-file_location|sra-file_\d*)", col)]
			for col in sra_cols:
				attribute_value = row[col]
				if pd.notnull(attribute_value) and attribute_value.strip() != "":
					attribute = etree.SubElement(addfiles, "Attribute", name=col.replace("sra-",""))
					attribute.text = row[col]
			if pd.notnull(row["bioproject"]) and row["bioproject"].strip() != "":
				attribute_ref_id = etree.SubElement(addfiles, "AttributeRefId", name="BioProject")
				refid = etree.SubElement(attribute_ref_id, "RefId")
				primaryid = etree.SubElement(refid, "PrimaryId")
				primaryid.text = row["bioproject"]
			if config_dict["Link_Sample_Between_NCBI_Databases"] and metadata.columns.str.contains("bs-sample_name").any():
				attribute_ref_id = etree.SubElement(addfiles, "AttributeRefId", name="BioSample")
				refid = etree.SubElement(attribute_ref_id, "RefId")
				spuid = etree.SubElement(refid, "SPUID", spuid_namespace=config_dict["Spuid_Namespace"])
				spuid.text = metadata.loc[metadata["sra-sample_name"] == row["sra-sample_name"], "bs-sample_name"].iloc[0]
			identifier = etree.SubElement(addfiles, "Identifier")
			spuid = etree.SubElement(identifier, "SPUID", spuid_namespace=config_dict["Spuid_Namespace"])
			spuid.text = row["sra-sample_name"]
	# Pretty print xml
	xml_str = etree.tostring(root, encoding="utf-8", pretty_print=True, xml_declaration=True)
	return xml_str

# Create list of raw read paths inside sra submission folder
def create_raw_reads_list(submission_dir: str, raw_files_list: Set[str]) -> None:
	with open(os.path.join(submission_dir, "raw_reads_location.txt"), "w+") as file:
		for line in raw_files_list:
			file.write(line + "\n")

# Main create function for BioSample/SRA
def create_biosample_sra_submission(organism: str, database: str, submission_name: str, submission_dir: str, config_dict: Dict[str, Any], metadata: pd.DataFrame):
	if database == "SRA":
		# Validate and write raw reads location
		raw_files_list = check_raw_read_files(submission_name=submission_name, submission_dir=submission_dir, metadata=metadata)
		create_raw_reads_list(submission_dir=submission_dir, raw_files_list=raw_files_list)
	manual_df = metadata.copy()
	create_manual_submission_files(database=database, submission_dir=submission_dir, metadata=manual_df, config_dict=config_dict)
	xml_str = create_submission_xml(organism=organism, database=database, submission_name=submission_name, metadata=metadata, config_dict=config_dict, failed_seqs_auto_removed=True)
	file_handler.save_xml(xml_str, submission_dir)

# Read xml report and get status of the submission
def process_biosample_sra_report(report_file: str, database: str, submission_dir: str) -> Tuple[str, str]:
	report_dict, submission_status, submission_id = ncbi_handler.process_report_header(report_file=report_file)
	sample_name_prefix = {"BIOSAMPLE":"bs-", "SRA":"sra-"}
	sample_info = []
	if "Action" not in report_dict["SubmissionStatus"]:
		return submission_status, submission_id
	try:
		for action_dict in report_dict["SubmissionStatus"]["Action"]:
			# Skip if incorrect database
			if "@target_db" not in action_dict or action_dict["@target_db"].lower() != database.lower():
				continue
			# Skip if missing sample info
			elif "Response" not in action_dict:
				continue
			elif isinstance(action_dict["Response"], list) and any("Object" in response_dict for response_dict in action_dict["Response"]) and "@status" in action_dict:
				# print(action_dict["Response"])
				for index in range(len(action_dict["Response"])):
					response_dict = action_dict["Response"][index]
					if "Object" in response_dict and all(field in response_dict["Object"] for field in ["@accession", "@spuid"]):
						sample_name = response_dict["Object"]["@spuid"]
						accession = response_dict["Object"]["@accession"]
						sample_name_col = f"{sample_name_prefix[database]}sample_name"
						column_prefix = database.lower()
						sample_info.append({sample_name_col:sample_name, f"{column_prefix}_status":response_dict["@status"], f"{column_prefix}_accession":accession, f"{column_prefix}_message":""})
						break
			elif isinstance(action_dict["Response"], dict) and "Object" in action_dict["Response"] and "@status" in action_dict:
				if "@accession" in action_dict["Response"]["Object"] and "@spuid" in action_dict["Response"]["Object"]:
					sample_name = action_dict["Response"]["Object"]["@spuid"]
					accession = action_dict["Response"]["Object"]["@accession"]
					sample_name_col = f"{sample_name_prefix[database]}sample_name"
					column_prefix = database.lower()
					sample_info.append({sample_name_col:sample_name, f"{column_prefix}_status":action_dict["@status"], f"{column_prefix}_accession":accession, f"{column_prefix}_message":""})
	except:
		pass
	if sample_info:
		update_df = pd.DataFrame(sample_info)
		upload_log.update_submission_status_csv(submission_dir=submission_dir, update_database=database, update_df=update_df)
	return submission_status, submission_id

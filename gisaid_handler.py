#!/usr/bin/env python3

###########################    Description    ##################################
# Functions for GISAID input/output handling
################################################################################

import shutil
import subprocess
from typing import Dict, Any, List, Optional, Match, Any
import os
import pandas as pd
import file_handler
import sys
import time
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
import warnings
warnings.filterwarnings("ignore", 'This pattern has match groups')
import re

import upload_log
import tools
from settings import GISAID_REGEX

# Create directory and files for GISAID submission
def create_gisaid_files(organism: str, database: str, submission_name: str, submission_dir: str, config_dict: Dict[str, Any], metadata: pd.DataFrame) -> None:
	# Get column names for gisaid submission only
	gisaid_df = metadata.filter(regex=GISAID_REGEX).copy()
	gisaid_df.columns = gisaid_df.columns.str.replace("gs-","").str.strip()
	# Add required GISAID fields
	# covCLI returns an error when authors or collection_date are capitalized
	if organism in ["COV", "POX", "ARBO", "RSV"]:
		if organism == "COV":
			sample_name_column = "covv_virus_name"
			sex_field_sanitized = "covv_sex"
			sex_field_required = "covv_gender"
		else:
			sample_name_column = organism.lower() + "_virus_name"
			sex_field_sanitized = organism.lower() + "_sex"
			sex_field_required = organism.lower() + "_gender"
		gisaid_df = gisaid_df.rename(columns = {"sample_name": sample_name_column, sex_field_sanitized: sex_field_required})
		gisaid_df["submitter"] = config_dict["Username"]
		# fn field is for fasta file name
		gisaid_df["fn"] = "sequence.fsa"
		first_cols = ["submitter", "fn", sample_name_column]
	elif "FLU" in organism:
		gisaid_df = gisaid_df.rename(columns = {"authors": "Authors", "Host_Sex": "Host_Gender"})
		# Parse out dates into respective columns
		gisaid_df[["Collection_Date", "Collection_Year", "Collection_Month"]] = gisaid_df["collection_date"].apply(process_flu_dates)
		gisaid_df["Isolate_Id"] = ""
		gisaid_df["Segment_Ids"] = ""
		# Pivot FLU segment names from long form to wide form
		gisaid_df["segment"] = "Seq_Id (" + gisaid_df["segment"].astype(str) + ")"
		group_df = gisaid_df.pivot(index="Isolate_Name", columns="segment", values="sample_name").reset_index()
		gisaid_df = gisaid_df.drop(columns=["sample_name", "segment", "collection_date"])
		gisaid_df = gisaid_df.drop_duplicates(keep="first")
		gisaid_df = gisaid_df.merge(group_df, on="Isolate_Name", how="inner", validate="1:1")
		first_cols = ["Isolate_Id","Segment_Ids","Isolate_Name"]
	# Restructure column order
	last_cols = [col for col in gisaid_df.columns if col not in first_cols]
	gisaid_df = gisaid_df[first_cols + last_cols]
	# Create submission files
	file_handler.save_csv(df=gisaid_df, file_path=submission_dir, file_name="metadata.csv")
	shutil.copy(os.path.join(submission_dir, "metadata.csv"), os.path.join(submission_dir, "orig_metadata.csv"))
	file_handler.create_fasta(database="GISAID", metadata=metadata, submission_dir=submission_dir, config_dict = config_dict)
	shutil.copy(os.path.join(submission_dir, "sequence.fsa"), os.path.join(submission_dir, "orig_sequence.fsa"))

# Flu collection dates require partial dates to use different columns
def process_flu_dates(row: Any) -> pd.Series:
	sections = row.strip().split("-")
	if len(sections) == 1:
		full_date = ""
		year = sections[0]
		month = ""
	elif len(sections) == 2:
		full_date = ""
		year = sections[0]
		month = sections[1]
	elif len(sections) == 3:
		full_date = row.strip()
		year = ""
		month = ""
	else:
		print(f"Error: Unable to process 'Collection_Date' column for FLU GISAID submission. The field should be in format 'YYYY-MM-DD'. Value unable to process: {row.strip()}", file=sys.stderr)
		sys.exit(1)
	return pd.Series([full_date, year, month])


# Read output log from gisaid submission script
def process_gisaid_log(log_file: str, submission_dir: str) -> pd.DataFrame:
	file_handler.validate_file(file_type="GISAID log", file_path=log_file)
	# Read in log file
	gisaid_isolate_log = []
	gisaid_segment_log = []
	with open(log_file, "r") as file:
		line = file.readline().strip()
		while line:
			# If accession generated record it
			# Pattern options: "msg:": "<Sample Name>; <EPI_ISL/EPI_ID>_<Accession Numbers>" OR <epi_isl_id/epi_id>: <Sample Name>; <EPI_ISL/EPI_ID>_<Accession Numbers>
			if re.search(r"(?i)(\W|^)(\"msg\":\s*\"\S+.*;\s*(EPI_ISL|EPI_ID)_\d*\"|(epi_id|epi_isl_id):\s*\S.*;\s*(EPI_ISL_|EPI)\d+)(\W|$)", line):
				gisaid_string_search = re.findall(r'(?:[a-zA-Z0-9_-]+(?:/[a-zA-Z0-9_-]+)+|EPI_\w*)', line)
				gisaid_string = ' '.join(gisaid_string_search)
				gisaid_string_list: List[str] = gisaid_string.split(' ')
				sample_name = gisaid_string_list[0].strip()
				accession_string = gisaid_string_list[1].strip()
				if re.match(r"EPI_ISL_\d+", accession_string):
					gisaid_isolate_log.append({"gs-sample_name":sample_name, "gisaid_accession_epi_isl_id":accession_string, "gisaid_message":"Successfully uploaded"})
				elif re.match(r"EPI\d+", accession_string):
					gisaid_segment_log.append({"gs-segment_name":sample_name, "gisaid_accession_epi_id":accession_string, "gisaid_message":"Successfully uploaded"})
			elif re.search(r"(?i)SUCCESS;\s+\S.*assigned to accession:\s*EPI_", line):
				success_search = re.search(r"(?i)SUCCESS;\s+(\S.*?)\s+assigned to accession:\s*(EPI_(?:ISL_)?\d+)", line)
				if success_search:
					sample_name = success_search.group(1).strip()
					accession_string = success_search.group(2).strip()
					if re.match(r"EPI_ISL_\d+", accession_string):
						gisaid_isolate_log.append({"gs-sample_name":sample_name, "gisaid_accession_epi_isl_id":accession_string, "gisaid_message":"Successfully uploaded"})
					elif re.match(r"EPI_\d+", accession_string):
						gisaid_segment_log.append({"gs-segment_name":sample_name, "gisaid_accession_epi_id":accession_string, "gisaid_message":"Successfully uploaded"})
			# Handling if submitting samples have already been registered in GISAID (WARNING plain-text format)
			elif re.search(r"(?i)WARNING;\s+\S.*;\s*validation_error;.*already exists.*existing_ids:", line):
				sample_name_search = re.search(r"(?i)WARNING;\s+(\S[^;]+);", line)
				if sample_name_search:
					sample_name = sample_name_search.group(1).strip()
					if re.search(r"\['(EPI_ISL_\d+)'\]", line):
						accession_search = re.search(r"\['(EPI_ISL_\d+)'\]", line)
						accession = accession_search.group(1) if accession_search else ""
						gisaid_isolate_log.append({"gs-sample_name":sample_name, "gisaid_accession_epi_isl_id":accession, "gisaid_message":"Already exists in GISAID"})
					elif re.search(r"\['(EPI_\d+)'\]", line):
						accession_search = re.search(r"\['(EPI_\d+)'\]", line)
						accession = accession_search.group(1) if accession_search else ""
						gisaid_segment_log.append({"gs-segment_name":sample_name, "gisaid_accession_epi_id":accession, "gisaid_message":"Already exists in GISAID"})
			# Handling if submitting samples have already been registered in GISAID (legacy JSON format)
			elif re.search(r'"code":\s*"validation_error".*?already exists;\s*existing_virus_name:', line):
				sample_name_search = re.search(r"(hCoV[^;]+);", line)
				if sample_name_search:
					sample_name = sample_name_search.group(1)
					if re.search(r"\['(EPI_ISL_\d+)'\]", line):
						accession_search = re.search(r"\['(EPI_ISL_\d+)'\]", line)
						if accession_search:
							accession = accession_search.group(1)
						else:
							accession = ""
						gisaid_isolate_log.append({"gs-sample_name":sample_name, "gisaid_accession_epi_isl_id":accession, "gisaid_message":"Already exists in GISAID"})
					elif re.search(r"\['(EPI_\d+)'\]", line):
						accession_search = re.search(r"\['(EPI_\d+)'\]", line)
						if accession_search:
							accession = accession_search.group(1)
						else:
							accession = ""
						gisaid_segment_log.append({"gs-segment_name":sample_name, "gisaid_accession_epi_id":accession, "gisaid_message":"Already exists in GISAID"})
			else:
				print("Finished reading GISAID log. If workflow has failed here, it's likely no GISAID IDs were returned. Check results in GISAID upload log.")
			line = file.readline().strip()
	gisaid_isolate_df = pd.DataFrame(gisaid_isolate_log)
	gisaid_segment_df = pd.DataFrame(gisaid_segment_log)
	# Update GISAID submission status
	if not gisaid_isolate_df.empty and not gisaid_segment_df.empty:
		print("GISAID isolates and GISAID segments found.")
		upload_log.update_submission_status_csv(submission_dir=submission_dir, update_database="GISAID", update_df=gisaid_isolate_df)
		upload_log.update_submission_status_csv(submission_dir=submission_dir, update_database="GISAID", update_df=gisaid_segment_df)
	elif not gisaid_isolate_df.empty:
		print("GISAID isolates found.")
		upload_log.update_submission_status_csv(submission_dir=submission_dir, update_database="GISAID", update_df=gisaid_isolate_df)
	elif not gisaid_segment_df.empty:
		print("GISAID segments found.")
		upload_log.update_submission_status_csv(submission_dir=submission_dir, update_database="GISAID", update_df=gisaid_segment_df)
	else:
		print("Warning: no GISAID isolates or segments found")
	if not gisaid_isolate_df.empty:
		gisaid_isolate_df = gisaid_isolate_df[~gisaid_isolate_df["gisaid_accession_epi_isl_id"].str.contains(r"EPI_ISL_\d*", regex = True, na = False)].copy()
	return gisaid_isolate_df[["gs-sample_name"]] if not gisaid_isolate_df.empty else pd.DataFrame(columns=["gs-sample_name"])

# Submit to GISAID
def submit_gisaid(organism: str, submission_dir: str, submission_name: str, config_dict: Dict[str, Any], submission_type: str) -> str:
	# Gather all required files
	metadata = os.path.join(submission_dir, "metadata.csv")
	orig_metadata = os.path.join(submission_dir, "orig_metadata.csv")
	fasta = os.path.join(submission_dir, "sequence.fsa")
	orig_fasta = os.path.join(submission_dir, "orig_sequence.fsa")
	submission_status_file = os.path.join(os.path.dirname(submission_dir), "submission_status_report.csv")
	# Extract user credentials (e.g. username, password, client-id)
	tools.check_credentials(config_dict=config_dict, database="GISAID")
	gisaid_cli = file_handler.validate_gisaid_installer(submission_dir=submission_dir, organism=organism, config_dict=config_dict)
	print(f"Uploading sample files to GISAID-{organism}, as a '{submission_type}' submission. If this is not intended, interrupt immediately.", file=sys.stdout)
	time.sleep(5)
	# Set number of attempt to 3 if erroring out occurs
	attempts = 0
	# Submit to GISAID
	while attempts <= 3:
		attempts += 1
		print("\n"+"Submission attempt: " + str(attempts), file=sys.stdout)
		# Create a log submission for each attempt
		log_file = os.path.join(submission_dir, "gisaid_upload_log_" + str(attempts) + ".txt")
		# If log file exists, removes it
		if os.path.isfile(log_file) == True:
			os.remove(log_file)
		# Upload submission
		command = subprocess.run([gisaid_cli, "upload", "--username", config_dict["Username"], "--password", config_dict["Password"], "--clientid", config_dict["Client-Id"], "--metadata", metadata, "--fasta", fasta, "--log", log_file, "--debug"],
			cwd=submission_dir, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
		# Check if uploading is successful
		if command.returncode != 0:
			print("Error: upload command error", file=sys.stderr)
			print(command.stdout)
			print(command.stderr)
			sys.exit(1)
		# Check if log file exists
		while not os.path.exists(log_file):
			time.sleep(10)
		# Check submission log to see if all samples are uploaded successfully
		process_gisaid_log(log_file=log_file, submission_dir=submission_dir)
		# Read in the submission status report
		status_df = pd.read_csv(submission_status_file, header = 0, dtype = str, engine = "python", encoding="utf-8", index_col=False)
		# Gather all required files
		metadata = os.path.join(submission_dir, "metadata.csv")
		fasta = os.path.join(submission_dir, "sequence.fsa")
		# Filter out samples with accession
		if "FLU" in organism:
			metadata_column_name = "Isolate_Name"
			fasta_column_name = "gs-segment_name"
			gisaid_status_df = status_df[~status_df["gisaid_accession_epi_isl_id"].str.contains("EPI_ISL_", na=False) & ~status_df["gisaid_accession_epi_id"].str.contains("EPI", na=False)].copy()
			gisaid_status_df = gisaid_status_df[["gs-sample_name", "gs-segment_name"]]
		elif "COV" in organism:
			metadata_column_name = "covv_virus_name"
			fasta_column_name = "gs-sample_name"
			gisaid_status_df = status_df[~status_df["gisaid_accession_epi_isl_id"].str.contains("EPI_ISL_", na=False)].copy()
			gisaid_status_df = gisaid_status_df[["gs-sample_name"]]
		else:
			metadata_column_name = organism.lower() + "_virus_name"
			fasta_column_name = "gs-sample_name"
			gisaid_status_df = status_df[~status_df["gisaid_accession_epi_isl_id"].str.contains("EPI_ISL_", na=False)].copy()
			gisaid_status_df = gisaid_status_df[["gs-sample_name"]]
		# Identify remaining samples
		metadata_df = pd.read_csv(orig_metadata, header = 0, dtype = str, engine = "python", encoding="utf-8", index_col=False)
		metadata_df = metadata_df.merge(gisaid_status_df, how="inner", left_on=metadata_column_name, right_on="gs-sample_name")
		if metadata_df.empty:
			print("Uploading successfully", file=sys.stdout)
			print("Log file is stored at: " + submission_dir + "/gisaid_upload_log_attempt_" + str(attempts) +  ".txt", file=sys.stdout)
			return "PROCESSED"
		# Update metadata file
		fasta_names = gisaid_status_df[fasta_column_name].tolist()
		metadata_df = metadata_df.drop(columns=["gs-sample_name", "gs-segment_name"], errors="ignore")
		metadata_df.to_csv(metadata, header = True, index = False)
		# Update fasta file
		fasta_dict = []
		with open(orig_fasta, "r") as fsa:
			records = SeqIO.parse(fsa, "fasta")
			for record in records:
				if record.id in fasta_names:
					fasta_dict.append(record)
		with open(fasta, "w+") as fasta_file:
			SeqIO.write(fasta_dict, fasta_file, "fasta")
	if not metadata_df.empty:
		print("Error: " + str(len(metadata_df.index)) + " sample(s) failed to upload to GISAID", file=sys.stderr)
		print("Please check log file at: " + submission_dir + "/gisaid_upload_log_attempt_{1,2,3}.txt", file=sys.stderr)
		return "ERROR"
	else:
		return "PROCESSED"

def update_gisaid_files(organism: str, submission_dir: str, submission_status_file: str) -> None:
	# Read in the submission status report
	status_df = pd.read_csv(submission_status_file, header = 0, dtype = str, engine = "python", encoding="utf-8", index_col=False)
	# Gather all required files
	metadata = os.path.join(submission_dir, "metadata.csv")
	orig_metadata = os.path.join(submission_dir, "orig_metadata.csv")
	fasta = os.path.join(submission_dir, "sequence.fsa")
	orig_fasta = os.path.join(submission_dir, "orig_sequence.fsa")
	# Filter out genbank that has accession number
	genbank_status_df = status_df[status_df["genbank-status"].str.contains("processed-ok", na=False)].copy()
	# Add required gisaid fields
	metadata_df = pd.read_csv(metadata, header = 0, dtype = str, engine = "python", encoding="utf-8", index_col=False)
	if "FLU" in organism:
		metadata_column_name = "Isolate_Name"
		fasta_column_name = "gs-segment_name"
		gisaid_status_df = genbank_status_df[["gs-sample_name", "gs-segment_name"]]
	elif "COV" in organism:
		metadata_column_name = "covv_virus_name"
		fasta_column_name = "gs-sample_name"
		gisaid_status_df = genbank_status_df[["gs-sample_name"]]
	else:
		metadata_column_name = organism.lower() + "_virus_name"
		fasta_column_name = "gs-sample_name"
		gisaid_status_df = genbank_status_df[["gs-sample_name"]]
	metadata_df = metadata_df.merge(gisaid_status_df, how="inner", left_on=metadata_column_name, right_on="gs-sample_name")
	fasta_names = metadata_df[fasta_column_name].tolist()
	metadata_df = metadata_df.drop(columns=["gs-sample_name", "gs-segment_name"], errors="ignore")
	metadata_df.to_csv(orig_metadata, header = True, index = False)
	fasta_dict = []
	with open(orig_fasta, "r") as fsa:
		records = SeqIO.parse(fsa, "fasta")
		for record in records:
			if record.id in fasta_names:
				fasta_dict.append(record)
	with open(fasta, "w+") as fasta_file:
		SeqIO.write(fasta_dict, fasta_file, "fasta")
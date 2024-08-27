#!/usr/bin/env python3

###########################    Description    ##################################
# Functions for GISAID input/output handling
################################################################################

import shutil
import subprocess
from typing import Dict, Any, List
import os
import pandas as pd
import file_handler
import sys
import time
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
import re
import json

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
	if organism in ["COV", "POX", "ARBO"]:
		gisaid_df = gisaid_df.rename(columns = {"sample_name": "virus_name"})
		if organism in ["POX", "ARBO"]:
			gisaid_df = gisaid_df.rename(columns = {"authors": "Authors", "collection_date": "Collection_Date"})
		if organism == "COV":
			prefix_name = "covv_"
		else:
			prefix_name = organism.lower() + "_"
		gisaid_df = gisaid_df.add_prefix(prefix_name)
		gisaid_df["submitter"] = config_dict["Username"]
		# Copy FASTA headers from sequence_name column to fn column
		gisaid_df["fn"] = gisaid_df["covv_sequence_name"]
		# Dropping the original column so it doesn't cause issues later.
		gisaid_df.drop(columns=['covv_sequence_name'], inplace=True)
		first_cols = ["submitter", "fn", (prefix_name + "virus_name")]
	elif "FLU" in organism:
		gisaid_df = gisaid_df.rename(columns = {"authors": "Authors", "collection_date": "Collection_Date"})
		gisaid_df["Isolate_Id"] = ""
		gisaid_df["Segment_Ids"] = ""
		# Pivot FLU segment names from long form to wide form
		gisaid_df["segment"] = "Seq_Id (" + gisaid_df["segment"].astype(str) + ")"
		group_df = gisaid_df.pivot(index="Isolate_Name", columns="segment", values="sample_name").reset_index()
		gisaid_df = gisaid_df.drop(columns=["sample_name", "segment"])
		gisaid_df = gisaid_df.drop_duplicates(keep="first")
		gisaid_df = gisaid_df.merge(group_df, on="Isolate_Name", how="inner", validate="1:1")
		first_cols = ["Isolate_Id","Segment_Ids","Isolate_Name"]
	# Restructure column order
	last_cols = [col for col in gisaid_df.columns if col not in first_cols]
	gisaid_df = gisaid_df[first_cols + last_cols]
	# If a submission fails and a new one is initiated without cleaning up metadata.csv/orig_metadata.csv, the column prefix doubles up.
	# This works for now, might set up some housekeeping in the future
	gisaid_df.columns = gisaid_df.columns.str.replace('^covv_covv_', 'covv_', regex=True)
	# Create submission files
	file_handler.save_csv(df=gisaid_df, file_path=submission_dir, file_name="metadata.csv")
	shutil.copy(os.path.join(submission_dir, "metadata.csv"), os.path.join(submission_dir, "orig_metadata.csv"))
	file_handler.create_fasta(database="GISAID", metadata=metadata, submission_dir=submission_dir)
	shutil.copy(os.path.join(submission_dir, "sequence.fsa"), os.path.join(submission_dir, "orig_sequence.fsa"))

# Read output log from gisaid submission script
def process_gisaid_log(log_file: str, submission_dir: str) -> pd.DataFrame:
	file_handler.validate_file(file_type="GISAID log", file_path=log_file)
	# Read in log file
	gisaid_isolate_log = []
	gisaid_segment_log = []
	with open(log_file, "r") as file:
		line = file.readline().strip()
		while line:
			line2 = json.loads(line)
			# If accession generated record it
			# Pattern options: "msg:": "<Sample Name>; <EPI_ISL/EPI_ID>_<Accession Numbers>" OR <epi_isl_id/epi_id>: <Sample Name>; <EPI_ISL/EPI_ID>_<Accession Numbers>
			# Changed re.match to re.search to return tokens that aren't at the beginning. Changed the regex to capture an arbitrary amount of numbers after EPI_
			if re.search("(?i)(\W|^)(\"msg\":\s*\"\S+.*;\s*(EPI_ISL|EPI_ID)_\d*\"|(epi_id|epi_isl_id):\s*\S.*;\s*(EPI_ISL_|EPI)\d+)(\W|$)", line):
				gisaid_string = re.findall(r'(?:[a-zA-Z0-9_-]+(?:/[a-zA-Z0-9_-]+)+|EPI_\w*)', line)
				gisaid_string = ' '.join(gisaid_string)
				gisaid_string_list: List[str] = gisaid_string.split(' ')
				sample_name = gisaid_string_list[0].strip()
				accession = gisaid_string_list[1].strip()
				if re.match("EPI_ISL_\d+", accession):
					gisaid_isolate_log.append({"gs-sample_name":sample_name, "gisaid_accession_epi_isl_id":accession})
				elif re.match("EPI\d+", accession):
					gisaid_segment_log.append({"gs-segment_name":sample_name, "gisaid_accession_epi_id":accession})
			# handling if submitting samples that have already been registered in GISAID
			elif re.search(r'"code":\s*"validation_error".*?already exists;\s*existing_virus_name:', line):
				sample_name = re.search(r"(hCoV[^;]+);", line)
				sample_name = sample_name.group(1) if sample_name else None
				if re.search(r"\['(EPI_ISL_\d+)'\]", line):
					accession =  re.search(r"\['(EPI_ISL_\d+)'\]", line)
					accession = accession.group(1)
					gisaid_isolate_log.append({"gs-sample_name":sample_name, "gisaid_accession_epi_isl_id":accession})
				elif re.search(r"\['(EPI_\d+)'\]", line):
					accession = re.search(r"\['(EPI_\d+)'\]", line)
					accession = accession.group(1)
					gisaid_segment_log.append({"gs-segment_name":sample_name, "gisaid_accession_epi_id":accession})
			else:
				print("Finished reading GISAID log. If workflow has failed here, it's likely no GISAID IDs were returned. Check results in GISAID upload log.")
			line = file.readline().strip()
	gisaid_isolate_df = pd.DataFrame(gisaid_isolate_log)
	gisaid_segment_df = pd.DataFrame(gisaid_segment_log)
	
	# Avoid passing an empty data frame to update_submission_status_csv since it was causing issues.
	# Output types of GISAID entities found
	if not gisaid_isolate_df.empty and not gisaid_segment_df.empty:
		print("GISAID isolates and GISAID segments found.")
		upload_log.update_submission_status_csv(submission_dir=submission_dir.replace("/GISAID", "/"), update_database="GISAID", update_df=gisaid_isolate_df)
		upload_log.update_submission_status_csv(submission_dir=submission_dir.replace("/GISAID", "/"), update_database="GISAID", update_df=gisaid_segment_df)
	elif not gisaid_isolate_df.empty:
		print("GISAID isolates found.")
		upload_log.update_submission_status_csv(submission_dir=submission_dir.replace("/GISAID", "/"), update_database="GISAID", update_df=gisaid_isolate_df)
	elif not gisaid_segment_df.empty:
		print("GISAID segments found.")
		upload_log.update_submission_status_csv(submission_dir=submission_dir.replace("/GISAID", "/"), update_database="GISAID", update_df=gisaid_segment_df)
	else:
		print("Warning: no GISAID isolates or segments found")
	gisaid_isolate_df = gisaid_isolate_df[~gisaid_isolate_df["gisaid_accession_epi_isl_id"].str.contains("EPI_ISL_\d*", regex = True, na = False)].copy()
	return gisaid_isolate_df[["gs-sample_name"]]

# Submit to GISAID
def submit_gisaid(organism: str, submission_dir: str, submission_name: str, config_dict: Dict[str, Any], submission_type: str) -> str:
	# Gather all required files
	metadata = os.path.join(submission_dir, "metadata.csv")
	orig_metadata = os.path.join(submission_dir, "orig_metadata.csv")
	fasta = os.path.join(submission_dir, "sequence.fsa")
	orig_fasta = os.path.join(submission_dir, "orig_sequence.fsa")
	# Extract user credentials (e.g. username, password, client-id)
	tools.check_credentials(config_dict=config_dict, database="GISAID")
	gisaid_cli = file_handler.validate_gisaid_installer(submission_dir=submission_dir, organism=organism)
	print(f"Uploading sample files to GISAID-{organism}, as a '{submission_type}' submission. If this is not intended, interrupt immediately.", file=sys.stdout)
	time.sleep(5)
	# Set number of attempt to 3 if erroring out occurs
	attempts = 1
	# Submit to GISAID
	while attempts <= 3:
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
		# Check if not_submitted_df has contents before trying to process it
		not_submitted_df = process_gisaid_log(log_file=log_file, submission_dir=submission_dir)
		# If submission completed, no more attempts
		if not_submitted_df.empty:
			print("Uploading successfully", file=sys.stdout)
			print("Log file is stored at: " + submission_dir + "/gisaid_upload_log_attempt_" + str(attempts) +  ".txt", file=sys.stdout)
			return "PROCESSED"
		else:
			# If submission is not completed, try again
			metadata_df = pd.read_csv(metadata, header = 0, dtype = str, engine = "python", encoding="utf-8", index_col=False)
			if "FLU" in organism:
				column_name = "Isolate_Name"
			elif "COV" in organism:
				column_name = "covv_virus_name"
			metadata_df = metadata_df.merge(not_submitted_df, how="inner", left_on=column_name, right_on="gs-sample_name")
			fasta_names = metadata_df["sequence_name"].tolist()
			if "sequence_name" in metadata_df.columns:
				metadata_df.rename(columns={"sequence_name": "fn"}, inplace=True)
			metadata_df = metadata_df.drop(columns=["gs-sample_name"])
			metadata_df = metadata_df.columns.str.replace('^covv_covv_', 'covv_', regex=True)
			metadata_df.to_csv(orig_metadata, header = True, index = False)
			fasta_dict = []
			with open(orig_fasta, "r") as fsa:
				records = SeqIO.parse(fsa, "fasta")
				for record in records:
					if record.id in fasta_names:
						fasta_dict.append(record)
			with open(fasta, "w+") as fasta_file:
				SeqIO.write(fasta_dict, fasta_file, "fasta")
			attempts += 1
	if not not_submitted_df.empty:
		print("Error: " + str(len(not_submitted_df.index)) + " sample(s) failed to upload to GISAID", file=sys.stderr)
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
	gisaid_status_df = genbank_status_df[["gs-sample_name", "sequence_name"]]
	# Add required gisaid fields
	metadata_df = pd.read_csv(metadata, header = 0, dtype = str, engine = "python", encoding="utf-8", index_col=False)
	if "FLU" in organism:
		column_name = "Isolate_Name"
	elif "COV" in organism:
		column_name = "virus_name"
	metadata_df = metadata_df.merge(gisaid_status_df, how="inner", left_on=column_name, right_on="gs-sample_name")
	fasta_names = metadata_df["sequence_name"].tolist()
	metadata_df = metadata_df.drop(columns=["gs-sample_name", "sequence_name"])
	metadata_df.to_csv(orig_metadata, header = True, index = False)
	fasta_dict = []
	with open(orig_fasta, "r") as fsa:
		records = SeqIO.parse(fsa, "fasta")
		for record in records:
			if record.id in fasta_names:
				fasta_dict.append(record)
	with open(fasta, "w+") as fasta_file:
		SeqIO.write(fasta_dict, fasta_file, "fasta")

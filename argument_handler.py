#!/usr/bin/env python3

###########################    Description    ##################################
# Parsers for handling SeqSender input
################################################################################

import argparse
from typing import List
from settings import ORGANISM_CHOICES

def args_parser():
	"""
	Argument parser setup and build.
	"""
	parser = argparse.ArgumentParser(formatter_class=argparse.ArgumentDefaultsHelpFormatter,
									description="Genomic tool to simplify/automate the process of submitting organism samples to public repositories. With built-in tools to create/submit/link/log organism samples for the databases: BioSample, SRA, GenBank, and GISAID.")
	database_parser = argparse.ArgumentParser(add_help=False)
	organism_parser = argparse.ArgumentParser(add_help=False)
	validate_parser = argparse.ArgumentParser(add_help=False)
	submission_name_parser = argparse.ArgumentParser(add_help=False)
	submission_dir_parser = argparse.ArgumentParser(add_help=False)
	upload_log_submission_name_parser = argparse.ArgumentParser(add_help=False)
	config_file_parser = argparse.ArgumentParser(add_help=False)
	file_parser = argparse.ArgumentParser(add_help=False)
	test_parser = argparse.ArgumentParser(add_help=False)

	database_parser.add_argument("--biosample", "-b",
		help="Create/Submit BioSample data.",
		action="store_const",
		const="BIOSAMPLE",
		default="")
	database_parser.add_argument("--sra", "-s",
		help="Create/Submit SRA data.",
		action="store_const",
		const="SRA",
		default="")
	database_parser.add_argument("--genbank", "-n",
		help="Create/Submit GenBank data. (requires --fasta_file)",
		action="store_const",
		const="GENBANK",
		default="")
	database_parser.add_argument("--gisaid", "-g",
		help="Create/Submit GISAID data. (requires --fasta_file)",
		action="store_const",
		const="GISAID",
		default="")
	organism_parser.add_argument("--organism",
		help="Type of organism data. Listed organism options have unique submissions options/processes, if your specific organism is not listed, use 'OTHER' for options available to all organisms.",
		choices=ORGANISM_CHOICES,
		default="",
		required=True)
	validate_parser.add_argument("--skip_validation",
		help="Skip initial validation for metadata file. Validation will still occur for the 'config_file' and for any subsequent submissions made via 'submission_status'. Warning, this can cause unexpected errors using SeqSender if required columns are missing.",
		required=False,
		action="store_const",
		default=False,
		const=True)
	submission_name_parser.add_argument("--submission_name",
		help="Unique name for the submission of your data. Reusing the same name can cause issues during the submission process. A folder will be created at: 'submission_dir/submission_name'.",
		required=True)
	upload_log_submission_name_parser.add_argument("--submission_name",
		help="Unique name for the submission of your data. This is an optional field if you want Seqsender to only update the specified submission in the 'submission_log.csv'.",
		required=False)
	submission_dir_parser.add_argument("--submission_dir",
		help="Output directory where all files for your submission will be stored. A folder will be created at '<submission_dir>/<submission_name>'; this is the location where: all of the submission files will be created, SeqSender will stage each step of the submission process automatically, and where SeqSender will generate all the output from your submission.",
		required=True)
	config_file_parser.add_argument("--config_file",
		help="Config file to be used in the creation/submission of your samples. SeqSender will store this file location in your 'submission_log.csv' where it will use it to manage your submission, be careful when modifying and ensure SeqSender maintains access to this file. Input either full file path or if just file name it must be stored at '<submission_dir>/<submission_name>/<config_file>'.",
		required=True)
	file_parser.add_argument("--metadata_file",
		help="Metadata file to be used in the creation/submission of your samples. Input either full file path or if just file name it must be stored at '<submission_dir>/<submission_name>/<metadata_file>'.",
		required=True)
	file_parser.add_argument("--fasta_file",
		help="Fasta file used to generate submission files; fasta header should match the column 'sequence_name' stored in your metadata. Input either full file path or if just file name it must be stored at '<submission_dir>/<submission_name>/<fasta_file>'.",
		default = None)
	file_parser.add_argument("--table2asn",
		help="Perform a table2asn submission instead of GenBank FTP submission for organism choices 'FLU' or 'COV'.",
		required=False,
		action="store_const",
		default=False,
		const=True)
	file_parser.add_argument("--gff_file",
		help="Annotation file only available for table2asn submissions. (requires '--table2asn' for organism choices 'FLU', or 'COV').",
		default=None)
	test_parser.add_argument("--test",
		help="Perform a test submission.",
		action="store_const",
		default=False,
		const=True)

	# Create the submodule commands
	subparser_modules = parser.add_subparsers(dest="command")

	# prep command
	prep_module = subparser_modules.add_parser(
		"prep",
		formatter_class=argparse.ArgumentDefaultsHelpFormatter,
		description="Generate all files required to submit to databases selected.",
		parents=[database_parser, organism_parser, submission_name_parser, submission_dir_parser, config_file_parser, file_parser, validate_parser]
	)

	# submit command
	submit_module = subparser_modules.add_parser(
		"submit",
		formatter_class=argparse.ArgumentDefaultsHelpFormatter,
		description="Generate all files required and begin the submission process to databases selected.",
		parents=[database_parser, organism_parser, submission_name_parser, submission_dir_parser, config_file_parser, file_parser, test_parser, validate_parser]
	)

	# check_submission_status command
	update_module = subparser_modules.add_parser(
		"submission_status",
		formatter_class=argparse.ArgumentDefaultsHelpFormatter,
		description="Checks the submission status for (all/specified <submission_name>) submission('s) which will for each database: update the status of the submission('s), download output file('s), submit to subsequent specified databases if linking information requires output of previous database('s).",
		parents=[submission_dir_parser, upload_log_submission_name_parser]
	)

	# Generate test data command
	test_output_module = subparser_modules.add_parser(
		"test_data",
		formatter_class=argparse.ArgumentDefaultsHelpFormatter,
		description="Returns a set of test data examples (config_file, metadata_file, fasta_file, etc.) to the specified 'submission_dir' based on organism and database selections.",
		parents=[database_parser, organism_parser, submission_dir_parser]
	)

	# biosample xml download command
	biosample_xml_module = subparser_modules.add_parser(
		"update_biosample",
		formatter_class=argparse.ArgumentDefaultsHelpFormatter,
		description="Downloads the BioSample Package XML from NCBI and updates SeqSender's metadata schema options for the BioSample database."
	)

	# version command
	version_module = subparser_modules.add_parser(
		"version",
		formatter_class=argparse.ArgumentDefaultsHelpFormatter,
		description="Print version."
	)

	return parser

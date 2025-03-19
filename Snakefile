"""
MinION Metabarcoding Pipeline
Schyler Nunziata
"""
shell.executable('bash')
import os
import glob
import pandas as pd
from shutil import copyfile
from pathlib import Path

#Get configuration parameters from config file
configfile: "config.yaml"
R1_extension = config['R1_extension']
directory_name=config["directory_name"]
samples_file = Path(config["samples"])

#Read sample sheet if it exists, if not create from file names in data directory
if samples_file.is_file():
	#Get the sample list
	samples_df = pd.read_table(config["samples"]).set_index("sample", drop=False)
	samples_df.index = samples_df.index.map(str)
	barcodes = list(samples_df['sample'])
else:
	input_directory = os.getcwd() + "/data/"
	filenames = glob.glob(input_directory + "*" + R1_extension)
	barcodes = [os.path.split(this)[1].replace(R1_extension, "")
							for this in filenames]

#Run rules
localrules: all, align_ref_clusters, format_noscaffold_consensus_seqs, blast_consensus_seqs, summarize_blast_output, count_mapped_reads, blastn_excel_summary, blastn_excel_summary_all, summarize_read_mapping

rule all:
	input:
		expand("summary_output/{run_name}_readqc_report.html", run_name=directory_name),
		expand("summary_output/{run_name}_mapping_summary.xlsx", run_name=directory_name)

'''
onsuccess:
        print("Workflow finished, no error")
        shell("mail -s 'VirusDetection workflow finished, no error' {config[email][email]} < {log}")

onerror:
        print("An error occurred")
        shell("mail -s 'VirusDetection workflow error occurred' {config[email][email]} < {log}")
'''


##### load rules #####
include: "rules/QC_trim.smk"
include: "rules/mapping.smk"
include: "rules/consensus_seq_cluster_based.smk"
include: "rules/create_db.smk"

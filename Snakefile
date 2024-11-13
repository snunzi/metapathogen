"""
MinION Metabarcoding Pipeline
Schyler Nunziata
"""
shell.executable('bash')
import os
import glob
import pandas as pd
from shutil import copyfile

configfile: "config.yaml"

R1_extension = config['R1_extension']
directory_name=config["directory_name"]

#Get the sample list
samples_df = pd.read_table(config["samples"]).set_index("sample", drop=False)
samples_df.index = samples_df.index.map(str)
barcodes = list(samples_df['sample'])




rule all:
	input:
		#expand("demultiplexed/{barcode}_raw.fastq.gz", barcode=barcodes)
		expand("summary_output/{run_name}_readqc_report.html", run_name=directory_name),
		expand("summary_output/{run_name}_mapping_summary.xlsx", run_name=directory_name)


""""
onsuccess:
	print("Workflow finished, no error")
	shell("mail -s 'VirusDetection workflow finished, no error' {config[email][email]} < {log}")

onerror:
	print("An error occurred")
	shell("mail -s 'VirusDetection workflow error occurred' {config[email][email]} < {log}")
"""


##### load rules #####
include: "rules/QC_trim.smk"
include: "rules/mapping.smk"
include: "rules/consensus_seq_cluster_blased.smk"

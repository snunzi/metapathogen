#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Jan 16 14:06:56 2024

@author: sonunziata

"""

import sys
import pandas as pd
from Bio.SeqIO.FastaIO import SimpleFastaParser
from functools import reduce
from Bio import SeqIO
from itertools import groupby
import re
import numpy as np

#sys.argv = ['script.py', 'barcode02_COI_blast.out', 'barcode02_read_number.txt', 'barcode02_consensus.fasta', 'barcode02.fasta.clstr', 'barcode02', 'output.csv']

blast_file = sys.argv[1]
cluster1_read_file = sys.argv[2]
final_consensus_file = sys.argv[3]
cluster_file = sys.argv[4]
bar_name = sys.argv[5]
out_file = sys.argv[6]

#Open blast output and delete the _ref in cluster name
blastn_out = pd.read_csv(blast_file, delimiter= "\t", names = ['Cluster_Name', 'pident', 'Alignment_Length', 'mismatch', 'gapopen', 'evalue', 'bitscore', 'Ref_Seqeunce'])
blastn_out['Cluster_Name'] = blastn_out.Cluster_Name.str.replace('_ref,?' , '')

#Open read counts per cluster
cluster1_read_counts = pd.read_csv(cluster1_read_file, delimiter=":", names = ['Cluster_Name', 'Number Reads'])

cluster1_read_counts["Cluster2 Name"] = np.nan

#Get recluster info
with open(cluster_file, 'r') as inputFile:
	for line in inputFile:	# read in each line
		if ">Cluster " in line.strip():	 # grab number of amplier
			clus = line.replace(">Cluster ", "").strip()
		else:	# grab sequence identifier
			id = re.findall(r'\>([^]]*)\_ref', line)
			cluster1_read_counts['Cluster2 Name']=(cluster1_read_counts['Cluster2 Name'].mask(cluster1_read_counts['Cluster_Name'].eq(id[0]), clus))

#Open consensus seqeunce file
with open(sys.argv[3]) as fasta_file:
    ids = []
    seqs = []
    for title, sequence in SimpleFastaParser(fasta_file):
        ids.append(title)  # First word is ID
        seqs.append(str(sequence))

#Create pandas df from contig sequences
df = pd.DataFrame(list(zip(ids, seqs)),
               columns =['Cluster_Name', 'Sequence'])
df['Cluster_Name'] = df.Cluster_Name.str.replace('_ref,?' , '')

#Merge all dataframes on cluster name
data_frames = [blastn_out, cluster1_read_counts, df]
summary_out = reduce(lambda  left,right: pd.merge(left,right,on=['Cluster_Name'],
                                            how='outer'), data_frames)

#summary_out = summary_out[summary_out['pident'].notna()]
summary_out['Cluster_Name'] = summary_out['Cluster_Name'].str.split(pat='/').str[2]
summary_out['Barcode'] = bar_name

summary_out = summary_out[['Barcode', 'Cluster_Name', 'Cluster2 Name', 'pident', 'Alignment_Length', 'mismatch', 'gapopen', 'evalue', 'bitscore', 'Ref_Seqeunce', 'Number Reads', 'Sequence']]

#Output summary file
summary_out.to_csv(out_file, index=False)

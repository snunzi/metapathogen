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

#sys.argv = ['script.py', 'barcode04_blast.out', 'barcode04_read_number.txt', 'barcode04_consensus.fasta', 'barcode04', 'output.csv', 'output_filtered.csv']

blast_file = sys.argv[1]
cluster1_read_file = sys.argv[2]
final_consensus_file = sys.argv[3]
bar_name = sys.argv[4]
out_file = sys.argv[5]
out_file_filtered = sys.argv[6]

#Open blast output and delete the _ref in cluster name
blastn_out = pd.read_csv(blast_file, delimiter= "\t", names = ['Cluster_Name', 'pident', 'Acc_Length', 'Query_Length', 'Alignment_Length', 'Query Cov', 'mismatch', 'gapopen', 'evalue', 'bitscore', 'Ref_Sequence'])
blastn_out['Cluster_Name'] = blastn_out.Cluster_Name.str.replace('_ref,?' , '')


#Open read counts per cluster
cluster1_read_counts = pd.read_csv(cluster1_read_file, delimiter=":", names = ['Cluster_Name', 'Number Reads'])
cluster1_read_counts['Cluster_Name'] = cluster1_read_counts.Cluster_Name.str.replace('.fa,?' , '')


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
#summary_out['Cluster_Name'] = summary_out['Cluster_Name'].str.split(pat='/').str[1]
summary_out['Barcode'] = bar_name


#Filter by parameters
#summary_out['Number Reads'] = summary_out[summary_out['Number Reads'] > int(5)]['Number Reads']
summary_out = summary_out.dropna()

summary_out['%_Ref_Cov'] = 100* (summary_out['Alignment_Length']/summary_out['Acc_Length'])

#Calculate Grade
#summary_out['Grade_Cov'] = 50 * summary_out['Query Cov']
#summary_out['Grade_evalue'] = 25 * (np.maximum(0,(1 - summary_out['evalue']/10**-20)))
#summary_out['Grade_pident'] = 25 * (np.maximum(0,(summary_out['pident'] - 50)/(100-50)))
#summary_out['Grade'] = 50 * summary_out['Query Cov'] + 25 * (np.maximum(0,(1 - summary_out['evalue']/10**-20))) + 25 * (np.maximum(0,(summary_out['pident'] - 50)/(100-50)))

summary_out = summary_out[['Barcode', 'Cluster_Name', 'Number Reads', 'Ref_Sequence', 'pident', 'Alignment_Length', 'Acc_Length', 'mismatch', 'gapopen', 'evalue', 'bitscore', '%_Ref_Cov', 'Sequence']]


#Output summary file
summary_out.to_csv(out_file, index=False)

#Filter by parameters
summary_out= summary_out[summary_out['%_Ref_Cov'] > 50]
summary_out = summary_out[summary_out['pident'] > 95]
summary_out = summary_out[summary_out['mismatch'] < 11 ]
summary_out = summary_out.dropna()

groups = summary_out.groupby(by=['Cluster_Name'])
summary_filtered = groups.apply(lambda g: g[g['pident'] == g['pident'].max()])
summary_filtered.to_csv(out_file_filtered, index=False)

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


#sys.argv = ['script.py', 'barcode23_blast.out', 'barcode23_scaffolded_read_count', 'barcode23_consensus.fasta', 'barcode23', 'output.csv', 'output_filtered.csv', 'erictest']

blast_file = sys.argv[1]
cluster1_read_file = sys.argv[2]
final_consensus_file = sys.argv[3]
bar_name = sys.argv[4]
out_file = sys.argv[5]
out_file_filtered = sys.argv[6]
run_name = sys.argv[7]
ref_cov_cutoff = float(sys.argv[8])
pident_cutoff = float(sys.argv[9])
mismatch_cutoff = float(sys.argv[10])
con_len_cutoff = float(sys.argv[11])

#Open blast output and delete the _ref in cluster name
blastn_out = pd.read_csv(blast_file, delimiter= "\t", names = ['Cluster_Name', 'pident', 'Acc_Length', 'Query_Length', 'Alignment_Length', 'Query Cov', 'mismatch', 'gapopen', 'evalue', 'bitscore', 'Ref_Sequence'])
blastn_out['Cluster_Name'] = blastn_out.Cluster_Name.str.replace('_ref,?' , '')


#Open read counts per cluster
cluster1_read_counts = pd.read_csv(cluster1_read_file, delimiter="\t", names = ['Cluster_Name', 'read_count'])

#Open consensus seqeunce file
with open(sys.argv[3]) as fasta_file:
    ids = []
    seqs = []
    length = []
    for title, sequence in SimpleFastaParser(fasta_file):
        ids.append(title)  # First word is ID
        seqs.append(str(sequence))
        length.append(len(sequence))

#Create pandas df from contig sequences
df = pd.DataFrame(list(zip(ids, seqs, length)),
               columns =['Cluster_Name', 'Sequence', 'Consensus_Length'])
df['Cluster_Name'] = df.Cluster_Name.str.replace('_ref,?' , '')

#Merge all dataframes on cluster name
data_frames = [blastn_out, cluster1_read_counts, df]
summary_out = reduce(lambda  left,right: pd.merge(left,right,on=['Cluster_Name'], how='outer'), data_frames)

summary_out['Barcode'] = bar_name


#Filter by parameters
summary_out['read_count'] = summary_out[summary_out['read_count'] > int(5)]['read_count']
summary_out = summary_out.dropna()

summary_out['%_Ref_Cov'] = 100* (summary_out['Alignment_Length']/summary_out['Acc_Length'])


summary_out['%_Ref_Cov'] = summary_out['%_Ref_Cov'].apply(lambda x: round(x, 2))
summary_out['pident'] = summary_out['pident'].apply(lambda x: round(x, 2))

summary_out = summary_out[['Barcode', 'Cluster_Name', 'read_count', 'Ref_Sequence', 'pident', 'Alignment_Length', 'Consensus_Length', 'Acc_Length', 'mismatch', 'gapopen', 'evalue', 'bitscore', '%_Ref_Cov', 'Sequence']]

summary_out.insert(0, "Run", run_name)
#Output summary file
summary_out.to_csv(out_file, index=False)

#Filter by parameters
summary_out= summary_out[summary_out['%_Ref_Cov'] > ref_cov_cutoff]
summary_out = summary_out[summary_out['pident'] > pident_cutoff]
summary_out = summary_out[summary_out['mismatch'] < mismatch_cutoff ]
summary_out = summary_out[summary_out['Consensus_Length'] < con_len_cutoff]
summary_out = summary_out.dropna()

groups = summary_out.groupby(by=['Cluster_Name'], as_index=False, sort=False)
summary_filtered = groups.apply(lambda g: g[g['pident'] == g['pident'].max()])
summary_filtered.to_csv(out_file_filtered, index=False)

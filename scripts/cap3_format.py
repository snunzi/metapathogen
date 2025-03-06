#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Dec 13 13:01:36 2024

@author: sonunziata
"""

from Bio.Sequencing import Ace
import pandas as pd
import sys

counts_file = sys.argv[1]
ace_file = sys.argv[2]
count_summary = sys.argv[3]

#Get number of reads used to build each contig
read_count = pd.read_csv(counts_file, names=['contig'])
#remove contigs with no reads
read_count = read_count[read_count['contig'].str.contains('_')]

#Get Read count
read_count['read_count']=[x.split('_')[1] for x in read_count['contig']]
read_count['read_count']=read_count['read_count'].astype(int)


#Get which contigs contributed to each scaffold
scaffolds = {}
for record in Ace.parse(ace_file):
    for i in range(len(record.reads)):
        scaffold_name = record.reads[i].rd.name
        scaffolds[scaffold_name] = record.name 
        #print(record.name," ",record.reads[i].rd.name)


read_count['Cluster_Name'] = read_count['contig'].map(scaffolds)
read_count['Cluster_Name'] = read_count['Cluster_Name'].fillna(read_count['contig'])

#Sum number of reads used to build each scffold
read_count = read_count.groupby(['Cluster_Name'])['read_count'].sum()

read_count.to_csv(count_summary, sep="\t", header=False)
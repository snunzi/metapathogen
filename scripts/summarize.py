#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Apr 26 13:29:47 2022

@author: sonunziata
"""
import pandas as pd
import sys
import os
from pathlib import Path

summary_all = pd.DataFrame()
read_all = pd.DataFrame()

raw_read_count = pd.read_csv(sys.argv[4], delimiter="\t", header=None, names=['Sample_Name', 'Raw_Reads'])

my_file = Path(sys.argv[2])
if os.path.exists(my_file) and os.path.getsize(my_file) > 0:
    blastn_summary = pd.read_csv(sys.argv[2])
else:
    blastn_summary = None

my_file = Path(sys.argv[3])
if os.path.exists(my_file) and os.path.getsize(my_file) > 0:
    blastn_summary_unfiltered = pd.read_csv(sys.argv[3])
else:
    blastn_summary_unfiltered = None


for file in sys.argv[5:]:
    file_name = file.replace('/','.')
    name = (file_name.split(".")[1])
    if not os.stat(file).st_size == 0:
        species_list = pd.read_csv(file, delimiter= "\t", header=0)
        if not species_list.empty:
            read_count = pd.read_csv("mapping/"+name+"_summary.txt", delimiter="\t", skiprows=2, nrows=1, header=None, names=['parameter', 'Cleaned_Reads'])
            species_list['Mapped_Reads'] = species_list['Plus_reads'] + species_list['Minus_reads']
            species_list = species_list.sort_values('Mapped_Reads', ascending=False)
            species_list = species_list.nlargest(5,'Mapped_Reads')
            species_list['Proportion_Mapped'] = round(100 * species_list['Mapped_Reads'].astype('float')/float(read_count['Cleaned_Reads'][0]),3)
            species_list['Sample_Name'] = name
            read_count['Sample_Name'] = name
            #summary_all = summary_all.append(species_list[['Sample_Name', 'ID', 'Mapped_Reads', 'Proportion_Mapped']])
            summary_all = pd.concat([summary_all, species_list[['Sample_Name', 'ID', 'Mapped_Reads', 'Proportion_Mapped']]], ignore_index=True)
            #read_all = read_all.append(read_count[['Sample_Name', 'Cleaned_Reads']])
            read_all = pd.concat([read_all, read_count[['Sample_Name', 'Cleaned_Reads']]], ignore_index=True)
            #print(name)

read_all_raw_trimmed = pd.merge(raw_read_count, read_all, on='Sample_Name')

# create excel writer
writer = pd.ExcelWriter(sys.argv[1])
# write dataframe to excel sheet
read_all_raw_trimmed.to_excel(writer, 'read_summary', index=False)
summary_all.to_excel(writer, 'mapping_summary', index=False)
if blastn_summary is not None:
	blastn_summary.to_excel(writer, 'blastn_summary', index=False)
if blastn_summary_unfiltered is not None:
	blastn_summary_unfiltered.to_excel(writer, 'blastn_unfiltered', index=False)

# save the excel file
#writer.save()
writer.close()

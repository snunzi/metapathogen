# -*- coding: utf-8 -*-
"""
Spyder Editor

This is a temporary script file.
"""

import sys
from Bio import SeqIO
import pandas as pd

#sys.argv = ['script', 'Phplcontrol_positive_control.fq.fa', 'ref.fa']


with open(sys.argv[1]) as fasta_file:  # Will close handle cleanly
    identifiers = []
    lengths = []
    for seq_record in SeqIO.parse(fasta_file, 'fasta'):  # (generator)
        identifiers.append(seq_record.id)
        lengths.append(len(seq_record.seq))

fasta_len = pd.DataFrame(list(zip(identifiers, lengths)),
              columns=['ID','length']).sort_values('length', ascending=False).reset_index(drop = True)
        


drop_longest_n = int(round(len(fasta_len) * 0.10))

ref_read = fasta_len.iloc[drop_longest_n]['ID']

                   
for record in SeqIO.parse(str(sys.argv[1]), "fasta"):
    if record.id in ref_read:
         SeqIO.write(record, sys.argv[2], "fasta")

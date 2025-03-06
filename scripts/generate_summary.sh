#!/bin/bash

list=("Run" "Barcode" "Cluster_Name" "Read Count" "Ref_Sequence" "pident" "Alignment_Length" "Consensus_Length" "Acc_Length" "mismatch" "gapopen" "evalue" "bitscore" "%_Ref_Cov" "Sequence")

printf "%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s\n" "${list[@]}"
# PPCDL MinION Metabarcoding Bioinformatic Pipeline


The purpose of this pipeline is for the analysis amplicon sequencing of marker genes from Nanopore reads. The outlined bioinformatics pipeline accepts raw sequencing reads and proceeds through a reproducible workflow of multiple processing steps from read filtering, generation of consensus sequences, and taxonomic classification. 

# Workflow

![schema of workflow](doc/minmetabar_scheme.png?raw=true)

## Tools in pipeline

1. Read preprocessing
	1. Read quality assessment
		1. Read quality before and after filtering (['nanoplot'](https://github.com/wdecoster/NanoPlot))
	2. Read trimming and filtering
		1. Split chimeric reads (['duplex_tools'](https://github.com/nanoporetech/duplex-tools))
		2. Remove adapters([`porechop`](https://github.com/rrwick/Porechop)) 
		3. Length and quality filter (['chopper'](https://github.com/wdecoster/chopper))
2. Read clustering
	1. Cluster reads by shared kmer with reference database ([`seal.sh`](https://github.com/BioInfoTools/BBMap/blob/master/sh/seal.sh))
	2. BP cutoff for making consensus sequence (['filtlong'](https://github.com/rrwick/Filtlong))
3. Build Consensus Sequence
	1. Align reads within clusters ([`minimap2`](https://github.com/lh3/minimap2)
	2. Build consensus resuence ([`racon`](https://github.com/isovic/racon))
	3. Polish consensus sequence ([`medaka`](https://github.com/nanoporetech/medaka))
4.	Taxonomic assignment
	1. Taxonomic classification (['BLASTn'](https://blast.ncbi.nlm.nih.gov/Blast.cgi))
5.	Summary consensus sequences and assignments in excel workbook (python script)


## Pipeline Description

This pipeline accepts raw reads in “fastq.gz” format. Chimeric reads are identified and split by duplex_tools (v0.2.17). Sequencing adapters are identified and trimmed using porechop (0.2.4) and resulting reads quality metrics and read counts are calculated using NanoPlot (v1.41.6). Reads are afterwards filtered on a minimum quality score of 10, with a read length filter between 500 and 1500 using chopper. Quality metrics and read counts of filtered reads are calculated using Nanoplot (v1.41.6) and summarized along with unfiltered reads using MultiQC (v1.17). Filtered reads are then clustered against a custom reference sequence database (custom gene DB) using Seal.sh (v38.76), which bins sequences based on which reference sequences sharing the most kmers with the query. Ambiguous sequences with an equal number of kmer matches to multiple reference sequences are assigned to all best-matching sequences and reads that which do not hit a sequence in the reference database are excluded from downstream analysis. Consensus sequences are created from each cluster larger than 5 reads, using the highest quality reads up to 750,000bp. First, the clustered reads are aligned, using the longest read (after excluding top 10% of longest reads) within the cluster as a reference with Minimap2 (v2.26). The alignments are used to build an initial draft consensus sequence of each cluster using Racon (v.1.5.0), which is then polished by Medaka (v.1.11.3). The resulting consensus sequence(s) are then classified using a BLASTn search (NCBI, version 2.15.0) against the custom gene DB. The 5 top BLAST hits are reported in an unfiltered summary report table. The top hits are then filtered at a minimum reference coverage of 50%, 95% identity, and ≤11 mismatches, and a single top hit with the highest percent identity is reported for each cluster passing filtering parameters in a filtered summary table report.

# Run the Pipeline

Reads in ".fastq.gz" format must be uploaded to the server following directions here. 

## Parameter Choices

1. Reference Database
	* Phytophthora
	* Phytoplasma
2. Read Filtering Parameters
	* Minimum Read Length
	* Maximum Read Length
	* Minimum Read Quality

# Results Analyses

a. Description of DADA2 output and interpretation

	1.	readqc_summary.html

	This file contains the read counts and quality metrics of the “raw reads” after chimera and adapter removal, and the “clean reads” after quality and length filtering.

	2.	mapping_summary.xlsx
	The main output of the minion_metabarcode pipeline is the mapping_summary.xlsx file. The first tab of this spreadsheet, ‘read_summary’, includes the total number of reads passing quality filtering for each barcode. The next tab, ‘blastn_summary’, contains the filtered blastn top hit for each consensus sequence. The first column in the spreadsheet identifies the barcode of the sample, and then identifies the reference sequence the reads were clustered with. Subsequent columns report the blast hit alignment statistics for the consensus sequence. Detailed information regarding these columns is explained in Table 1. 
	
	The final tab in the workbook, ‘blastn_unfiltered’, contains the same columns as the ‘blastn_summary’ but with no filtering and with 3 tops hits per consensus sequence reported. 

	The taxonomic annotations provided in the ‘blastn_summary’ spreadsheet are not meant to provide a final determination as to the identity of the consensus sequences in the sample, but as a means to organize and sort the sequence data for further analysis. 


Column | Description 
--- | ---
Barcode | Barcode name as assigned from the directory name in the data folder.
Cluster_Name | Reference sequence reads were clustered with.
Number Reads | Number of reads used to create consensus.
Ref_Sequence | Name of blast hit from the reference database. 
pident | Percentage of identical matches between reference and consensus.
Alignment_Length | Length of alignment.
Acc_Length | Length of reference.
mismatch | Number of mismatches.
gapopen | Number of gap openings.
evalue | Expect value.
bitscore | Bit score.
%_Ref_Cov | Percent of reference covered, as calculated by Alignment_length/Acc_Length*100.
Sequence | Nucleotide sequence of the consensus.

		
b. Sample assessment

	Occasionally, reads from a particular taxon may be mis-assigned to a sample after demultiplexing as the result of index-hopping. Therefore, caution should be used when determining the presence of species with low sequencing depth (< ~50 reads) to ensure that the OTUs detected in a sample correspond to true biological sequences and are not experimental artifacts. This is primarily of concern when multiplexed libraries include high-titer samples. A specific cut-off has not been established for this method and may depend on the pathogen, testing matrix, and sequencing depth, among other factors. The cut-off guidelines are given below as a rule of thumb and should be assessed on a case-by-case basis. Consult a subject matter expert in the case of inconclusive results.

	1.	Quality Control

	Before proceeding with analysis of test samples, examine the ‘readqc_summary.html’ output file. For a high-quality dataset, expect at least 50% of the initial reads to remain after quality filtering. If less than 50% of the initial reads passed quality filtering, then taxa present in the dataset may not be detected. Consult a bioinformatician as the pipeline parameters may be adjusted to optimize the number of reads that pass quality filtering prior to downstream analysis.

	2.	Sequence analysis

	In the ‘blastn_summary’ tab of the ‘mapping_summary.xlsx’ output, select the column headers in the spreadsheet. Then, under the “sort and filter” tab, click the filter button. Parse the spreadsheet to examine each barcoded sample individually. Begin by sorting the table by descending ‘pident’, and then look at the ‘Number Reads’ used to build the high identity sequences. Consensus sequences built from a larger number of reads (≥ 50) are more likely to contain true biological sequences. 

	Generally, a consensus sequence should display ≥ 97% identity and cover ≥ 80% of the reference sequence for species level identification. However, consensus sequences built with fewer than 50 reads, but with high sequence identity (≥ 97% identity) and reference coverage (≥ 98%), and low mismatches (≤ 5) should be evaluated further. The identity of any suspect sequences should subsequently be verified through BLASTn analysis of NCBI GenBank (https://blast.ncbi.nlm.nih.gov/Blast.cgi) or other sequence databases.

	The ‘blastn_summary’ tab of the ‘mapping_summary.xlsx’ output should also be assessed for missed taxonomic assignments. The spreadsheet should be filtered as above and sorted by ‘Number Reads’. Sequences not reported in ‘blastn_summary’ with ≥ 50 number of reads should be flagged for further analysis.

c. Phytophthora Sample Assessment

	•	POSITIVE: A sample was considered positive for a particular Phytophthora species if the pathogen consensus sequence should have a sequencing depth of approximately 50 reads, ≥ 97% identity with the reference match, cover ≥ 80% of the reference sequence and ≤ 1 mismatch. Or in low concentration and mixed samples, if the depth consensus is low (10-50 reads) but maintains high sequence identity (≥ 97% identity) and reference coverage (≥ 98%), and low mismatches (≤ 5). In these cases, the identity of any suspect sequences was subsequently verified through BLASTn analysis of NCBI GenBank (Sayers et al., 2022).
	•	INCONCLUSIVE: A sample was considered inconclusive if the consensus sequence had a sequencing depth of approximately 50 reads, a 90 – 96 % identity with the reference match, and a cover of the reference sequence ranging from 70 – 79%. 
	•	NEGATIVE: A sample was considered negative for a particular Phytophthora species when its consensus sequence does not match with this species (that is ≤ 90% identity with the reference match, cover ≤ 70% of the reference sequence, and ≥ 5 mismatches). 
	•	FAILLED: Since the Phy-MinION assay is a confirmatory assay to be used with samples that previously tested positive for the Phytophthora genus with screening -more conventional- methods (such as conventional or  real time PCR), the complete lack of results (lack of consensus sequence or consensus sequencing not matching with any Phytophthora species in the curated database) was not considered as an indication of pathogen absence but as a problem with the sample while conducting the assay, and the sample was assessed as fail sample. Those samples were set to be re-analyzed repeating the assay from the first step (PCR amplification).

d. Phytoplasma Sample Assessment

COMING SOON
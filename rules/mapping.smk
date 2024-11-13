
rule map_reads_database:
	input:
		fastq_clean = "demultiplexed/{barcode}_clean.fastq.gz"
	params:
		database = config[config['database']]
	output:
		mapping = "mapping/{barcode}.sam"
	conda:
		"envs/bwa.yaml"
	shell:
		"""
		if [ -s {input.fastq_clean} ]
		then
			bwa mem -x ont2d {params.database} {input.fastq_clean}  > {output.mapping}
		else
			touch {output}
		fi
		"""

rule count_mapped_reads:
	input:
		mapping = "mapping/{barcode}.sam"
	output:
		mapping_stats = "mapping/{barcode}.txt",
		read_count = "mapping/{barcode}_summary.txt"
	conda:
		"envs/bbmap.yaml"
	shell:
		"""
		if [ -s {input.mapping} ]
		then
			pileup.sh in={input.mapping} out={output.mapping_stats} secondary=false nzo=true headerpound=false 2>&1 | tee {output.read_count}
		else
			touch {output}
		fi
		"""

rule summarize_read_mapping:
	input:
		mapping_stats = expand("mapping/{barcode}.txt", barcode=barcodes),
		blast = "summary_output/blastn_summary.csv",
		blast_unfiltered = "summary_output/blastn_summary_unfiltered.csv",
	output:
		summary_file = "summary_output/{run_name}_mapping_summary.xlsx"
	conda:
		"envs/biopy.yaml"
	shell:
		"""
		python {config[snakefile_dir]}/scripts/summarize.py {output} {input.blast} {input.blast_unfiltered} {input.mapping_stats}
		rm -r mapping
		rm -r demultiplexed
		#rm -r ref_consensus_seq
		rm {input.blast}
		rm {input.blast_unfiltered}
		"""

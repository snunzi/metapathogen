
checkpoint align_ref_clusters:
	"""
	Align reads with minimap2
	"""
	input:
		fastq = "demultiplexed/{barcode}_clean.fastq.gz",
		db = config[config['database']]
	output:
		outdir = directory("ref_consensus_{barcode}"),
	params:
		sam = "ref_consensus_{barcode}/%.fq"
	conda:
		"envs/bbmap.yaml"
	threads: 20
	shell:
		"""
		if [ ! -z $(gzip -cd {input.fastq} | head -c1) ];
		then
			seal.sh in={input.fastq} ref={input.db} pattern={params.sam} ambig=all --ignorejunk qin=33
			mkdir -p {output}
		else
			mkdir {output}
		fi
		"""

def get_ref_cluster_names(wildcards):
	cp_output = checkpoints.align_ref_clusters.get(**wildcards).output[0]
	CLUS, = glob_wildcards(os.path.join(cp_output, "{cluster}.fq"))
	print(CLUS)
	return expand(os.path.join(cp_output, "{CLUSTER}.fq"), CLUSTER=CLUS)


rule consensus_from_ref_clusters:
	"""
	Align reads with minimap2
	"""
	input:
		clusters = get_ref_cluster_names
	output:
		consensus = "ref_consensus_seq/{barcode}_consensus_v1.fasta",
		read_count = "ref_consensus_seq/{barcode}_read_number.txt"
	params:
		"ref_consensus_{barcode}"
	conda:
		"envs/minimap2.yaml"
	threads: 20
	log: "log/{barcode}_map.txt"
	shell:
		"""
		if [ -z {input} ]; then
			touch {output}
		else
			for file in {input} ;  do
				if ! [ "${{file}}" == .fq] ; then
					total_reads=$(awk '{{s++}}END{{print s/4}}' ${{file}})
					if (( $(echo "$total_reads > 5" | bc -l) )); then
						filtlong --target_bases 750000 ${{file}} | gzip > ${{file}}_hq.fq.gz
						seqtk seq -a ${{file}}_hq.fq.gz  > ${{file}}.fa
						python {config[snakefile_dir]}/scripts/get_ref_read.py ${{file}}.fa ${{file}}_ref.fasta
						minimap2 -ax map-ont -k 15 ${{file}}_ref.fasta ${{file}}.fa -t {threads} > ${{file}}.sam
						racon -m 8 -x -6 -g -8 -w 500 -t {threads} ${{file}}.fa ${{file}}.sam ${{file}}_ref.fasta > ${{file}}_racon.fasta
						medaka_consensus -i ${{file}}.fa -d ${{file}}_racon.fasta -o ${{file}}_medaka -t 20
						grep --with-filename -c ">" ${{file}}.fa >> {output.read_count}
					else
						awk '{{s++}}END{{print s/4}}' ${{file}} >> {output.read_count}
						mkdir ref_consensus_barcode10/null_medaka
						touch ref_consensus_barcode10/null_medaka/consensus.fasta
					fi
				fi
			done
			cat {params}/*medaka/consensus.fasta 2>/dev/null > {output.consensus}
		fi || true
		rm -r {params}
		"""

rule scaffold_consensus_seqs:
	input:
		consensus = "ref_consensus_seq/{barcode}_consensus_v1.fasta"
	output:
		consensus = "ref_consensus_seq/{barcode}_consensus.fasta"
	conda:
		"envs/cap3.yaml"
	params:
		contigs = "ref_consensus_seq/{barcode}_consensus_v1.fasta.cap.contigs",
		singlets = "ref_consensus_seq/{barcode}_consensus_v1.fasta.cap.singlets",
		junk1= "ref_consensus_seq/{barcode}_consensus_v1.fasta.cap.contigs.links",
		junk2= "ref_consensus_seq/{barcode}_consensus_v1.fasta.cap.contigs.qual",
		junk3="ref_consensus_seq/{barcode}_consensus_v1.fasta.cap.ace",
		junk4="ref_consensus_seq/{barcode}_consensus_v1.fasta.cap.info",
	shell:
		"""
		if [ -s {input.consensus} ]
		then
			cap3 {input.consensus} -p 99
			cat {params.contigs} {params.singlets} > {output}
			#rm {params}
		else
			touch {output}
		fi
		"""

rule blast_consensus_seqs:
	input:
		consensus = "ref_consensus_seq/{barcode}_consensus.fasta",
	params:
		database = config[config['database']]
	output:
		"ref_consensus_seq/{barcode}_blast.out"
	conda:
		"envs/blast.yaml"
	shell:
		"""
		blastn -query {input} -db {params.database} -perc_identity 90 -outfmt '6 qseqid pident slen qlen length qcovs mismatch gapopen evalue bitscore salltitles' -max_target_seqs 5 -max_hsps 5 -num_threads {threads} > {output}
		"""

rule summarize_blast_output:
	input:
		blast = "ref_consensus_seq/{barcode}_blast.out",
		read_count = "ref_consensus_seq/{barcode}_read_number.txt",
		consensus = "ref_consensus_seq/{barcode}_consensus.fasta",
	params:
		barcode="{barcode}",
		run = directory_name
	output:
		all = "ref_consensus_seq/{barcode}_summary.csv",
		filtered = "ref_consensus_seq/{barcode}_summary_filtered.csv",
	conda:
		"envs/biopy.yaml"
	shell:
		"""
		if [ -s {input.blast} ]
		then
			python {config[snakefile_dir]}/scripts/summarize_blast_refv2.py {input.blast} {input.read_count} {input.consensus} {params.barcode} {output.all} {output.filtered} {params.run}
		else
			touch {output}
		fi
		"""


rule blastn_excel_summary:
		input:
				expand("ref_consensus_seq/{barcode}_summary_filtered.csv", barcode=barcodes)
		output:
				"summary_output/blastn_summary.csv"
		shell:
				"""
				head -1 {input[0]} > {output}; tail -n +2 -q {input} >> {output}
				"""

rule blastn_excel_summary_all:
		input:
				expand("ref_consensus_seq/{barcode}_summary.csv", barcode=barcodes)
		output:
				"summary_output/blastn_summary_unfiltered.csv"
		shell:
				"""
				head -1 {input[0]} > {output}; tail -n +2 -q {input} >> {output}
				#rm -r ref_consensus_seq
				"""

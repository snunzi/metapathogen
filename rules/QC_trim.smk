rule trim_barcodes:
	"""
	Trim adapters with porechop
	"""
	input:
		reads = expand("data/{{barcode}}{exten}", exten=config['R1_extension']),
	output:
		fastq = "demultiplexed/{barcode}_raw.fastq.gz",
		read_count = "mapping/{barcode}_raw_reads.txt"
	message:
		"""Concatenating reads from {wildcards.barcode}."""
	conda:
		"envs/porechop.yaml"
	threads: 20
	shell:
		"""
		porechop --discard_middle -i {input.reads} -o {output.fastq} --threads {threads}
		printf "%s\t%s\n" {wildcards.barcode} \ $(echo $(zcat {input.reads}|wc -l)/4|bc) > {output.read_count}
		"""

rule nanoplot_minion_raw:
	"""
	Run QC stats on raw reads using nanoplot
	"""
	input:
		fastq = "demultiplexed/{barcode}_raw.fastq.gz"
	output:
		out1 = "read_QC/{barcode}_raw_NanoStats.txt",
		out2 = "read_QC/{barcode}_raw_NanoPlot-report.html"
	message:
		"""QC on raw {wildcards.barcode} reads with NanoPlot."""
	params:
		out_dir = "{barcode}_raw_QC",
		out_prefix = "{barcode}_raw_",
		out1_temp = "{barcode}_raw_QC/{barcode}_raw_NanoStats.txt",
		out2_temp = "{barcode}_raw_QC/{barcode}_raw_NanoPlot-report.html"
	threads: 15
	conda:
		"envs/nanoplot.yaml"
	shell:
		"""
		NanoPlot --fastq {input} -o {params.out_dir} -p {params.out_prefix} -t {threads}
		mv {params.out1_temp} {output.out1}
		mv {params.out2_temp} {output.out2}
		rm -r {params.out_dir}
		"""

rule trim_minion:
	"""
	Discard short and low quality with chopper
	"""
	input:
		fastq = "demultiplexed/{barcode}_raw.fastq.gz"
	output:
		fastq_clean = "demultiplexed/{barcode}_clean.fastq.gz"
	message:
		"""Quality filter {wildcards.barcode} with chopper."""
	threads: 15
	log:
		"log/chopper/{barcode}.log"
	conda:
		"envs/chopper.yaml"
	shell:
		"gunzip -c {input} | chopper --maxlength {config[maxlength]} --minlength {config[minlength]} --quality {config[minqual]} | gzip > {output} 2>&1 | tee {log}"

rule nanoplot_minion_trimmed:
	"""
	Run QC stats on cleaned reads using nanoplot
	"""
	input:
		fastq_clean = "demultiplexed/{barcode}_clean.fastq.gz"
	output:
		out1 = "read_QC/{barcode}_clean_NanoStats.txt",
		out2 = "read_QC/{barcode}_clean_NanoPlot-report.html"
	message:
		"""QC on clean {wildcards.barcode} reads with NanoPlot."""
	params:
		out_dir = "{barcode}_clean_QC",
		out_prefix = "{barcode}_clean_",
		out1_temp = "{barcode}_clean_QC/{barcode}_clean_NanoStats.txt",
		out2_temp = "{barcode}_clean_QC/{barcode}_clean_NanoPlot-report.html"
	threads: 15
	conda:
		"envs/nanoplot.yaml"
	shell:
		"""
		if [ ! -z $(gzip -cd {input} | head -c1) ];
		then
			NanoPlot --fastq {input} -o {params.out_dir} -p {params.out_prefix} -t {threads}
			mv {params.out1_temp} {output.out1}
			mv {params.out2_temp} {output.out2}
			rm -r {params.out_dir}
		else
			touch {output.out1}
			touch {output.out2}
		fi
		"""

rule multiqc_minion:
	"""
	Summarize QC reports with MultiQC.
	"""
	input:
		raw = expand("read_QC/{barcode}_raw_NanoStats.txt",barcode=barcodes),
		clean = expand("read_QC/{barcode}_clean_NanoStats.txt",barcode=barcodes),
	output:
		"summary_output/{run_name}_readqc_report.html",
	message:
		"""Summarize QC into final QC report."""
	params:
		summary = "read_QC/multiqc_report.html",
		out_dir = directory("read_QC")
	conda:
		"envs/multiqc.yaml"
	shell:
		"""
		multiqc -o {params.out_dir} {input}
		mv {params.summary} {output}
		rm -r read_QC
		"""


rule create_blast_database:
	input:
		db_fasta = config[config['database']]
	output:
		db_done = "log/blast_db_done.txt"
	conda:
		"envs/blast.yaml"
	shell:
		"""
		makeblastdb -in {input.db_fasta} -dbtype nucl
		touch {output}
		"""

rule create_bwa_database:
	input:
		db_fasta = config[config['database']]
	output:
		db_done = "log/bwa_db_done.txt"
	conda:
		"envs/bwa.yaml"
	shell:
		"""
		bwa index {input} {input}
		touch {output}
		"""

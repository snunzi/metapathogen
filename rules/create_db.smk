
rule create_blast_database:
	input:
		db_fasta = config[config['database']]
	output:
		db_done = "log/blast_db_done.txt"
	conda:
		"envs/blast.yaml"
	shell:
		"""
		sed -i '/^>/s/\W/_/2g' {input.db_fasta}
		makeblastdb -in {input.db_fasta} -dbtype nucl
		touch {output}
		"""

rule create_bwa_database:
	input:
		db_fasta = config[config['database']],
		db_done = "log/blast_db_done.txt"
	output:
		db_done = "log/bwa_db_done.txt"
	conda:
		"envs/bwa.yaml"
	shell:
		"""
		bwa index {input} {input}
		touch {output}
		"""

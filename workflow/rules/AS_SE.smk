import os
import yaml 

rule all:
	input:
		qc_raw_report = "results/intermediate_files/multiqc/AS/raw/multiqc_report.html",
		qc_trim_report = "results/intermediate_files/multiqc/AS/trim/multiqc_report.html",
		silvadb = "results/intermediate_files/silva_db/silva_nr99_v138.1_train_set.fa.gz", 
		top_taxa_and_host_names = "results/intermediate_files/top_taxa_host.txt",
		as_abundance_plot = "results/final/AS/as_abundance_plot.svg",
		as_based_database = "results/final/prot_db/as_based_database.fa",
		diff_abun_as = "results/final/diff_abun/taxa-maaslin2-AS/maaslin2.log"

rule fastqc_raw_se:
	input:
		r1=expand("data/AS/raw/{sample}_1.fastq.gz", sample=config["AS_samples"])
	output:
		directory("results/intermediate_files/fastqc/raw/AS")
	conda:
		srcdir("../envs/fastqc.yaml")
	threads: 16
	shell:
		"""
		mkdir -p {output}
		fastqc -t {threads} {input.r1} -o {output}
		"""

rule multiqc_raw:
	input:
		fastqc_html = "results/intermediate_files/fastqc/raw/AS"
	output:
		output = "results/intermediate_files/multiqc/AS/raw/multiqc_report.html"
	params:
		output_dir =  "results/intermediate_files/multiqc/AS/raw/"
	conda:
		srcdir("../envs/multiqc.yaml")
	shell:
		"multiqc {input} -o {params.output_dir}"

rule trim_se:
	input:
		r1="data/AS/raw/{sample}_1.fastq.gz"
	output:
		o1="results/intermediate_files/trimmed/AS/{sample}_1.fastq.gz"
	params:
		params = config["parameters"]["trimmomatic"]
	conda:
		srcdir("../envs/trimmomatic.yaml")
	shell:
		"""
		trimmomatic SE -threads {threads} {input.r1} {output.o1} {params.params}
		"""

rule fastqc_trim_se:
	input:
		r1=expand("results/intermediate_files/trimmed/AS/{sample}_1.fastq.gz", sample=config["AS_samples"])
	output:
		directory("results/intermediate_files/fastqc/trim/AS")
	params:
		out="results/intermediate_files/fastqc/trim/AS"
	conda:
		srcdir("../envs/fastqc.yaml")
	threads: 16
	shell:
		"""
		mkdir -p {output}
		fastqc -t {threads} {input.r1} -o {params.out}
		"""

rule multiqc_trim:
	input:
		fastqc_html = "results/intermediate_files/fastqc/trim/AS"
	output:
		output = "results/intermediate_files/multiqc/AS/trim/multiqc_report.html"
	params:
		output_dir =  "results/intermediate_files/multiqc/AS/trim/"
	conda:
		srcdir("../envs/multiqc.yaml")
	shell:
		"multiqc {input} -o {params.output_dir}"

rule gunzip_se:
	input:
		input="results/intermediate_files/trimmed/AS/{sample}_1.fastq.gz"
	output:
		f1="results/intermediate_files/merged/AS/{sample}.extendedFrags.fastq"
	shell:
		"gunzip -c {input.input} > {output.f1}"

rule download_silvadb:
	output:
		"results/intermediate_files/silva_db/silva_nr99_v138.1_train_set.fa.gz"
	shell:
		"""
		mkdir -p results/intermediate_files/silva_db
		cd results/intermediate_files/silva_db
		echo 6b41db7139834c71171f8ce5b5918fc6	silva_nr99_v138.1_train_set.fa.gz > silva_nr99_v138.1_train_set.fa.gz.md5
		wget "https://zenodo.org/record/4587955/files/silva_nr99_v138.1_train_set.fa.gz"
		# Recommended checksum to confirm the download was successful:
		md5sum -c silva_nr99_v138.1_train_set.fa.gz.md5
		# Must return silva_nr99_v138.1_train_set.fa.gz: OK*
		# If not - try downloading the file again as it may be a corrupted copy.
		"""

rule dada2:
	input:
		silva_db = "results/intermediate_files/silva_db/silva_nr99_v138.1_train_set.fa.gz",
		samples = expand("results/intermediate_files/merged/AS/{sample}.extendedFrags.fastq", sample=config["AS_samples"])
	output:
		d1 = "results/intermediate_files/top_taxa.txt",
		d2 = "results/final/AS/as_abundance_plot.svg",
		d3 = "results/final/AS/as_physeq.rds",
		d4 = "results/final/diff_abun/taxa-maaslin2-AS/maaslin2.log"
	params:
		param1 = config["parameters"]["group"],
		param2 = config["parameters"]["taxa_rank"],
		param3 = config["parameters"]["top_taxa"]	
	conda:
		srcdir("../envs/R.yaml")
	shell:
		"Rscript workflow/scripts/seq_process_as.R -g {params.param1} -t {params.param2} -n {params.param3}"

rule add_host:
	input:
		"results/intermediate_files/top_taxa.txt"
	output:
		"results/intermediate_files/top_taxa_host.txt"
	params:
		params = config["host"]
	shell:
		"""
		echo {params.params} > results/intermediate_files/host.txt 
		cat {input} results/intermediate_files/host.txt > {output}
		"""

checkpoint downloadgenomes:
	input:
		"results/intermediate_files/top_taxa_host.txt"
	output:
		directory("results/intermediate_files/genomes/")
	conda:
		srcdir("../envs/ncbi.yaml")
	shell:
		"ncbi-genome-download -g {input} -l complete -o results/intermediate_files/genomes -F protein-fasta -r 3 --flat-output all"

def aggregate_input(wildcards):
	checkpoint_output = checkpoints.downloadgenomes.get(**wildcards).output[0]
	return expand("results/intermediate_files/genomes/{i}.faa.gz",
		   i=glob_wildcards(os.path.join(checkpoint_output, "{i}.faa.gz")).i)

rule prep_proteins:
	input:
		aggregate_input
	output:
		temp("results/intermediate_files/prot_db/db.faa.gz")
	shell:
		"cat {input} > {output}"

rule gunzip_db:
	input:
		"results/intermediate_files/prot_db/db.faa.gz"
	output:
		temp("results/intermediate_files/prot_db/db.fa")
	shell:
		"gunzip -c {input} > {output}"

rule prep_headers:
	input:
		"results/intermediate_files/prot_db/db.fa"
	output:
		temp("results/intermediate_files/prot_db/db_v2.fa")
	shell:
		"sed -e 's/\[//g;s/\]//g;s/\.//g;s/-/_/g;s/://g;s/ /_/g' {input} > {output}"

rule rm_dup:
	input:
		"results/intermediate_files/prot_db/db_v2.fa"
	output:
		"results/final/prot_db/as_based_database.fa"
	conda:
		srcdir("../envs/seqkit.yaml")
	shell:
		"seqkit rmdup -i -s < {input} > {output}"
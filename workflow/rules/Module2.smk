# 2:  Amplicon sequencing + Metaproteomics

import os
import yaml 

# Config file to store initial input data
configfile: "config.yaml"

ruleorder: trimPE > trimSE

rule all:
	input:
		silvadb = "results/intermediate_files/silva_db/silva_nr99_v138.1_train_set.fa.gz", 
		top_taxa_and_host_names = "results/intermediate_files/top_taxa_host.txt",
		as_abundance_plot = "results/final/AS/as_abundance_plot.svg",
		as_based_database = "results/final/prot_db/as_based_database.fa",
		diff_abun_as = "results/final/diff_abun/taxa-maaslin2-AS/maaslin2.log",
		pep_abundance_table = "results/final/MP/pep_abundance_table.txt",
		peptide_list = "results/final/MP/unipept_list.txt",
		unipept_results = "results/final/MP/unipept_results.csv",
		diff_abun_mp = "results/final/diff_abun/taxa-maaslin2-MP/maaslin2.log",
		combi_results = "results/final/integrated/combi_plot.png"

rule trimPE:
	input:
		r1="data/AS/raw/{sample}_1.fastq.gz",
		r2="data/AS/raw/{sample}_2.fastq.gz"
	output:
		o1="results/intermediate_files/trimmed/AS/{sample}_1.fastq.gz",
		o2="results/intermediate_files/trimmed/AS/{sample}_2.fastq.gz",
		o1un="results/intermediate_files/trimmed/AS/{sample}_1un.trim.fastq.gz",
		o2un="results/intermediate_files/trimmed/AS/{sample}_2un.trim.fastq.gz"
	params:
		params = config["parameters"]["trimmomatic"]
	conda:
		"../envs/trimmomatic.yaml"
	shell:
		"trimmomatic PE {input.r1} {input.r2} {output.o1} {output.o1un} {output.o2} {output.o2un} {params.params}"

rule trimSE:
	input:
		r1="data/AS/raw/{sample}.fastq.gz"
	output:
		o1="results/intermediate_files/trimmed/AS/{sample}.fastq.gz"
	params:
		params = config["parameters"]["trimmomatic"]
	conda:
		srcdir("../envs/trimmomatic.yaml")
	shell:
		"trimmomatic SE {input.r1} {output.o1} {params.params}"

rule gunzipSE:
	input:
		trimmed_se = "results/intermediate_files/trimmed/AS/{sample}.fastq.gz"
	output:
		output = "results/intermediate_files/merged/AS/{sample}.extendedFrags.fastq"
	shell:
		"gunzip -c {input} > {output}"

rule mergePE:
	input:
		o1="results/intermediate_files/trimmed/AS/{sample}_1.fastq.gz",
		o2="results/intermediate_files/trimmed/AS/{sample}_2.fastq.gz"
	output:
		f1="results/intermediate_files/merged/AS/{sample}.extendedFrags.fastq",
		f2="results/intermediate_files/merged/AS/{sample}.notCombined_1.fastq",
		f3="results/intermediate_files/merged/AS/{sample}.notCombined_2.fastq"
	params:
		dir = "results/intermediate_files/merged/AS/",
		gid = "{sample}"
	conda:
		srcdir("../envs/flash2.yaml")
	shell:
		"flash2 {input.o1} {input.o2} -d {params.dir} -o {params.gid} --allow-outies"

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
		echo {params.params} > resources/host.txt 
		cat {input} resources/host.txt > {output}
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

rule cat_AS:
	input:
		aggregate_input
	output:
		temp("results/intermediate_files/prot_db/db.faa.gz")
	shell:
		"cat {input} > {output}"

rule gunzip_AS:
	input:
		"results/intermediate_files/prot_db/db.faa.gz"
	output:
		temp("results/intermediate_files/prot_db/db.fa")
	shell:
		"gunzip -c {input} > {output}"

rule prep_headers_AS:
	input:
		"results/intermediate_files/prot_db/db.fa"
	output:
		temp("results/intermediate_files/prot_db/db_v2.fa")
	shell:
		"sed -e 's/\[//g;s/\]//g;s/\.//g;s/-/_/g;s/://g;s/ /_/g' {input} > {output}"

rule seqkit_AS:
	input:
		"results/intermediate_files/prot_db/db_v2.fa"
	output:
		"results/final/prot_db/as_based_database.fa"
	conda:
		srcdir("../envs/seqkit.yaml")
	shell:
		"seqkit rmdup -i -s < {input} > {output}"

rule build_msgf_db:
	input:
		"results/final/prot_db/as_based_database.fa"
	output:
		"results/final/prot_db/as_based_database.canno"
	conda:
		srcdir("../envs/msgfplus.yaml")
	shell:
		"msgf_plus -Xmx16000M edu.ucsd.msjava.msdbsearch.BuildSA -d {input}"

rule run_msgf:
	input:
		input = "data/MP/spectra/{spectra}.mgf",
		prot_db = "results/final/prot_db/as_based_database.fa",
		auxiliary = "results/final/prot_db/as_based_database.canno"
	output:
		output="results/intermediate_files/msgf_plus/{spectra}.mzid"
	conda:
		srcdir("../envs/msgfplus.yaml")
	params:
		params = config["parameters"]["msgf"]
	shell:
		"""
		msgf_plus -Xmx200g -s {input.input} -d {input.prot_db} {params.params} -o {output.output}
		"""

rule mzid_to_tsv:
	input:
		input = "results/intermediate_files/msgf_plus/{spectra}.mzid"
	output:
		output = "results/intermediate_files/msgf_plus/{spectra}.tsv"
	conda:
		srcdir("../envs/msgfplus.yaml") 
	shell:
		"msgf_plus edu.ucsd.msjava.ui.MzIDToTsv -Xmx20g -i {input.input} -o {output.output}"

rule process_pept:
	input:
		input = "results/intermediate_files/msgf_plus/{spectra}.tsv"
	output:
		output = "results/intermediate_files/msgf_plus/peptides_{spectra}.txt"
	shell:
		"""
		awk 'BEGIN{{FS="\t"; OFS=FS;}} $16<0.05 {{ print; }}' {input} | cut -f10 | sed 's/^.\.//' | sed 's/[0-9\+,]//g' | sed 's/\..$//' | sed 's/\.//g'  > {output}
		"""

rule unipept_prep:
	input:
		input = expand("results/intermediate_files/msgf_plus/peptides_{spectra}.txt", spectra=config["MP_samples"])
	output:
		output1 = "results/final/MP/unipept_list.txt",
		output2 = "results/final/MP/pep_abundance_table.txt"
	conda:
		srcdir("../envs/R.yaml")
	shell:
		"Rscript workflow/scripts/pep_process_mp.R"

rule run_unipept:
	input:
		input="results/final/MP/unipept_list.txt"
	output:
		output="results/final/MP/unipept_results.csv"
	conda:
		srcdir("../envs/pyteomics.yaml")
	shell:
		"python workflow/scripts/unipept_get_peptinfo.py -i {input.input} -o {output}"

rule diff_abund_MP:
	input:
		abundance = "results/final/MP/pep_abundance_table.txt",
		metadata = "resources/metadata.txt",
		taxonomy = "results/final/MP/unipept_results.csv"
	output:
		"results/final/diff_abun/taxa-maaslin2-MP/maaslin2.log"
	params:
		outdir = "results/diff_abun/taxa-maaslin2-MP/",
		param1 = config["parameters"]["group"],
		param2 = config["parameters"]["taxa_rank"]
	conda:
		srcdir("../envs/R.yaml")
	shell:
		"Rscript workflow/scripts/taxa_diff_abun_mp.R --group {params.param1} --taxa_rank {params.param2}"

rule combi:
	input:
		abundance_mp = "results/final/MP/pep_abundance_table.txt",
		metadata = "resources/metadata.txt",
		abundance_as = "results/final/AS/as_physeq.rds"
	output:
		plot = "results/final/integrated/combi_plot.png"
	params:
		param1 = config["parameters"]["group"]
	conda:
		srcdir("../envs/R.yaml")
	shell:
		"Rscript workflow/scripts/visual_integration.R --group {params.param1}"
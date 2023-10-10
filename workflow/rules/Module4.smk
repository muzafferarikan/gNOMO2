# 4:  Metagenomics + Metatranscriptomics

import os
import yaml 

ruleorder: trimPE > trimSE

rule all:
	input:
		interproscan_db = "results/final/prot_db/database_MP.fasta",
		interproscan_setup = "results/intermediate_files/interproscan/interproscan_setup.txt",
		diff_abun_results_mg = "results/final/diff_abun/taxa-maaslin2-MG/maaslin2.log",
		diff_abun_results_mt = "results/final/diff_abun/taxa-maaslin2-MT/maaslin2.log",
		diff_abun_tigrfam_results_mg = "results/final/diff_abun/tigrfam-maaslin2-MG/maaslin2.log",
		diff_abun_tigrfam_results_mt = "results/final/diff_abun/tigrfam-maaslin2-MT/maaslin2.log",
		pathview_results = "results/final/integrated/pathview/log.txt",
		combi_results = "results/final/integrated/combi_plot.png"

rule trimPE:
	input:
		r1="data/{omics}/raw/{sample}_1.fastq.gz",
		r2="data/{omics}/raw/{sample}_2.fastq.gz"
	output:
		o1="results/intermediate_files/trimmed/{omics}/{omics}_{sample}_1.fastq.gz",
		o2="results/intermediate_files/trimmed/{omics}/{omics}_{sample}_2.fastq.gz",
		o1un="results/intermediate_files/trimmed/{omics}/{omics}_{sample}_1un.trim.fastq.gz",
		o2un="results/intermediate_files/trimmed/{omics}/{omics}_{sample}_2un.trim.fastq.gz"
	params:
		params = config["parameters"]["trimmomatic"]
	conda:
		srcdir("../envs/trimmomatic.yaml")
	shell:
		"trimmomatic PE {input.r1} {input.r2} {output.o1} {output.o1un} {output.o2} {output.o2un} {params.params}"

rule trimSE:
	input:
		r1="data/{omics}/raw/{sample}.fastq.gz"
	output:
		o1="results/intermediate_files/trimmed/{omics}/{sample}.fastq.gz"
	params:
		params = config["parameters"]["trimmomatic"]
	conda:
		srcdir("../envs/trimmomatic.yaml")
	shell:
		"trimmomatic SE {input.r1} {output.o1} {params.params}"

rule gunzipSE:
	input:
		input="results/intermediate_files/trimmed/{omics}/{omics}_{sample}.fastq.gz"
	output:
		output="results/intermediate_files/merged/{omics}/{omics}_{sample}.extendedFrags.fastq"
	shell:
		"gunzip -c {input} > {output}"

rule mergePE:
	input:
		o1="results/intermediate_files/trimmed/{omics}/{omics}_{sample}_1.fastq.gz",
		o2="results/intermediate_files/trimmed/{omics}/{omics}_{sample}_2.fastq.gz"
	output:
		f1="results/intermediate_files/merged/{omics}/{omics}_{sample}.extendedFrags.fastq",
		f2="results/intermediate_files/merged/{omics}/{omics}_{sample}.notCombined_1.fastq",
		f3="results/intermediate_files/merged/{omics}/{omics}_{sample}.notCombined_2.fastq"
	params:
		dir = "results/intermediate_files/merged/{omics}/",
		gid = "{omics}_{sample}"
	conda:
		srcdir("../envs/flash2.yaml")
	shell:
		"flash2 {input.o1} {input.o2} -d {params.dir} -o {params.gid} --allow-outies"
	
rule cat:
	input:
		f1="results/intermediate_files/merged/{omics}/{omics}_{sample}.extendedFrags.fastq",
		f2="results/intermediate_files/merged/{omics}/{omics}_{sample}.notCombined_1.fastq",
		f3="results/intermediate_files/merged/{omics}/{omics}_{sample}.notCombined_2.fastq"
	output:
		output="results/intermediate_files/merged/{omics}/{omics}_{sample}/all.fastq"
	shell:
		"cat {input} > {output}"

rule sed:
	input:
		MG="results/intermediate_files/merged/{omics}/{omics}_{sample}/all.fastq"
	output:
		output="results/intermediate_files/kaiju/kaiju_input/{omics}_{sample}_ns.fastq"
	shell:
		"sed -r 's/\s+//g' {input} > {output}"

rule kaiju_db:
	output:
		nr = "results/intermediate_files/kaiju/kaiju_db/nr/kaiju_db_nr.fmi",
		nodes = "results/intermediate_files/kaiju/kaiju_db/nodes.dmp",
		names = "results/intermediate_files/kaiju/kaiju_db/names.dmp"
	conda:
		srcdir("../envs/kaiju.yaml")
	shell:
		"""
		cd results/intermediate_files/kaiju/kaiju_db/
		kaiju-makedb -s nr -t 2
		"""

rule kaiju_classify:
	input:
		input="results/intermediate_files/kaiju/kaiju_input/{omics}_{sample}_ns.fastq",
		nodes="results/intermediate_files/kaiju/kaiju_db/nodes.dmp",
		fmi="results/intermediate_files/kaiju/kaiju_db/nr/kaiju_db_nr.fmi"
	output:
		output="results/intermediate_files/kaiju/kaiju_output/{omics}/{omics}_{sample}_ns_kaiju.out"
	conda:
		srcdir("../envs/kaiju.yaml")
	shell:
		"""
		kaiju -t {input.nodes} -f {input.fmi} -i {input.input} -o {output}
		"""

rule kaiju_summarize:
	input:
		in1="results/intermediate_files/kaiju/kaiju_output/{omics}/{omics}_{sample}_ns_kaiju.out",
		in2="results/intermediate_files/kaiju/kaiju_db/nodes.dmp",
		in3="results/intermediate_files/kaiju/kaiju_db/names.dmp"
	output:
		output="results/intermediate_files/kaiju/kaiju_output/{omics}/{omics}_{sample}_kaiju_summary.tsv"
	conda:
		srcdir("../envs/kaiju.yaml")
	shell:
		"""
		kaiju2table -t {input.in2} -n {input.in3} -r genus -o {output} {input.in1} -u
		"""		

rule diff_abund_MG:
	input:
		input=expand("results/intermediate_files/kaiju/kaiju_output/MG/MG_{sample}_kaiju_summary.tsv", omics=config["omics"], sample=config["MG_samples"])
	output:
		log = "results/final/diff_abun/taxa-maaslin2-MG/maaslin2.log",
		abundance_mg = "results/final/MG/mg_taxa_abundance.txt"
	params:
		param1 = config["parameters"]["group"]
	conda:
		srcdir("../envs/R.yaml")
	shell:
		"Rscript workflow/scripts/taxa_diff_abun_mg.R --group {params.param1}"

rule diff_abund_MT:
	input:
		input=expand("results/intermediate_files/kaiju/kaiju_output/MT/MT_{sample}_kaiju_summary.tsv", omics=config["omics"], sample=config["MG_samples"])
	output:
		log = "results/final/diff_abun/taxa-maaslin2-MT/maaslin2.log",
		abundance_mt = "results/final/MT/mt_taxa_abundance.txt"
	params:
		param1 = config["parameters"]["group"]
	conda:
		srcdir("../envs/R.yaml")
	shell:
		"Rscript workflow/scripts/taxa_diff_abun_mt.R --group {params.param1}"

rule metaspades:
	input:
		i1="results/intermediate_files/trimmed/MG/MG_{sample}_1.fastq.gz",
		i2="results/intermediate_files/trimmed/MG/MG_{sample}_2.fastq.gz"
	params:
		outdir = "results/intermediate_files/spades/MG_{sample}/"
	output:
		output = "results/intermediate_files/spades/MG_{sample}/contigs.fasta"
	conda:
		srcdir("../envs/spades.yaml")
	threads: 16
	shell:
		"metaspades.py -1 {input.i1} -2 {input.i2} -o {params.outdir}"

rule rnaspades:
	input:
		i1="results/intermediate_files/trimmed/MT/MT_{sample}_1.fastq.gz",
		i2="results/intermediate_files/trimmed/MT/MT_{sample}_2.fastq.gz"
	params:
		outdir = "results/intermediate_files/spades/MT_{sample}/"
	output:
		output = "results/intermediate_files/spades/MT_{sample}/transcripts.fasta"
	conda:
		srcdir("../envs/spades.yaml")
	threads: 16
	shell:
		"rnaspades.py -1 {input.i1} -2 {input.i2} -o {params.outdir}"

rule augustus:
	input:
		input1 = "results/intermediate_files/spades/MG_{sample}/contigs.fasta",
		input2 = "results/intermediate_files/spades/MT_{sample}/transcripts.fasta"
	output:
		output1 = "results/intermediate_files/augustus/MG_{sample}/augustus_output.gff",
		output2 = "results/intermediate_files/augustus/MT_{sample}/augustus_output.gff"
	conda:
		srcdir("../envs/augustus.yaml")
	params:
		params = config["parameters"]["augustus"]
	shell:
		"""
		augustus {input.input1} {params.params} > {output.output1}
		augustus {input.input2} {params.params} > {output.output2}
		"""

rule getAnnoFasta:
	input:
		sample = "results/intermediate_files/augustus/{omics}_{sample}/augustus_output.gff"
	output:
		output = "results/intermediate_files/augustus/{omics}_{sample}/augustus_output.aa"
	conda:
		srcdir("../envs/perl.yaml")
	shell:
		"perl ./workflow/scripts/get_anno_fasta.pl {input}"

rule seqkit_host:
	input:
		sample = "results/intermediate_files/augustus/{omics}_{sample}/augustus_output.aa"
	output:
		unique = "results/intermediate_files/augustus/{omics}_{sample}/augustus_output_uniq.fasta"
	conda:
		srcdir("../envs/seqkit.yaml")
	shell:
		"seqkit rename -1 -O uniq < {input} > {output}"

rule mark_host_protein_seqs:
	input:
		input = "results/intermediate_files/augustus/{omics}_{sample}/augustus_output_uniq.fasta"
	output:
		output = "results/intermediate_files/aug_prod/{omics}_{sample}/proteins_host.fasta"
	shell:
		"sed 's/^>/>host_/' {input} > {output}"

rule prodigal:
	input:
		input1 = "results/intermediate_files/spades/MG_{sample}/contigs.fasta",
		input2 = "results/intermediate_files/spades/MT_{sample}/transcripts.fasta"
	output:
		proteins1 = "results/intermediate_files/prodigal/MG_{sample}/Proteinas.fasta",
		proteins2 = "results/intermediate_files/prodigal/MT_{sample}/Proteinas.fasta"
	conda:
		srcdir("../envs/prodigal.yaml")
	shell:
		"""
		prodigal -i {input.input1} -p meta -a {output.proteins1}
		prodigal -i {input.input2} -p meta -a {output.proteins2}
		"""

rule seqkit_microbial:
	input:
		sample = "results/intermediate_files/prodigal/{omics}_{sample}/Proteinas.fasta"
	output:
		unique = "results/intermediate_files/prodigal/{omics}_{sample}/Proteinas_uniq.fasta"
	conda:
		srcdir("../envs/seqkit.yaml")
	shell:
		"seqkit rename -1 -O uniq < {input} > {output}"

rule mark_microbial_protein_seqs:
	input:
		input = "results/intermediate_files/prodigal/{omics}_{sample}/Proteinas_uniq.fasta"
	output:
		output = "results/intermediate_files/aug_prod/{omics}_{sample}/proteins_microbial.fasta"
	shell:
		"sed 's/^>/>microbial_/' {input} > {output}"

rule aug_prod:
	input:
		host = "results/intermediate_files/aug_prod/{omics}_{sample}/proteins_host.fasta",
		microbial = "results/intermediate_files/aug_prod/{omics}_{sample}/proteins_microbial.fasta"
	output:
		output = "results/intermediate_files/aug_prod/{omics}_{sample}/all_proteins.fasta"
	shell:
		"cat {input.host} {input.microbial} > {output}"

rule sed_prot:
	input:
		sample = "results/intermediate_files/aug_prod/{omics}_{sample}/all_proteins.fasta"
	output:
		output = "results/intermediate_files/aug_prod/{omics}_{sample}/all_proteins_clean.fasta"
	shell:
		"sed 's/*//g' {input} > {output}"

rule merge_all_samples:
	input:
		sample = expand("results/intermediate_files/aug_prod/{omics}_{sample}/all_proteins_clean.fasta", omics=config["omics"], sample=config["MG_samples"]),
	output:
		output = "results/final/prot_db/database_MP.fasta"
	shell:
		"cat {input} > {output}"

rule interproscan_setup:
	input:
		database = "results/final/prot_db/database_MP.fasta"
	output:
		output = "results/intermediate_files/interproscan/interproscan_setup.txt"
	shell:
		"""
		mkdir -p results/intermediate_files/interproscan 
		cd results/intermediate_files/interproscan
		wget ftp://ftp.ebi.ac.uk/pub/software/unix/iprscan/5/5.62-94.0/interproscan-5.62-94.0-64-bit.tar.gz >> interproscan_setup.txt
		wget ftp://ftp.ebi.ac.uk/pub/software/unix/iprscan/5/5.62-94.0/interproscan-5.62-94.0-64-bit.tar.gz.md5 >> interproscan_setup.txt

		# Recommended checksum to confirm the download was successful:
		md5sum -c interproscan-**.tar.gz.md5 >> interproscan_setup.txt
		# Must return *interproscan-**.tar.gz: OK*
		# If not - try downloading the file again as it may be a corrupted copy.
		tar -pxvzf interproscan-**.tar.gz >> interproscan_setup.txt
		cd interproscan-5.62-94.0
		python3 setup.py -f interproscan.properties
	"""

rule interproscan:
	input:
		sample = "results/intermediate_files/aug_prod/{omics}_{sample}/all_proteins_clean.fasta",
		setup = "results/intermediate_files/interproscan/interproscan_setup.txt"
	output:
		output = "results/intermediate_files/interproscan/{omics}_{sample}/TIGRFAM.tsv"
	params:
		prefix = "results/intermediate_files/interproscan/{omics}_{sample}/TIGRFAM"
	conda:
		srcdir("../envs/interproscan.yaml")
	shell:
		"./results/intermediate_files/interproscan/interproscan*/interproscan.sh -appl TIGRFAM -pa -dra -b {params} -i {input}"

rule cut_tigrfam:
	input:
		sample = "results/intermediate_files/interproscan/{omics}_{sample}/TIGRFAM.tsv"
	output:
		interproscan_files = "results/intermediate_files/tigrfam/{omics}_{sample}_TIGRFAM_v2.tsv"
	shell:
		"cut -f 1,5 {input} > {output}"

rule seqkit_rmdup:
	input:
		"results/final/prot_db/database_MP.fasta"
	output:
		"results/final/prot_db/db_ready.fasta"
	conda:
		srcdir("../envs/seqkit.yaml")
	shell:
		"seqkit rmdup -i -s < {input} > {output}"

rule tigrfam_diff:
	input:
		input=expand("results/intermediate_files/tigrfam/{omics}_{sample}_TIGRFAM_v2.tsv", omics=config["omics"], sample=config["MG_samples"])
	output:
		diff_abun_tigrfam_results_mg = "results/final/diff_abun/tigrfam-maaslin2-MG/maaslin2.log",
		diff_abun_tigrfam_results_mt = "results/final/diff_abun/tigrfam-maaslin2-MT/maaslin2.log"
	params:
		outdir = "results/final/diff_abun/",
		param1 = config["parameters"]["group"]
	conda:
		srcdir("../envs/R.yaml")
	shell:
		"Rscript workflow/scripts/tigrfam_diff_abun.R --group {params.param1}"

rule eggnog_db:
	input: 
		sample = expand("results/intermediate_files/aug_prod/{omics}_{sample}/all_proteins_clean.fasta",  omics=config["omics"], sample=config["MG_samples"])
	output:
		eggdb = "results/intermediate_files/eggnog/eggnog.db",
		diamond = "results/intermediate_files/eggnog/eggnog_proteins.dmnd"
	conda:
		srcdir("../envs/eggnog.yaml")
	shell:
		"download_eggnog_data.py --data_dir results/intermediate_files/eggnog/  -y bact euk > ./results/intermediate_files/eggnog/eggnog_download.log 2>&1"

rule eggnog:
	input:
		sample = "results/intermediate_files/aug_prod/{omics}_{sample}/all_proteins_clean.fasta",
		eggdb = "results/intermediate_files/eggnog/eggnog.db",
		diamond = "results/intermediate_files/eggnog/eggnog_proteins.dmnd"
	output:
		output = "results/intermediate_files/eggnog/{omics}_{sample}_eggnog.emapper.annotations"
	params:
		prefix = "results/intermediate_files/eggnog/{omics}_{sample}_eggnog"
	conda:
		srcdir("../envs/eggnog.yaml")
	threads: 16
	shell:
		"emapper.py -o {params.prefix} -i {input.sample} --data_dir results/intermediate_files/eggnog/ -m diamond --dmnd_db {input.diamond} --cpu {threads} --resume"

rule pathway_integration:
	input:
		input=expand("results/intermediate_files/eggnog/{omics}_{sample}_eggnog.emapper.annotations", omics=config["omics"], sample=config["MG_samples"])
	output:
		output = "results/final/integrated/pathview/log.txt"
	params:
		param1 = config["parameters"]["group"]
	conda:
		srcdir("../envs/R.yaml")
	shell:
		"Rscript workflow/scripts/pathway_integration.R --group {params.param1}"

rule combi:
	input:
		metadata = "resources/metadata.txt",
		abundance_mg = "results/final/MG/mg_taxa_abundance.txt",
		abundance_mt = "results/final/MT/mt_taxa_abundance.txt",
	output:
		plot = "results/final/integrated/combi_plot.png"
	params:
		param1 = config["parameters"]["group"]
	conda:
		srcdir("../envs/R.yaml")
	shell:
		"Rscript workflow/scripts/visual_integration.R --group {params.param1}"
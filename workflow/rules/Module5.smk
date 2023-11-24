# Module:		5
# Description:	This module processes metagenomics, metatranscriptomics and metaproteomics data

import os
import yaml 

ruleorder: trimPE > trimSE
ruleorder: fastqc_raw_pe > fastqc_raw_se
ruleorder: fastqc_trim_pe > fastqc_trim_se
ruleorder: mergePE > gunzipSE
ruleorder: cat_pe > gunzipSE
ruleorder: rnaspades_pe > rnaspades_se

rule all:
	input:
		qc_raw_report = expand("results/intermediate_files/multiqc/{omics}/multiqc_raw_report.html", omics=config["omics"]),
		qc_trim_report = expand("results/intermediate_files/multiqc/{omics}/multiqc_trim_report.html", omics=config["omics"]),
		interproscan_db = "results/final/prot_db/database_MP.fasta",
		interproscan_setup = "results/intermediate_files/interproscan/interproscan_setup.txt",
		diff_abun_results_mg = "results/final/diff_abun/taxa-maaslin2-MG/maaslin2.log",
		diff_abun_results_mt = "results/final/diff_abun/taxa-maaslin2-MT/maaslin2.log",
		diff_abun_tigrfam_results_mg = "results/final/diff_abun/tigrfam-maaslin2-MG/maaslin2.log",
		diff_abun_tigrfam_results_mt = "results/final/diff_abun/tigrfam-maaslin2-MT/maaslin2.log",
		pathview_results = "results/final/integrated/pathview/log.txt",
		combi_results = "results/final/integrated/combi_plot.svg",
		msgf_db = "results/final/prot_db/db_ready.fasta",
		pep_abundance_table = "results/final/MP/pep_abundance_table.txt",
		peptide_list = "results/final/MP/unipept_list.txt",
		unipept_results = "results/final/MP/unipept_results.csv",
		diff_abun_mp = "results/final/diff_abun/taxa-maaslin2-MP/maaslin2.log"
		
rule fastqc_raw_pe:
	input:
		r1="data/{omics}/raw/{sample}_1.fastq.gz",
		r2="data/{omics}/raw/{sample}_2.fastq.gz"
	output:
		fastqc_html = "results/intermediate_files/fastqc/raw/{omics}/{sample}_fastqc.html"
	params:
		output_folder = "results/intermediate_files/fastqc/raw/{omics}/"
	conda:
		srcdir("../envs/fastqc.yaml")
	threads: 4
	shell:
		"fastqc -t {threads} {input.r1} {input.r2} -o {params} --extract"

rule fastqc_raw_se:
	input:
		r1="data/{omics}/raw/{sample}_1.fastq.gz"
	output:
		fastqc_html = "results/intermediate_files/fastqc/raw/{omics}/{sample}_fastqc.html"
	params:
		output_folder = "results/intermediate_files/fastqc/raw/{omics}/"
	conda:
		srcdir("../envs/fastqc.yaml")
	threads: 4
	shell:
		"fastqc -t {threads} {input.r1} -o {params} --extract"

rule multiqc_raw:
	input:
		sample = expand("results/intermediate_files/fastqc/raw/{omics}/{sample}_fastqc.html", omics=config["omics"], sample=config["MG_samples"])
	output:
		output = "results/intermediate_files/multiqc/{omics}/multiqc_raw_report.html"
	params:
		input_dir = "results/intermediate_files/fastqc/raw/{omics}/",
		output_dir =  "results/intermediate_files/multiqc/{omics}/"
	conda:
		srcdir("../envs/multiqc.yaml")
	shell:
		"multiqc {params.input_dir} -o {params.output_dir}"

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
		r1="data/MT/raw/{sample}_1.fastq.gz"
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
		"trimmomatic SE {input.r1} {output.o1} {params.params}"

rule fastqc_trim_pe:
	input:
		r1="results/intermediate_files/trimmed/{omics}/{omics}_{sample}_1.fastq.gz",
		r2="results/intermediate_files/trimmed/{omics}/{omics}_{sample}_2.fastq.gz",
	output:
		fastqc_html = "results/intermediate_files/fastqc/trim/{omics}/{sample}_fastqc.html"
	params:
		output_folder = "results/intermediate_files/fastqc/trim/{omics}/"
	conda:
		srcdir("../envs/fastqc.yaml")
	threads: 4
	shell:
		"fastqc -t {threads} {input.r1} {input.r2} -o {params} --extract"

rule fastqc_trim_se:
	input:
		r1="results/intermediate_files/trimmed/{omics}/{omics}_{sample}_1.fastq.gz"
	output:
		fastqc_html = "results/intermediate_files/fastqc/trim/{omics}/{sample}_fastqc.html"
	params:
		output_folder = "results/intermediate_files/fastqc/trim/{omics}/"
	conda:
		srcdir("../envs/fastqc.yaml")
	threads: 4
	shell:
		"fastqc -t {threads} {input.r1} -o {params} --extract"

rule multiqc_trim:
	input:
		sample = expand("results/intermediate_files/fastqc/trim/{omics}/{sample}_fastqc.html", omics=config["omics"], sample=config["MG_samples"])
	output:
		output = "results/intermediate_files/multiqc/{omics}/multiqc_trim_report.html"
	params:
		input_dir = "results/intermediate_files/fastqc/trim/{omics}/",
		output_dir =  "results/intermediate_files/multiqc/{omics}"
	conda:
		srcdir("../envs/multiqc.yaml")
	shell:
		"multiqc {params.input_dir} -o {params.output_dir}"

rule gunzipSE:
	input:
		input="results/intermediate_files/trimmed/{omics}/{omics}_{sample}_1.fastq.gz"
	output:
		output="results/intermediate_files/merged/{omics}/{omics}_{sample}/all.fastq"
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
	
rule cat_pe:
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
		nr = "results/intermediate_files/kaiju/kaiju_db/kaiju_db_nr.fmi",
		nodes = "results/intermediate_files/kaiju/kaiju_db/nodes.dmp",
		names = "results/intermediate_files/kaiju/kaiju_db/names.dmp"
	conda:
		srcdir("../envs/kaiju.yaml")
	shell:
		"""
		cd results/intermediate_files/kaiju/kaiju_db/
		wget https://kaiju.binf.ku.dk/database/kaiju_db_nr_2020-05-25.tgz
		tar xzf kaiju_db_nr_2020-05-25.tgz
		"""

rule kaiju_classify:
	input:
		sample="results/intermediate_files/kaiju/kaiju_input/{omics}_{sample}_ns.fastq",
		nodes="results/intermediate_files/kaiju/kaiju_db/nodes.dmp",
		fmi="results/intermediate_files/kaiju/kaiju_db/kaiju_db_nr.fmi"
	output:
		output="results/intermediate_files/kaiju/kaiju_output/{omics}/{omics}_{sample}_ns_kaiju.out"
	conda:
		srcdir("../envs/kaiju.yaml")
	resources:
		gpu=1
	threads: 16
	shell:
		"""
		kaiju -t {input.nodes} -f {input.fmi} -i {input.sample} -o {output} -z {threads}
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
		param1 = config["parameters"]["group"],
		param2 = config["parameters"]["taxa_rank"],
		param3 = config["parameters"]["top_taxa"]
	conda:
		srcdir("../envs/R.yaml")
	shell:
		"Rscript workflow/scripts/taxa_diff_abun_mg.R -g {params.param1} -t {params.param2} -n {params.param3}"

rule diff_abund_MT:
	input:
		input=expand("results/intermediate_files/kaiju/kaiju_output/MT/MT_{sample}_kaiju_summary.tsv", omics=config["omics"], sample=config["MG_samples"])
	output:
		log = "results/final/diff_abun/taxa-maaslin2-MT/maaslin2.log",
		abundance_mt = "results/final/MT/mt_taxa_abundance.txt"
	params:
		param1 = config["parameters"]["group"],
		param2 = config["parameters"]["taxa_rank"],
		param3 = config["parameters"]["top_taxa"]
	conda:
		srcdir("../envs/R.yaml")
	shell:
		"Rscript workflow/scripts/taxa_diff_abun_mt.R -g {params.param1} -t {params.param2} -n {params.param3}"

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

rule rnaspades_pe:
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

rule rnaspades_se:
	input:
		i1="results/intermediate_files/trimmed/MT/MT_{sample}_1.fastq.gz"
	params:
		outdir = "results/intermediate_files/spades/MT_{sample}/"
	output:
		output = "results/intermediate_files/spades/MT_{sample}/transcripts.fasta"
	conda:
		srcdir("../envs/spades.yaml")
	threads: 16
	shell:
		"rnaspades.py -s {input.i1} -o {params.outdir}"

rule seqkit_contigs_mg:
	input:
		contigs = "results/intermediate_files/spades/MG_{sample}/contigs.fasta"
	output:
		output = "results/intermediate_files/seqkit/MG_{sample}/contigs_filtered.fasta"
	conda:
		srcdir("../envs/seqkit.yaml")
	shell:
		"seqkit seq -m 1000 {input} > {output.output}"

rule seqkit_contigs_mt:
	input:
		contigs = "results/intermediate_files/spades/MT_{sample}/transcripts.fasta"
	output:
		output = "results/intermediate_files/seqkit/MT_{sample}/contigs_filtered.fasta"
	conda:
		srcdir("../envs/seqkit.yaml")
	shell:
		"seqkit seq -m 200 {input} > {output.output}"

rule eukrep:
	input:
		input1 = "results/intermediate_files/seqkit/{omics}_{sample}/contigs_filtered.fasta"
	output:
		output1 = "results/intermediate_files/eukrep/{omics}_{sample}/prok_contigs.fasta",
		output2 = "results/intermediate_files/eukrep/{omics}_{sample}/euk_contigs.fasta"
	conda:
		srcdir("../envs/eukrep.yaml")
	shell:
		"EukRep -i {input.input1} -o {output.output2} --prokarya {output.output1}"

rule augustus_mg:
	input:
		input1 = "results/intermediate_files/eukrep/{omics}_{sample}/euk_contigs.fasta"
	output:
		output1 = "results/intermediate_files/augustus/{omics}_{sample}/augustus_output.gff"
	resources:
		gpu=1
	conda:
		srcdir("../envs/augustus.yaml")
	params:
		params = config["parameters"]["augustus"]
	shell:
		"augustus {input.input1} {params.params} > {output.output1}"

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
		contigs = "results/intermediate_files/eukrep/{omics}_{sample}/prok_contigs.fasta"
	output:
		proteins = "results/intermediate_files/prodigal/{omics}_{sample}/Proteinas.fasta",
		genes = "results/intermediate_files/prodigal/{omics}_{sample}/Genes.fasta",
		potential_genes = "results/intermediate_files/prodigal/{omics}_{sample}/PotentialGenes.fasta",
		gbk = "results/intermediate_files/prodigal/{omics}_{sample}/prodigal_output.gbk"
	conda:
		srcdir("../envs/prodigal.yaml")
	shell:
		"prodigal -i {input} -d {output.genes} -s {output.potential_genes} -p meta -o {output.gbk} -a {output.proteins}"

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
		"./results/intermediate_files/interproscan/interproscan*/interproscan.sh -appl NCBIFAM -pa -dra -b {params} -i {input}"

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

rule build_msgf_db:
	input:
		"results/final/prot_db/db_ready.fasta"
	output:
		"results/final/prot_db/db_ready.canno"
	conda:
		srcdir("../envs/msgfplus.yaml")
	shell:
		"msgf_plus -Xmx16000M edu.ucsd.msjava.msdbsearch.BuildSA -d {input}"

rule run_msgf:
	input:
		input = "data/MP/spectra/{spectra}.mgf",
		prot_db = "results/final/prot_db/db_ready.fasta",
		auxiliary = "results/final/prot_db/db_ready.canno"
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
		mp_maaslin2_results = "results/final/diff_abun/taxa-maaslin2-MP/maaslin2.log",
		mp_abundance_plot = "results/final/MP/mp_abundance_plot.svg",
		ec_abundance = "results/final/MP/ec_abundance.txt"
	params:
		outdir = "results/final/diff_abun/taxa-maaslin2-MP/",
		param1 = config["parameters"]["group"],
		param2 = config["parameters"]["taxa_rank"],
		param3 = config["parameters"]["top_taxa"]
	conda:
		srcdir("../envs/R.yaml")
	shell:
		"Rscript workflow/scripts/taxa_diff_abun_mp.R --group {params.param1} --taxa_rank {params.param2} --top_taxa {params.param3}"

rule seqkit_microbial_genes:
	input:
		genes = "results/intermediate_files/prodigal/{omics}_{sample}/Genes.fasta"
	output:
		unique = "results/intermediate_files/prodigal/{omics}_{sample}/Genes_uniq.fasta"
	conda:
		srcdir("../envs/seqkit.yaml")
	shell:
		"seqkit rename -1 -O uniq < {input} > {output}"

rule mark_microbial_gene_seqs:
	input:
		input = "results/intermediate_files/prodigal/{omics}_{sample}/Genes_uniq.fasta"
	output:
		output = "results/intermediate_files/aug_prod/{omics}_{sample}/genes_microbial.fasta"
	shell:
		"sed 's/^>/>microbial_/' {input} > {output}"

rule seqkit_host_genes:
	input:
		sample = "results/intermediate_files/augustus/{omics}_{sample}/augustus_output.codingseq"
	output:
		unique = "results/intermediate_files/augustus/{omics}_{sample}/augustus_host_uniq.fasta"
	conda:
		srcdir("../envs/seqkit.yaml")
	shell:
		"seqkit rename -1 -O uniq < {input} > {output}"

rule mark_host_gene_seqs:
	input:
		input = "results/intermediate_files/augustus/{omics}_{sample}/augustus_host_uniq.fasta"
	output:
		output = "results/intermediate_files/aug_prod/{omics}_{sample}/genes_host.fasta"
	shell:
		"sed 's/^>/>host_/' {input} > {output}"

rule aug_prod_genes:
	input:
		host = "results/intermediate_files/aug_prod/{omics}_{sample}/genes_host.fasta",
		microbial = "results/intermediate_files/aug_prod/{omics}_{sample}/genes_microbial.fasta"
	output:
		output = "results/intermediate_files/aug_prod/{omics}_{sample}/all_genes.fasta"
	shell:
		"cat {input.host} {input.microbial} > {output}"

rule sed_genes:
	input:
		sample = "results/intermediate_files/aug_prod/{omics}_{sample}/all_genes.fasta"
	output:
		output = "results/intermediate_files/aug_prod/{omics}_{sample}/all_genes_clean.fasta"
	shell:
		"sed 's/*//g' {input} > {output}"

rule calculate_coverage:
	input:
		in1 = "results/intermediate_files/trimmed/{omics}/{omics}_{sample}_1.fastq.gz",
		in2 = "results/intermediate_files/trimmed/{omics}/{omics}_{sample}_2.fastq.gz",
		genes = "results/intermediate_files/aug_prod/{omics}_{sample}/all_genes_clean.fasta"
	output:
		coverage = "results/intermediate_files/aug_prod/{omics}_{sample}/stats.txt",
		histogram = "results/intermediate_files/aug_prod/{omics}_{sample}/histogram.txt"
	conda:
		srcdir("../envs/bbmap.yaml")
	shell:
		"bbmap.sh in1={input.in1} in2={input.in2} ref={input.genes} nodisk covstats={output.coverage} covhist={output.histogram}"


rule pathway_integration:
	input:
		input1 = expand("results/intermediate_files/eggnog/{omics}_{sample}_eggnog.emapper.annotations", omics=config["omics"], sample=config["MG_samples"]),
		input2 = expand("results/intermediate_files/aug_prod/{omics}_{sample}/stats.txt", omics=config["omics"], sample=config["MG_samples"]),
		input3 = "results/final/MP/ec_abundance.txt"
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
		abundance_mp = "results/final/MP/pep_abundance_table.txt"
	output:
		plot = "results/final/integrated/combi_plot.svg"
	params:
		param1 = config["parameters"]["group"]
	conda:
		srcdir("../envs/R.yaml")
	shell:
		"Rscript workflow/scripts/visual_integration.R --group {params.param1}"
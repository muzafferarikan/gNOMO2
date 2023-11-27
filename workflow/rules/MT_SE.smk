import os
import yaml 

rule all:
	input:
		qc_raw_report = "results/intermediate_files/multiqc/MT/raw/multiqc_report.html",
		qc_trim_report = "results/intermediate_files/multiqc/MT/trim/multiqc_report.html",
		diff_abun_results_mt = "results/final/diff_abun/taxa-maaslin2-MT/maaslin2.log",
		diff_abun_tigrfam_results_mt = "results/final/diff_abun/tigrfam-maaslin2-MT/maaslin2.log"

rule fastqc_raw_se:
	input:
		r1=expand("data/MT/raw/{sample}_1.fastq.gz", sample=config["MT_samples"])
	output:
		directory("results/intermediate_files/fastqc/raw/MT")
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
		fastqc_html = "results/intermediate_files/fastqc/raw/MT"
	output:
		output = "results/intermediate_files/multiqc/MT/raw/multiqc_report.html"
	params:
		output_dir =  "results/intermediate_files/multiqc/MT/raw/"
	conda:
		srcdir("../envs/multiqc.yaml")
	shell:
		"multiqc {input} -o {params.output_dir}"

rule trim_se:
	input:
		r1="data/MT/raw/{sample}_1.fastq.gz"
	output:
		o1="results/intermediate_files/trimmed/MT/MT_{sample}_1.fastq.gz",
		o2="results/intermediate_files/trimmed/MT/MT_{sample}_2.fastq.gz",
		o1un="results/intermediate_files/trimmed/MT/MT_{sample}_1un.trim.fastq.gz",
		o2un="results/intermediate_files/trimmed/MT/MT_{sample}_2un.trim.fastq.gz"
	params:
		params = config["parameters"]["trimmomatic"]
	conda:
		srcdir("../envs/trimmomatic.yaml")
	threads: 16
	shell:
		"""
		trimmomatic SE -threads {threads} {input.r1} {output.o1} {params.params}
		touch {output.o2} {output.o1un} {output.o2un}
		"""

rule fastqc_trim_se:
	input:
		r1=expand("results/intermediate_files/trimmed/MT/MT_{sample}_1.fastq.gz", sample=config["MT_samples"])
	output:
		directory("results/intermediate_files/fastqc/trim/MT")
	params:
		out="results/intermediate_files/fastqc/trim/MT"
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
		fastqc_html = "results/intermediate_files/fastqc/trim/MT"
	output:
		output = "results/intermediate_files/multiqc/MT/trim/multiqc_report.html"
	params:
		output_dir =  "results/intermediate_files/multiqc/MT/trim/"
	conda:
		srcdir("../envs/multiqc.yaml")
	shell:
		"multiqc {input} -o {params.output_dir}"

rule gunzip_se:
	input:
		input="results/intermediate_files/trimmed/MT/MT_{sample}_1.fastq.gz"
	output:
		f1="results/intermediate_files/merged/MT/MT_{sample}/all.fastq"
	shell:
		"gunzip -c {input.input} > {output.f1}"

rule sed_seqs:
	input:
		input="results/intermediate_files/merged/MT/MT_{sample}/all.fastq"
	output:
		output="results/intermediate_files/kaiju/kaiju_input/MT_{sample}_ns.fastq"
	shell:
		"sed -r 's/\s+//g' {input} > {output}"

rule kaiju_classify:
	input:
		sample="results/intermediate_files/kaiju/kaiju_input/MT_{sample}_ns.fastq",
		nodes="results/intermediate_files/kaiju/kaiju_db/nodes.dmp",
		fmi="results/intermediate_files/kaiju/kaiju_db/kaiju_db_nr.fmi"
	output:
		output="results/intermediate_files/kaiju/kaiju_output/MT/MT_{sample}_ns_kaiju.out"
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
		in1="results/intermediate_files/kaiju/kaiju_output/MT/MT_{sample}_ns_kaiju.out",
		in2="results/intermediate_files/kaiju/kaiju_db/nodes.dmp",
		in3="results/intermediate_files/kaiju/kaiju_db/names.dmp"
	output:
		output="results/intermediate_files/kaiju/kaiju_output/MT/MT_{sample}_kaiju_summary.tsv"
	conda:
		srcdir("../envs/kaiju.yaml")
	shell:
		"""
		kaiju2table -t {input.in2} -n {input.in3} -r genus -o {output} {input.in1} -u
		"""		

rule diff_abund:
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

rule lenfilt_contigs:
	input:
		contigs = "results/intermediate_files/spades/MT_{sample}/transcripts.fasta"
	output:
		output = "results/intermediate_files/seqkit/MT_{sample}/contigs_filtered.fasta"
	conda:
		srcdir("../envs/seqkit.yaml")
	shell:
		"seqkit seq -m 200 {input} > {output.output}"

rule classify_contigs:
	input:
		input1 = "results/intermediate_files/seqkit/MT_{sample}/contigs_filtered.fasta"
	output:
		output1 = "results/intermediate_files/eukrep/MT_{sample}/prok_contigs.fasta",
		output2 = "results/intermediate_files/eukrep/MT_{sample}/euk_contigs.fasta"
	conda:
		srcdir("../envs/eukrep.yaml")
	shell:
		"EukRep -i {input.input1} -o {output.output2} --prokarya {output.output1}"

rule augustus:
	input:
		input1 = "results/intermediate_files/eukrep/MT_{sample}/euk_contigs.fasta"
	output:
		output1 = "results/intermediate_files/augustus/MT_{sample}/augustus_output.gff"
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
		sample = "results/intermediate_files/augustus/MT_{sample}/augustus_output.gff"
	output:
		output = "results/intermediate_files/augustus/MT_{sample}/augustus_output.aa",
		output2 = "results/intermediate_files/augustus/MT_{sample}/augustus_output.codingseq"
	conda:
		srcdir("../envs/perl.yaml")
	shell:
		"perl ./workflow/scripts/get_anno_fasta.pl {input}"

rule rename_host_proteins:
	input:
		sample = "results/intermediate_files/augustus/MT_{sample}/augustus_output.aa"
	output:
		unique = "results/intermediate_files/augustus/MT_{sample}/augustus_output_uniq.fasta"
	conda:
		srcdir("../envs/seqkit.yaml")
	shell:
		"seqkit rename -1 -O uniq < {input} > {output}"

rule mark_host_proteins:
	input:
		input = "results/intermediate_files/augustus/MT_{sample}/augustus_output_uniq.fasta"
	output:
		output = "results/intermediate_files/aug_prod/MT_{sample}/proteins_host.fasta"
	shell:
		"sed 's/^>/>host_/' {input} > {output}"

rule prodigal:
	input:
		contigs = "results/intermediate_files/eukrep/MT_{sample}/prok_contigs.fasta"
	output:
		proteins = "results/intermediate_files/prodigal/MT_{sample}/Proteinas.fasta",
		genes = "results/intermediate_files/prodigal/MT_{sample}/Genes.fasta",
		potential_genes = "results/intermediate_files/prodigal/MT_{sample}/PotentialGenes.fasta",
		gbk = "results/intermediate_files/prodigal/MT_{sample}/prodigal_output.gbk"
	conda:
		srcdir("../envs/prodigal.yaml")
	shell:
		"prodigal -i {input} -d {output.genes} -s {output.potential_genes} -p meta -o {output.gbk} -a {output.proteins}"

rule rename_microbial_proteins:
	input:
		sample = "results/intermediate_files/prodigal/MT_{sample}/Proteinas.fasta"
	output:
		unique = "results/intermediate_files/prodigal/MT_{sample}/Proteinas_uniq.fasta"
	conda:
		srcdir("../envs/seqkit.yaml")
	shell:
		"seqkit rename -1 -O uniq < {input} > {output}"

rule mark_microbial_proteins:
	input:
		input = "results/intermediate_files/prodigal/MT_{sample}/Proteinas_uniq.fasta"
	output:
		output = "results/intermediate_files/aug_prod/MT_{sample}/proteins_microbial.fasta"
	shell:
		"sed 's/^>/>microbial_/' {input} > {output}"

rule merge_aug_prod_proteins:
	input:
		host = "results/intermediate_files/aug_prod/MT_{sample}/proteins_host.fasta",
		microbial = "results/intermediate_files/aug_prod/MT_{sample}/proteins_microbial.fasta"
	output:
		output = "results/intermediate_files/aug_prod/MT_{sample}/all_proteins.fasta"
	shell:
		"cat {input.host} {input.microbial} > {output}"

rule clean_proteins:
	input:
		sample = "results/intermediate_files/aug_prod/MT_{sample}/all_proteins.fasta"
	output:
		output = "results/intermediate_files/aug_prod/MT_{sample}/all_proteins_clean.fasta"
	shell:
		"sed 's/*//g' {input} > {output}"

rule merge_all_samples:
	input:
		sample = expand("results/intermediate_files/aug_prod/MT_{sample}/all_proteins_clean.fasta", omics=config["omics"], sample=config["MG_samples"]),
	output:
		output = "results/final/prot_db/mt_based_database.fa"
	shell:
		"cat {input} > {output}"

rule interproscan:
	input:
		sample = "results/intermediate_files/aug_prod/MT_{sample}/all_proteins_clean.fasta",
		setup = "results/intermediate_files/interproscan/interproscan_setup.txt"
	output:
		output = "results/intermediate_files/interproscan/MT_{sample}/TIGRFAM.tsv"
	params:
		prefix = "results/intermediate_files/interproscan/MT_{sample}/TIGRFAM"
	conda:
		srcdir("../envs/interproscan.yaml")
	shell:
		"./results/intermediate_files/interproscan/interproscan*/interproscan.sh -appl NCBIFAM -pa -dra -b {params} -i {input}"

rule cut_tigrfam:
	input:
		sample = "results/intermediate_files/interproscan/MT_{sample}/TIGRFAM.tsv"
	output:
		interproscan_files = "results/intermediate_files/tigrfam/MT_{sample}_TIGRFAM_v2.tsv"
	shell:
		"cut -f 1,5 {input} > {output}"

rule tigrfam_diff:
	input:
		input=expand("results/intermediate_files/tigrfam/MT_{sample}_TIGRFAM_v2.tsv", omics=config["omics"], sample=config["MG_samples"])
	output:
		diff_abun_tigrfam_results_mt = "results/final/diff_abun/tigrfam-maaslin2-MT/maaslin2.log"
	params:
		outdir = "results/final/diff_abun/",
		param1 = config["parameters"]["group"]
	conda:
		srcdir("../envs/R.yaml")
	shell:
		"Rscript workflow/scripts/tigrfam_diff_abun.R --group {params.param1}"

rule eggnog:
	input:
		sample = "results/intermediate_files/aug_prod/MT_{sample}/all_proteins_clean.fasta",
		eggdb = "results/intermediate_files/eggnog/eggnog.db",
		diamond = "results/intermediate_files/eggnog/eggnog_proteins.dmnd"
	output:
		output = "results/intermediate_files/eggnog/MT_{sample}_eggnog.emapper.annotations"
	params:
		prefix = "results/intermediate_files/eggnog/MT_{sample}_eggnog"
	conda:
		srcdir("../envs/eggnog.yaml")
	threads: 16
	shell:
		"emapper.py -o {params.prefix} -i {input.sample} --data_dir results/intermediate_files/eggnog/ -m diamond --dmnd_db {input.diamond} --cpu {threads} --resume"

rule rename_microbial_genes:
	input:
		genes = "results/intermediate_files/prodigal/MT_{sample}/Genes.fasta"
	output:
		unique = "results/intermediate_files/prodigal/MT_{sample}/Genes_uniq.fasta"
	conda:
		srcdir("../envs/seqkit.yaml")
	shell:
		"seqkit rename -1 -O uniq < {input} > {output}"

rule mark_microbial_genes:
	input:
		input = "results/intermediate_files/prodigal/MT_{sample}/Genes_uniq.fasta"
	output:
		output = "results/intermediate_files/aug_prod/MT_{sample}/genes_microbial.fasta"
	shell:
		"sed 's/^>/>microbial_/' {input} > {output}"

rule rename_host_genes:
	input:
		sample = "results/intermediate_files/augustus/MT_{sample}/augustus_output.codingseq"
	output:
		unique = "results/intermediate_files/augustus/MT_{sample}/augustus_host_uniq.fasta"
	conda:
		srcdir("../envs/seqkit.yaml")
	shell:
		"seqkit rename -1 -O uniq < {input} > {output}"

rule mark_host_genes:
	input:
		input = "results/intermediate_files/augustus/MT_{sample}/augustus_host_uniq.fasta"
	output:
		output = "results/intermediate_files/aug_prod/MT_{sample}/genes_host.fasta"
	shell:
		"sed 's/^>/>host_/' {input} > {output}"

rule merge_aug_prod_genes:
	input:
		host = "results/intermediate_files/aug_prod/MT_{sample}/genes_host.fasta",
		microbial = "results/intermediate_files/aug_prod/MT_{sample}/genes_microbial.fasta"
	output:
		output = "results/intermediate_files/aug_prod/MT_{sample}/all_genes.fasta"
	shell:
		"cat {input.host} {input.microbial} > {output}"

rule clean_genes:
	input:
		sample = "results/intermediate_files/aug_prod/MT_{sample}/all_genes.fasta"
	output:
		output = "results/intermediate_files/aug_prod/MT_{sample}/all_genes_clean.fasta"
	shell:
		"sed 's/*//g' {input} > {output}"

rule calculate_coverage_se:
	input:
		in1 = "results/intermediate_files/trimmed/MT/MT_{sample}_1.fastq.gz",
		genes = "results/intermediate_files/aug_prod/MT_{sample}/all_genes_clean.fasta"
	output:
		coverage = "results/intermediate_files/aug_prod/MT_{sample}/stats.txt",
		histogram = "results/intermediate_files/aug_prod/MT_{sample}/histogram.txt"
	conda:
		srcdir("../envs/bbmap.yaml")
	shell:
		"bbmap.sh in1={input.in1} ref={input.genes} nodisk covstats={output.coverage} covhist={output.histogram}"
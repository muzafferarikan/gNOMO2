import os
import yaml 

rule all:
	input:
		pep_abundance_table = "results/final/MP/pep_abundance_table.txt",
		peptide_list = "results/final/MP/unipept_list.txt",
		unipept_results = "results/final/MP/unipept_results.csv",
		diff_abun_results_mp = "results/final/diff_abun/taxa-maaslin2-MP/maaslin2.log"

rule build_msgf_db:
	input:
		"results/final/prot_db/db_ready.fasta"
	output:
		"results/final/prot_db/db_ready.canno"
	conda:
		srcdir("../envs/msgfplus.yaml")
	shell:
		"msgf_plus -Xmx200G edu.ucsd.msjava.msdbsearch.BuildSA -d {input}"

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

rule diff_abund:
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
		param3 = config["parameters"]["top_taxa"],
		param4 = config["parameters"]["covariates"],
		param5 = config["parameters"]["transformation"],
		param6 = config["parameters"]["normalization"]
	conda:
		srcdir("../envs/R.yaml")
	shell:
		"Rscript workflow/scripts/taxa_diff_abun_mp.R -g {params.param1} -t {params.param2} -n {params.param3} -c {params.param4} -f {params.param5} -m {params.param6}"

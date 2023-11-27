# Module:		4
# Description:	This module processes metagenomics and metatranscriptomics data


from snakemake.utils import min_version
min_version("6.0")

configfile: "config/config.yaml"

layout = config['MT_layout']

def define_rules(layout_type, MT_snakefile):
	module MT:
		snakefile: MT_snakefile
		config: config
	use rule * from MT as MT_*

	module MG:
		snakefile: "MG.smk"
		config: config
	use rule * from MG as MG_*

	rule merge_db:
		input:
			mg_based_db = "results/final/prot_db/mg_based_database.fa",
			mt_based_db = "results/final/prot_db/mt_based_database.fa"
		output:
			merged_db = "results/final/prot_db/merged_db.fa"
		shell:
			"cat {input} > {output}"

	rule seqkit_rmdup:
		input:
			"results/final/prot_db/merged_db.fa"
		output:
			"results/final/prot_db/db_ready.fasta"
		conda:
			srcdir("../envs/seqkit.yaml")
		shell:
			"seqkit rmdup -i -s < {input} > {output}"
			
	rule pathway_integration:
		input:
			input1 = expand("results/intermediate_files/eggnog/{omics}_{sample}_eggnog.emapper.annotations", omics=config["omics"], sample=config["MG_samples"]),
			input2 = expand("results/intermediate_files/aug_prod/{omics}_{sample}/stats.txt", omics=config["omics"], sample=config["MG_samples"])
		output:
			output = "results/final/integrated/pathview/pathview_analysis_summary.txt"
		params:
			param1 = config["parameters"]["group"]
		conda:
			srcdir("../envs/R.yaml")
		shell:
			"Rscript workflow/scripts/pathway_integration.R --group {params.param1}"

	rule visual_integration:
		input:
			metadata = "resources/metadata.txt",
			abundance_mg = "results/final/MG/mg_taxa_abundance.txt",
			abundance_mt = "results/final/MT/mt_taxa_abundance.txt"
		output:
			plot = "results/final/integrated/combi/plot.svg"
		params:
			param1 = config["parameters"]["group"]
		conda:
			srcdir("../envs/R.yaml")
		shell:
			"Rscript workflow/scripts/visual_integration.R --group {params.param1}"

	rule all:
		input:
			rules.MG_all.input,
			rules.MT_all.input,
			"results/final/integrated/combi/plot.svg",
			"results/final/integrated/pathview/pathview_analysis_summary.txt"
		default_target: True

if layout == 'PE':
	define_rules('PE', "MT_PE.smk")

elif layout == 'SE':
	define_rules('SE', "MT_SE.smk")

else:
	print("Please define MT_layout in config file!")
# Module:		6
# Description:	Amplicon Sequencing + Metagenomics + Metatranscriptomics + Metaproteomics

from snakemake.utils import min_version

min_version("6.0")

configfile: "config/config.yaml"

AS_layout = config['AS_layout']
MT_layout = config['MT_layout']

def define_rules(layout1, layout2, AS_snakefile, MT_snakefile):
	module AS:
		snakefile: AS_snakefile
		config: config
	use rule * from AS as AS_*

	module MG:
		snakefile: "MG.smk"
		config: config
	use rule * from MG as MG_*

	module MT:
		snakefile: MT_snakefile
		config: config
	use rule * from MT as MT_*

	rule merge_db:
		input:
			as_based_db = "results/final/prot_db/as_based_database.fa",
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
			input2 = expand("results/intermediate_files/aug_prod/{omics}_{sample}/stats.txt", omics=config["omics"], sample=config["MG_samples"]),
			input3 = "results/final/MP/ec_abundance.txt"
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
			abundance_mt = "results/final/MT/mt_taxa_abundance.txt",
			abundance_mp = "results/final/MP/pep_abundance_table.txt"
		output:
			plot = "results/final/integrated/combi/plot.svg"
		params:
			param1 = config["parameters"]["group"]
		conda:
			srcdir("../envs/R.yaml")
		shell:
			"Rscript workflow/scripts/visual_integration.R --group {params.param1}"

	module MP:
		snakefile: "MP.smk"
		config: config
	use rule * from MP as MP_*

	rule all:
		input:
			rules.AS_all.input,
			rules.MG_all.input,
			rules.MT_all.input,
			rules.MP_all.input,
			"results/final/integrated/combi/plot.svg",
			"results/final/integrated/pathview/pathview_analysis_summary.txt"
		default_target: True

if AS_layout == 'PE' and MT_layout == 'PE':
	define_rules('PE', 'PE', "AS_PE.smk", "MT_PE.smk")

elif AS_layout == 'SE' and MT_layout == 'SE':
	define_rules('SE', 'SE', "AS_SE.smk", "MT_SE.smk")

elif AS_layout == 'SE' and MT_layout == 'PE':
	define_rules('SE', 'PE', "AS_SE.smk", "MT_PE.smk")

elif AS_layout == 'PE' and MT_layout == 'SE':
	define_rules('PE', 'SE', "AS_PE.smk", "MT_SE.smk")

else:
	print("Please define AS_layout and MT_layout in the config file!")

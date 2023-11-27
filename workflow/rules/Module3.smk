# Module:		3
# Description:	Metagenomics + Metaproteomics

from snakemake.utils import min_version
min_version("6.0")

configfile: "config/config.yaml"

module MG:
	snakefile: "MG.smk"
	config: config
use rule * from MG as MG_*

rule seqkit_rmdup:
	input:
		mg_based_db = "results/final/prot_db/mg_based_database.fa"
	output:
		"results/final/prot_db/db_ready.fasta"
	conda:
		srcdir("../envs/seqkit.yaml")
	shell:
		"seqkit rmdup -i -s < {input} > {output}"

module MP:
	snakefile: "MP.smk"
	config: config
use rule * from MP as MP_*

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
		abundance_mp = "results/final/MP/pep_abundance_table.txt",
		metadata = "resources/metadata.txt",
		abundance_mg = "results/final/MG/mg_taxa_abundance.txt"
	output:
		plot = "results/final/integrated/combi_plot.svg"
	conda:
		srcdir("../envs/R.yaml")
	shell:
		"Rscript workflow/scripts/visual_integration.R"

rule all:
	input:
		rules.MG_all.input,
		rules.MP_all.input,
		"results/final/integrated/combi_plot.svg",
		"results/final/integrated/pathview/pathview_analysis_summary.txt"
	default_target: True
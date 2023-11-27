# Module:		2
# Description:	Amplicon Sequencing + Metaproteomics

from snakemake.utils import min_version
min_version("6.0")

configfile: "config/config.yaml"

layout = config['AS_layout']

def define_rules(layout_type, AS_snakefile):
    module AS:
        snakefile: AS_snakefile
        config: config
    use rule * from AS as AS_*

    module MP:
        snakefile: "MP.smk"
        config: config
    use rule * from MP as MP_*

	rule visual_integration:
		input:
			abundance_mp = "results/final/MP/pep_abundance_table.txt",
			metadata = "resources/metadata.txt",
			abundance_as = "results/final/AS/as_physeq.rds"
		output:
			plot = "results/final/integrated/combi_plot.svg"
		params:
			param1 = config["parameters"]["group"]
		conda:
			srcdir("../envs/R.yaml")
		shell:
			"Rscript workflow/scripts/visual_integration.R --group {params.param1}"

	rule all:
		input:
			rules.AS_all.input,
			rules.MP_all.input,
			"results/final/integrated/combi_plot.svg"
		default_target: True

if layout == 'PE':
	define_rules('PE', "AS_PE.smk")

elif layout == 'SE':
	define_rules('SE', "AS_SE.smk")

else:
	print("Please define AS_layout in config file!")
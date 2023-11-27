# Module:		1
# Description:	This module processes amplicon sequencing data 
# and generates a protein database for metaproteomics

from snakemake.utils import min_version
min_version("6.0")

configfile: "config/config.yaml"

layout = config['AS_layout']

if layout == 'PE':
	module AS_PE:
		snakefile: "AS_PE.smk"
		config: config
	use rule * from AS_PE

elif layout == 'SE':
	module AS_SE:
		snakefile: "AS_SE.smk"
		config: config
	use rule * from AS_SE

else:
	print("Please define AS_layout in config file!")
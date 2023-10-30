# gNOMO2
A comprehensive and modular pipeline for integrated multi-omics analysis of microbiomes

If you use this tool, please cite the preprint:
Arikan M, Muth T. (2023) gNOMO2: a comprehensive and modular pipeline for integrated 
multi-omics analysis of microbiomes. bioRxiv. doi: Link

# Installation
1. If you do not have conda installed: Install [conda](https://conda.io/projects/conda/en/latest/user-guide/install/index.html).
2. Create a Snakemake environment in conda:
```
conda create -n snakemake -c bioconda snakemake
```
3. Clone gNOMO2 repository:
```
git clone --recursive https://github.com/muzafferarikan/gNOMO2.git
```

# Setup
## Data
Copy your raw data to the relevant `raw` folder in `data`. For example if you have amplicon sequencing data, copy your files to `data/AS/raw`.  
**Important**: Please make sure you sample names have following structure: 
For AS, MG and MT samples: "samplename"_1.fastq.gz, "samplename"_2.fastq.gz
For MP: "samplename".mgf

## Metadata
gNOMO2 require a metadata file to perform sample group comparisons. Create a tab delimited metadata file containig information about samples and copy it to `resources` folder.  
**Important**: The name of first column in metadata file must be "SampleID".

## Config
After copying your data and metadata files to the relevant folders, run the following script (from your main project folder) to create a config file for your run: 
```
bash workflow/scripts/create_config.sh
```
Now you will find a config.yaml file created inside `config` folder. Check sample names and if you need modify analysis parameters.

## Running
After completing the steps described above, gNOMO2 is ready to run:
1. Activate you snakemake environment in conda:
```
conda activate snakemake
```
2. Run gNOMO2 from your project folder. The value for --cores should be adjusted to reflect the number of cores available. 
```
snakemake -s workflow/Snakefile --cores 2 --use-conda
```

## Running on a cluster
1. Fill in the "gnomo_slurm_template.sh" file in main gNOMO2 folder according to your cluster settings.
2. From your home directory in cluster environemnt run gNOMO2:
```
sbatch path/to/gNOMO2/gnomo_slurm_template.sh
```

# Outputs
A main `results` directory will be created for all gNOMO2 outputs which includes `final` and `intermediate` folders. 
## Final outputs
`final` folder includes final outputs of the pipeline: differential abundance results for each omics dataset, integrated analysis results (joint-visualization - `combi` and pathway level integration results - `pathview`), proteogenomic datadase - `prot_db` and results for each omics datasets (abundance tables, phyloseq objects and plots) which can be used for downstream analyses using other microbiome analysis tools. 

## Intermediate outputs
`intermediate` folder contains results folders for each step executed by gNOMO2. 

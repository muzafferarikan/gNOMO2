<h1 align="center">gNOMO2</h1>
<h2 align="center">A comprehensive and modular pipeline for integrated multi-omics analyses of microbiomes</h2>

If you use this tool, please cite the preprint:  
  
Arikan M, Muth T. (2023) gNOMO2: a comprehensive and modular pipeline for integrated 
multi-omics analyses of microbiomes. bioRxiv. doi: Link

# Table of contents
- [Installation](#installation)
- [Setup](#setup)
    - [Data](#data)
    - [Metadata](#metadata)
    - [Config](#config)
- [Running](#running)
    - [Running locally](#running-locally)
    - [Running on a cluster](#running-on-a-cluster)
- [Outputs](#outputs)
    - [Final outputs](#final-outputs)
    - [Intermediate outputs](#intermediate-outputs)

# Installation
To use gNOMO2, ensure you have `conda` and  `snakemake` installed:  
1. **Install conda**: If you do not have conda installed, [install conda](https://conda.io/projects/conda/en/latest/user-guide/install/index.html).
2. **Create a Snakemake environment in conda**:
```
conda create -n snakemake -c bioconda snakemake
```
3. **Clone gNOMO2 repository**:
```
git clone --recursive https://github.com/muzafferarikan/gNOMO2.git
```
**Note**: Once conda and snakemake are set up, gNOMO2 manages the installation of all other tools and dependencies automatically in their respective environments during the first run. 

# Setup
## Data
Copy your raw data to the relevant subfolder in `data` directory. For example:  
* If you have amplicon sequencing data, copy your files to `data/AS/raw`  
* If you have metaproteomics data copy your files to `data/MP/spectra`.   
**Important**: Please make sure sample names have following format: 
For AS, MG and MT samples: "samplename"_1.fastq.gz, "samplename"_2.fastq.gz  
For MP: "samplename".mgf

## Metadata
gNOMO2 requires a metadata file to perform sample group comparisons. Create a tab delimited metadata file containig information about samples and place it in the `resources` folder.  

**Important**: The name of first column in metadata file must be "SampleID".

## Config
After copying your data and metadata, run the following script from your main project folder to generate a config file:   
```
bash workflow/scripts/prepare_config.sh
```
This script generates a config.yaml file based on contents of "data" directory within `config` folder. Review and modify analysis parameters in this file according to your analysis and metadata.

# Running
## Running locally
Once setup is complete, follow these steps to run gNOMO2:  
1. **Activate your snakemake environment in conda**:
```
conda activate snakemake
```
2. **Run gNOMO2**:  
Execute the following command from your project folder:
```
snakemake -s workflow/Snakefile --cores 2 --use-conda
```
**Note**: Adjust the `--cores` value to reflect the number of cores available.  

## Running on a cluster
To run gNOMO2 on a cluster:  
1. **Configure the cluster settings**:  
Edit the provided "gnomo_slurm_template.sh" file in the main gNOMO2 folder according to your cluster settings.
2. **Run gNOMO2**:  
Execute the following command from your home directory in the cluster environemnt:
```
sbatch path/to/gNOMO2/gnomo_slurm_template.sh
```

# Outputs
When gNOMO2 pipeline starts, it generates a `results` folder within your project directory, containing both `final` and `intermediate` outputs.

## Final outputs
The `final` folder includes:  
* Differential abundance analysis results for each omics dataset (`diff_abun`)  
* Integrated multi-omics analysis results (joint-visualization results (`combi`) and pathway level integration results (`pathview`))  
* A proteogenomic database (`prot_db`)  
* Results for each omics datasets (abundance tables, phyloseq objects and plots) within subfolders named accordingly (`AS`,`MG`,`MT`,`MP`). These files are suitable for further analyses using other microbiome analysis tools. 

## Intermediate outputs
`intermediate` folder contains outputs of each step executed by the gNOMO2 pipeline. 

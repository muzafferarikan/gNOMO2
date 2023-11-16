#!/bin/bash

# Title: Prepare Snakemake Configuration File
# Author: Muzaffer Arikan
# Date: Sep 2023
# Description:
#   This script generates a config.yaml file based on contents of "data" directory.
#   It assigns relevant omics combination and writes sample file names and parameters.
#   User should check if assigned omics combination is correct and all sample files are provided
#   Parameters are listed with default values and can be changed by the user before running gNOMO2. 

# Function to format the files for AS, MG, and MT
format_files() {
    folder="$1"
    sample_name="$2"

    r1_file="data/$folder/raw/${sample_name}_1.fastq.gz"
    r2_file="data/$folder/raw/${sample_name}_2.fastq.gz"

    if [ -f "$r1_file" ] && [ -f "$r2_file" ]; then
        echo "         \"$sample_name\": [r1:\"$r1_file\", r2:\"$r2_file\"]" >> "$output_file"
    fi
}

# Function to format the files for MP
format_mp_files() {
    folder="MP"
    sample_name="$1"

    mgf_file="data/$folder/spectra/${sample_name}.mgf"

    if [ -f "$mgf_file" ]; then
        echo "         \"$sample_name\": $mgf_file" >> "$output_file"
    fi
}

# Check if the 'data' directory exists
if [ -d "data" ]; then
    # Create a 'config' directory if it doesn't exist
    mkdir -p config
    # Create a text file to store the filenames
    output_file="config/config.yaml"
    > "$output_file" # Clears the file or creates a new one

    # Add the required text at the beginning of the file
    cat <<EOF >> "$output_file"
# Please check if corresponding number for the "omics_combination" parameter is correct:
# 1:  Amplicon sequencing
# 2:  Amplicon sequencing + Metaproteomics
# 3:  Metagenomics + Metaproteomics
# 4:  Metagenomics + Metatranscriptomics
# 5:  Metagenomics + Metatranscriptomics + Metaproteomics
# 6:  Amplicon sequencing + Metagenomics + Metatranscriptomics + Metaproteomics

EOF

    # Initialize omics_combination variable
    omics_combination=""

    # Check for non-empty folders
    AS_non_empty=$(ls -A "data/AS/raw")
    MG_non_empty=$(ls -A "data/MG/raw")
    MT_non_empty=$(ls -A "data/MT/raw")
    MP_non_empty=$(ls -A "data/MP/spectra")

    # Determine the omics_combination number based on conditions
    if [ -n "$AS_non_empty" ] && [ -z "$MG_non_empty" ] && [ -z "$MT_non_empty" ] && [ -z "$MP_non_empty" ]; then
        omics_combination="1"
    elif [ -n "$AS_non_empty" ] && [ -z "$MG_non_empty" ] && [ -z "$MT_non_empty" ] && [ -n "$MP_non_empty" ]; then
        omics_combination="2"
    elif [ -z "$AS_non_empty" ] && [ -n "$MG_non_empty" ] && [ -z "$MT_non_empty" ] && [ -n "$MP_non_empty" ]; then
        omics_combination="3"
    elif [ -z "$AS_non_empty" ] && [ -n "$MG_non_empty" ] && [ -n "$MT_non_empty" ] && [ -z "$MP_non_empty" ]; then
        omics_combination="5"
    elif [ -n "$AS_non_empty" ] && [ -n "$MG_non_empty" ] && [ -n "$MT_non_empty" ] && [ -n "$MP_non_empty" ]; then
        omics_combination="6"
    fi


    if [ -n "$omics_combination" ]; then
    	echo "omics_combination: \"$omics_combination\"" >> "$output_file"
    	echo >> "$output_file" # Empty line
   		echo "Module${omics_combination}:" >> "$output_file" # Include the desired text

        case $omics_combination in
            "1")
                cat <<EOF >> "$output_file"
   host: ""
   parameters:
      trimmomatic: "SLIDINGWINDOW:4:20 MINLEN:25 ILLUMINACLIP:resources/adapters.fa:2:30:10"
      group: "Group"
      taxa_rank: "Genus"
      top_taxa: "11"
EOF
                ;;
            "2")
                cat <<EOF >> "$output_file"
   host: ""
   parameters:
      trimmomatic: "SLIDINGWINDOW:4:20 MINLEN:25 ILLUMINACLIP:resources/adapters.fa:2:30:10"
      group: "Group"
      taxa_rank: "Genus"
      top_taxa: "11"
      msgf: "-t 6ppm -tda 1 -m 1 -inst 1 -e 1"
EOF
                ;;
            "3")
                cat <<EOF >> "$output_file"
   omics: ["MG"]
   parameters:
      trimmomatic: "SLIDINGWINDOW:4:20 MINLEN:25 ILLUMINACLIP:resources/adapters.fa:2:30:10"
      augustus: "--species=human --protein=on --gff3=on --uniqueGeneId=true --codingseq=on"
      group: "Group"
      taxa_rank: "Genus"
      top_taxa: "11"
      msgf: "-t 6ppm -tda 1 -m 1 -inst 1 -e 1 -maxMissedCleavages 2"
EOF
                ;;
            "4")
                cat <<EOF >> "$output_file"
   omics: ["MG","MT"]
   parameters:
      trimmomatic: "SLIDINGWINDOW:4:20 MINLEN:25 ILLUMINACLIP:resources/adapters.fa:2:30:10"
      augustus: "--species=human --protein=on --gff3=on --uniqueGeneId=true --codingseq=on"
      group: "Group"
      taxa_rank: "Genus"
      top_taxa: "11"
EOF
                ;;
            "5")
                cat <<EOF >> "$output_file"
   omics: ["MG","MT"]
   parameters:
      trimmomatic: "SLIDINGWINDOW:4:20 MINLEN:25 ILLUMINACLIP:resources/adapters.fa:2:30:10"
      augustus: "--species=human --protein=on --gff3=on --uniqueGeneId=true --codingseq=on"
      group: "Group"
      taxa_rank: "Genus"
      top_taxa: "11"
      msgf: "-t 6ppm -tda 1 -m 1 -inst 1 -e 1 -maxMissedCleavages 2"
EOF
                ;;
            "6")
                cat <<EOF >> "$output_file"
   omics: ["MG","MT"]
   parameters:
      trimmomatic: "SLIDINGWINDOW:4:20 MINLEN:25 ILLUMINACLIP:resources/adapters.fa:2:30:10"
      augustus: "--species=human --protein=on --gff3=on --uniqueGeneId=true --codingseq=on"
      group: "Group"
      taxa_rank: "Genus"
      top_taxa: "11"
      msgf: "-t 6ppm -tda 1 -m 1 -inst 1 -e 1 -maxMissedCleavages 2"
EOF

                ;;
            *)
                echo "Unknown omics_combination: $omics_combination"
                ;;
        esac

        echo "Configuration file successfully generated. See $output_file for details."
    else
        echo "No valid omics combination found."
    fi
else
    echo "'data' directory not found."
fi

        # List files for AS folder
        if [ -n "$AS_non_empty" ]; then
        	echo >> "$output_file" # Empty line
            echo "   AS_Samples:" >> "$output_file"
            for file in data/AS/raw/*_1.fastq.gz; do
                base=$(basename "$file" _1.fastq.gz)
                format_files "AS" "$base"
            done
        fi

        # List files for MG folder
        if [ -n "$MG_non_empty" ]; then
            echo >> "$output_file" # Empty line
            echo "   MG_Samples:" >> "$output_file"
            for file in data/MG/raw/*_1.fastq.gz; do
                base=$(basename "$file" _1.fastq.gz)
                format_files "MG" "$base"
            done
        fi

        # List files for MT folder
        if [ -n "$MT_non_empty" ]; then
            echo >> "$output_file" # Empty line
            echo "   MT_Samples:" >> "$output_file"
            for file in data/MT/raw/*_1.fastq.gz; do
                base=$(basename "$file" _1.fastq.gz)
                format_files "MT" "$base"
            done
        fi

        # List files for MP folder
        if [ -n "$MP_non_empty" ]; then
            echo >> "$output_file" # Empty line
            echo "   MP_Samples:" >> "$output_file"
            for file in data/MP/spectra/*.mgf; do
                base=$(basename "$file" .mgf)
                format_mp_files "$base"
            done
        fi


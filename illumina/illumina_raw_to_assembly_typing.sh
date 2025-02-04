#!/bin/bash
source ~/anaconda3/etc/profile.d/conda.sh

# ensures that conda is initialized in the script's shell session + sets up the necessary environment variables + functions for conda to operate		

conda_envs='~/mambaforge/envs'

	# to execute the script:
# 	a) the whole script = all steps:
# bash name_of_the_script.sh all

# 	b) only specific steps (if you want to execute from 3 onward, you have to specify all the subsequent steps)
#		1: qc_raw = fastqc+multiqc on raw reads
#		2: trim = trim raw reads w/trimmomatic  
#		3: qc_trimmed = fastqc+multiqc on trimmed reads
#		4: organism = determine the organism w/kmerfinder
#		5: assembly = de novo assembly w/SPAdes
#		6: quast = quality check of the assemblies w/quast
#		7: busco = quality check of the assemblies w/busco
#		8: sam_bam_coverage = map reads to ref genome - generate sam and bam files and calculate the covrage of all assemblies
#		9. mlst = ST determination w/mlst
#		10. sistr = serotyping for Salmonella w/sistr
#		11: report = generate a merged report for QC (kmerfinder, avg cov, QUAST, ...)
# bash name_of_the_script.sh qc_raw trim qc_trimmed organism assembly qc_assembly sam_bam_coverage

	
	# CHANGE this:
	# raw data:
raw_reads_dir='path/to/fastq/dir'
	# inside the fastq dir, the .fastq files are in subfolders (organism specific)	

output_dir='path/to/desired/output/dir' && mkdir -p "$output_dir"

run_name='user_defined_run_name'
samplesheet='path/to/samplesheet.tsv'

trimmomatic_q_score="4:20"
	# 4:20 = average qual.score of 4 bases above q20 (99% accuracy) --> we ususally use this
	# 4:30 = average qual.score of 4 bases above q30 (99,9% accuracy)
	# 1:20 = each base q20 (99% accuracy)
	# 1:30 = each base q30	(99,9% accuracy)

	# reference genomes:
ref_genomes_dir="location/of/reference/genomes/for/QUAST"


	# kmerfinder databases:
kmerfinder_db='location/of/kmerfinder/databases/bacteria/bacteria.ATG'
kmerfinder_tax='location/of/kmerfinder/databases/bacteria/bacteria.tax'


	
	
	# DON'T CHANGE
		
	# function to get the current date and time = the date and the start time of each process (and writes it into a log file)
log_start_time() {
	local log_text="$1"
	echo "[$(date)] PROCESS: $log_text" >> $output_dir/${run_name}_assembly_pipeline_report.log
}


update_conda_envs(){
	
	echo -e "\nConda envs available to update: \nqc (fastqc, multiqc), trimmomatic, kmerfinder, spades, quast, busco, minimap2 (minimap2, samtools), mlst, sistr."
	echo -e "\nEnter which env you want to update (if multiple, separate them with a space):"
	read conda_env_user
	
	
	for env in $conda_env_user; do
		
		echo -e "\nUpdating env $env."
		conda activate $conda_envs/$env
		conda env export --from-history
			# if you only want to see which tool was installed (not all packages)
			# the output:
		# name: ~/mambaforge/envs/qc
		# channels:
		#  - conda-forge
		#  - bioconda
		#  - defaults
		# dependencies:
		#  - multiqc
		#  - fastqc
		#  - ca-certificates
		#  - certifi
		#  - openssl
		# prefix: ~/mambaforge/envs/qc
		
		echo -e "\nEnter which tool to update (copy the name from the field - dependencies:"
		read tool_user
		
		for tool in $tool_user; do
		
			echo -e "Updating tool $tool"
			 
			conda update -y $tool_user
			
			echo -e "\n$tool finished updating."
		
		done
			
	done 


}

	# quality check - raw reads
quality_check_raw() {	
	
	echo -e "\nStarting quality check of raw reads (R1 & R2) w/fastqc and multiqc. Location: $raw_reads_dir"
	log_start_time "	QUALITY CHECK OF RAW READS (R1 & R2) w/FastQC AND MultiQC. Location: $raw_reads_dir"
	
	
	conda activate $conda_envs/qc

		# get the versions of the tools and save the output of the command into a variable
	version_fastqc=$(fastqc -v | head -n 1)
		# this command outputs: FastQC v0.12.1
	version_multiqc=$(multiqc --version | head -n 1 | awk '{print "MultiQC v"$3}')
		# the output: MultiQC v1.14
	log_start_time "Tool version: 1. $version_fastqc; 2. $version_multiqc"
	
		
	for organism in $raw_reads_dir/*/; do 
		# $raw_reads_dir/*/ - ensurest it only takes directories (not other files)
		
		organism_name=$(basename $organism)
		mkdir -p $output_dir/01_QC_reads/01_raw_QC/01_fastQC/$organism_name $output_dir/01_QC_reads/01_raw_QC/02_multiQC/$organism_name
		
		for i in $organism/*; do

			fastqc $i -o $output_dir/01_QC_reads/01_raw_QC/01_fastQC/$organism_name -t 32
		
		done
		
		multiqc $output_dir/01_QC_reads/01_raw_QC/01_fastQC/$organism_name -o $output_dir/01_QC_reads/01_raw_QC/02_multiQC/$organism_name
		
		log_start_time "QC finished for $organism_name."
	done
	
	echo -e "\nQuality check finished --> results in $output_dir/01_QC_reads/01_raw_QC"
	log_start_time "QUALITY CHECK OF RAW READS FINISHED --> RESULTS IN: $output_dir/01_QC_reads/01_raw_QC"
}


	# low quality base trimming (trimmomatic):
trim_raw_reads() {
	echo -e "\nStarting trimmomatic (low quality base + adapter trimming)."
	
	IFS=":" read -r sliding_window qscore <<< $trimmomatic_q_score

	log_start_time "	LOW QUALITY BASE + ADAPTER TRIMMING w/trimmomatic (SLIDINGWINDOW:$trimmomatic_q_score - ${sliding_window}bp quality average above Q$qscore + MINLEN:36) ."
	
	conda activate $conda_envs/trimmomatic
	
		# get the versions of the tool and save the output of the command into a variable
	version_trimmomatic=$(trimmomatic -version | head -n 1)
		# this command outputs: 0.39
	log_start_time "Tool version: trimmomatic v$version_trimmomatic"	

	for organism in $raw_reads_dir/*/; do
		organism_name=$(basename $organism)
		mkdir -p $output_dir/02_trimmomatic/$organism_name/paired $output_dir/02_trimmomatic/$organism_name/unpaired 
 
		for fwd in $organism/*R1_001.fastq.gz; do
			rev="${fwd/R1/R2}"
			# fwd and rev = paths to files
			fwd_name=$(basename $fwd _001.fastq.gz)
			# name of the file	
			trimmomatic PE $fwd $rev $output_dir/02_trimmomatic/$organism_name/paired/${fwd_name}_trimmomatic_paired.fastq.gz $output_dir/02_trimmomatic/$organism_name/unpaired/${fwd_name}_trimmomatic_unpaired.fastq.gz $output_dir/02_trimmomatic/$organism_name/paired/${fwd_name/_R1/_R2_trimmomatic_paired.fastq.gz} $output_dir/02_trimmomatic/$organism_name/unpaired/${fwd_name/_R1/_R2_trimmomatic_unpaired.fastq.gz} -phred33 -threads 32 LEADING:3 TRAILING:3 SLIDINGWINDOW:$trimmomatic_q_score MINLEN:36 ILLUMINACLIP:/home/martin/mambaforge/envs/trimmomatic/share/trimmomatic-0.39-2/adapters/NexteraPE-PE.fa:2:30:10
		done
		# if you do not want to cut adapters, remove ILLUMINACLIP
		
		echo -e "\nTrimming finished for $organism_name."
		log_start_time "Trimming finished for $organism_name."
	
	done
		
	echo -e "\nLow quality base trimming finished --> results in $output_dir/02_trimmomatic"
	log_start_time "TRIMMING FINISHED --> RESULTS IN: $output_dir/02_trimmomatic"
}


	# quality check - trimmed reads	
quality_check_trimmed() {
	echo -e "\nStarting quality check of trimmed reads (R1 & R2) w/fastqc and multiqc."
	log_start_time "	QUALITY CHECK OF TRIMMED READS (R1 & R2) w/FastQC and MultiQC."

	conda activate $conda_envs/qc
	
		# get the versions of the tools and save the output of the command into a variable
	version_fastqc=$(fastqc -v | head -n 1)
		# this command outputs: FastQC v0.12.1
	version_multiqc=$(multiqc --version | head -n 1)
		# the output: multiqc, version 1.14
	log_start_time "Tools versions: 1. $version_fastqc &  2. $version_multiqc"

	
	for organism in $output_dir/02_trimmomatic/*; do
		organism_name=$(basename $organism)
		mkdir -p $output_dir/01_QC_reads/02_trimmomatic_QC/01_fastQC/$organism_name $output_dir/01_QC_reads/02_trimmomatic_QC/02_multiQC/$organism_name
 	
		for i in $organism/paired/*; do
		
			fastqc $i -o $output_dir/01_QC_reads/02_trimmomatic_QC/01_fastQC/$organism_name/ -t 32

		done 
	
		multiqc $output_dir/01_QC_reads/02_trimmomatic_QC/01_fastQC/$organism_name/ -o $output_dir/01_QC_reads/02_trimmomatic_QC/02_multiQC/$organism_name/ --force

	done

	echo -e "\nQuality check of trimmed reads finished --> results in $output_dir/01_fastQC/02_trimmomatic_QC"
	log_start_time "FINISHED QUALITY CHECK OF TRIMMED READS --> RESULTS IN: $output_dir/01_fastQC/02_trimmomatic_QC"	
}
	
	
	
	# kmerfinder - determine the organism 
determine_organism() {

	echo -e "\nStarting kmerfinder (determine the organism)."
	log_start_time "	DETERMINING THE ORGANISM w/kmerfinder (bacteria.ATG database)."

	mkdir -p $output_dir/03_kmerfinder
	
	conda activate $conda_envs/kmerfinder	
	
		# get the versions of the tool and save the output of the command into a variable
	version_kmerfinder=$(mamba list kmerfinder | sed -n '4p' | awk '{print $1 " " "v"$2}')
 		# this command outputs: kmerfinder v3.0.2
 		
	log_start_time "Tool version: $version_kmerfinder"
	
	
		# fq	
	for organism in $output_dir/02_trimmomatic/*; do
		# fasta
#	for organism in $output_dir/assembly-all_ONT-illumina/*; do	
	
			# 1. fq	
		for fwd in $organism/paired/*R1*; do
			rev="${fwd/R1/R2}"
			# fwd and rev = paths to files
			fwd_name=$(basename $fwd)
			# name of the file

			time kmerfinder.py -i $fwd $rev -o $output_dir/03_kmerfinder/${fwd_name/trimmomatic_paired.fastq.gz/R2} -db $kmerfinder_db -tax $kmerfinder_tax -x
			
			# 2. fasta:
#		for assembly in $organism/*.fasta; do 
#			ass_name=$(basename $assembly .fasta)
			
#			kmerfinder.py -i $assembly -o $output_dir/03_kmerfinder/$ass_name -db $kmerfinder_db -tax $kmerfinder_tax -x
		
		done
	done

	echo -e "\nGenerating summary of k-mer detected organims ..."
		# -e : enables backslash interpretation --> \n = new line
	log_start_time "Generating summary of k-mer detected organims ..."


		# copy all results.txt (detected organisms) to 1 folder 
	mkdir -p $output_dir/03_kmerfinder/01_results_all
	
	for folder in $output_dir/03_kmerfinder/*; do 
		folder_name=$(basename $folder)
		cp "$folder/results.txt" "$output_dir/03_kmerfinder/01_results_all/${folder_name/_R1_R2/}_kmerfinder.txt"
	done

		
		# create a new txt file with only species names and query cov + pval
		# header:
	echo -e "Sample\tSpecies\tNCBI assembly accession\tQuery Coverage\tp value" > $output_dir/03_kmerfinder/${run_name}_kmerfinder_species_summary.txt


	for i in $output_dir/03_kmerfinder/01_results_all/*_kmerfinder.txt; do	
		
		awk -v samplename=$(basename $i _kmerfinder.txt) 'NR>1 {print samplename "\t" $(NF-1) " " $NF "\t" $1 "\t" $6 "\t" $13}' "$i" >> $output_dir/03_kmerfinder/${run_name}_kmerfinder_species_summary.txt
			
	done

	echo -e "\nKmerfinder process finished --> results in: $output_dir/03_kmerfinder"
	log_start_time "FINISHED Kmerfinder PROCESS --> RESULTS IN: $output_dir/03_kmerfinder"


	
	# using INLINE python code
	
	inputfile=$output_dir/03_kmerfinder/${run_name}_kmerfinder_species_summary.txt
	outputfile=$output_dir/03_kmerfinder/${run_name}_kmerfinder_species_summary_merged.txt	
	
python3 <<END
			# everything between << END and END is treated as Python code (block)
		
		# importing Python libraries 	
import csv
	# to work with csv,txt,tsv files 
	
from collections import defaultdict
	# Python Defaultdict is a container-like dictionaries

# DEFINE input and output files:
input_file = "$inputfile"
output_file = "$outputfile"	


	# DICTIONARY to store merged data
	# Create a defaultdict and passing lambda as default_factory argument:
kmerfinder_dict = defaultdict(lambda: {
    "ncbi": [],
    "coverage": 0
    })


# 1. READ the kmerfinder_species_summary.txt file
with open(input_file, 'r') as infile:
	
	reader = csv.reader(infile, delimiter='\t') 
		# read the file - each line is split into columns based on a delimiter
	next(reader)
		# skips the header

	# 2. PROCESS each row in the file (row[] = columns)	
	for row in reader:
		sample_name = row[0]
		species = row[1]
		ncbi_acc = row[2]
		query_coverage = float(row[3])
		p_value = row[4]
		
		# 3. COMBINE data by sample and species
		key = (sample_name, species)	
			# uniquely identify rows by this key
		kmerfinder_dict[key]["ncbi"].append(ncbi_acc)
		kmerfinder_dict[key]["coverage"] += query_coverage 
		# kmerfinder_dict[key]["pvalue"] = p_value	
		

	
# 4. WRITE the merged data to the output file
with open(output_file, 'w') as outfile:
	
	writer = csv.writer(outfile, delimiter='\t')
	# write header:
	writer.writerow(["Filename", "Species (kmerfinder)", "NCBI assembly accession numbers", "Sum of query coverages"])
		# writer.writerow() expects exactly one argument --> if multiple values, store them into a list: []	
	
	for key, values in kmerfinder_dict.items():
		# values: ncbi & coverage
		sample_name, species = key
		ncbi_merged = ", ".join(values["ncbi"])
		query_coverage_sum = values["coverage"]
		# pvalue_last = values("pvalue")
		writer.writerow([sample_name, species, ncbi_merged, query_coverage_sum])



print(f"Processed file written to: {output_file}")
	
END



}





	# assembly w/SPAdes for each organism in a separate folder (= $i) inside the trimmomatic folder
denovo_assembly() {
	echo -e "\nStarting SPAdes (de-novo assembly)."
	log_start_time "	GENERATING DE-NOVO ASSEMBLY w/SPAdes."

	conda activate $conda_envs/spades

	version_spades=$(spades.py -v | head -n 1)
		# the output: SPAdes v3.13.1
	log_start_time "Tool version: $version_spades"
	
	
	for organism in $output_dir/02_trimmomatic/*; do 

		organism_name=$(basename $organism)	
		mkdir -p $output_dir/04_spades/$organism_name $output_dir/04_spades/01_scaffolds_all/$organism_name	
		
		
		for fwd in $output_dir/02_trimmomatic/$organism_name/paired/*R1*; do
			rev="${fwd/_R1_trimmomatic/_R2_trimmomatic}"
			fwd_name=$(basename $fwd)
			
			time spades.py -1 $fwd -2 $rev --isolate --cov-cutoff auto -o $output_dir/04_spades/$organism_name/${fwd_name/_R1_trimmomatic_paired.fastq.gz/_SPAdes} -t 32
				# --isolate mode does not include read error correction by default, BUT IF your data has high uniform coverage depth, it is strongly recommended to use --isolate option.
			
			cp "$output_dir/04_spades/$organism_name/${fwd_name/_R1_trimmomatic_paired.fastq.gz/_SPAdes}/scaffolds.fasta" "$output_dir/04_spades/01_scaffolds_all/$organism_name/${fwd_name/_R1_trimmomatic_paired.fastq.gz/_scaffolds.fasta}" 			
			
		done
	
		log_start_time "Assembly process finished for $organism_name."

	done
	
	echo -e "\nSPAdes process finished --> assemblies in: $output_dir/04_spades/01_scaffolds_all"
	log_start_time "ASSEMBLY GENERATION FINISHED --> ASSEMBLIES IN: $output_dir/04_spades/01_scaffolds_all"
}


	# quality assesment - assemblies
assembly_quality_check_QUAST() {
	echo -e "\nStarting quality check of the assemblies w/QUAST."
	log_start_time "	QUALITY CHECK OF THE ASSEMBLIES w/QUAST."
	
	conda activate $conda_envs/quast
	
	version_quast=$(quast -v | head -n 1)
		# the output: QUAST v5.2.0
	log_start_time "Tool version: $version_quast"

	for i in $output_dir/04_spades/01_scaffolds_all/*; do

		organism_name=$(basename $i)
		mkdir -p $output_dir/05_QC_assembly/01_QUAST/$organism_name
		
	# IF YOU WANT TO MAKE IT EASIER FOR YOU, just input your ref genome paths here and uncomment this section:
#		echo -e "\nPlease enter the path to the most suitable reference genome for $organism_name --> some of the genomes in: $ref_genomes_dir: 
#				\nBordetella_holmesii_ASM962735v1/GCF_009627355.1_ASM962735v1_genomic.fna
#				\nBordetella_pertussis_ASM400897v1/GCF_004008975.1_ASM400897v1_genomic.fna
#				\nCampylobacter_coli_ASM973039v1/GCF_009730395.1_ASM973039v1_genomic.fna
#				\nCampylobacter_jejuni_ASM908v1/GCF_000009085.1_ASM908v1_genomic.fna
#				\nEnterobacter_cloacae_AI2999v1_cpp/GCF_905331265.2_AI2999v1_cpp_genomic.fna
#				\nEscherichia_coli_ASM886v2/GCF_000008865.2_ASM886v2_genomic.fna
#				\nHaemophillus_influenzae_ASM93157v1/GCF_000931575.1_ASM93157v1_genomic.fna
#				\nKlebsiella_pneumoniae_ASM24018v2/GCF_000240185.1_ASM24018v2_genomic.fna
#				\nListeria_innocua_ASM964857v1/GCF_009648575.1_ASM964857v1_genomic.fna
#				\nListeria_moncytogenes_ASM19603v1/GCF_000196035.1_ASM19603v1_genomic.fna
#				\nNeisseria_meningitidis_ASM833080v1/GCF_008330805.1_ASM833080v1_genomic.fna
#				\nSalmonella_enterica_ASM694v2/GCF_000006945.2_ASM694v2_genomic.fna
#				\nShigella_boydii_ASM229048v1/GCF_002290485.1_ASM229048v1_genomic.fna
#				\nShigella_flexneri_ASM692v2/GCF_000006925.2_ASM692v2_genomic.fna
#				\nShigella_sonnei_ASM1337481v1/GCF_013374815.1_ASM1337481v1_genomic.fna
#				\nStreptococcus_pneumoniae_NCTC7465/GCF_001457635.1_NCTC7465_genomic.fna
#				\nStreptococcus_pyogenes_41965_F01/GCF_900475035.1_41965_F01_genomic.fna
#				\nYersinia_enterocolitica_ASM2575863v1/GCF_025758635.1_ASM2575863v1_genomic.fna
#				\npath to $organism_name:"
		
	# IF YOU UNCOMMENTED THE SECTION ABOVE, THEN "COMMENT THE NEXT LINE":	
		echo -e "\nPlease enter the path to the most suitable reference genome for $organism_name"
		
		read ref_genome_user 

		echo -e "\nYou have entered: $ref_genome_user\n"
		log_start_time "User entered: $ref_genome_user as the most suitable reference genome for quality check"
	

		quast.py $output_dir/04_spades/01_scaffolds_all/$organism_name/*.fasta -R $ref_genomes_dir/$ref_genome_user -o $output_dir/05_QC_assembly/01_QUAST/$organism_name
		
		cp $output_dir/05_QC_assembly/01_QUAST/$organism_name/transposed_report.tsv $output_dir/05_QC_assembly/01_QUAST/QUAST_report_transposed_${run_name}_$organism_name.tsv
		cp $output_dir/05_QC_assembly/01_QUAST/$organism_name/report.html $output_dir/05_QC_assembly/01_QUAST/QUAST_report_${run_name}_$organism_name.html	

		echo -e "\nFinished QUAST quality check of assemblies --> results in $output_dir/05_QC_assembly/01_QUAST/$organism_name"
		log_start_time "FINISHED QUAST QUALITY CHECK OF THE $organism_name ASSEMBLIES"		

	done

	log_start_time "FINISHED QUAST QUALITY CHECK --> RESULTS IN $output_dir/05_QC_assembly/01_QUAST"

}




assembly_quality_check_BUSCO() {
	echo -e "\nStarting quality check of the assemblies w/BUSCO."
	log_start_time "	QUALITY CHECK OF THE ASSEMBLIES w/BUSCO."
			
	conda activate $conda_envs/busco
	version_busco=$(busco -v | head -n 1 | awk '{print $1 " " "v"$2}')
		# BUSCO 5.8.2
	log_start_time "Tool version: $version_busco"
	
	
	mkdir -p $output_dir/05_QC_assembly/02_BUSCO
	cd $output_dir/05_QC_assembly/02_BUSCO
		# BUSCO uses pwd as a work environment - it outputs the files here, even though you provide the path to the desired output 

	for organism in $output_dir/04_spades/01_scaffolds_all/*; do
		organism_name=$(basename $organism)

			
		for assembly in $organism/*.fasta; do
			ass_name=$(basename $assembly _scaffolds.fasta)
			
			busco -i $assembly -m genome --auto-lineage --cpu 32 -o "$organism_name/${ass_name}_BUSCO"
		# -m genome --> The genome mode is used to assess the completeness of genome assemblies. It uses the --mode option set to genome. The input file should be a nucleotide fasta file. There are distinct pipelines for eukaryotic and prokaryotic genomes. BUSCO decides which pipeline to use depending on the dataset selected. Each pipeline uses a different gene predictor tool, explained below. The predicted genes are then passed to HMMER, which scores them against the HMM profiles for the BUSCO marker genes in a given dataset.
		#  Recommended parameters:
		# -l or --lineage_dataset --> Specify the name of the BUSCO lineage dataset to be used, e.g. kitasatospora_odb12. A full list of available datasets can be viewed by entering busco --list-datasets --> you need all from -bacteria_odb12. You should always select the dataset that is most closely related to the assembly or gene set you are assessing. If you are unsure, you can use the --auto-lineage option to automatically select the most appropriate dataset. BUSCO will automatically download the requested dataset if it is not already present in the download folder. 
	
		done

	
		echo -e "\nFinished BUSCO quality check of assemblies --> results in $output_dir/05_QC_assembly/02_BUSCO/$organism_name"
		log_start_time "Finished BUSCO quality check of the $organism_name assemblies." 		

	
		# create a merged report
		log_start_time "Creating merged reports for $organism_name."
		merged_report=$output_dir/05_QC_assembly/02_BUSCO/${organism_name}_merged_BUSCO_report.xlsx
		merged_report_filtered=$output_dir/05_QC_assembly/02_BUSCO/${organism_name}_merged_BUSCO_report_filtered.xlsx
		# for python code (a temp file to hold merged data)

python3 <<PYTHON

import json
import pandas as pd
import os
import glob # for iterating through files


# initialize the list to store merged data
merged_data = [] # list to store dfs
merged_data_filtered = []

# go through all .json files
# 1. list all files including subfolders
json_files = glob.glob("$output_dir/05_QC_assembly/02_BUSCO/$organism_name/**/*specific*.json", recursive=True)


# read and process    
for json_f in json_files:
    with open(json_f, 'r') as f:
        data = json.load(f)

    # extract sample name from "parameters" + extract "results" and "metrics" sections
    parameters = data.get("parameters", {})
    # the name of the sample is in the section: "out": "Escherichia_coli/8941-EC_S259_BUSCO"
    out_field = parameters.get("out", "")
    sample_name = out_field.split("/")[-1].replace("_BUSCO", "")
	# split() function with : as delimiter to split the given string into smaller strings
	# string_to_be_split.split(":")
	# [-1] --> take the last part of the split list (=8941-EC_S259_BUSCO)
    
    results = data.get("results", {})
    # define specific keys to be extracted (filter) + rename them (if you want)
    # results_keys{old_keys: new_keys}
    results_keys = {
        "one_line_summary": "Summary",
        "Complete percentage": "Complete BUSCOs [%]",
        "Complete BUSCOs": "Complete BUSCOs (C)",
        "Single copy BUSCOs": "Complete and single-copy BUSCOs (S)",
        "Multi copy BUSCOs": "Complete and duplicated BUSCOs (D)",
        "Fragmented BUSCOs": "Fragmented BUSCOs (F)",
        "Missing BUSCOs": "Missing BUSCOs (M)",
        "n_markers": "Total BUSCO groups searched",
        "Number of scaffolds": "Number of scaffolds",
        "Number of contigs": "Number of contigs",
        "Total length": "Total length",
        "Percent gaps": "Gaps [%]",
        "Scaffold N50": "Scaffold N50",
        "Contigs N50": "Contigs N50"
    }
    
    # EXTRACT + RENAME the old_keys from results object 
    results_filtered = {new_key: results.get(old_key) for old_key, new_key in results_keys.items()}

    
    # REMOVE % from value in json field "Percent gaps"
    if "Gaps [%]" in results_filtered:
        results_filtered["Gaps [%]"] = float(results_filtered["Gaps [%]"].strip('%'))
    
    # 1. combine unfiltered data into a DataFrame
    df = pd.DataFrame([{"Sample name": sample_name, **results}])
    	# append dataframe to the global list
    merged_data.append(df)
    # 2. combine filtered data:
    df_filtered = pd.DataFrame([{"Sample name": sample_name, **results_filtered}])
    merged_data_filtered.append(df_filtered)


    # combine all dfs + save merged dfs to excel file
final_merged_df = pd.concat(merged_data, ignore_index=True)
final_merged_df.to_excel("$merged_report", index=False)
final_merged_df_filtered = pd.concat(merged_data_filtered, ignore_index=True)
final_merged_df_filtered.to_excel("$merged_report_filtered", index=False)
print(f"Saved: $merged_report & $merged_report_filtered")


PYTHON
	log_start_time "Finished merging reports for $organism_name."
	done
	
	rm -R $output_dir/05_QC_assembly/02_BUSCO/busco_downloads
	rm $output_dir/05_QC_assembly/02_BUSCO/*.log


	log_start_time "FINISHED BUSCO QUALITY CHECK --> RESULTS IN $output_dir/05_QC_assembly/02_BUSCO"
}


	# for assembly coverage -->
	# you need to map reads to assembly (use minimap2)
mapping_reads_to_assembly_bam_sam_coverage() {
	echo -e "\nStarting process: mapping raw reads to assemblies (generation of .sam and .bam files + calculating coverage of assemblies)."
	log_start_time "	MAPPING TRIMMED READS TO THE ASSEMBLIES w/minimap2 & samtools (generation of .sam and .bam files + calculating coverage of assemblies)."
	
	conda activate $conda_envs/minimap2
	
	version_minimap2=$(minimap2 --version | head -n 1)
		# the output: 2.26-r1175
	version_samtools=$(samtools --version | head -n 1 | awk '{print $1 " " "v"$2}')
		# the output: samtools v1.18
	log_start_time "Tools versions: 1. minimap2 v$version_minimap2; 2. $version_samtools"

	
	for organism in $output_dir/04_spades/01_scaffolds_all/*; do
		organism_name=$(basename $organism)
		mkdir -p $output_dir/06_minimap2/$organism_name
	
		for assembly in $organism/*; do
			assembly_name=$(basename $assembly)
			samplename=${assembly_name/_scaffolds.fasta/}
			R1=$output_dir/02_trimmomatic/$organism_name/paired/${samplename}_R1_trimmomatic_paired.fastq.gz
			R2=$output_dir/02_trimmomatic/$organism_name/paired/${samplename}_R2_trimmomatic_paired.fastq.gz
	
				# map the reads to assembly
			minimap2 -ax sr $assembly $R1 $R2 > $output_dir/06_minimap2/$organism_name/${samplename}_assembly.sam
				# -ax sr : the preset for short paired-end reads (Illumina)

				# convert to BAM + sort + index 
			samtools view -S -b $output_dir/06_minimap2/$organism_name/${samplename}_assembly.sam > $output_dir/06_minimap2/$organism_name/${samplename}_assembly.bam
			samtools sort $output_dir/06_minimap2/$organism_name/${samplename}_assembly.bam -o $output_dir/06_minimap2/$organism_name/${samplename}_assembly_sorted.bam
			samtools index $output_dir/06_minimap2/$organism_name/${samplename}_assembly_sorted.bam
	
				# calculate coverage
			samtools depth $output_dir/06_minimap2/$organism_name/${samplename}_assembly_sorted.bam > $output_dir/06_minimap2/$organism_name/${samplename}_assembly_coverage.txt
		done

		echo -e "\nFinished mapping trimmed reads to the assemblies --> results in $output_dir/06_minimap2/$organism_name"
		log_start_time "Mapping of the trimmed reads to the $organism_name assemblies finished --> results in $output_dir/06_minimap2/$organism_name"

		
		echo -e "\nCalculating coverage of the assemblies ..."
		log_start_time "Calculating coverage of the assemblies ..."
		
			# create coverage report
			# header:
		echo -e "sample\taverage_coverage_of_assembly\torganism" > $output_dir/06_minimap2/$organism_name/avg_cov_assembly_"$organism_name"_summary.txt	
	

		for i in $output_dir/06_minimap2/$organism_name/*_coverage.txt; do
				
			i_name=$(basename $i)
			# calculate the average coverage using awk and save it to a variable
			average=$(awk '{sum+=$3; count++} END {print sum/count}' $i)

			# print the filename (without extension) + the average coverage to the output file (multiple assembly coverages in 1 file) + the organism name
			echo "${i_name/_coverage.txt/}	$average	$organism_name" >> $output_dir/06_minimap2/$organism_name/avg_cov_assembly_"$organism_name"_summary.txt
		done

			# sort in descending order 
		(head -n 1 $output_dir/06_minimap2/$organism_name/avg_cov_assembly_"$organism_name"_summary.txt && tail -n +2 $output_dir/06_minimap2/$organism_name/avg_cov_assembly_"$organism_name"_summary.txt | sort -k2,2gr) > $output_dir/06_minimap2/$organism_name/avg_cov_assembly_"$organism_name"_summary_sorted.txt
			
		cp $output_dir/06_minimap2/$organism_name/avg_cov_assembly_"$organism_name"_summary.txt $output_dir/06_minimap2/avg_cov_assembly_"$organism_name"_summary.txt
		
		
		
		echo -e "\nFinished calculating coverage of the $organism_name assemblies --> results in $output_dir/06_minimap2/$organism_name"
		log_start_time "Finished calculating coverage of the $organism_name assemblies."
	done
	
		# merge all coverages from all of the species into a new file
		# header:
	echo -e "sample\taverage_coverage_of_assembly\torganism" > $output_dir/06_minimap2/${run_name}_avg_cov_assembly_all.txt

	for txt in $output_dir/06_minimap2/*_summary.txt; do
		
		tail -n +2 $txt >> $output_dir/06_minimap2/${run_name}_avg_cov_assembly_all.txt
		rm $txt
	
	done
	
		# sort in descending order
	(head -n 1 $output_dir/06_minimap2/${run_name}_avg_cov_assembly_all.txt && tail -n +2 $output_dir/06_minimap2/${run_name}_avg_cov_assembly_all.txt | sort -k2,2gr) > $output_dir/06_minimap2/${run_name}_avg_cov_assembly_all_sorted.txt
	
	log_start_time "FINISHED MAPPING READS TO THE ASSEMBLIES --> results in $output_dir/06_minimap2"
	

}

mlst_all(){

	echo -e "\nStarting sequence typing (determine the ST) w/mlst."
	log_start_time "	DETERMINING SEQUENCE TYPE (ST) w/mlst."

	mkdir -p $output_dir/07_MLST

		# determine ST:
	conda activate $conda_envs/mlst	
	
		# get the versions of the tool and save the output of the command into a variable
	version_mlst=$(mlst -v | head -n 1 | awk '{print $1 " " "v"$2}')
		# the output: mlst v2.23.0
	log_start_time "Tool version: $version_mlst"

	
	mlst_report_merged="$output_dir/07_MLST/${run_name}_MLST_report.tsv"
	echo -e "filename\tmatching PubMLST scheme name\tST (sequence type)\tallele ID\tallele ID\tallele ID\tallele ID\tallele ID\tallele ID\tallele ID" > $mlst_report_merged
	
	for organism in $output_dir/04_spades/01_scaffolds_all/*; do
		organism_name=$(basename $organism)
		mlst_report="$output_dir/07_MLST/${organism_name}_MLST_report.tsv"

		
		echo -e "filename\tmatching PubMLST scheme name\tST (sequence type)\tallele ID\tallele ID\tallele ID\tallele ID\tallele ID\tallele ID\tallele ID" > $mlst_report 		
		
		for assembly in $organism/*.fasta; do
			samplename=$(basename $assembly _scaffolds.fasta)
			
			echo -e "\nSeq typing sample:\n$assembly"
		
			mlst --label $samplename $assembly | tee -a $mlst_report $mlst_report_merged	
			
		done	
	done


	echo -e "\nFinished sequence typing (ST) w/mlst."
	log_start_time "FINISHED SEQUENCE TYPING (ST) w/mlst --> RESULTS IN: $output_dir/07_MLST"
	
}




sistr_serotyping(){


	echo -e "\nStarting serotyping w/SISTR.\n"
	log_start_time "	SEROTYPING w/SISTR STARTED."


	conda activate $conda_envs/sistr
	
	version_sistr=$(sistr -V | head -n 1 | awk '{print "sistr v"$2}')
		# the output: sistr_cmd v1.1.2
	log_start_time "Tool version: $version_sistr"
	
	# find if there is a folder with Salmonella in its name:
	salmonella_dir=$(find $output_dir/04_spades/01_scaffolds_all -type d -name "*Salmonella*")

		
	if [[ -n $salmonella_dir ]] ; then 
	
		mkdir -p "$output_dir/08_serotyping/01_SISTR/reports"
		
		for organism in $salmonella_dir; do				
		
			for assembly in $organism/*.fasta; do

				
				sistr --qc -f tab -o "$output_dir/08_serotyping/01_SISTR/reports/$(basename $assembly .fasta)_sistr-output.tab" $assembly
				
				echo "$(basename $assembly) processed."
				log_start_time "$(basename $assembly) processed."

			done
		
		log_start_time "Sistr process finished."		
		
		done


			#create report:
		
		log_start_time "Creating a merged SISTER report."	
			
		merged_report="$output_dir/08_serotyping/01_SISTR/serotyping_report_SISTR_all_${run_name}_Salmonella.txt"

		echo -e "Serotyping results\n" >  $merged_report

		for report in $output_dir/08_serotyping/01_SISTR/reports/*.tab; do 
		
			assembly_name=$(awk -F'\t' 'NR==2 {print $8}' $report)
			serovar=$(awk -F'\t' 'NR==2 {print $15}' $report)
			serogroup=$(awk -F'\t' 'NR==2 {print $14}' $report)
			ser_antigen=$(awk -F'\t' 'NR==2 {print $16}' $report)
			h1=$(awk -F'\t' 'NR==2 {print $9}' $report)
			h2=$(awk -F'\t' 'NR==2 {print $10}' $report)
			o_antigen=$(awk -F'\t' 'NR==2 {print $11}' $report)
			qc_stat=$(awk -F'\t' 'NR==2 {print $13}' $report)
			qc_info=$(awk -F'\t' 'NR==2 {print $12}' $report)
		
			echo -e "\nSample:	$assembly_name" >> $merged_report	
			echo -e "Serovar:	$serovar (serogroup $serogroup)" >> $merged_report
			echo -e "Serovar antigen:	$ser_antigen" >> $merged_report
			echo -e "Antigen formula" >> $merged_report
			echo -e " h1:	$h1" >> $merged_report
			echo -e " h2:	$h2" >> $merged_report
			echo -e " o_antigen:	$o_antigen" >> $merged_report
		
				
			if [ $qc_stat == "WARNING" ]; then
				echo -e "quality check message:	$qc_info" >> $merged_report
			fi	
		# *Serovar:* Choleraesuis (*serogroup* C1)  

		# *Serovar antigen:* Hissar | Choleraesuis | Paratyphi C | Typhisuis | Chiredzi  

		# *Antigen formula:*

		# - h1: c
		# - h2: 1,5
		# - o_antigen: 6,7

		done


		declare -A dictionary

		dictionary=(
			[8]="Sample:"
			[15]="Serovar:"
			[14]="Serogroup:"
			[16]="Serovar antigen:"
			[100]="Antigen formula"
			[9]=" h1:"
			[10]=" h2:"
			[11]=" o_antigen:"
			[13]="Quality check:"
			[12]="Quality check message:"
			)
			
			
		column_numbers="8 15 14 16 100 9 10 11 13 12"	
			
		output=""

			
		for column in $column_numbers; do
				
		new_line_in_output_file="${dictionary[$column]}\t"
			
			for report in $output_dir/08_serotyping/01_SISTR/reports/*.tab; do
				
				if [ $column = "100" ]; then
					new_line_in_output_file+="\t"
				
				else

					value=$(awk -F'\t' -v col=$column 'NR==2 {print $col}' $report)
					new_line_in_output_file+="$value\t"	
				fi
				
			done
			
			output+="${new_line_in_output_file}\n"
			
		done

		echo -e $output >> "$output_dir/08_serotyping/01_SISTR/serotyping_report_SISTR_all_for_excel_${run_name}_Salmonella.txt"


#convert .txt to .xlsx		
python3 <<END

input = "$output_dir/08_serotyping/01_SISTR/serotyping_report_SISTR_all_for_excel_${run_name}_Salmonella.txt"
output = "$output_dir/08_serotyping/01_SISTR/serotyping_report_SISTR_all_for_excel_${run_name}_Salmonella.xlsx"

import pandas as pd

df_report = pd.read_csv(input, sep="\t")

df_report.to_excel(output, index=False)



END

		rm "$output_dir/08_serotyping/01_SISTR/serotyping_report_SISTR_all_for_excel_${run_name}_Salmonella.txt"
		
		echo -e "Merged report in: $output_dir/08_serotyping/01_SISTR"	
		log_start_time "FINISHED SEROTYPING w/SISTR --> RESULTS IN: $output_dir/08_serotyping/01_SISTR"
	
	else 
	
		echo -e "No directory named Salmonella (serotyping w/SISTR impossible)"
		log_start_time "No directory named Salmonella (serotyping w/SISTR impossible)"
	fi	
		

		
}

merge_all_reports() {

	
	log_start_time "	CREATING MERGED QC REPORT."
	
	# python code
		
python3 <<END

# python script for merging all reports
# run script:
# python3 name_of_the_script.py

import os
import re
import glob
import pandas as pd
from collections import defaultdict
	# for left align all cells + set the same font:
from openpyxl import load_workbook
from openpyxl.styles import Alignment, Font
	# for checking if a file exists
from pathlib import Path


output_dir = "$output_dir"	#take the var from the main bash script
run_name = "$run_name"
samples_S = "$samplesheet"

# 1.1
#  initialize paths and directories (define and set up file paths and dirs)
def init_paths():
    global input_files, merged_output_QC       # global vars -> can be accessed throughout the script (all functions); input_files a DICT w/all the reports + merged_output_QC
    # all reports:
    input_files = {                           
        "samples_S": samples_S,    
        "kmerfinder_report": os.path.join(output_dir, "03_kmerfinder", f"{run_name}_kmerfinder_species_summary_merged.txt"),
        	# os.path.join automatically handles the correct path separator based on the OS    
        "MLST_report": os.path.join(output_dir, "07_MLST", f"{run_name}_MLST_report.tsv"),
        "serotype_report": os.path.join(output_dir, "08_serotyping", "01_SISTR", f"serotyping_report_SISTR_all_for_excel_{run_name}_Salmonella.xlsx"),
        "coverages_report": os.path.join(output_dir, "06_minimap2", f"{run_name}_avg_cov_assembly_all.txt"),
        "QUAST_dir": os.path.join(output_dir, "05_QC_assembly", "01_QUAST"),
        "BUSCO_dir": os.path.join(output_dir, "05_QC_assembly", "02_BUSCO")
    }
    # merged output report
    merged_output_QC = os.path.join(output_dir, f"{run_name}_merged_QC.xlsx")



# 1.2
# initialize data structure
def init_data_structure():
    global sample_merged_QC_dict 
    sample_merged_QC_dict = defaultdict(lambda: {       # defaultdict = a type of dict, that provides a default value for a missing key --> when you try to access a key that doesn't exist it creates it automatically and assigns a default value 
    							  # a function(lambda) defines the default value for any key --> EFFICIENT: you don't need to define all possible sample keys bc the structure is created on the fly when you access a new key
            # samplenames from S:
        "sample_notes": "",       
            # kmer report:
        "seq_samplename": "failed species identification",    
        "species": "N/A",
        "ncbi_acc": "N/A",
        "sum_coverage_ncbi": 0.0,
            # MLST report
        "ST_mlst": "failed mlst",
        "MLST_scheme": "N/A",
            # serotyping report
        "serovar": "N/A",
        "serogroup": "N/A",
            # coverages report
        "avg_cov_assembly": "failed assembly generation",
        "organism_folder": "N/A",
            # QUAST report
        "quast_num_contigs": 0,
        "quast_ref_length": 0,
        "quast_total_length": 0,
        "quast_ref_GC": 0.0,
        "quast_GC": 0.0,
        "quast_N50": 0,
        "quast_N_per_100": 0.0,
        "quast_mismatch_per_100": 0.0,
        "quast_indels_per_100": 0.0,
        	# BUSCO report
        "one_line_summary_BUSCO": "failed assembly generation",
        "Complete_percentage_BUSCO": 0.0,
        "Complete_BUSCOs": 0,
        "Single_copy_BUSCOs": 0,
        "Multi_copy_BUSCOs": 0,
        "Fragmented_BUSCOs": 0,
        "Missing_BUSCOs": 0,
        "n_markers_BUSCO": 0,
    })



# 2.
# Helper function for matching orig sample names and sample names in the reports

def match_sample_name(report_names, orig_name):
    for report_name in report_names:
        if report_name.startswith(orig_name):
            remaining_string = report_name[len(orig_name)]
            if not remaining_string or remaining_string[0] in "_":		# if there isn't a remaining string in report name OR 
            									# if the remaining string starts with _ (illumina has _S000 after each sample)
                return orig_name, report_name, "match"
    return orig_name, None, "notmatched"
         
    

# 3.
# process sample names from S and get ORIG.names + concentration + seq run
def process_report_S():

	# check if the file exists (if it doesn't just proceed to the next def/function):
    if not os.path.exists(input_files['samples_S']):
        print(f"{input_files['samples_S']} file not found.")
        return

    global orig_names
    samples_S_df = pd.read_csv(input_files['samples_S'], sep="\t")
    orig_names = []
    for _, row in samples_S_df.iterrows():
        
        sample = ''.join(char for char in row.iloc[0].strip() if char.isprintable())	# .join with "" (=empty/no char) each char in the sample name that .isprintable()
        
        #sample = row.iloc[0].strip()			
        # build the sample list (orig_names)
        if sample not in orig_names:                
            orig_names.append(sample)
            sample_merged_QC_dict[sample].update({
                "sample_notes": row.iloc[1].strip()
            })            
            
    
    
# 4.
# process kmerfinder data
def process_kmerfinder():

	# check if the file exists (if it doesn't just proceed to the next def/function):
    if not os.path.exists(input_files['kmerfinder_report']):
        print(f"{input_files['kmerfinder_report']} file not found.")
        return
        
    kmerfinder_df = pd.read_csv(input_files['kmerfinder_report'], sep="\t")      
    report_names = set(kmerfinder_df.iloc[:, 0].str.strip())			# create a set of report names from the kmerfinder report (iloc.[:, 0] --> : = "select all rows")  
    

    for orig_name in orig_names:
        matched_name, report_name, matched_status = match_sample_name(report_names, orig_name)  

        if matched_status == "match":
            row = kmerfinder_df[kmerfinder_df.iloc[:, 0].str.strip() == report_name].iloc[0]
            sample_merged_QC_dict[matched_name].update({
                "seq_samplename": row.iloc[0].strip(),
                "species": row.iloc[1].strip(),
                "ncbi_acc": row.iloc[2].strip(),
                "sum_coverage_ncbi": float(row.iloc[3])
            })
        else:
            sample_merged_QC_dict[matched_name].update({
                "seq_samplename": "match not found"
            })




# 6.
# process MLST report
def process_mlst():

	# check if the file exists (if it doesn't just proceed to the next def/function):
    if not os.path.exists(input_files['MLST_report']):
        print(f"{input_files['MLST_report']} file not found.")
        return

    mlst_df = pd.read_csv(input_files['MLST_report'], sep="\t")
    report_names = set(mlst_df.iloc[:, 0].str.strip())
    
    for orig_name in orig_names:
        matched_name, report_name, matched_status = match_sample_name(report_names, orig_name)

        if matched_status == "match":
            row = mlst_df[mlst_df.iloc[:, 0].str.strip() == report_name].iloc[0]
            sample_merged_QC_dict[matched_name].update({
                "ST_mlst": row.iloc[2],
                "MLST_scheme": row.iloc[1].strip()
            })
        else:
            sample_merged_QC_dict[matched_name].update({
                "ST_mlst": "match not found"
            })        



# 7.
# process coverages report
def process_coverages():

	# check if the file exists (if it doesn't just proceed to the next def/function):
    if not os.path.exists(input_files['coverages_report']):
        print(f"{input_files['coverages_report']} file not found.")
        return

    coverages_df = pd.read_csv(input_files['coverages_report'], sep="\t")
    report_names = set(coverages_df.iloc[:, 0].str.strip())
    
    for orig_name in orig_names:
        matched_name, report_name, matched_status = match_sample_name(report_names, orig_name)
        
        if matched_status == "match":
            row = coverages_df[coverages_df.iloc[:, 0].str.strip() == report_name].iloc[0]
            sample_merged_QC_dict[matched_name].update({
                "avg_cov_assembly": row.iloc[1],
                "organism_folder": row.iloc[2].strip()
            })
        else:
            sample_merged_QC_dict[matched_name].update({
                "avg_cov_assembly": "match not found"
            }) 



# 8. 
# process QUAST reports
def process_quast():
    
    quast_df = pd.DataFrame() 	# initialize quast_df (you will then recursively store all values into it - for each quast report)
    
    for file in glob.glob(f"{input_files['QUAST_dir']}/*report_transposed*.tsv"):
    		# glob.glob already filters for existing files --> IF NOT FOUND, the loop doesn't run'     
    
        quast_species_report_df = pd.read_csv(file, sep="\t")
        quast_df = pd.concat([quast_df, quast_species_report_df], ignore_index=True)

    report_names = set(quast_df.iloc[:, 0].str.strip())
        
    for orig_name in orig_names:
        matched_name, report_name, matched_status = match_sample_name(report_names, orig_name)
            
        if matched_status == "match":
            row = quast_df[quast_df.iloc[:, 0].str.strip() == report_name].iloc[0]  
            sample_merged_QC_dict[matched_name].update({
                "quast_num_contigs": row.iloc[13],
                "quast_total_length": row.iloc[15],
                "quast_ref_length": row.iloc[16],
                "quast_GC": row.iloc[17],
                "quast_ref_GC": row.iloc[18],
                "quast_N50": row.iloc[19],
                "quast_N_per_100": row.iloc[40],
                "quast_mismatch_per_100": row.iloc[41],
                "quast_indels_per_100": row.iloc[42]
            })
        else:
            sample_merged_QC_dict[matched_name].update({
                "quast_num_contigs": "match not found"
            })


# 9.
# process BUSCO reports
def process_busco():

    busco_df = pd.DataFrame() 	# initialize busco_df (you will then recursively store all values into it - for each busco report)
    
    for file in glob.glob(f"{input_files['BUSCO_dir']}/*filtered.xlsx"):
    		# glob.glob already filters for existing files --> IF NOT FOUND, the loop doesn't run'     
    
        busco_species_report_df = pd.read_excel(file)
        busco_df = pd.concat([busco_df, busco_species_report_df], ignore_index=True)

    report_names = set(busco_df.iloc[:, 0].str.strip())

    for orig_name in orig_names:
        matched_name, report_name, matched_status = match_sample_name(report_names, orig_name) 
        
        if matched_status == "match":
            row = busco_df[busco_df.iloc[:, 0].str.strip() == report_name].iloc[0]
            sample_merged_QC_dict[matched_name].update({
                "one_line_summary_BUSCO": row.iloc[1].strip(),
                "Complete_percentage_BUSCO": row.iloc[2],
                "Complete_BUSCOs": row.iloc[3],
                "Single_copy_BUSCOs": row.iloc[4],
                "Multi_copy_BUSCOs": row.iloc[5],
                "Fragmented_BUSCOs": row.iloc[6],
                "Missing_BUSCOs": row.iloc[7],
                "n_markers_BUSCO": row.iloc[8],
            })
        else:
            sample_merged_QC_dict[matched_name].update({
                "one_line_summary_BUSCO": "match not found"
            })        


# 10.
# process serotyping report
def process_serotype():
    if Path(input_files['serotype_report']).exists():									# check if serotype report exists
        serotype_df = pd.read_excel(input_files['serotype_report'], sheet_name="Sheet1", header=None).transpose()
        
        report_names = set(serotype_df.iloc[:, 0].str.strip())
    
        for orig_name in orig_names:
            matched_name, report_name, matched_status = match_sample_name(report_names, orig_name)

            if matched_status == "match":
                row = serotype_df[serotype_df.iloc[:, 0].str.strip() == report_name].iloc[0]
                
                if row.iloc[8] == "PASS":								
                    sample_merged_QC_dict[matched_name].update({
                        "serovar": row.iloc[1].strip(),
                        "serogroup": row.iloc[2].strip()
                    })
                else:
                     sample_merged_QC_dict[matched_name].update({
                        "serovar": row.iloc[9].strip()
                    })
            else:
                sample_merged_QC_dict[matched_name].update({
                    "serovar": "not analyzed"
                })




# 11. 
# generate merged report
def write_report():
    header_report = [
    "Sample Notes",
    "Sequencing sample name",  
    "Species (kmerfinder)", 
    "NCBI species assembly (kmerfinder)",
    "Species coverage (kmerfinder) [%]", 
    "ST (mlst)",
    "Matching PubMLST scheme",
    "Serovar",
    "Serogroup",    
    "Assembly average coverage",
    "Organism folder",
    "Q:Number of contigs",
    "Q:Reference length",
    "Q:Assembly total length",
    "Q:Reference GC [%]",
    "Q:Assembly GC [%]",
    "Q:Assembly N50",
    "Q:N's per 100 kbp",
    "Q:Mismatches per 100 kbp",
    "Q:Indels per 100 kbp",  
    "Summary BUSCO",
    "Complete BUSCOs [%]",
    "Complete BUSCOs (C)",
    "Complete and single-copy BUSCOs (S)",
    "Complete and duplicated BUSCOs (D)",
    "Fragmented BUSCOs (F)",
    "Missing BUSCOs (M)",
    "Total BUSCO groups searched"
    ]
    pd.DataFrame.from_dict(sample_merged_QC_dict, orient="index").to_excel(merged_output_QC, index_label="Sample name", header=header_report, sheet_name=run_name, float_format="%.2f")


# APPLY VISUAL CHANGES TO xlsx    
    xlsx_file = load_workbook(merged_output_QC)    	
    xlsx_worksheet = xlsx_file[run_name]    		
    
    alignment_style = Alignment(horizontal="left")			# left align each cell
    header_font_style = Font(name="Calibri", size=11, bold=True)	# set the same font and BOLD for header cell + 1st column (sample names)
    body_font_style = Font(name="Calibri", size=11)			# set the same font for body cells
    
    for row_id, row in enumerate(xlsx_worksheet.iter_rows(), start=1): 	
        for cell in row:
            if row_id == 1:						# header
                cell.alignment = alignment_style			
                cell.font = header_font_style
            elif cell.column == 1:					# 1st column    	
                cell.alignment = alignment_style
                cell.font = header_font_style
            else:
                cell.alignment = alignment_style
                cell.font = body_font_style

    xlsx_file.save(merged_output_QC)			# save realigned file
    
    print(f"Merged QC report saved to: {merged_output_QC}")



# 10. 
# main function
def main():     # if in the .sh script: def main(output_dir, run_name):
    init_paths()
    init_data_structure()
    process_report_S()
    process_kmerfinder()
    process_mlst()
    process_coverages()
    process_quast()
    process_busco()
    process_serotype()
    write_report()

main()


END

log_start_time "MERGED QC REPORT in $output_dir"

}






	# WORKFLOW
	
	# flags to control the workflow	
run_update_env=false
run_quality_check_raw=false
run_trim_raw_reads=false
run_quality_check_trimmed=false
run_determine_organism=false
run_denovo_assembly=false
run_assembly_quality_check_Q=false
run_assembly_quality_check_B=false
run_mapping_reads_to_assembly_bam_sam_coverage=false
run_mlst=false
run_sistr=false
run_generate_report=false


	# check command line arguments (loop to iterate through all CLI arguments)    
	
for arg in "$@"; do
	case "$arg" in
		update)
			run_update_env=true
			;; 
		qc_raw)
			run_quality_check_raw=true
			;;
		trim)
			run_trim_raw_reads=true
			;;
		qc_trimmed)
			run_quality_check_trimmed=true
			;;
		organism)
			run_determine_organism=true
			;;
		assembly)
			run_denovo_assembly=true
			;;
		quast)
			run_assembly_quality_check_Q=true
			;;
		busco)
			run_assembly_quality_check_B=true
			;;	
		sam_bam_coverage)
			run_mapping_reads_to_assembly_bam_sam_coverage=true
			;;
		mlst)
			run_mlst=true
			;;
		sistr)
			run_sistr=true
			;;	
		report)
			run_generate_report=true
			;;
		all)
			run_quality_check_raw=true
			run_trim_raw_reads=true
			run_quality_check_trimmed=true
			run_determine_organism=true
			run_denovo_assembly=true
			run_assembly_quality_check_Q=true
			run_assembly_quality_check_B=true
			run_mapping_reads_to_assembly_bam_sam_coverage=true
			run_mlst=true
			run_sistr=true			
			run_generate_report=true
			;;
	esac
done


	# run the functions/steps based on the flags
$run_update_env && update_conda_envs
$run_quality_check_raw && quality_check_raw
$run_trim_raw_reads && trim_raw_reads
$run_quality_check_trimmed && quality_check_trimmed
$run_determine_organism && determine_organism
$run_denovo_assembly && denovo_assembly
$run_assembly_quality_check_Q && assembly_quality_check_QUAST
$run_assembly_quality_check_B && assembly_quality_check_BUSCO
$run_mapping_reads_to_assembly_bam_sam_coverage && mapping_reads_to_assembly_bam_sam_coverage
$run_mlst && mlst_all
$run_sistr && sistr_serotyping
$run_generate_report && merge_all_reports











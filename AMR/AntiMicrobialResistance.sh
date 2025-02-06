#!bin/bash

source ~/anaconda3/etc/profile.d/conda.sh
# ensures that conda is initialized in the script's shell session + sets up the necessary environment variables + functions for conda to operate		

conda_envs='~/mambaforge/envs/'

	
	
	# to execute the script:
# 	a) the whole script = all steps:
# bash AMR.sh all

# 	b) only specific steps
#		1: amrfinder
#		2: resfinder
#		3: rgi	
# bash AMR.sh rgi / bash AMR.sh amrfinder rgi


	
	# CHANGE
	
input_dir="path/to/assembly/dir"
	# input_dir should have subsequent folders based on the organism (input_folder/Salmonella_enterica, input_folder/Escherichia_coli, ...)

output_dir_user="path/to/desired/output/dir"

run_name='user_defined_run_name'

	# DON'T CHANGE
	
output_dir="$output_dir_user/AMR" && mkdir -p $output_dir


	# function to get the current date and time = the date and the start time of each process (and writes it into a log file)
log_start_time() {
	local log_text="$1"		
	echo "[$(date)] PROCESS: $log_text" >> "$output_dir/${run_name}_analysis_pipeline_report.log"
	#echo "[$(date)] PROCESS: $1" >> $output_dir_user/analysis_pipeline_report.log
}


amrfinderplus_fun() {

	echo -e "\nStarting AMRFinderPlus.\n"

	conda activate $conda_envs/amrfinderplus

	log_start_time "	STARTING AMRfinderplus version:"
	amrfinder --database_version >> $output_dir/${run_name}_analysis_pipeline_report.log
	# or -V: Print out complete version information of both the database and software.


	for organism in $input_dir/*/; do
		organism_name=$(basename $organism)

		echo -e "\nStarting AMRFinder."	
		log_start_time "Starting AMRfinder for organism: $organism_name"
		
		mkdir -p $output_dir/01_amrfinderplus/$organism_name
				
		amrfinder -l	
		echo -e "\nPlease choose the correct scheme from the list of organisms for: $organism_name\n"
		
		read organism_name_amrfinder 
		
		log_start_time "User entered: $organism_name_amrfinder as the name of the organism."
	
		for assembly in $organism/*.fasta; do
			assembly_name=$(basename $assembly .fasta) 
			
			amrfinder -n $assembly -O $organism_name_amrfinder --name $assembly_name --threads 20 -o $output_dir/01_amrfinderplus/$organism_name/"${assembly_name}_amrfinder.txt"

		done
		
	echo -e "\nProcess finished for $organism_name"		
	log_start_time "Process finished for: $organism_name"

	done
	
	log_start_time "Finished AMRFinderplus process --> results in $output_dir/01_amrfinderplus"
	
}



resfinder_fun() {

	echo -e "\nStarting ResFinder.\n"

	conda activate $conda_envs/resfinder
	log_start_time "	STARTING ResFinder version:"
	run_resfinder.py --version >> $output_dir/${run_name}_analysis_pipeline_report.log
	
		
	for organism in $input_dir/*/; do
		organism_name=$(basename $organism)
		
		mkdir -p $output_dir/02_resfinder/$organism_name
	
		echo -e "\nPlease choose the correct name from the list (to input into resfinder command) for your organism: $organism_name\n
		Available species: Campylobacter, Campylobacter jejuni, Campylobacter coli, Enterococcus faecalis, Enterococcus faecium, Escherichia coli, Helicobacter pylori, Klebsiella, Mycobacterium tuberculosis, Neisseria gonorrhoeae, Plasmodium falciparum, Salmonella, Salmonella enterica, Staphylococcus aureus; \"Other\" can be used for metagenomic samples"		
		read organism_name_resfinder 
		log_start_time "User entered $organism_name_resfinder as the name of the organism."
	
		for assembly in $organism/*.fasta; do
			assembly_name=$(basename $assembly .fasta)  
			
			run_resfinder.py -o $output_dir/02_resfinder/$organism_name/$assembly_name -s "${organism_name_resfinder}" -l 0.6 -t 0.8 --point --acquired -ifa $assembly
		
		done
		
		echo -e "\nProcess finished for $organism_name"
		log_start_time "Process finished for $organism_name"
	
	done
	
	log_start_time "Finished ResFinder process --> results in $output_dir/02_resfinder"

}


rgi_fun() {
	
	echo -e "\nStarting rgi.\n"

	conda activate $conda_envs/rgi
	log_start_time "	STARTING RGI version:"
	rgi main -v >> $output_dir/${run_name}_analysis_pipeline_report.log
	
		
	for organism in $input_dir/*/; do
		organism_name=$(basename $organism)
		
		mkdir -p $output_dir/03_rgi/$organism_name
		
		for assembly in $organism/*.fasta; do 
			assembly_name=$(basename $assembly .fasta)
			
			rgi main --input_sequence $assembly --output_file $output_dir/03_rgi/$organism_name/${assembly_name}_rgi --clean 
	
		done
		
		echo -e "\nProcess finished for $organism_name"
		log_start_time "Process finished for $organism_name"			

	done

	log_start_time "Finished Rgi process --> results in $output_dir/03_rgi"

}




# WORKFLOW
	
	# flags to control the workflow (initialize = false --> if user inouts into command line the arg, then it changes into "true" in the next step)	
run_amrfinderplus=false
run_resfinder=false
run_rgi=false
  
	
for arg in "$@"; do
	case "$arg" in 
		amrfinder)
			run_amrfinderplus=true
			;;
		resfinder)
			run_resfinder=true
			;;
		rgi)
			run_rgi=true
			;;
		all)
			run_amrfinderplus=true
			run_resfinder=true
			run_rgi=true
			;;
	esac
done

	# run the steps based on the flags
$run_amrfinderplus && amrfinderplus_fun
$run_resfinder && resfinder_fun
$run_rgi && rgi_fun







# **AMR Detection Pipeline**

## **Pipeline Overview**

***AntiMicrobialResistance.sh***

This pipeline automates the identification of antimicrobial resistance (AMR) genes for assembled genomes.

### **Features**

- **AMR detection**  (`AMRFinderPlus`, `ResFinder`, `RGI`)

## **Installation & Dependencies**

To run this pipeline, you need:
- **Linux OS** (`Ubuntu` or `HPC cluster`)
- **Conda** (for managing dependencies)
- **Required tools installed in Conda environments:**
  - `AMRFinderPlus`
  - `ResFinder`
  - `RGI`

### **Conda Environment Setup**

To install required environments (the env names listed below are directly referenced in the script - if the env names differ, the user must alter the script accordingly):

```bash
    # Activate conda base
source ~/anaconda3/etc/profile.d/conda.sh

    # Create environments (if not installed)
conda create -n amrfinderplus amrfinderplus -c bioconda -c conda-forge -y
conda create -n resfinder resfinder -c bioconda -c conda-forge -y
conda create -n rgi rgi -c bioconda -c conda-forge -y
```


## **Pipeline Execution**

### Pre-Execution Configuration

The input files (assemblies) are expected to be in the organism specific directories:
```
input_dir/
│── Name_of_organism1/
│   ├── assembly1.fasta
│   ├── assembly2.fasta  
├── Name_of_organism2/
│   ├── assembly3.fasta
│   ├── assembly4.fasta
...
```

#### 1. Update paths in the script
Before executing the pipeline, the user **must edit** the **.sh script** to update variables with paths specific to your environment:

- **Conda Environments Directory**
```bash
conda_envs='~/mambaforge/envs'
```
- **Input Files Directory (.fastq)**
```bash
raw_reads_dir='path/to/fastq/dir'
```

- **input_dir='path/to/assembly/dir'**
<br > *--> inside the input dir, the .fasta files are in subfolders:*

- **output_dir='path/to/desired/output/dir'**
<br > *--> Set the location where all pipeline results will be saved.*
- **run_name='user_defined_run_name'**

After saving the **.sh script**, the pipeline can be executed.

### Execution
#### 1. Run Entire Pipeline:
```bash
bash AntiMicrobialResistance.sh all
```
#### 2. Run Specific Steps
```bash
bash AntiMicrobialResistance.sh <step1 step2 ...>
```

#### **Available Steps**:

| step | Description |
|------|-------------|
| amrfinder | AMR genes detection |
| resfinder | AMR genes detection |
| rgi | AMR genes detection|

### File Structure 
The pipeline generates the following directory structure:

```
output_dir/
│── AMR/
│   ├── 01_amrfinderplus/
│   ├── 02_resfinder/
│   ├── 03_rgi/
```

### Logging
All pipeline steps log output to:
```
output_dir/assembly_pipeline_report.log 
```

## Credits & Contact
- <big>Tea Janko, <small>MSc Bioinformatics
- oxytotia@gmail.com
- National Laboratory of Health, Environment and Food\
Department of Public Health Microbiology\
Slovenia  
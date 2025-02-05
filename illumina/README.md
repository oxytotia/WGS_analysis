# **Illumina WGS Processing Pipeline**

## **Pipeline Overview**
***illumina_raw_to_assembly_typing.sh***

A comprehensive Whole Genome Sequencing (WGS) pipeline for processing Illumina sequencing data - designed for automated high-throughput bacterial genome analysis, from raw reads (FASTQ) to quality control, trimming, assembly, typing, and final QC report generation.


### **Features**

- **Quality control** before and after trimming (`FastQC`, `MultiQC`)
- **Read trimming** (`Trimmomatic`)
- **Organism detection** (`KmerFinder`)
- **Genome assembly** (`SPAdes`)
- **Assembly quality assessment** (`QUAST`, `BUSCO`)
- **Read mapping & coverage analysis** (`Minimap2`, `Samtools`)
- **Sequence typing** (`MLST`)
- **Salmonella serotyping** (`SISTR`)
- **Final QC report merging** (`Python script`)
- **Updating tools** (`Bash script`)

## **Installation & Dependencies**

To run this pipeline, you need:
- **Linux OS** (`Ubuntu` or `HPC cluster`)
- **Conda** (for managing dependencies)
- **Python 3** (for report merging)

- **Required tools installed in Conda environments:**
  - `FastQC`, `MultiQC`
  - `Trimmomatic`
  - `KmerFinder`
  - `SPAdes`
  - `QUAST`
  - `BUSCO`
  - `Minimap2`, `Samtools`
  - `MLST`
  - `SISTR`

### **Conda Environment Setup**

To install required environments:

```bash
    # Activate conda base
source ~/anaconda3/etc/profile.d/conda.sh

    # Create environments (if not installed)
conda create -n qc fastqc multiqc -c bioconda -c conda-forge -y
conda create -n trimmomatic trimmomatic -c bioconda -c conda-forge -y
conda create -n kmerfinder kmerfinder -c bioconda -c conda-forge -y
conda create -n spades spades -c bioconda -c conda-forge -y
conda create -n quast quast -c bioconda -c conda-forge -y
conda create -n busco busco -c bioconda -c conda-forge -y
conda create -n minimap2 minimap2 samtools -c bioconda -c conda-forge -y
conda create -n mlst mlst -c bioconda -c conda-forge -y
conda create -n sistr sistr -c bioconda -c conda-forge -y
```

## **Pipeline Execution**

### Pre-Execution Configuration

This pipeline is currently designed to process **Illumina NextSeq sequencing output**.\
The input files (**raw reads**) are expected to be **organized into organism-specific folders** within the input directory and to have **a specific appendix** to the sample name:
```
input_dir/
│── Name_of_organism1/
│   ├── samplename_1_R1_001.fastq.gz
│   ├── samplename_1_R2_001.fastq.gz  
├── Name_of_organism2/
│   ├── samplename_2_R1_001.fastq.gz
│   ├── samplename_2_R2_001.fastq.gz
```

> [!NOTE]  
> If your files follow a different naming scheme, you can adjust the pipeline to match the unique format. However, be **extremely careful** to modify **all occurences** of the file naming pattern throughout the script to avoid errors.


#### 1. Prepare the sample sheet
The pipeline requires a **sample sheet**, which should be a **tab-separated file (.tsv)** with the following structure:
```
Sample name         Sample notes
19-402-2024-Yers    Isolation - conc.: 34.8
19-231-2024-STm     Isolation - conc.: 50.6
...
```
**Sample name** must be a prefix of the filenames of raw reads (.fastq), or a complete match.\
**Sample notes** field is optional but can contain additional info about the sample.

#### 2. Update paths in the script
Before executing the pipeline, the user **must edit** the **.sh script** to update variables with paths specific to your environment:

- **Conda Environments Directory**
```bash
conda_envs='~/mambaforge/envs'
```

- **Input Files Directory (.fastq)**
```bash
raw_reads_dir='path/to/fastq/dir'
```

- **Pipeline Results Directory**
```bash
output_dir='path/to/desired/output/dir'
```

- **Run Name**
```bash
run_name='user_defined_run_name'
```

- **Trimming Quality Score**
```bash
trimmomatic_q_score="4:20"
```
> (default) 4:20 = average quality score of 4 bases above q20 [ 99% accuracy ]
<br > - 4:30 = average q score of 4 bases above q30 [ 99,9% accuracy ]
<br > - 1:20 = each base q20 [ 99% accuracy ]
<br > - 1:30 = each base q30 [ 99,9% accuracy ]
<br > ...

- **Reference Genomes Directory**
```bash
ref_genomes_dir="location/of/reference/genomes/for/QUAST"
```
> [!NOTE]  
> During execution, the script will **prompt the user in CLI** to enter reference genome path.  
> Since the path to the directory with ref genomes is specified here, the user only needs to enter the name of the ref genome in the CLI and the script will automatically append the ful path and continue.

- **KmerFinder Database Directory**
```bash
kmerfinder_db='location/of/kmerfinder/databases/bacteria/*bacteria.ATG*'
kmerfinder_tax='location/of/kmerfinder/databases/bacteria/*bacteria.tax*'
```

After saving the **.sh script**, the pipeline can be executed.

### Execution
#### 1. Run Entire Pipeline:
```bash
bash illumina_processing_raw_pipeline_linux.sh all
```
#### 2. Run Specific Steps
```bash
bash illumina_processing_raw_pipeline_linux.sh <step1 step2 step3 ...>
```

#### **Available Steps**:

| step | Description |
|------|-------------|
| qc_raw | Quality control of raw reads (FastQC, MultiQC) |
| trim | Read trimming (Trimmomatic) |
| qc_trimmed | QC on trimmed reads (FastQC, MultiQC)|
| organism | Organism identification (KmerFinder) |
| assembly | Genome assembly (SPAdes) |
| quast | Assembly quality check (QUAST) |
| busco | Assembly completeness check (BUSCO) |
| sam_bam_coverage | Read mapping & coverage analysis (Minimap2, Samtools) |
| mlst | Sequence typing (MLST) |
| sistr | Salmonella serotyping (SISTR) |
| report | Generate final QC report |
| update | Update specific tools within their respective conda environments |


### File Structure 
The pipeline generates the following directory structure:

```
/output_dir/
│── 01_QC_reads/
│   ├── 01_raw_QC/
│       ├── 01_fastQC/
│       ├── 02_multiQC/
│   ├── 02_trimmomatic_QC/
│       ├── 01_fastQC/
│       ├── 02_multiQC/
│── 02_trimmomatic/
│── 03_kmerfinder/
|   ├── 01_results_all/
|   ├── kmerfinder_species_summary.txt
|   ├── kmerfinder_species_summary_merged.txt
│── 04_spades/
|   ├── 01_scaffolds_all/
│── 05_QC_assembly/
│   ├── 01_QUAST/
│       ├── QUAST_report.html
│       ├── QUAST_report_transposed.tsv
│   ├── 02_BUSCO/
│       ├── merged_BUSCO_report.xlsx
│       ├── merged_BUSCO_report_filtered.xlsx
│── 06_minimap2/
│   ├── avg_cov_assembly_all_sorted.txt
│── 07_MLST/
│   ├── MLST_report.tsv
│── 08_serotyping/
│   ├──01_SISTR
│       ├── reports/
│       ├── serotyping_report_SISTR_Salmonella.xlsx
│── merged_QC_report.xlsx
```

### Logging
All pipeline steps log output to:
```
output_dir/assembly_pipeline_report.log 
```
### Final Merged QC Report

The report provides a sumary of sequencing quality.

Some samples may have **"No Match Found"** in certain fields. This indicates that no corresponding entry was found for that sample in a specific report. 
#### Possible Reasons: 
1. **The sample was not processed** in one of the pipeline steps (e.g., missing sequencing data, failed species identification, failed assembly generation). 
2. **Naming discrepancies** between the sample sheet and the reports (e.g., extra characters, different formatting). 
3. **The corresponding report was missing** when generating the merged file. 
4. **The sample did not pass a quality check**, and no results were recorded. 
#### What to Do: 
- **Double-check sample names** in the sample sheet and reports to ensure consistency. 
- **Verify that all required reports exist** in the expected directories.
- **Resequence** the failed sample. 


| Sample name | Notes | Sequencing sample name |Species (kmerfinder) | NCBI species assembly (kmerfinder) | Species coverage (kmerfinder) [%] | ST (mlst) | Matching PubMLST scheme | Serovar | Serogroup | Assembly average coverage | Organism folder | Q:Number of contigs | Q:Reference length | Q:Assembly total length | Q:Reference GC [%] | Q:Assembly GC [%] |	Q:Assembly N50 | Q:N's per 100 kbp | Q:Mismatches per 100 kbp |	Q:Indels per 100 kbp | Summary BUSCO | Complete BUSCOs [%] | Complete BUSCOs (C) | Complete and single-copy BUSCOs (S) | Complete and duplicated BUSCOs (D) | Fragmented BUSCOs (F) | Missing BUSCOs (M) | Total BUSCO groups searched |
|-|-|-|-|-|-|-|-|-|-|-|-|-|-|-|-|-|-|-|-|-|-|-|-|-|-|-|-|-|
| 19-401-2019-Shsonn | Isolation - conc.: 56.8 | 19-401-2019-Shsonn_S50 | Shigella sonnei | GCF_002224625.1 | 92.48 | 152 | ecoli_achtman_4 | not analyzed | N/A | 97.26 | Shigella_sonnei | 383 | 4762774 | 4523572 | 50.73 | 50.85 | 25752 | 11.05 | 1112.79 |	26.3 | C:99.2%[S:99.1%,D:0.1%],F:0.6%,M:0.2%,n:874 | 99.2 | 867 | 866 | 1 | 5 | 2 | 874 |
| 19-402-2024-Yers | Isolation - conc.: 34.8 | match not found | N/A | N/A | 0 | match not found | N/A | match not found | N/A | match not found | N/A | match not found | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | match not found | 0 | 0 | 0 | 0 | 0 | 0 | 0 |
| 19-231-2024-STm | Isolation - conc.: 50.6 | 19-231-2024-STm_S30 | Salmonella enterica | GCF_002761055.1 | 98.65 | 323 | senterica_achtman_2 | Typhimurium | B | 127.65 | Salmonella_enterica | 44 | 4951383 | 4856091 | 52.24 | 52.19 | 293037 | 4.12 | 20.19 |	2.5 | C:99.3%[S:99.3%,D:0.0%],F:0.1%,M:0.6%,n:874 | 99.3 | 868 | 868 | 0 | 1 | 5 | 874 |
| 19-411-2024-Sdiarizonae | Isolation - conc.: 95.6 | 19-411-2024-Sdiarizonae_S39 | Salmonella enterica | GCF_002794415.1, GCF_019222725.1, GCF_003324755.1 | 92.65 | 152 | senterica_achtman_2 | INFO: Number of cgMLST330 loci found (n=329) WARNING: Only matched 274 cgMLST330 loci. Min threshold for confident serovar prediction from cgMLST is 297.0 | N/A | 100.13 | Salmonella_enterica | 41 | 4951383 | 4820138 | 52.24 | 51.56 | 407194 | 8.3 | 3676.57 | 85.24 | C:99.4%[S:99.4%,D:0.0%],F:0.1%,M:0.5%,n:874 | 99.4 | 869 | 869 | 0 | 1 | 4 | 874 |



## Credits & Contact
- <big>Tea Janko, <small>MSc Bioinformatics
- oxytotia@gmail.com
- National Laboratory of Health, Environment and Food\
Department of Public Health Microbiology\
Slovenia  
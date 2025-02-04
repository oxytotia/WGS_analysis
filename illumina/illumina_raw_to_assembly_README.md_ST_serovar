This repository contains a comprehensive Whole Genome Sequencing (WGS) pipeline for processing Illumina sequencing data, from raw reads (FASTQ) to quality control, trimming, assembly, typing, and final QC report generation.

The pipeline is designed for automated high-throughput bacterial genome analysis.

  1. Quality Control & Preprocessing
• FastQC & MultiQC: Assess raw read quality
• Trimmomatic: Trim low-quality bases and adapters (sliding window can be adjusted - default: 4:20)
• FastQC & MultiQC (Trimmed Reads): Assess quality after trimming

  2. Organism Identification
• KmerFinder: Identify species from raw reads

  3. Genome Assembly
• SPAdes: De novo genome assembly

  4. Assembly Quality Assessment
• QUAST: Evaluate assembly quality
• BUSCO: Assess genome completeness

  5. Read Mapping & Coverage
• Minimap2 & Samtools: Map reads to assemblies and calculate coverage

  6. Typing & Serotyping
• MLST: Multi-locus sequence typing (ST determination)
• SISTR: Serovar identification for Salmonella

  7. Merging Reports & Summary Generation
• Python script: Merges QC and typing results into a structured Excel (.xlsx) summary

Installation & Dependencies



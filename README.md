# Shewanella Genome Analysis Pipeline 🧬

[![Bash Script](https://img.shields.io/badge/bash-script-blue?logo=gnu-bash)](https://www.gnu.org/software/bash/)
[![VirSorter2](https://img.shields.io/badge/VirSorter2-ready-success)](https://github.com/jiarong/VirSorter2)
[![Ruby](https://img.shields.io/badge/Ruby-%E2%89%A52.4.1-red?logo=ruby)](https://www.ruby-lang.org/)
[![Python](https://img.shields.io/badge/Python-%E2%89%A53.7-blue?logo=python)](https://www.python.org/)
[![AutoFlow](https://img.shields.io/badge/AutoFlow-integrated-orange)](https://github.com/linlabcode/AutoFlow)
[![BLAST+](https://img.shields.io/badge/BLAST+-≥2.15.0-lightgrey?logo=ncbi)](https://blast.ncbi.nlm.nih.gov/)
[![PHASTEST](https://img.shields.io/badge/PHASTEST-prophage--detection-success)](https://phastest.ca/)
[![TarSynFlow](https://img.shields.io/badge/TarSynFlow-functional--annotation-brightgreen)](https://github.com/RouxLab/TarSynFlow)
[![EMBOSS](https://img.shields.io/badge/EMBOSS-installed-blue)](http://emboss.open-bio.org/)
[![BBMap](https://img.shields.io/badge/BBMap-installed-lightgrey)](https://sourceforge.net/projects/bbmap/)
[![FastQC](https://img.shields.io/badge/FastQC-installed-green)](https://www.bioinformatics.babraham.ac.uk/projects/fastqc/)
[![DFAST](https://img.shields.io/badge/DFAST-functional--annotation-brightgreen)](https://dfast.nig.ac.jp/)  
[![Sibelia](https://img.shields.io/badge/Sibelia-genome--synteny-lightblue)](https://github.com/medvedevgroup/Sibelia)  
[![Circos](https://img.shields.io/badge/Circos-visualization-purple)](http://circos.ca/)  
[![pyANI](https://img.shields.io/badge/pyANI-ANI--calculation-yellow)](https://github.com/widdowquinn/pyani)


This repository provides bash scripts for automated genome analysis of Shewanella and related bacteria. It integrates genome download, gene identification, transposon analysis, genomic islands, prophage detection, PCA analysis, and report generation.

# 🚀 Features

1- Automatic download of complete genomes and plasmids from NCBI.

2- Gene identification and functional annotation using TarSynFlow.

3- Transposon analysis and matrix generation.

4- Detection and analysis of genomic islands.

5- Prophage detection using PHASTEST.

6- PCA analysis for genomic features.

7- Fully automated HTML report generation.

# 📂 Repository Structure
project/

├── genome_daemon.sh      # Main bash script to run analyses

├── data/                 # Input data (RefSeq.tsv, genome files, etc.)

├── script/               # Scripts for genome parsing and analysis

├── sh_file/              # sbatch scripts for HPC execution

├── templates/            # AutoFlow templates and report templates

└── results/              # Folder where final results and reports are saved

# 🛠 Requirements

Bash (Linux/macOS)

Ruby ≥ 2.4.1

Python 3.7+ (required for PCA and AutoFlow integration)

AutoFlow

BLAST+ (≥ 2.15.0)

Genome files in FASTA format

Optional: Circos for visualization

# ⚡ Usage

Run the main script with a specific mode:

`bash genome_daemon.sh <mode>`

## Mode	Description
### 🧩 down	
This module automatically retrieves complete Shewanella genomes and plasmids from the NCBI RefSeq database.
It requires an input file named RefSeq.tsv located in the data/ directory.

The file can be obtained directly from the NCBI Datasets portal using the following link:
🔗 NCBI Datasets – Shewanella genomes (https://www.ncbi.nlm.nih.gov/datasets/genome/?taxon=22&annotated_only=true&refseq_annotation=true&assembly_level=2:3)

This table contains metadata for each genome to be downloaded, including:
- Organism name
- Taxonomic ID
- Assembly accession (GCF/GCA)
- Strain name
- FTP download path
  
During execution, the script reads the RefSeq.tsv file, downloads each listed genome and plasmid from NCBI, and organizes them into structured directories for downstream analysis.
Download and separate genomes and plasmids from NCBI.

### 🧩 ab1_clean--Clean AB1 sequencing files.
This module processes raw **AB1 Sanger sequencing files** to generate high-quality FASTQ and FASTA sequences suitable for downstream 16S or genomic analyses.

**Required inputs:**
1. `$data_path/ab1/` — folder containing `.ab1` files (one per sample).  
2. `$data_path/ab1_name` — plain text file listing AB1 filenames without extension
   **Workflow steps:**
- Convert AB1 → FASTQ → FASTA using **EMBOSS seqret**  
- Perform quality trimming with **BBDuk** (`qtrim=rl`, `trimq=20`)  
- Generate FASTA files with proper headers  
- Run **FastQC** to create quality control reports in an `analysis/` folder

**Outputs:**
- Cleaned FASTQ: `$data_path/16S_file/<sample>/<sample>clean.fastq`  
- FASTA: `$data_path/16S_file/<sample>/<sample>.fasta`  
- FastQC reports: `$data_path/16S_file/<sample>/analysis/`

**Dependencies / Software required:**

[![EMBOSS](https://img.shields.io/badge/EMBOSS-installed-blue)](http://emboss.open-bio.org/)  
[![BBMap](https://img.shields.io/badge/BBMap-installed-lightgrey)](https://sourceforge.net/projects/bbmap/)  
[![FastQC](https://img.shields.io/badge/FastQC-installed-green)](https://www.bioinformatics.babraham.ac.uk/projects/fastqc/)

### 🧩 char--Reassign genomes using 16S, ANI, and annotation comparison.
This module performs a comprehensive characterization of genomes, including:

1. **16S rRNA BLAST** comparison against the NCBI 16S database  
2. **Genome annotation** using DFAST (CDS, genome length, COG/TIGR counts)  
3. **Synteny analysis** with Sibelia and Circos visualization  
4. **Average Nucleotide Identity (ANI)** calculation using pyANI

**Required Inputs:**
- `$data_path/total_genomes/*.fasta` → genomes to analyze  
- `$data_path/genomes_problem/*.fasta` → any problematic genomes  
- `$data_path/16S_file/` → 16S sequences  
- `$data_path/16S_NCBI/16S_ribosomal_RNA.fasta` → 16S reference database  

**Outputs:**
- 16S BLAST results: `blast_16S/blast_<sample>`  
- Annotated genomes and DFAST parser tables: `genome_annotation/results_dfast_parser/`  
- Synteny analysis results and Circos plots  
- ANI matrices and visualizations: `genome_pyani_anim/`  

**Dependencies:**
[![DFAST](https://img.shields.io/badge/DFAST-functional--annotation-brightgreen)](https://dfast.nig.ac.jp/)  
[![Sibelia](https://img.shields.io/badge/Sibelia-genome--synteny-lightblue)](https://github.com/medvedevgroup/Sibelia)  
[![Circos](https://img.shields.io/badge/Circos-visualization-purple)](http://circos.ca/)  
[![pyANI](https://img.shields.io/badge/pyANI-ANI--calculation-yellow)](https://github.com/widdowquinn/pyani)

### 🧩 protein_db	
Build a protein database for transposon analysis.
### 🧩 tp_case	
Analyze transposons per genome.
### 🧩 tp_matrix	
Generate transposon matrices.
### 🧩 genes	
Identify genes with TarSynFlow.
### 🧩 genes_comps	
Compare gene presence across genomes.
### 🧩 genes_results_protein
Extract protein results.
### 🧩 genes_results_annotation	
Extract annotation results.
### 🧩 seqs	
Extract specific sequences.
### 🧩 GI_up	
Detect genomic islands.
### 🧩 GI_result	
Summarize genomic island results.
### 🧩 GI_clean	
Clean genomic island outputs.
### 🧩 phage_visorter	
Identify prophages with PHASTEST.
### 🧩 results
Integrate all analysis results.
### 🧩 PCA	
Run PCA on genomic features.
### 🧩 PCA_filter	
Filter PCA results.
### 🧩 report	
Generate final HTML report with figures.

# Example Workflow
### 1️⃣ Download genomes
bash genome_daemon.sh down

### 2️⃣ Clean AB1 sequencing files
bash genome_daemon.sh ab1_clean

### 3️⃣ Perform genome characterization and 16S analysis
bash genome_daemon.sh char

### 4️⃣ Identify genes and compare across genomes
bash genome_daemon.sh genes

bash genome_daemon.sh genes_comps

### 5️⃣ Analyze transposons
bash genome_daemon.sh tp_case

bash genome_daemon.sh tp_matrix

### 6️⃣ Detect genomic islands and prophages
bash genome_daemon.sh GI_up

bash genome_daemon.sh GI_result

bash genome_daemon.sh phage_visorter

### 7️⃣ Integrate results and generate report
bash genome_daemon.sh results

bash genome_daemon.sh report

# 🔗 References

NCBI Datasets API

PHASTEST

Roux et al., 2021. TarSynFlow: Genome-wide gene and mobile element analysis.

# Shewanella Genome Analysis Pipeline 🧬

[![Bash Script](https://img.shields.io/badge/bash-script-blue?logo=gnu-bash)](https://www.gnu.org/software/bash/)
[![VirSorter2](https://img.shields.io/badge/VirSorter2-ready-success)](https://github.com/jiarong/VirSorter2)
[![Ruby](https://img.shields.io/badge/Ruby-%E2%89%A52.4.1-red?logo=ruby)](https://www.ruby-lang.org/)
[![Python](https://img.shields.io/badge/Python-%E2%89%A53.7-blue?logo=python)](https://www.python.org/)
[![AutoFlow](https://img.shields.io/badge/AutoFlow-integrated-orange)](https://github.com/linlabcode/AutoFlow)
[![BLAST+](https://img.shields.io/badge/BLAST+-≥2.15.0-lightgrey?logo=ncbi)](https://blast.ncbi.nlm.nih.gov/)
[![PHASTEST](https://img.shields.io/badge/PHASTEST-prophage--detection-success)](https://phastest.ca/)
[![TarSynFlow](https://img.shields.io/badge/TarSynFlow-functional--annotation-brightgreen)](https://github.com/RouxLab/TarSynFlow)
[![Circos](https://img.shields.io/badge/Circos-optional-lightblue)](http://circos.ca/)


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
Download and separate genomes and plasmids from NCBI.
### 🧩 ab1_clean	
Clean AB1 sequencing files.
### 🧩 char	
Reassign genomes using 16S, ANI, and annotation comparison.
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

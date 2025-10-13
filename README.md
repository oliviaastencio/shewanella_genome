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
[![curl](https://img.shields.io/badge/curl-download-orange)](https://curl.se/)
[![make_user_db.rb](https://img.shields.io/badge/make_user_db.rb-local_database-blueviolet)](https://github.com/)
[![ISEScan](https://img.shields.io/badge/ISEScan-transposase--detection-green)](https://github.com/xiezhq/ISEScan)
[![CD-HIT](https://img.shields.io/badge/CD--HIT-clustering-lightblue)](http://weizhongli-lab.org/cd-hit/)
[![merge_tabular.rb](https://img.shields.io/badge/merge_tabular.rb-local--Ruby--script-blueviolet)]
[![R](https://img.shields.io/badge/R-data--visualization-blue)](https://www.r-project.org/)
[![UniProt API](https://img.shields.io/badge/UniProt-annotations-purple)](https://www.uniprot.org/help/api)
[![jq](https://img.shields.io/badge/jq-json--parsing-lightblue)](https://stedolan.github.io/jq/)
[![Ruby](https://img.shields.io/badge/Ruby-scripting-red)](https://www.ruby-lang.org/)
[![Python](https://img.shields.io/badge/Python-script-blue)]
[![transpose_table](https://img.shields.io/badge/transpose_table-tool-lightgrey)]
[![Slurm](https://img.shields.io/badge/Slurm-job_scheduler-orange)]


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

### 🧩 ab1_clean -- Clean AB1 sequencing files.
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

### 🧩 char -- Reassign genomes using 16S, ANI, and annotation comparison.
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

**Dependencies / Software required:**

[![DFAST](https://img.shields.io/badge/DFAST-functional--annotation-brightgreen)](https://dfast.nig.ac.jp/)  
[![Sibelia](https://img.shields.io/badge/Sibelia-genome--synteny-lightblue)](https://github.com/medvedevgroup/Sibelia)  
[![Circos](https://img.shields.io/badge/Circos-visualization-purple)](http://circos.ca/)  
[![pyANI](https://img.shields.io/badge/pyANI-ANI--calculation-yellow)](https://github.com/widdowquinn/pyani)

### 🧩 protein_db	-- Build a protein database for transposon analysis.

This module constructs a custom protein database from UniProt for Shewanella genomes. The database is later used for transposon identification and analysis in subsequent modules (tp_case and tp_matrix).

**Required Inputs:**

- $script_path → folder containing scripts like make_user_db.rb
- $data_path → project data folder
- $transposon_analysis_path → folder to store transposon analysis results
- keyword → UniProt search term (e.g., Shewanella)


**Outputs:**

$data_path/tp_data/total_prots.fasta → downloaded protein sequences

$data_path/tp_data/local_database/ → local protein database ready for transposon analysis

**Dependencies / Software required:**

[![curl](https://img.shields.io/badge/curl-download-orange)](https://curl.se/)  
[![make_user_db.rb](https://img.shields.io/badge/make_user_db.rb-local_database-blueviolet)](https://github.com/) 

### 🧩 tp_case -- Analyze transposons per genome.

#### This module:
- Performs comprehensive transposon analysis for each genome using the protein database from protein_db.
- Identifies transposases, analyzes their genomic environment, filters candidates, and generates summary tables for downstream analysis.

**Required Inputs:**

- $data_path/all_genome_list → list of genome FASTA files to analyze
- $data_path/total_genomes/ → genome FASTA files
- $data_path/tp_data/ → protein database from protein_db

**Outputs:**

$genome_analysis_path/transposon/executions/<genome>/ → per-genome transposon analysis results

tp_case/ subfolders for each transposon, containing:

- proteins.fasta → clustered protein sequences
- tp.fasta → transposon sequences
- blast_summary → BLASTX results
- result_summary → final transposon summary table
  
tp_case/Total_tp → total number of transposons per genome

**Dependencies / Software required:**

[![ISEScan](https://img.shields.io/badge/ISEScan-transposase--detection-green)](https://github.com/xiezhq/ISEScan)  
[![CD-HIT](https://img.shields.io/badge/CD--HIT-clustering-lightblue)](http://weizhongli-lab.org/cd-hit/)  

### 🧩 tp_matrix -- Generate transposon matrices.

#### This module: 
- Compiles the transposon data from tp_case into presence/absence and count matrices for each genome.
- Computes both absolute and relative transposon counts normalized by the total number of CDSs in each genome.
- Prepares per-genome and combined tables for downstream comparative analyses.

### Workflow / Steps:

1. Remove previous results and create the output directory.
2. For each transposon type (interrupt_names, transposase_names):
   - Merge all per-genome transposon identifiers into all_<tab> files.
   - For each genome:
       - Extract genome-specific transposon entries.
       - Compute absolute and relative counts normalized by CDSs.
       - Generate presence/absence table (tab_<tab>) and count table (<tab>_total).

3. Merge all per-genome tables into global matrices using merge_tabular.rb.
4. Concatenate absolute and relative transposon counts across all genomes into Tp_absolute and Tp_relative.

**Required Inputs:**

- $genome_analysis_path/transposon/executions/ → per-genome transposon results from tp_case
- $data_path/all_genome_list → list of genomes
- $genome_analysis_path → folder containing genome annotation results (CDSs)

**Outputs:**

- $genome_analysis_path/transposon/results/tab_interrupt_names → merged presence/absence table for interrupted genes
- $genome_analysis_path/transposon/results/tab_transposase_names → merged presence/absence table for transposases
- $genome_analysis_path/transposon/results/Tp_absolute → absolute transposon counts
- $genome_analysis_path/transposon/results/Tp_relative → relative transposon counts (normalized by CDSs)
- Per-genome subfolders containing detailed tables: tab_<tab>, <tab>_total, Tp_absolute, Tp_relative

**Dependencies / Software required:**

[![merge_tabular.rb](https://img.shields.io/badge/merge_tabular.rb-local--Ruby--script-blueviolet)]()

### 🧩 genes -- Identify genes with TarSynFlow.

#### This module:

- Prepares the protein dataset for TarSynFlow by removing redundancy using CD-HIT.
- Generates a clean proteome ready for gene identification and functional annotation.

#### Workflow / Steps:

1) Loads the CD-HIT module.
2) Counts the number of protein sequences in the original FASTA (total_prots.fasta).
3)  Runs CD-HIT to cluster proteins at 60% identity and remove redundant sequences:
   
  `cd-hit -i total_prots.fasta -o prots_clean.fasta -c 0.6 -M 0 -n 4`

4) Reports statistics: total proteins before, after, and percentage retained.
5) Outputs the cleaned proteome in:
   
`$genome_analysis_path/genes_identification/Tarsynflow/proteome/prots_clean.fasta`

**Required Inputs:**

- $data_path/tp_data/total_prots.fasta → raw protein sequences from UniProt
- $genome_analysis_path → base folder for gene identification outputs

**Outputs:**

- $genome_analysis_path/genes_identification/Tarsynflow/proteome/prots_clean.fasta → reduced non-redundant protein database

**Dependencies / Software required:**

[![CD-HIT](https://img.shields.io/badge/CD--HIT-clustering-lightblue)](http://weizhongli-lab.org/cd-hit/)  
[![TarSynFlow](https://img.shields.io/badge/TarSynFlow-gene--annotation-green)](https://github.com/jiarong/TarSynFlow)  
[![SLURM](https://img.shields.io/badge/SLURM-job--scheduler-blue)](https://slurm.schedmd.com/)

### 🧩 genes_comps -- Compare gene presence across genomes.

#### This module:

- Performs pairwise genome comparisons to evaluate gene presence/absence across all genomes.
- Uses the non-redundant proteome from genes and performs synteny and coverage-based comparisons with AutoFlow workflows.
- Prepares results for downstream extraction and analysis (genes_results_protein and genes_results_annotation).

#### Workflow / Steps:

1) Loads AutoFlow environment.
2) Reads query genomes from $data_path/gen_queries and reference genomes from $data_path/gen_refs.
3) Prepares genome FASTA files by prefixing sequence IDs (>Q_ for query, >R_ for reference).
4) Runs AutoFlow with the template workflow_genome_sinteny_with_proteins, specifying:
   
- Proteome dataset (prots_clean.fasta)
- Coverage and identity thresholds (85+85)
- Output folder per query–reference pair

5) Results are stored per query–reference genome pair for further analysis.

**Required Inputs:**

- $data_path/genomes_problem/ → genome FASTA files
- $data_path/gen_queries → list of query genomes
- $data_path/gen_refs → list of reference genomes
- $genome_analysis_path/genes_identification/Tarsynflow/proteome/prots_clean.fasta → cleaned proteome dataset
- $template_path/workflow_genome_sinteny_with_proteins → AutoFlow workflow template

**Outputs:**

- $genome_analysis_path/genes_identification/Tarsynflow/comps/<QUERY>/<REF>/ → per query–reference comparison results
- Presence/absence and coverage tables for each genome pair

**Dependencies / Software required:**

[![AutoFlow](https://img.shields.io/badge/AutoFlow-genome--synteny-orange)](https://github.com/linlabcode/AutoFlow)

### 🧩 genes_results_protein -- Extract protein results.

#### This module:

- Consolidates and analyzes protein-level results from gene identification and genome comparisons.
- Extracts protein matches from comparative analyses and generates frequency tables, filtered lists, and annotated datasets.
- Produces visualizations like heatmaps for genome–protein relationships.

#### Workflow / Steps:

1) Loads the AutoFlow environment.
2) Iterates over query and reference genomes (gen_queries and gen_refs).
3) Checks completed AutoFlow workflows (flow_logger) and copies Circos images.
4) Generates lists of:

- Query-specific proteins
- Reference-specific proteins
- Shared proteins (sure and putative)
- Not matching proteins (sure and putative)

5) Computes protein frequency tables across genomes.
6) Optionally, maps proteins to UniProt to retrieve annotations (get_annotations) using the REST API.
7) Generates a matrix of genome–protein relationships (pairs2matrix.rb) and visualizes it as a heatmap (plot_heatmap.R).

**Required Inputs:**

- $genome_analysis_path/genes_identification/Tarsynflow/comps/ → genome comparison outputs
- $data_path/ → project data folder with gen_queries and gen_refs
- $genome_analysis_path/genes_identification/Tarsynflow/results/ → TarSynFlow results
- $script_path/ → supporting scripts (pairs2matrix.rb, plot_heatmap.R)

**Outputs:**

- $output/matches_analysis/ → protein frequency and presence/absence tables
- $output/protein_lists/ → filtered protein lists (query, reference, shared, non-matching)
- Annotated protein tables (via UniProt API)
- Heatmaps: $output/matches_analysis/heatmap.pdf

**Dependencies / Software required:**

[![AutoFlow](https://img.shields.io/badge/AutoFlow-genome--synteny-orange)](https://github.com/linlabcode/AutoFlow)  
[![R](https://img.shields.io/badge/R-data--visualization-blue)](https://www.r-project.org/)  
[![UniProt API](https://img.shields.io/badge/UniProt-annotations-purple)](https://www.uniprot.org/help/api)  
[![jq](https://img.shields.io/badge/jq-json--parsing-lightblue)](https://stedolan.github.io/jq/)  
[![Ruby](https://img.shields.io/badge/Ruby-scripting-red)](https://www.ruby-lang.org/)

### 🧩 genes_results_annotation	-- Extract annotation results.

#### This module:

- Retrieves functional annotations for proteins identified in the gene analysis workflow.
- Filters proteins by frequency and maps them to UniProt entries to get:
    - Protein names
    - Gene names
    - Organism
    - Sequence length
    - Reviewed/unreviewed status

#### Workflow / Steps:

1) Calls get_all_results.sh with get_annotations mode.
2) Iterates over query and reference genomes.
3) Filters proteins by a minimum occurrence threshold (min_count).
4) Maps filtered proteins to UniProt via the REST API and waits for job completion.
5) Generates annotated TSV files for each protein set.
6) Produces matrices and visualizations (heatmaps) for genome–protein annotation relationships.

**Required Inputs:**

- $genome_analysis_path/genes_identification/Tarsynflow/results/ → gene annotation results
- $data_path/ → project data folder (gen_queries and gen_refs)
- $genome_analysis_path/genes_identification/Tarsynflow/comps/ → genome comparison outputs
- $script_path/ → supporting scripts (pairs2matrix.rb, plot_heatmap.R)

**Outputs:**

- $output/Query_specific_filtered_annotated_prots → annotated query proteins
- $output/Ref_specific_filtered_annotated_prots → annotated reference proteins
- Filtered protein lists for shared and non-matching proteins
- Matrices and heatmaps summarizing annotations

**Dependencies / Software required:**

[![AutoFlow](https://img.shields.io/badge/AutoFlow-genome--synteny-orange)](https://github.com/linlabcode/AutoFlow)  
[![R](https://img.shields.io/badge/R-data--visualization-blue)](https://www.r-project.org/)  
[![UniProt API](https://img.shields.io/badge/UniProt-annotations-purple)](https://www.uniprot.org/help/api)  
[![jq](https://img.shields.io/badge/jq-json--parsing-lightblue)](https://stedolan.github.io/jq/)  
[![Ruby](https://img.shields.io/badge/Ruby-scripting-red)](https://www.ruby-lang.org/)

### 🧩 seqs	-- Extract specific sequences.

#### This module:

- Extracts specific genomic sequences based on gene/protein annotations.
- Works on both query and reference genomes, using coordinate tables generated in the gene comparison workflow.
- Adds flanking regions (150 bp) around each feature for downstream analyses.

#### Workflow / Steps:

1) Prepare output folder: $output/extracted_seqs.
2) Generate lists of query-specific and reference-specific proteins from filtered annotation results (Query_specific_filtered_annotated_prots and Ref_specific_filtered_annotated_prots).
3) Extract the corresponding coordinates from the coordinate tables (coord_table_with_strand).
4) Use get_seqs.rb to retrieve sequences from genome FASTA files, including 150 bp flanking regions.
5) Repeat for each reference genome listed in $data/gen_refs.
6) Output is organized per genome and per query/reference type.

**Required Inputs:**

- $coord_table → coordinates for query genome sequences
- $genome_path → folder with genome FASTA files
- $output → results folder for extracted sequences
- $data → project data folder (gen_refs)
- $script_path → supporting scripts (get_seqs.rb)

**Outputs:**

- Query sequences: $output/extracted_seqs/seqs_query.fasta
- Reference sequences: $output/extracted_seqs/seqs_ref_<REF>.fasta
- Coordinate tables used for extraction: $output/extracted_seqs/seqs_coords_*.txt
- Lists of protein IDs used: $output/extracted_seqs/seqs_list_*.txt and $output/extracted_seqs/seqs_*.lst

**Dependencies / Software required:**

[![Ruby](https://img.shields.io/badge/Ruby-scripting-red)](https://www.ruby-lang.org/)  
[![get_seqs.rb](https://img.shields.io/badge/get__seqs.rb-sequence--extraction-lightgrey)](path_to_scripts)  

### 🧩 GI_up -- Detect genomic islands.

This module uploads annotated genome files (GenBank format) to IslandViewer4 for genomic island prediction. For each genome, it generates a job token that can later be used to download the results.

**Required Inputs:**

- $data_path/ref_genome_GI → list of reference genomes and accession numbers
- $genome_analysis_path/char/dfast_0000/genome_annotation/<genome>_file/genome.gbk → annotated genome files

**Outputs:**

- $genome_analysis_path/genomic_island/Island_Viewer4_results/<genome>_submit → submission response
- $genome_analysis_path/genomic_island/Island_Viewer4_results/<genome>_jobtoken → job tokens for each genome

 **Dependencies / Software required:** 
[![curl](https://img.shields.io/badge/curl-command--line--tool-blue)](https://curl.se/)

### 🧩 GI_result -- Summarize genomic island results.

This module downloads the results of Genomic Island predictions from Island Viewer 4 for all submitted genomes. Using previously obtained job tokens, it queries the Island Viewer REST API and saves the CSV files containing genomic island coordinates and annotations.

**Required Inputs:**

- $data_path/all_genome_list → list of genomes to download
- Job tokens from the GI_up module
- $genome_analysis_path/genomic_island/Island_Viewer4_results/ → submission results folder

**Outputs:**

- $genome_analysis_path/genomic_island/Island_Viewer4_results/file_scv/<genome>.csv → raw CSV results per genome
  
**Dependencies / Software required:** 

[![curl](https://img.shields.io/badge/curl-command--line--tool-blue)](https://curl.se/)

### 🧩 GI_clean -- Clean genomic island outputs.

This module processes the raw Genomic Island (GI) results downloaded from IslandViewer4. It performs:

1) Cleaning GI data

- Filters and formats GI tables (_Integrated, _LOCUS)
- Counts total GI regions per genome
- Generates GI locus lists for further analysis

2) Enrichment analysis

- Calculates relative GI sizes using genome length
- Creates LOCUS lists for COG categories
- Performs functional enrichment using clusters_to_enrichment.R
- Generates relative enrichment tables per genome
- Cleans temporary or intermediate files

**Required Inputs:**

- $genome_analysis_path/genomic_island/Island_Viewer4_results/file_scv/ → raw GI CSV files
- $data_path/COG_categories → list of COG categories
- $genome_analysis_path/char/dfast_0000/genome_annotation/results_dfast_parser/<genome>/length → genome lengths

**Outputs:**

- $genome_analysis_path/genomic_island/genomic_island_results/<genome>/* → cleaned GI tables, locus lists, and enrichment tables
- Merged tables across all genomes: enrichment_GI_category, enrichment_GI_category_relative, Total_GI, Total_GI_relative

**Dependencies / Software required:** 

[![R](https://img.shields.io/badge/R-statistics--environment-blue)](https://www.r-project.org/)
[![Ruby](https://img.shields.io/badge/Ruby-scripting-red)](https://www.ruby-lang.org/)

### 🧩 phage_visorter -- Identify prophages with PHASTEST.

This module identifies and quantifies prophages in all genomes using VirSorter. For each genome, it calculates the total number of prophages and generates both absolute and relative counts based on genome length.

**Required Inputs:**

- $data_path/total_genomes/*.fasta → genomes to analyze
- $data_path/all_genome_list → list of genome filenames
- 
**Outputs:**

- Per-genome prophage counts: prophage/<genome>/final_number, final_number_tab
- Absolute and relative prophage counts: final_number_absolute, final_number_relative
- Total across all genomes: prophage/Total_phage, Total_phage_relative

**Dependencies / Software required:** 

[![VirSorter](https://img.shields.io/badge/VirSorter-prophage--detection-purple)](https://github.com/simroux/VirSorter)

### 🧩 results -- Integrate all analysis results.
The results module centralizes and compiles the final outputs from the genomic workflow. It integrates multiple analyses, including 16S rRNA BLAST, whole-genome comparisons, synteny visualization, gene annotation, mobile elements detection (transposons and prophages), and genomic island enrichment. The module standardizes strain and gene names and generates visualization-ready tables and figures.

### Workflow Steps:

1) Directory Setup:
Cleans and creates the results_path and report_img directories for organized outputs.

2) Strain Name Standardization:
Reformats genome and strain identifiers for consistent visualization. Replaces uncultured or incomplete labels with standardized names.

3) 16S rRNA BLAST Processing:
Filters hits by identity and coverage thresholds and prepares tables for downstream comparison.

4) Whole Genome Comparison:
Processes pyANI results for identity and coverage matrices, selects reference genomes for synteny analysis, and extracts Sibelia synteny images for the highest-identity genome pairs.

5) Genome Annotation and Functional Analysis:
Processes DFAST outputs, generating COG functional annotation tables and genomic island enrichment data.

6) Mobile Element Analysis:
Handles transposons and prophages. Aggregates absolute and relative counts, identifies interrupted proteins, and generates unified tables for comparison across strains.

7) Integration and Report Generation:
Combines Genomic Islands, prophages, and transposons into integrated tables. Generates strain-specific and functional summary tables in integration/ and visualization-ready outputs in report_img/.

**Required Inputs:**

- Genome sequences
- DFAST annotation results
- Transposon and prophage analysis results
- Genomic island enrichment results
- pyANI and Sibelia synteny outputs
- External Dependencies:

#### Notes:

- Outputs are standardized for cross-strain comparison.
- Produces both absolute counts and normalized percentages for mobile elements and genomic islands.
- Visualization-ready tables and figures are generated for reporting or downstream analysis.
  
### 🧩 PCA -- Run PCA on genomic features.

The PCA module performs a principal component analysis (PCA) on genomic features for all strains. It consolidates genomic feature tables, transposes them for proper input formatting, and runs a multivariate analysis workflow to identify patterns of similarity or divergence among strains.

#### Workflow Steps:

1) Data Preparation:

- Collects integrated genomic features from previous modules:
- GI_Functions – COG functions of genomic islands
- Whole_Genome_Functions – COG functions across the whole genome
- Relative_Mobilome_Number – counts of transposons, genomic islands, and prophages
- Transposes tables to have strains as rows and features as columns.

2) Analysis Execution:

- Submits the preprocessed tables to the PCA workflow via mfa_analysis_all.sh.
- PCA is run to identify major axes of variation and relationships among strains.

3) Visualization and Interpretation:
   
- Produces PCA plots for visual comparison of strains.
- Can help highlight clustering based on genomic content or mobile element patterns.

**Required Inputs:**

- GI_Functions
- Whole_Genome_Functions
- Relative_Mobilome_Number

**Dependencies / Software required:** 

[![Python](https://img.shields.io/badge/Python-script-blue)]
[![transpose_table](https://img.shields.io/badge/transpose_table-tool-lightgrey)]
[![Slurm](https://img.shields.io/badge/Slurm-job_scheduler-orange)]

#### Notes:

- Tables are preprocessed and transposed for PCA.
- Output includes both PCA scores and plots for strain comparison.
  
### 🧩 PCA_filter -- Filter and run PCA on genomic features

The PCA_filter module performs filtering of genomic features based on clustering results and prepares the dataset for PCA (Principal Component Analysis) or MFA (Multiple Factor Analysis). This step ensures that only relevant samples and features are included in the downstream analysis.

#### Workflow

1) Filtering Samples

- Reads the clustering results from results/mfa/PCA_results/mfa_hcpc.txt.
- Removes samples assigned to cluster 4 ($23 != 4).
- Creates a list of samples to exclude (eliminar.tsv).

2) Filtering Features

- Removes the excluded samples from the following datasets:
    - GI_Functions
    - Whole_Genome_Functions
    - Relative_Mobilome_Number
- Stores filtered datasets in $data_path/mfa/filter.

3) Run PCA/MFA

- Sources Python environment (~soft_bio_267/initializes/init_python).
- Submits mfa_analysis_all.sh using sbatch on the filtered datasets.

**Dependencies / Software required:** 

[![Python](https://img.shields.io/badge/Python-script-blue)]
[![transpose_table](https://img.shields.io/badge/transpose_table-tool-lightgrey)]
[![Slurm](https://img.shields.io/badge/Slurm-job_scheduler-orange)]

### 🧩 report -- Generate final HTML report with figures.

This module generates a comprehensive genome analysis report in HTML format using the report_html Ruby tool and the template report_template.erb. The report integrates all the results produced by previous modules, visualizing the genomic features, annotations, and comparisons of Shewanella strains.

#### Functionality

1) Input data:

- Precomputed results from previous modules:
    - 16S RNA BLAST results
    - Complete genome annotation (DFAST)
    - Pyani genome comparisons (ANI and coverage)
    - Synteny blocks (Sibelia)
    - Transposon identification (Tp)
    - Genomic islands (GIs) and enrichment analysis
    - Prophage identification
    - Tarsynflow-specific gene analysis for probiotic/pathogenic strains

2) Report structure:

- Header: Project title and general experiment description.
- 16S RNA Analysis: Table of top 16S RNA hits.
- Genome Comparisons:
    - Pyani heatmaps (identity and coverage)
    - Sibelia synteny block visualizations
    - Custom comparisons between reference genomes and target strains
- Genome Annotation:
    - Heatmaps for COG functional categories (absolute and relative values)
- Transposon Analysis:
  - Absolute and relative TE counts
  - Heatmaps for transposases and interrupted proteins
- Genomic Islands:
  - Total GIs per strain
  - GI functional enrichment (absolute and relative)
- Prophage Analysis:
  - Total and relative counts per strain
- Pathogenic vs Probiotic Comparison:
  - Table of specific genes identified by Tarsynflow

3) Output:

- HTML report located at: $results_path/project_report
- Embedded heatmaps, barplots, and tables for interactive exploration.

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

- **NCBI Datasets API** – For downloading and managing genomic datasets.  
- **PHASTEST** – For prophage identification in bacterial genomes.  
- **Roux et al., 2021** – *TarSynFlow: Genome-wide gene and mobile element analysis*.  
- **DFAST** – Genome annotation and COG functional assignment.  
- **Pyani** – Average Nucleotide Identity (ANI) calculations for genome comparisons. [https://github.com/widdowquinn/pyani](https://github.com/widdowquinn/pyani)  
- **Sibelia** – Synteny block identification for genome structure comparisons.  
- **ISEScan** – Detection of transposable elements. [https://github.com/xiezhq/ISEScan](https://github.com/xiezhq/ISEScan)  
- **IslandViewer4** – Genomic island identification. [https://www.pathogenomics.sfu.ca/islandviewer/](https://www.pathogenomics.sfu.ca/islandviewer/)  
- **ExpHunterSuite** – Functional enrichment analysis for genomic islands. [https://github.com/seoanezonjic/ExpHunterSuite](https://github.com/seoanezonjic/ExpHunterSuite)  
- **VirSorter2** – Prophage detection in bacterial genomes. [https://github.com/jiarong/VirSorter2](https://github.com/jiarong/VirSorter2)  

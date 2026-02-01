# APG-OHEC-Retro-M2: Plasmid Analysis and Data Visualization Repository

This repository contains the data, scripts, and supplementary materials for the manuscript "The Escherichia coli F plasmidome: an Australian perspective" by Cummins, Watt and Donato et al., 2025. The study examines plasmid distribution across different E. coli sources and phylogroups, with a focus on F plasmids and their virulence and antimicrobial resistance gene carriage.

## Repository Structure
```
APG-OHEC-Retro-M2/
├── README.md                                    # This readme file
├── analysis/                                    # Primary analysis data files
│   ├── RapidNJ_5471.nwk                         # Phylogenetic tree file in Newick format
│   ├── cgMLST_dists_5471.txt.gz                 # Compressed core genome MLST distance matrix
│   ├── genotype.txt.gz                          # Compressed genotype data from ABRicate
│   ├── mobsuite_concatenated_mge_report.txt     # MOB-suite mobile genetic element report
│   ├── mobsuite_ref_plasmid_pipelord/           # MOB-suite plasmid reference analyses
│   │   ├── mobsuite_refs_IncF_RST.txt           # IncF replicon sequence typing data
│   │   ├── mobsuite_refs_IncF_RST_full.txt      # Full IncF RST results
│   │   ├── mobsuite_refs_abritamr_partials.txt  # Partial AMR genes from ABRicate
│   │   ├── mobsuite_refs_abritamr_resistance.txt # Resistance gene data from ABRicate
│   │   ├── mobsuite_refs_abritamr_virulence.txt # Virulence gene data from ABRicate
│   │   ├── mobsuite_refs_amrfinder_raw.txt      # Raw AMRFinder results
│   │   ├── mobsuite_refs_genotype.txt           # Plasmid genotype data
│   │   └── mobsuite_refs_pMLST.txt              # Plasmid MLST typing results
│   └── mobtyper_results_concatenated.txt        # Consolidated MOB-typer results
├── delims/                                      # Delimiter and metadata files
│   ├── mobsuite/
│   │   └── clusters.txt                         # MOB-suite plasmid cluster definitions
│   ├── plasmid_cluster_data.txt.gz              # Compressed plasmid cluster metadata
├── markdown/                                    # R Markdown files for analysis
│   └── Figure_Generation_Data_Processing.Rmd    # Main analysis workflow for figure generation
├── scripts/                                     # R scripts for data processing
│   └── plasmid_data_aggregate.R                 # Script to aggregate plasmid data
└── supplementary/                               # Supplementary tables for manuscript
    ├── Supplementary_Table_1.txt                # Metadata and plasmid presence for all strains
    ├── Supplementary_Table_1_base.txt           # Base data for Supplementary Table 1
    └── Supplementary_Table_2.txt                # Plasmid cluster characteristics
```

## Study Overview

This study analyzes the distribution and characteristics of plasmids in a large collection of E. coli isolates (n>5,000) from diverse sources including humans, companion animals, livestock, wild animals, food, and environmental samples. The analysis focuses on:

1. Plasmid distribution across different E. coli sources and phylogroups
2. Identification of major plasmid clusters and their associated replicon types
3. Virulence and antimicrobial resistance gene carriage on plasmids
4. Comparative analysis of ColV and senB-associated plasmids
5. Plasmid-mediated horizontal gene transfer across different sources

## Key Data Files

- RapidNJ_5471.nwk: Phylogenetic tree file for the 5,471 E. coli genomes analyzed
- cgMLST_dists_5471.txt.gz: Distance matrix of core genome MLST differences between isolates
- mobtyper_results_concatenated.txt: MOB-typer results for plasmid typing and classification
- plasmid_cluster_data.txt.gz: Processed data on plasmid clusters including replicon types, AMR/VF gene carriage, and source distribution

## Analysis Workflow
The main analysis workflow is contained in Figure_Generation_Data_Processing.Rmd, which:

Processes genomic and plasmid data from multiple sources
Generates tables and figures for manuscript publication
Performs statistical analyses on plasmid distribution by source and phylogroup
Creates visualizations for plasmid clusters, their gene content, and distribution patterns
The plasmid_data_aggregate.R script performs data aggregation of plasmid metadata, including:

- Replicon typing
- AMR gene detection
- Virulence gene detection
- Plasmid clustering
- Supplementary Materials

The repository includes supplementary tables referenced in the manuscript:

Supplementary_Table_1.txt: Comprehensive metadata for all strains, including plasmid carriage
Supplementary_Table_2.txt: Detailed characteristics of plasmid clusters identified in the study

# Usage
To reproduce the analysis:

- Clone this repository
- Install required R packages (tidyverse, vroom, ggplot2, ComplexHeatmap, ggpubr, etc.)
- Run the R Markdown file Figure_Generation_Data_Processing.Rmd

## Citation
If you use data or code from this repository, please cite:

Cummins ML, Watt A., Celese D. et al. (2025) The Escherichia coli F plasmidome: an Australian perspective. [Journal details pending]

## Contact
For questions or further information, please contact the corresponding authors.

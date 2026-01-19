# 🧫 Amplicon Sequencing Analysis of Microbial Community Shifts Following Biofilm Disruption

## 🧠 Project Summary

This project involved a statistical analysis of microbial communities to assess the effectiveness of a chemical treatment in disrupting complex biofilms in drainage systems. The primary goal was to evaluate shifts in microbial composition before and after treatment using taxonomic profiling techniques.

Due to confidentiality restrictions, raw data and detailed outputs are not included. However, the workflow, tools, and overall findings are shared below to demonstrate technical approach and insight.

---

## 🔬 Key Findings

- The chemical product effectively **eliminated several dominant bacterial genera** known to contribute to biofilm formation.
- However, the treatment **enabled proliferation of other resistant genera**, which eventually re-established a new biofilm community.
- These findings suggest the product alters microbial ecosystems, potentially requiring combination treatments to prevent biofilm reformation.

---

## 🧰 Tools and Methods

This repository contains R scripts for analysing amplicon sequencing data of microbiome. The scripts process, analyse, and visualise microbial community data — starting from raw data processing (quality control, filtering, trimming, and merging) to taxonomic assignment of amplicons and plotting richness and diversity metrics.

The workflow adapts an existing method from the DADA2 tutorial, tailored to the specific context of biofilm disruption. This adaptation demonstrates proficiency in:

- R scripting for bioinformatics workflows

- Statistical analysis of microbial diversity

- Data visualisation using ggplot2 and phyloseq

- Reproducibility through structured and annotated code


Environment
- RStudio 2024.12.1+563

- R version 4.4.2


Sequencing 

Amplicon sequencing of 16S rDNA V3–V4 (50K-WBI)



Analysis Focus

- **Alpha diversity analysis** using Shannon and Simpson indices to assess within-sample microbial diversity across time points.
- **Beta diversity analysis** using Bray–Curtis dissimilarity and non-metric multidimensional scaling (NMDS) to evaluate differences in community composition before and after treatment.

---



\## 📂 Repository Structure


metagenomics-analysis/


├── Amplicon_seq_analysis_tax_profiling.R    # Main R script containing all analysis steps

├── results/                  # Example outputs (non-sensitive)

├── LICENSE                  # MIT License

└── README.md                 # This documentation file



---

## 🛠 Requirements

- **R** ≥ 4.4.2  
- **RStudio** ≥ 2024.12.1  

### Required R packages:

- `dada2`
- `Biostrings`
- `phyloseq`
- `ggplot2`

**Reference database** (download separately):  
`silva_nr99_v138.2_toGenus_trainset.fa.gz`

### Install dependencies in R:

```r
install.packages(c("phyloseq", "ggplot2"))

if (!requireNamespace("BiocManager", quietly = TRUE))
    install.packages("BiocManager")

BiocManager::install(c("Biostrings", "dada2"))





🚀 Usage

1- Clone the repository:


git clone https://github.com/RazanAbb/metagenomics-analysis.git

cd metagenomics-analysis



2- Open metagenomic_analysis.R in RStudio.

3- Run the script in order (sections are clearly marked in the script):

Quality control and filtering

Error frequency plotting

Taxonomic assignment

Diversity analysis and plotting

4- Adjust file paths in the script to point to your input data.

5- Output figures and tables will be saved in the results/ folder.





Licence



This project is licensed under the [MIT License](https://opensource.org/license/mit-0).






Author


Dr. Razan Abbara
PhD in Microbiology | Metagenomics & Bioinformatics Specialist









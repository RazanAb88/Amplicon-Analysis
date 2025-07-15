\# Metagenomic Analysis Scripts



This repository contains R scripts for performing metagenomic analysis of microbiome data. The scripts process, analyse, and visualise microbial community data — starting from raw data processing (quality control, filtering, trimming, and merging) to taxonomic assignment of amplicons and plotting richness and diversity metrics.



Developed in RStudio 2024.12.1+563 with R version 4.4.2.



The dataset involves amplicon sequencing of 16S rDNA V3-V4 (50K-WBI), to investigate changes in microbial communities:



\- \*\*Alpha diversity\*\*: within-sample diversity.  

\- \*\*Beta diversity\*\*: between-sample diversity.  



---



\## 📂 Repository Structure



metagenomics-analysis/

├── metagenomic\_analysis.R    # Main R script containing all analysis steps

├── results/                  # Output figures and tables

└── README.md                 # This documentation file













\## 🛠 Requirements



R >= 4.4.2



RStudio >= 2024.12.1



R packages:



dada2



Biostrings



phyloseq



ggplot2



Database file: silva\_nr99\_v138.2\_toGenus\_trainset.fa.gz (please download separately)







Install dependencies with:

install.packages(c("phyloseq", "ggplot2"))



if (!requireNamespace("BiocManager", quietly = TRUE))

&nbsp; install.packages("BiocManager")

BiocManager::install(c("Biostrings", "dada2"))









🚀 Usage

1- Clone the repository:



git clone https://github.com/RazanAbb/metagenomics-analysis.git

cd metagenomics-analysis



2- Open the script file metagenomic analysis.R in RStudio.



3- Run the script in order:



Quality control and filtering



Error frequency plotting



Taxonomic assignment



Diversity analysis and plotting



(These are marked as sections in the single script.)



4- Adjust file paths inside the script to point to your input data as needed.



5- Generated figures and tables will be saved in the results/ folder.





Licence



This project is licensed under the [MIT License](https://opensource.org/license/mit-0).







Author



Dr. Razan Abbara








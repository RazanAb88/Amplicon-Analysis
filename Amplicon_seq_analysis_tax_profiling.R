library(dada2); packageVersion("dada2")
path <- "the path to you working directory/Rstudio"

##pltoing the quality profiles of two reads of fowrward and reverse raw reads
fnFs <- sort(list.files(path, pattern=".raw_1.fastq.gz", full.names = TRUE))
fnRs <- sort(list.files(path, pattern=".raw_2.fastq.gz", full.names = TRUE))
sample.names <- sapply(strsplit(basename(fnFs), "_"), `[`, 1)

plotQualityProfile(fnFs[1:2])
plotQualityProfile(fnRs[1:2])


# Place filtered files in filtered/ subdirectory
filtFs <- file.path(path, "filtered", paste0(sample.names, "_F_filt.fastq.gz"))
filtRs <- file.path(path, "filtered", paste0(sample.names, "_R_filt.fastq.gz"))
names(filtFs) <- sample.names
names(filtRs) <- sample.names

print(filtFs)

out <- filterAndTrim(fnFs, filtFs, fnRs, filtRs, truncLen=c(250,250),
                     maxN=0, maxEE=c(2,2), truncQ=2, rm.phix=TRUE,
                     compress=TRUE, multithread= FALSE ) # On Windows set multithread=FALSE
head(out)

##calculating and plotting of the error rates 

errF <- learnErrors(filtFs, multithread=FALSE)
errR <- learnErrors(filtRs, multithread=FALSE)
plotErrors(errF, nominalQ=TRUE)
dadaFs <- dada(filtFs, err=errF, multithread=FALSE)
dadaRs <- dada(filtRs, err=errR, multithread=FALSE)
dadaFs[[1]]
mergers <- mergePairs(dadaFs, filtFs, dadaRs, filtRs, verbose=TRUE)

# Inspect the merger data.frame from the first sample
head(mergers[[1]])
seqtab <- makeSequenceTable(mergers)
dim(seqtab)
# Inspect distribution of sequence lengths
table(nchar(getSequences(seqtab)))
seqtab.nochim <- removeBimeraDenovo(seqtab, method="consensus", multithread=TRUE, verbose=TRUE)
dim(seqtab.nochim)
sum(seqtab.nochim)/sum(seqtab)
getN <- function(x) sum(getUniques(x))
track <- cbind(out, sapply(dadaFs, getN), sapply(dadaRs, getN), sapply(mergers, getN), rowSums(seqtab.nochim))
# If processing a single sample, remove the sapply calls: e.g. replace sapply(dadaFs, getN) with getN(dadaFs)
colnames(track) <- c("input", "filtered", "denoisedF", "denoisedR", "merged", "nonchim")
rownames(track) <- sample.names
head(track)


####Assign taxononmy using Silva v138.2
taxa <- assignTaxonomy(seqtab.nochim, "C:/Users/razon/Downloads/silva_nr99_v138.2_toGenus_trainset.fa.gz", multithread=FALSE)
taxa.print <- taxa # Removing sequence rownames for display only
rownames(taxa.print) <- NULL
head(taxa.print)




library(phyloseq); packageVersion("phyloseq")
library(Biostrings); packageVersion("Biostrings")
library(ggplot2); packageVersion("ggplot2")
theme_set(theme_bw())

##for Studying the disctubution of microbial genus and species in samples, some changes need to be done on dealing with samples names 

# Split sample names by "."
split_samples <- strsplit(samples.out, "\\.")

# Extract farm/site (first part)
farm <- sapply(split_samples, `[`, 1)

# Extract day (second part, removing 'D' and converting to integer)
day <- as.integer(sub("D", "", sapply(split_samples, `[`, 2)))

# Extract time (Before/After) only if Day == 1
time <- ifelse(day == 1, sapply(split_samples, `[`, 3), NA)


# Extract site (always second-to-last element)
site <- sapply(split_samples, function(x) x[length(x) - 1])
samples.out <- gsub("\\.fastq\\.gz$|\\.gz$", "", samples.out)




# Convert "Day" into a factor to maintain order (D1, D2, ..., D14)
samdf$Day <- factor(samdf$Day, levels = sort(unique(samdf$Day)))

# Split data by Farm
samdf_F1 <- subset(samdf, farm == "F1")
samdf_F2 <- subset(samdf, farm == "F2")

# Print the two separate farm datasets
print(samdf_F1)
print(samdf_F2)
rownames(samdf) <- samples.out
ls()

sample_data(ps)$Day <- factor(sample_data(ps)$Day)
sd <- as.data.frame(sample_data(ps))
sd$Day <- factor(sd$Day)
sample_data(ps) <- sample_data(sd)


ps <- phyloseq(otu_table(seqtab.nochim, taxa_are_rows=FALSE), 
               sample_data(samdf), 
               tax_table(taxa))
ps <- prune_samples(sample_names(ps), ps)


dna <- Biostrings::DNAStringSet(taxa_names(ps))
names(dna) <- taxa_names(ps)
ps <- merge_phyloseq(ps, dna)
taxa_names(ps) <- paste0("ASV", seq(ntaxa(ps)))
ps

sample_data(ps)

##richness plot
plot_richness(ps, x="Day", measures=c("Shannon", "Simpson"), color="Subject")

##NMDS PLOT
ps.prop <- transform_sample_counts(ps, function(otu) otu/sum(otu))
ord.nmds.bray <- ordinate(ps.prop, method="NMDS", distance="bray")

plot_ordination(ps.prop, ord.nmds.bray, color="Subject", title="Bray NMDS")

top20 <- names(sort(taxa_sums(ps), decreasing=TRUE))[1:20]
ps.top20 <- transform_sample_counts(ps, function(OTU) OTU/sum(OTU))
ps.top20 <- prune_taxa(top20, ps.top20)
plot_bar(ps.top20, x="Day", fill="Family") + facet_wrap(~Subject, scales="free_x")

##genus level
print(ps.genus) 
colnames(tax_table(ps.genus))
sample_data(ps.genus)$Day <- factor(sample_data(ps.genus)$Day)
top20_genus <- names(sort(taxa_sums(ps.genus), decreasing=TRUE))[1:20]
ps.top20.genus <- prune_taxa(top20_genus, ps.genus.rel)
plot_bar(ps.top20.genus, x="Day", fill="Genus") + 
  facet_wrap(~Subject, scales="free_x")
top50_genus <- names(sort(taxa_sums(ps.genus), decreasing=TRUE))[1:50]
ps.top50.genus <- prune_taxa(top50_genus, ps.genus.rel)

plot_bar(ps.top50.genus, x="Day", fill="Genus") + 
  facet_wrap(~Subject, scales="free_x")

##family level
# Aggregate at Family level
ps.family <- tax_glom(ps, taxrank = "Family")  

# Normalise relative abundances (optional)
ps.family <- transform_sample_counts(ps.family, function(x) x / sum(x))  

# Plot at Family level
ps.family <- tax_glom(ps, "Family")
ps.family.rel <- transform_sample_counts(ps.family, function(x) x / sum(x))

top20_family <- names(sort(taxa_sums(ps.family), decreasing=TRUE))[1:20]
ps.top20.family <- prune_taxa(top20_family, ps.family.rel)

plot_bar(ps.top20.family, x="Day", fill="Family") + 
  facet_wrap(~Subject, scales="free_x")

top50_family <- names(sort(taxa_sums(ps.family), decreasing=TRUE))[1:50]
ps.top50.family <- prune_taxa(top50_family, ps.family.rel)

plot_bar(ps.top50.family, x="Day", fill="Family") + 
  facet_wrap(~Subject, scales="free_x")


library(dada2)

install.packages("dada2")
install.packages("survival")

library("devtools")
if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")
BiocManager::install("dada2")

devtools::install_github("benjjneb/dada2")
path <- "/home/svenja/sratoolkit.3.0.0-ubuntu64/Fastq-files/fastq_trim"
list.files(path)

fnFs <- sort(list.files(path, pattern="_1.trimmed.fastq.gz", full.names = TRUE))
fnRs <- sort(list.files(path, pattern="_2.trimmed.fastq.gz", full.names = TRUE))

sample.names <- sapply(strsplit(basename(fnFs), "_"), `[`, 1)


####visualization of the quality profile and determination of the trimming points
plotQualityProfile(fnFs[4:5])
Q_fnFS <- plotQualityProfile(fnFs, n = 5e+05, aggregate = TRUE)
Q_fnFS
Q_fnRS <- plotQualityProfile(fnRs, n = 5e+05, aggregate = TRUE)
Q_fnRS

filtFs <- file.path(path, "filtered", paste0(sample.names, "_F_filt.fastq.gz"))
filtRs <- file.path(path, "filtered", paste0(sample.names, "_R_filt.fastq.gz"))
names(filtFs) <- sample.names
names(filtRs) <- sample.names
##decied based on the two figures for 240 and 200

###trimming
out <- filterAndTrim(fnFs, filtFs, fnRs, filtRs, truncLen=c(240,200),
                     maxN=0, maxEE=c(2,2), truncQ=2, rm.phix=TRUE,
                     compress=TRUE, multithread=TRUE) # On Windows set multithread=FALSE
head(out)
View(out)

### learn the error rates

errF <- learnErrors(filtFs, multithread=TRUE)
errR <- learnErrors(filtRs, multithread=TRUE)

###Sample Inference

dadaFs <- dada(filtFs, err=errF, multithread=TRUE)
dadaRs <- dada(filtRs, err=errR, multithread=TRUE)



###merge paired reads
mergers <- mergePairs(dadaFs, filtFs, dadaRs, filtRs, verbose=TRUE)

###construct sequence table
seqtab <- makeSequenceTable(mergers)
dim(seqtab)

# Inspect distribution of sequence lengths
table(nchar(getSequences(seqtab)))

###remove chimeras
seqtab.nochim <- removeBimeraDenovo(seqtab, method="consensus", multithread=TRUE, verbose=TRUE)
dim(seqtab.nochim)

sum(seqtab.nochim)/sum(seqtab)


###track reads through the pipeline
getN <- function(x) sum(getUniques(x))
track <- cbind(out, sapply(dadaFs, getN), sapply(dadaRs, getN), sapply(mergers, getN), rowSums(seqtab.nochim))
# If processing a single sample, remove the sapply calls: e.g. replace sapply(dadaFs, getN) with getN(dadaFs)
colnames(track) <- c("input", "filtered", "denoisedF", "denoisedR", "merged", "nonchim")
rownames(track) <- sample.names
head(track)
View(track)

###assign taxonomy
taxa <- assignTaxonomy(seqtab.nochim, "/home/svenja/sratoolkit.3.0.0-ubuntu64/Fastq-files/silva_nr_v132_train_set.fa.gz", multithread=TRUE)
View(taxa)
taxa <- data.table(taxa)
summary(taxa$Genus)

table <- count(taxa, "Genus")
sort(table$freq, decreasing = TRUE)

install.packages("Phyloseq")
top20 <- names(sort(taxa_sums(Genus), TRUE) [1:20])
new <- taxa[taxa$Genus %in% names(sort(table(taxa$Genus), decreasing = TRUE)[1:10]), ]
unique(new$Genus)

ggplot(table, aes( x=Genus, y=freq)) +
  geom_bar(stat = "identity")

###species level
taxa_spec <- addSpecies(taxa, "/home/svenja/sratoolkit.3.0.0-ubuntu64/Fastq-files/silva_species_assignment_v132.fa.gz")

##inspection of taxonomic assignments
taxa.print <- taxa # Removing sequence rownames for display only
rownames(taxa.print) <- NULL
head(taxa.print)

View(taxa.print)
####
library(data.table)
metadata_reimagine <- fread("/home/svenja/sratoolkit.3.0.0-ubuntu64/Fastq-files/SraRunTable_regions.csv")



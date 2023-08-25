library(data.table)
library(ggplot2)

# dtA <- fread("data/mgx_abundances/dzhk_unmapped.csv"); dtA$cohort <- "This study"
# dtB <- fread("data/mgx_abundances/troci_unmapped.csv"); dtB$cohort <- "Troci et al.\n(2022)"
# dtC <- fread("data/mgx_abundances/lifelines_unmapped.csv"); dtC$cohort <- "Chen et al.\n(2021)"
# 
# dt <- rbindlist(list(dtA,dtB,dtC))
# dt$cohort <- factor(dt$cohort, levels =c("This study","Troci et al.\n(2022)","Chen et al.\n(2021)"))
# 
# p <- ggplot(dt, aes(cohort, (100 - unmapped)/100)) +
#   geom_boxplot(outlier.shape = NA, fill = "#cccccc") +
#   coord_cartesian(ylim = c(0,1)) +
#   scale_y_continuous(labels = scales::percent, expand = c(0,0)) +
#   labs(x = "Cohort", y = "QC reads mapped to HRGM genomes") +
#   theme_bw() +
#   theme(panel.grid.major.x = element_blank(),
#         panel.grid.minor.y = element_blank(),
#         panel.border = element_blank(),
#         axis.line = element_line())
# 
# ggsave("output/plots/mapped_reads.pdf", plot = p, width = 3.2, height = 4.7)
# 
# dt[, round(quantile(100-unmapped, probs = c(.25,0.5,.75)), digits = 1),  by = cohort]
# 
# 

Metadata <- fread("data/REPR_Genomes_metadata.tsv")
rel_hrgm <- Metadata[`Completeness (%)` >= 85 & `Contamination (%)` <= 2, `HRGM name`]

out <- list()

# DZHK
dt <- fread("data/mgx_abundances/dzhk_unmapped.csv")
dtmat <- read.table("data/mgx_abundances/dzhk_abun.csv")
dtmat <- dtmat[dt$sample, ]
dtmat <- dt[, (100-unmapped)/100 * Reads_pe] * dtmat

out[[1]] <- data.table(cohort = "This study",
                       sample = rownames(dtmat),
                       mappedHRGM = apply(dtmat, 1, sum)/dt$Reads_pe,
                       mappedHRGMfilt = apply(dtmat[, colnames(dtmat) %in% rel_hrgm], 1, sum)/dt$Reads_pe)

# Troci
dt <- fread("data/mgx_abundances/troci_unmapped.csv")
dtmat <- read.table("data/mgx_abundances/troci_abun.csv")
dtmat <- dtmat[dt$sample, ]
dtmat <- dt[, (100-unmapped)/100 * Reads_pe] * dtmat

out[[2]] <- data.table(cohort = "Troci et al.\n(2022)",
                       sample = rownames(dtmat),
                       mappedHRGM = apply(dtmat, 1, sum)/dt$Reads_pe,
                       mappedHRGMfilt = apply(dtmat[, colnames(dtmat) %in% rel_hrgm], 1, sum)/dt$Reads_pe)

# Chen
dt <- fread("data/mgx_abundances/lifelines_unmapped.csv")
dtmat <- read.table("data/mgx_abundances/lifelines_abun.csv")
dtmat <- dtmat[dt$sample, ]
dtmat <- dt[, (100-unmapped)/100 * Reads_pe] * dtmat

out[[3]] <- data.table(cohort = "Chen et al.\n(2021)",
                       sample = rownames(dtmat),
                       mappedHRGM = apply(dtmat, 1, sum)/dt$Reads_pe,
                       mappedHRGMfilt = apply(dtmat[, colnames(dtmat) %in% rel_hrgm], 1, sum)/dt$Reads_pe)


out <- rbindlist(out)
out <- melt(out, id.vars = c("cohort","sample"), value.name = "mapped",
            variable.name = "filter")

out[, filter := ifelse(filter == "mappedHRGM",
                       "all HRGM genomes (n=5 414)",
                       "filtered HRGM genomes (n=3 687)")]


out$cohort <- factor(out$cohort, levels =c("This study","Troci et al.\n(2022)","Chen et al.\n(2021)"))

p <- ggplot(out, aes(cohort, mapped, fill = filter)) +
  geom_boxplot(outlier.shape = NA) +
  coord_cartesian(ylim = c(0,1)) +
  scale_y_continuous(labels = scales::percent, expand = c(0,0)) +
  labs(x = "Cohort", y = "QC reads mapped to HRGM genomes",
       fill = "Filter") +
  theme_bw() +
  theme(panel.grid.major.x = element_blank(),
        panel.grid.minor.y = element_blank(),
        panel.border = element_blank(),
        axis.line = element_line())
p


ggsave("output/plots/mapped_reads.pdf", plot = p, width = 6, height = 4.7)

out[, round(quantile(mapped*100, probs = c(.25,0.5,.75)), digits = 1),  by = .(cohort,filter)]
out[, round(quantile(mapped*100, probs = c(.25,0.5,.75)), digits = 1), by = filter]

# Combined Popgen & Lifeline stability analysis

library(MicrobiomeGS2)
library(vegan)
library(abdiv)
library(ape)
library(ggplot2)
library(ggpubr)
library(ggtext)
library(egg)

sybil::SYBIL_SETTINGS("SOLVER","cplexAPI")

# get QC for MAGs
completeness_cuttoff <- 85
contamination_cutoff <- 2

# load HRGM model data
Metadata <- fread("data/REPR_Genomes_metadata.tsv")

# abundance data
LL_hrgm_abun <- t(read.table("data/mgx_abundances/lifelines_abun.csv"))
PG_hrgm_abun <- t(read.table("data/mgx_abundances/troci_abun.csv"))
all(rownames(LL_hrgm_abun) == rownames(PG_hrgm_abun))

# # rarefy?
# rcts_LL <- fread("data/mgx_abundances/lifelines_unmapped.csv")
# rcts_PG <- fread("data/mgx_abundances/troci_unmapped.csv")
# 
# setkey(rcts_LL, "sample")
# setkey(rcts_PG, "sample")
# rcts_LL <- rcts_LL[colnames(LL_hrgm_abun)]
# rcts_PG <- rcts_PG[colnames(PG_hrgm_abun)]
# 
# rcts_LL[, read_mapped := round(Reads_pe * (100 - unmapped) / 100)]
# rcts_PG[, read_mapped := round(Reads_pe * (100 - unmapped) / 100)]
# 
# LL_hrgm_reads <- round(t(t(LL_hrgm_abun) * rcts_LL$read_mapped))
# PG_hrgm_reads <- round(t(t(PG_hrgm_abun) * rcts_PG$read_mapped))
# 
# HQmods <- Metadata[`Completeness (%)`>= completeness_cuttoff & `Contamination (%)` <= contamination_cutoff, `HRGM name`]
# LL_hrgm_reads <- LL_hrgm_reads[rownames(LL_hrgm_reads) %in% HQmods,]
# PG_hrgm_reads <- PG_hrgm_reads[rownames(PG_hrgm_reads) %in% HQmods,]
# 
# # LL_hrgm_reads <- t(as.matrix(rrarefy(t(LL_hrgm_reads), min(colSums(LL_hrgm_reads)))))
# # PG_hrgm_reads <- t(as.matrix(rrarefy(t(PG_hrgm_reads), min(colSums(PG_hrgm_reads)))))
# LL_hrgm_reads <- t(as.matrix(rrarefy(t(LL_hrgm_reads), 5000)))
# PG_hrgm_reads <- t(as.matrix(rrarefy(t(PG_hrgm_reads), 5000)))
# 
# LL_hrgm_abun <- t(t(LL_hrgm_reads)/colSums(LL_hrgm_reads))
# PG_hrgm_abun <- t(t(PG_hrgm_reads)/colSums(PG_hrgm_reads))

# combine abundance data
hrgm_abun <- cbind(LL_hrgm_abun, PG_hrgm_abun)

rel_models <- intersect(names(which(apply(hrgm_abun,1, function(x) any(x > 0)))),
                        Metadata[`Completeness (%)`>= completeness_cuttoff & `Contamination (%)` <= contamination_cutoff, `HRGM name`])

hrgm_abun <- hrgm_abun[rel_models,]

# Calculate realtive abundance table
relabun <- data.table(as.table(hrgm_abun))
setnames(relabun, c("model","sample","prop"))
relabun <- relabun[model %in% rel_models]
relabun[, prop := prop/sum(prop), by = sample]
relabun <- relabun[!is.na(prop)]

# load models and predict auxotrophies
models <- fetch_model_collection("/mnt/nuuk/2022/HRGM/models/", # /mnt/nuuk/2022/DZHK_MGX/models/
                                 IDs = rel_models)
source("Scripts/predict_auxos.R")

# get Meta info
LL_metainfo <- fread("data/meta/summary_longitudinalMGX_Lifelines.csv")
LL_metainfo <- LL_metainfo[,.(subject, gender, mg.T1 = mg.bam.T1, mg.T2 = mg.bam.T2,
                              age.T1, age.T2)]
LL_metainfo[, cohort := "Lifelines"]
PG_metainfo <- fread("data/mgx_abundances/atlas_samples.csv")
PG_metainfo[, subject := str_match(Full_Name, "-038_\\s*(.*?)\\s*_F(1|2)-L")[,2]]
PG_metainfo[, timepoint := str_match(Full_Name, "_\\s*(F(1|2)?)\\s*-L0")[,2]]
tmp <- fread("data/meta/1.metadata2.tsv")
tmp[, subject := substr(SampleID, 1, 8)]
tmp[, gender := ifelse(sex == 2,"female","male")]
tmp2 <- fread("data/meta/1.metadata1.tsv")
tmp2[, gebdat := gsub("/","-",gebdat)]
tmp2[, tmp := gsub("^..-..-","", gebdat)]
tmp2[, gebdat := paste0(gsub("-..$","",gebdat),"-19",tmp)]
tmp2[, age.T1 := round(-as.numeric(difftime(readr::parse_date(gebdat, "%d-%m-%Y"), "2014-03-1")) / 365.25)]
tmp2[, tmp := NULL]
tmp2[, subject := substr(SampleID, 1, 8)]
tmp3 <- list()
for(iper in PG_metainfo[, .N, by = subject][N==2, subject]) {
  t_age <- tmp2[]
  
  tmp3[[iper]] <- data.table(subject = iper,
                             gender = tmp[subject == iper, gender],
                             mg.T1 = PG_metainfo[subject == iper & timepoint == "F1", atlas_name],
                             mg.T2 = PG_metainfo[subject == iper & timepoint == "F2", atlas_name],
                             age.T1 = tmp2[subject == iper, age.T1],
                             BMI.T1 = tmp2[subject == iper, BMI])
}
tmp3 <- rbindlist(tmp3)
tmp3[, age.T2 := age.T1 + 3]
PG_metainfo <- copy(tmp3)
PG_metainfo[, cohort := "Troci"]
rm(tmp,tmp2, tmp3)

metainfo <- rbind(LL_metainfo, PG_metainfo, fill = TRUE)
metainfo[,quantile(age.T1, na.rm = TRUE, probs = c(0.25,0.5,0.75)), by = cohort]
metainfo[,quantile(BMI.T1, na.rm = TRUE, probs = c(0.25,0.5,0.75)), by = cohort]
metainfo[cohort == "Troci", table(gender)]
metainfo[cohort == "Lifelines", table(gender)]

# calculate Bray-Curtis distance
relabunMat <- dcast(relabun, model ~ sample, value.var = "prop")
tmpn <- relabunMat$model
relabunMat <- as.matrix(relabunMat[,-1])
rownames(relabunMat) <- tmpn
BCdist <- vegdist(t(relabunMat), method="bray")
BCdist <- as.matrix(BCdist)

# add BC distances between T1 and T2 to meta table
metainfo[, BCdist := BCdist[mg.T1, mg.T2], by = subject]

# UniFrac distances
hrgm_tree <- read.tree("data/hrgm_bac.rooted.tree")
relbac <- intersect(rownames(hrgm_abun),
                    Metadata[!grepl("d__Archaea", `GTDB Taxonomy`), `HRGM name`])
tmp_abunmat <- hrgm_abun[relbac,]
tmp_abunmat <- t(t(tmp_abunmat)/colSums(tmp_abunmat))

rmhrgm <- hrgm_tree$tip.label[!(hrgm_tree$tip.label %in% relbac)]
hrgm_tree <- drop.tip(hrgm_tree, rmhrgm)
tmp_abunmat <- tmp_abunmat[hrgm_tree$tip.label,]

for(i in 1:nrow(metainfo)) {
  cat("\r",i, "/",nrow(metainfo))
  aspl <- metainfo[i, mg.T1]
  bspl <- metainfo[i, mg.T2]
  
  metainfo$UF.unweighted[i] <- unweighted_unifrac(tmp_abunmat[,aspl], tmp_abunmat[,bspl], hrgm_tree)
  metainfo$UF.weighted[i] <- weighted_unifrac(tmp_abunmat[,aspl], tmp_abunmat[,bspl], hrgm_tree)
  metainfo$UF.weighted_normalized[i] <- weighted_normalized_unifrac(tmp_abunmat[,aspl], tmp_abunmat[,bspl], hrgm_tree)
  metainfo$UF.variance_adjusted[i] <- variance_adjusted_unifrac(tmp_abunmat[,aspl], tmp_abunmat[,bspl], hrgm_tree)
  metainfo$UF.information[i] <- information_unifrac(tmp_abunmat[,aspl], tmp_abunmat[,bspl], hrgm_tree)
  metainfo$UF.generalized[i] <- generalized_unifrac(tmp_abunmat[,aspl], tmp_abunmat[,bspl], hrgm_tree)
}

# calc auxo frequencies
Nr_auxos <- rowSums(Auxotrophy[,1:20] == 0)
names(Nr_auxos) <- Auxotrophy$Genomes
auxmat <- as.matrix(Auxotrophy[,1:20])
auxmat <- abs(auxmat-1) # now 1 is auxotrophy and 0 prototrophy
rownames(auxmat) <- Auxotrophy$Genomes

avgAuxNr <- Nr_auxos %*% relabunMat
avgAuxNr <- avgAuxNr[1,]

metainfo[, avgAuxNr.T1 := avgAuxNr[mg.T1]]
metainfo[, avgAuxNr.T2 := avgAuxNr[mg.T2]]

auxfreq <- t(t(auxmat) %*% relabunMat)

# calc Hamming distance
for(i in 1:nrow(metainfo)) {
  # T1
  relorgs <- names(which(relabunMat[,metainfo[i, mg.T1]] > 0))
  tmpcomb <- combn(relorgs,2)
  tmp_hamming <- mean(apply(tmpcomb, 2, function(x) {
    dist(auxmat[x,], method = "manhattan")[1]
  }))
  metainfo[i, avgHammingDist.T1 := tmp_hamming]
  
  # T2
  relorgs <- names(which(relabunMat[,metainfo[i, mg.T2]] > 0))
  tmpcomb <- combn(relorgs,2)
  tmp_hamming <- mean(apply(tmpcomb, 2, function(x) {
    dist(auxmat[x,], method = "manhattan")[1]
  }))
  metainfo[i, avgHammingDist.T2 := tmp_hamming]
}


# stability tests
stabi_stat <- list()
k <- 1
for(imetric in colnames(auxfreq)) {
  for(imethod in c("UF.unweighted","UF.weighted","UF.weighted_normalized",
                   "UF.variance_adjusted","UF.variance_adjusted",
                   "UF.information","UF.generalized","BCdist")) {
    
    y.troci <- 1- metainfo[cohort == "Troci"][[imethod]]
    x.troci <- auxfreq[metainfo[cohort == "Troci", mg.T1], imetric]
    y.LL <- 1- metainfo[cohort == "Lifelines"][[imethod]]
    x.LL <- auxfreq[metainfo[cohort == "Lifelines", mg.T1], imetric]
    
    
    sptest.troci <- cor.test(y.troci, x.troci, method = "spearman")
    sptest.LL <- cor.test(y.LL, x.LL, method = "spearman")
    
    
    stabi_stat[[k]] <- data.table(aa = imetric,
                                  method = imethod,
                                  cohort = c("Troci","Lifelines"),
                                  pval = c(sptest.troci$p.value,
                                           sptest.LL$p.value),
                                  rho = c(sptest.troci$estimate,
                                          sptest.LL$estimate))
    k <- k + 1
  }
  
  
}
stabi_stat <- rbindlist(stabi_stat)
stabi_stat[aa != "Gly", padj := p.adjust(pval, method = "fdr"), by = .(cohort,method)]
stabi_stat <- stabi_stat[!is.na(padj)]

cor.test(metainfo[cohort == "Lifelines", 1- UF.unweighted],
         metainfo[cohort == "Lifelines", avgAuxNr.T1], method = "spearman")
cor.test(metainfo[cohort == "Troci", 1- UF.unweighted],
         metainfo[cohort == "Troci", avgAuxNr.T1], method = "spearman")
cor.test(metainfo[cohort == "Lifelines", 1- UF.unweighted],
         metainfo[cohort == "Lifelines", avgHammingDist.T1], method = "spearman")
cor.test(metainfo[cohort == "Troci", 1- UF.unweighted],
         metainfo[cohort == "Troci", avgHammingDist.T1], method = "spearman")

# calc alpha diversity
adiv <- diversity(t(relabunMat), index = "shannon")
metainfo[, D.Shannon.T1 := adiv[mg.T1]]
metainfo[, D.Shannon.T2 := adiv[mg.T2]]

#––––––––––#
# Plotting #
#––––––––––#
metainfo[, cohort2 := ifelse(cohort == "Lifelines", "Chen *et al.* (2021)", "Troci *et al.* (2022)")]


# stability vs. abun.-weighted Average nr of auxotrophs per genotype
p_stabi_avgAux <- ggplot(metainfo, aes(avgAuxNr.T1, 1 - UF.unweighted, col = cohort2, fill = cohort2)) +
  geom_point(shape = 1) +
  geom_smooth(method = "lm") +
  labs(y = "Microbiome stability\n(1 – UniFrac)",
       x = "Abundance-weighted average of auxotrophies",
       fill = "Cohort", col = "Cohort") +
  stat_cor(method = "spearman", cor.coef.name = "rho",
           label.x = 1, label.y = c(0.87, 0.83)) +
  scale_fill_manual(values = c("#785ef0","#fe6100")) +
  scale_color_manual(values = c("#785ef0","#fe6100")) +
  theme_bw() +
  theme(panel.grid.minor = element_blank(),
        legend.text = element_markdown(),
        legend.position = "bottom")  +
  guides(col = guide_legend(title = "Cohort", override.aes = aes(label = "")))

# stability vs. avg Hamming distance
p_stabi_avgHam <- ggplot(metainfo, aes(avgHammingDist.T1, 1 - UF.unweighted, col = cohort2, fill = cohort2)) +
  geom_point(shape = 1) +
  geom_smooth(method = "lm") +
  labs(y = "Microbiome stability\n(1 – UniFrac)",
       x = "Average Hamming Distance",
       fill = "Cohort", col = "Cohort") +
  stat_cor(method = "spearman", cor.coef.name = "rho",
           label.x = 2, label.y = c(0.87, 0.83)) +
  scale_fill_manual(values = c("#785ef0","#fe6100")) +
  scale_color_manual(values = c("#785ef0","#fe6100")) +
  theme_bw() +
  theme(panel.grid.minor = element_blank(),
        legend.text = element_markdown(),
        legend.position = "none")

# stability vs rel. abundance of individual amino acids
stabi_stat[, cohort2 := ifelse(cohort == "Lifelines", "Chen *et al.*<br>(2021)", "Troci *et al.*<br>(2022)")]
stabi_stat[padj < 0.1, plab := "< 0.1"]
stabi_stat[padj < 0.05, plab := "< 0.05"]


p_stabi_heat <- ggplot(stabi_stat[method == "UF.unweighted"],
                       aes(cohort2, aa, col = rho, fill = rho)) +
  geom_tile() +
  geom_point(aes(shape = plab), color = "black") +
  scale_fill_gradient2(high = "#ca0020", mid = "#f7f7f7", low = "#0571b0") +
  scale_color_gradient2(high = "#ca0020", mid = "#f7f7f7", low = "#0571b0") +
  scale_x_discrete(expand = c(0,0)) +
  scale_y_discrete(expand = c(0,0), limits=rev) +
  scale_shape_manual(values = c(19,1), na.translate = F) +
  labs(x = "Cohort", y = "Auxotrophy",
       col = "Spearman's &rho;", fill = "Spearman's &rho;",
       shape = "p<sub>adj</sub>") +
  theme_bw() +
  theme(panel.grid.minor = element_blank(),
        axis.text.x = element_markdown(),
        legend.title = element_markdown())

p_stabi <- egg::ggarrange(p_stabi_avgAux, p_stabi_avgHam, p_stabi_heat,
                          nrow = 1, widths = c(0.4,0.4,0.2),
                          labels = c("A","B","C"),
                          label.args = list(gp = grid::gpar(fontface = "bold")),
                          draw = FALSE)

ggsave("output/plots/stability.pdf", plot = p_stabi,
       width = 10, height = 3.8, device = cairo_pdf)

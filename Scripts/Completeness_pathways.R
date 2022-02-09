
library(MicrobiomeGS2)
library(stringr)
library(dplyr)
library(ggplot2)
library(data.table)
library(ggpubr)
sybil::SYBIL_SETTINGS("SOLVER","cplexAPI")

Metadata <- fread("/mnt/nuuk/2021/HRGM/REPR_Genomes_metadata.tsv")

relGenomes <- Metadata[`Completeness (%)`>= 85 & `Contamination (%)` <=2 & !grepl("^d__Archaea", `GTDB Taxonomy`), `HRGM name`]


models <- fetch_model_collection("/mnt/nuuk/2021/HRGM/models/", IDs = relGenomes)


View(rxns[[1]])
rxns[[1]]
rxns <- fetch_model_collection("/mnt/nuuk/2021/HRGM/models", file.type = "reactions", IDs = relGenomes)
pwys <- fetch_model_collection("/mnt/nuuk/2021/HRGM/models", file.type = "pathways", IDs = relGenomes)
# get pathway coverage for all amino acids
pwys_cov_Arg <- get_pathway_coverage(models,rxns,pwys, pathways.of.interest = c("ARGSYNBSUB-PWY","PWY-5154"))
pwys_cov_Asn <- get_pathway_coverage(models,rxns,pwys, pathways.of.interest = c("ASPARAGINE-BIOSYNTHESIS","ASPARAGINESYN-PWY"))
pwys_cov_Chor <- get_pathway_coverage(models,rxns,pwys, pathways.of.interest = c("PWY-6163"))
pwys_cov_Cys <- get_pathway_coverage(models,rxns,pwys, pathways.of.interest = c("PWY-7870","CYSTSYN-PWY","HOMOCYSDEGR-PWY"))
pwys_cov_Gln <- get_pathway_coverage(models,rxns,pwys, pathways.of.interest = c("GLNSYN-PWY"))
pwys_cov_Gly <- get_pathway_coverage(models,rxns,pwys, pathways.of.interest = c("GLYSYN-ALA-PWY","GLYSYN-PWY","GLYSYN-THR-PWY"))
pwys_cov_His <- get_pathway_coverage(models,rxns,pwys, pathways.of.interest = c("HISTSYN-PWY"))
pwys_cov_Ile <- get_pathway_coverage(models,rxns,pwys, pathways.of.interest = c("ILEUSYN-PWY","PWY-5101","PWY-5103","PWY-5104","PWY-5108"))
pwys_cov_Leu <- get_pathway_coverage(models,rxns,pwys, pathways.of.interest = c("LEUSYN-PWY"))
pwys_cov_Lys <- get_pathway_coverage(models,rxns,pwys, pathways.of.interest = c("PWY-2941", "PWY-2942", "PWY-3081","PWY-5097","DAPLYSINESYN-PWY"))
pwys_cov_Met <- get_pathway_coverage(models,rxns,pwys, pathways.of.interest = c("PWY-7977","HOMOSER-METSYN-PWY"))
pwys_cov_Phe <- get_pathway_coverage(models,rxns,pwys, pathways.of.interest = c("PHESYN","PWY-3462"))
pwys_cov_Pro <- get_pathway_coverage(models,rxns,pwys, pathways.of.interest = c("PROSYN-PWY","PWY-4981"))
pwys_cov_Ser <- get_pathway_coverage(models,rxns,pwys, pathways.of.interest = c("SERSYN-PWY"))
pwys_cov_Thr <- get_pathway_coverage(models,rxns,pwys, pathways.of.interest = c("HOMOSER-THRESYN-PWY","THRBIOSYN-GS"))
pwys_cov_Trp <- get_pathway_coverage(models,rxns,pwys, pathways.of.interest = c("TRPSYN-PWY","TRPSYN-PWY2"))
pwys_cov_Tyr <- get_pathway_coverage(models,rxns,pwys, pathways.of.interest = c("PWY-3461","PWY-6120","TYRSYN"))
pwys_cov_Val <- get_pathway_coverage(models,rxns,pwys, pathways.of.interest = c("VALSYN-PWY"))


relrxns <- unique(pwys_cov_Chor$rxn.metacyc)

compl <- list()
k <- 1
for (rxni in relrxns) {
  rxnperc <- nrow(pwys_cov_Chor[pwys_cov_Chor$rxn.metacyc == rxni & pwys_cov_Chor$prediction == "FALSE"]) / nrow(pwys_cov_Chor[pwys_cov_Chor$rxn.metacyc == rxni])
  rxns <- data.frame(rxnperc)
  rxns$rxn <- rxni
  compl[[k]] <- rxns
  k <- k+1
}
compl


View(pwys_cov_Chor)







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

model.auxo <- lapply(models, FUN = predict_auxotrophies)

Auxotrophie <- data.frame(model.auxo)
head(Auxotrophie) 
summary(Auxotrophie)
str(Auxotrophie)
is.data.frame(Auxotrophie)

#column und rows tauschen
Auxotroph <- t(Auxotrophie)

is.matrix(Auxotroph)
#data frame erzeugen
Auxotrophy <- data.frame(Auxotroph)
is.data.frame(Auxotrophy)
str(Auxotrophy)
Auxotrophy
Genome <- rownames(Auxotrophy)
Auxotrophy$Genomes <- Genome
# ----
Auxotrophy <- as.data.table(Auxotrophy)

Auxo_info <- merge(Auxotrophy, Metadata, by.x = "Genomes",
                   by.y = "HRGM name")

# get pathway coverage for all amino acids
#Arginine
#get arginine auxotrophic genomes
relGenomes <- Auxo_info[Arg == 0, `Genomes`]
models <- fetch_model_collection("/mnt/nuuk/2021/HRGM/models/", IDs = relGenomes)
rxns <- fetch_model_collection("/mnt/nuuk/2021/HRGM/models", file.type = "reactions", IDs = relGenomes)
pwys <- fetch_model_collection("/mnt/nuuk/2021/HRGM/models", file.type = "pathways", IDs = relGenomes)
pwys_cov_Arg <- get_pathway_coverage(models,rxns,pwys, pathways.of.interest = c("ARGSYNBSUB-PWY","PWY-5154"))
head(rxns)
l <- rxns[[1]]
x <- l[,c(1,3)]
test <- distinct(x)
View(x)

relrxns <- unique(pwys_cov_Arg$rxn.metacyc)
relpathway <- unique(pwys_cov_Arg$pathway)
Arg <- list()
k <- 1
for (rxni in relrxns) {
  Argperc <- nrow(pwys_cov_Arg[pwys_cov_Arg$rxn.metacyc == rxni & pwys_cov_Arg$prediction == "FALSE"]) / nrow(pwys_cov_Arg[pwys_cov_Arg$rxn.metacyc == rxni])
  argperc<- data.frame(Argperc)
  colnames(argperc) <- "perc"
  argperc$rxn <- rxni
  argperc$AA <- "Arginine"
  Arg[[k]] <- argperc
  k <- k+1  
}
arg <- rbindlist(Arg)

#Asparagine
relGenomes <- Auxo_info[Asn == 0, `Genomes`]
models <- fetch_model_collection("/mnt/nuuk/2021/HRGM/models/", IDs = relGenomes)
rxns <- fetch_model_collection("/mnt/nuuk/2021/HRGM/models", file.type = "reactions", IDs = relGenomes)
pwys <- fetch_model_collection("/mnt/nuuk/2021/HRGM/models", file.type = "pathways", IDs = relGenomes)
pwys_cov_Asn <- get_pathway_coverage(models,rxns,pwys, pathways.of.interest = c("ASPARAGINE-BIOSYNTHESIS","ASPARAGINESYN-PWY"))
relrxns <- unique(pwys_cov_Asn$rxn.metacyc)
Asn <- list()
k <- 1
for (rxni in relrxns) {
  Asnperc <- nrow(pwys_cov_Asn[pwys_cov_Asn$rxn.metacyc == rxni & pwys_cov_Asn$prediction == "FALSE"]) / nrow(pwys_cov_Asn[pwys_cov_Asn$rxn.metacyc == rxni])
  asnperc <- data.frame(Asnperc)
  colnames(asnperc) <- "perc"
  asnperc$rxn <- rxni
  asnperc$AA <- "Asparagine"
  Asn[[k]] <- asnperc
  k <- k+1
}
asn <- rbindlist(Asn)

#Chorismate
relGenomes <- Auxo_info[Chor == 0, `Genomes`]
models <- fetch_model_collection("/mnt/nuuk/2021/HRGM/models/", IDs = relGenomes)
rxns <- fetch_model_collection("/mnt/nuuk/2021/HRGM/models", file.type = "reactions", IDs = relGenomes)
pwys <- fetch_model_collection("/mnt/nuuk/2021/HRGM/models", file.type = "pathways", IDs = relGenomes)
pwys_cov_Chor <- get_pathway_coverage(models,rxns,pwys, pathways.of.interest = c("PWY-6163"))
relrxns <- unique(pwys_cov_Chor$rxn.metacyc)
Chor <- list()
k <- 1
for (rxni in relrxns) {
  Chorperc <- nrow(pwys_cov_Chor[pwys_cov_Chor$rxn.metacyc == rxni & pwys_cov_Chor$prediction == "FALSE"]) / nrow(pwys_cov_Chor[pwys_cov_Chor$rxn.metacyc == rxni])
  chorperc <- data.frame(Chorperc)
  colnames(chorperc) <- "perc"
  chorperc$rxn <- rxni
  chorperc$AA <- "Chorismate"
  Chor[[k]] <- chorperc
  k <- k+1
}
chor <- rbindlist(Chor)

#Cysteine
relGenomes <- Auxo_info[Cys == 0, `Genomes`]
models <- fetch_model_collection("/mnt/nuuk/2021/HRGM/models/", IDs = relGenomes)
rxns <- fetch_model_collection("/mnt/nuuk/2021/HRGM/models", file.type = "reactions", IDs = relGenomes)
pwys <- fetch_model_collection("/mnt/nuuk/2021/HRGM/models", file.type = "pathways", IDs = relGenomes)
pwys_cov_Cys <- get_pathway_coverage(models,rxns,pwys, pathways.of.interest = c("PWY-7870","CYSTSYN-PWY","HOMOCYSDEGR-PWY"))
relrxns <- unique(pwys_cov_Cys$rxn.metacyc)
Cys <- list()
k <- 1
for (rxni in relrxns) {
  Cysperc <- nrow(pwys_cov_Cys[pwys_cov_Cys$rxn.metacyc == rxni & pwys_cov_Cys$prediction == "FALSE"]) / nrow(pwys_cov_Cys[pwys_cov_Cys$rxn.metacyc == rxni])
  cysperc <- data.frame(Cysperc)
  colnames(cysperc) <- "perc"
  cysperc$rxn <- rxni
  cysperc$AA <- "Cysteine"
  Cys[[k]] <- cysperc
  k <- k+1
}
cys <- rbindlist(Cys)


#Glutamine
relGenomes <- Auxo_info[Gln == 0, `Genomes`]
models <- fetch_model_collection("/mnt/nuuk/2021/HRGM/models/", IDs = relGenomes)
rxns <- fetch_model_collection("/mnt/nuuk/2021/HRGM/models", file.type = "reactions", IDs = relGenomes)
pwys <- fetch_model_collection("/mnt/nuuk/2021/HRGM/models", file.type = "pathways", IDs = relGenomes)
pwys_cov_Gln <- get_pathway_coverage(models,rxns,pwys, pathways.of.interest = c("GLNSYN-PWY"))
relrxns <- unique(pwys_cov_Gln$rxn.metacyc)
Gln <- list()
k <- 1
for (rxni in relrxns) {
  Glnperc <- nrow(pwys_cov_Gln[pwys_cov_Gln$rxn.metacyc == rxni & pwys_cov_Gln$prediction == "FALSE"]) / nrow(pwys_cov_Gln[pwys_cov_Gln$rxn.metacyc == rxni])
  glnperc <- data.frame(Glnperc)
  colnames(glnperc) <- "perc"
  glnperc$rxn <- rxni
  glnperc$AA <- "Glutamine"
  Gln[[k]] <- glnperc
  k <- k+1
}
gln <- rbindlist(Gln)

#Glycine
relGenomes <- Auxo_info[Gly == 0, `Genomes`]
models <- fetch_model_collection("/mnt/nuuk/2021/HRGM/models/", IDs = relGenomes)
rxns <- fetch_model_collection("/mnt/nuuk/2021/HRGM/models", file.type = "reactions", IDs = relGenomes)
pwys <- fetch_model_collection("/mnt/nuuk/2021/HRGM/models", file.type = "pathways", IDs = relGenomes)
pwys_cov_Gly <- get_pathway_coverage(models,rxns,pwys, pathways.of.interest = c("GLYSYN-ALA-PWY","GLYSYN-PWY","GLYSYN-THR-PWY"))
relrxns <- unique(pwys_cov_Gly$rxn.metacyc)
Gly <- list()
k <- 1
for (rxni in relrxns) {
  Glyperc <- nrow(pwys_cov_Gly[pwys_cov_Gly$rxn.metacyc == rxni & pwys_cov_Gly$prediction == "FALSE"]) / nrow(pwys_cov_Gly[pwys_cov_Gly$rxn.metacyc == rxni])
  glyperc <- data.frame(Glyperc)
  colnames(glyperc) <- "perc"
  glyperc$rxn <- rxni
  glyperc$AA <- "Glycine"
  Gly[[k]] <- glyperc
  k <- k+1
}
gly <- rbindlist(Gly)


#Histidine
relGenomes <- Auxo_info[His == 0, `Genomes`]
models <- fetch_model_collection("/mnt/nuuk/2021/HRGM/models/", IDs = relGenomes)
rxns <- fetch_model_collection("/mnt/nuuk/2021/HRGM/models", file.type = "reactions", IDs = relGenomes)
pwys <- fetch_model_collection("/mnt/nuuk/2021/HRGM/models", file.type = "pathways", IDs = relGenomes)
pwys_cov_His <- get_pathway_coverage(models,rxns,pwys, pathways.of.interest = c("HISTSYN-PWY"))
relrxns <- unique(pwys_cov_His$rxn.metacyc)
His <- list()
k <- 1
for (rxni in relrxns) {
  Hisperc <- nrow(pwys_cov_His[pwys_cov_His$rxn.metacyc == rxni & pwys_cov_His$prediction == "FALSE"]) / nrow(pwys_cov_His[pwys_cov_His$rxn.metacyc == rxni])
  hisperc <- data.frame(Hisperc)
  colnames(hisperc) <- "perc"
  hisperc$rxn <- rxni
  hisperc$AA <- "Histidine"
  His[[k]] <- hisperc
  k <- k+1
}
his <- rbindlist(His)

#Isoleucine
relGenomes <- Auxo_info[Ile == 0, `Genomes`]
models <- fetch_model_collection("/mnt/nuuk/2021/HRGM/models/", IDs = relGenomes)
rxns <- fetch_model_collection("/mnt/nuuk/2021/HRGM/models", file.type = "reactions", IDs = relGenomes)
pwys <- fetch_model_collection("/mnt/nuuk/2021/HRGM/models", file.type = "pathways", IDs = relGenomes)
pwys_cov_Ile <- get_pathway_coverage(models,rxns,pwys, pathways.of.interest = c("ILEUSYN-PWY","PWY-5101","PWY-5103","PWY-5104","PWY-5108"))
relrxns <- unique(pwys_cov_Ile$rxn.metacyc)
Ile <- list()
k <- 1
for (rxni in relrxns) {
  Ileperc <- nrow(pwys_cov_Ile[pwys_cov_Ile$rxn.metacyc == rxni & pwys_cov_Ile$prediction == "FALSE"]) / nrow(pwys_cov_Ile[pwys_cov_Ile$rxn.metacyc == rxni])
  ileperc <- data.frame(Ileperc)
  colnames(ileperc) <- "perc"
  ileperc$rxn <- rxni
  ileperc$AA <- "Isoleucine"
  Ile[[k]] <- ileperc
  k <- k+1
}
ile <- rbindlist(Ile)

#Leucine
relGenomes <- Auxo_info[Leu == 0, `Genomes`]
models <- fetch_model_collection("/mnt/nuuk/2021/HRGM/models/", IDs = relGenomes)
rxns <- fetch_model_collection("/mnt/nuuk/2021/HRGM/models", file.type = "reactions", IDs = relGenomes)
pwys <- fetch_model_collection("/mnt/nuuk/2021/HRGM/models", file.type = "pathways", IDs = relGenomes)
pwys_cov_Leu <- get_pathway_coverage(models,rxns,pwys, pathways.of.interest = c("LEUSYN-PWY"))
relrxns <- unique(pwys_cov_Leu$rxn.metacyc)
Leu<- list()
k <- 1
for (rxni in relrxns) {
  Leuperc <- nrow(pwys_cov_Leu[pwys_cov_Leu$rxn.metacyc == rxni & pwys_cov_Leu$prediction == "FALSE"]) / nrow(pwys_cov_Leu[pwys_cov_Leu$rxn.metacyc == rxni])
  leuperc <- data.frame(Leuperc)
  colnames(leuperc) <- "perc"
  leuperc$rxn <- rxni
  leuperc$AA <- "Leucine"
  Leu[[k]] <- leuperc
  k <- k+1
}
leu <- rbindlist(Leu)

#Lysine
relGenomes <- Auxo_info[Lys == 0, `Genomes`]
models <- fetch_model_collection("/mnt/nuuk/2021/HRGM/models/", IDs = relGenomes)
rxns <- fetch_model_collection("/mnt/nuuk/2021/HRGM/models", file.type = "reactions", IDs = relGenomes)
pwys <- fetch_model_collection("/mnt/nuuk/2021/HRGM/models", file.type = "pathways", IDs = relGenomes)
pwys_cov_Lys <- get_pathway_coverage(models,rxns,pwys, pathways.of.interest = c("PWY-2941", "PWY-2942", "PWY-3081","PWY-5097","DAPLYSINESYN-PWY"))
relrxns <- unique(pwys_cov_Lys$rxn.metacyc)
Lys <- list()
k <- 1
for (rxni in relrxns) {
  Lysperc <- nrow(pwys_cov_Lys[pwys_cov_Lys$rxn.metacyc == rxni & pwys_cov_Lys$prediction == "FALSE"]) / nrow(pwys_cov_Lys[pwys_cov_Lys$rxn.metacyc == rxni])
  lysperc <- data.frame(Lysperc)
  colnames(lysperc) <- "perc"
  lysperc$rxn <- rxni
  lysperc$AA <- "Lysine"
  Lys[[k]] <- lysperc
  k <- k+1
}
lys <- rbindlist(Lys)


#Methionine
relGenomes <- Auxo_info[Met == 0, `Genomes`]
models <- fetch_model_collection("/mnt/nuuk/2021/HRGM/models/", IDs = relGenomes)
rxns <- fetch_model_collection("/mnt/nuuk/2021/HRGM/models", file.type = "reactions", IDs = relGenomes)
pwys <- fetch_model_collection("/mnt/nuuk/2021/HRGM/models", file.type = "pathways", IDs = relGenomes)
pwys_cov_Met <- get_pathway_coverage(models,rxns,pwys, pathways.of.interest = c("PWY-7977","HOMOSER-METSYN-PWY"))
relrxns <- unique(pwys_cov_Met$rxn.metacyc)
Met <- list()
k <- 1
for (rxni in relrxns) {
  Metperc <- nrow(pwys_cov_Met[pwys_cov_Met$rxn.metacyc == rxni & pwys_cov_Met$prediction == "FALSE"]) / nrow(pwys_cov_Met[pwys_cov_Met$rxn.metacyc == rxni])
  metperc <- data.frame(Metperc)
  colnames(metperc) <- "perc"
  metperc$rxn <- rxni
  metperc$AA <- "Methionine"
  Met[[k]] <- metperc
  k <- k+1
}
met <- rbindlist(Met)

#Phenylalanine
relGenomes <- Auxo_info[Phe == 0, `Genomes`]
models <- fetch_model_collection("/mnt/nuuk/2021/HRGM/models/", IDs = relGenomes)
rxns <- fetch_model_collection("/mnt/nuuk/2021/HRGM/models", file.type = "reactions", IDs = relGenomes)
pwys <- fetch_model_collection("/mnt/nuuk/2021/HRGM/models", file.type = "pathways", IDs = relGenomes)
pwys_cov_Phe <- get_pathway_coverage(models,rxns,pwys, pathways.of.interest = c("PHESYN","PWY-3462"))
relrxns <- unique(pwys_cov_Phe$rxn.metacyc)
Phe <- list()
k <- 1
for (rxni in relrxns) {
  Pheperc <- nrow(pwys_cov_Phe[pwys_cov_Phe$rxn.metacyc == rxni & pwys_cov_Phe$prediction == "FALSE"]) / nrow(pwys_cov_Phe[pwys_cov_Phe$rxn.metacyc == rxni])
  pheperc <- data.frame(Pheperc)
  colnames(pheperc) <- "perc"
  pheperc$rxn <- rxni
  pheperc$AA <- "Phenylalanine"
  Phe[[k]] <- pheperc
  k <- k+1
}
phe <- rbindlist(Phe)

#Proline
relGenomes <- Auxo_info[Pro == 0, `Genomes`]
models <- fetch_model_collection("/mnt/nuuk/2021/HRGM/models/", IDs = relGenomes)
rxns <- fetch_model_collection("/mnt/nuuk/2021/HRGM/models", file.type = "reactions", IDs = relGenomes)
pwys <- fetch_model_collection("/mnt/nuuk/2021/HRGM/models", file.type = "pathways", IDs = relGenomes)
pwys_cov_Pro <- get_pathway_coverage(models,rxns,pwys, pathways.of.interest = c("PROSYN-PWY","PWY-4981"))
relrxns <- unique(pwys_cov_Pro$rxn.metacyc)
Pro <- list()
k <- 1
for (rxni in relrxns) {
  Properc <- nrow(pwys_cov_Pro[pwys_cov_Pro$rxn.metacyc == rxni & pwys_cov_Pro$prediction == "FALSE"]) / nrow(pwys_cov_Pro[pwys_cov_Pro$rxn.metacyc == rxni])
  properc <- data.frame(Properc)
  colnames(properc) <- "perc"
  properc$rxn <- rxni
  properc$AA <- "Proline"
  Pro[[k]] <- properc
  k <- k+1
}
pro <- rbindlist(Pro)

#Serine
relGenomes <- Auxo_info[Ser == 0, `Genomes`]
models <- fetch_model_collection("/mnt/nuuk/2021/HRGM/models/", IDs = relGenomes)
rxns <- fetch_model_collection("/mnt/nuuk/2021/HRGM/models", file.type = "reactions", IDs = relGenomes)
pwys <- fetch_model_collection("/mnt/nuuk/2021/HRGM/models", file.type = "pathways", IDs = relGenomes)
pwys_cov_Ser <- get_pathway_coverage(models,rxns,pwys, pathways.of.interest = c("SERSYN-PWY"))
relrxns <- unique(pwys_cov_Ser$rxn.metacyc)
Ser <- list()
k <- 1
for (rxni in relrxns) {
  Serperc <- nrow(pwys_cov_Ser[pwys_cov_Ser$rxn.metacyc == rxni & pwys_cov_Ser$prediction == "FALSE"]) / nrow(pwys_cov_Ser[pwys_cov_Ser$rxn.metacyc == rxni])
  serperc <- data.frame(Serperc)
  colnames(serperc) <- "perc"
  serperc$rxn <- rxni
  serperc$AA <- "Serine"
  Ser[[k]] <- serperc
  k <- k+1
}
ser <- rbindlist(Ser)

#Threonine
relGenomes <- Auxo_info[Thr == 0, `Genomes`]
models <- fetch_model_collection("/mnt/nuuk/2021/HRGM/models/", IDs = relGenomes)
rxns <- fetch_model_collection("/mnt/nuuk/2021/HRGM/models", file.type = "reactions", IDs = relGenomes)
pwys <- fetch_model_collection("/mnt/nuuk/2021/HRGM/models", file.type = "pathways", IDs = relGenomes)
pwys_cov_Thr <- get_pathway_coverage(models,rxns,pwys, pathways.of.interest = c("HOMOSER-THRESYN-PWY","THRBIOSYN-GS"))
relrxns <- unique(pwys_cov_Thr$rxn.metacyc)
Thr <- list()
k <- 1
for (rxni in relrxns) {
  Thrperc <- nrow(pwys_cov_Thr[pwys_cov_Thr$rxn.metacyc == rxni & pwys_cov_Thr$prediction == "FALSE"]) / nrow(pwys_cov_Thr[pwys_cov_Thr$rxn.metacyc == rxni])
  thrperc <- data.frame(Thrperc)
  colnames(thrperc) <- "perc"
  thrperc$rxn <- rxni
  thrperc$AA <- "Threonine"
  Thr[[k]] <- thrperc
  k <- k+1
}
thr <- rbindlist(Thr)

#Tryptophan
relGenomes <- Auxo_info[Trp == 0, `Genomes`]
models <- fetch_model_collection("/mnt/nuuk/2021/HRGM/models/", IDs = relGenomes)
rxns <- fetch_model_collection("/mnt/nuuk/2021/HRGM/models", file.type = "reactions", IDs = relGenomes)
pwys <- fetch_model_collection("/mnt/nuuk/2021/HRGM/models", file.type = "pathways", IDs = relGenomes)
pwys_cov_Trp <- get_pathway_coverage(models,rxns,pwys, pathways.of.interest = c("TRPSYN-PWY","TRPSYN-PWY2"))
relrxns <- unique(pwys_cov_Trp$rxn.metacyc)
Trp <- list()
k <- 1
for (rxni in relrxns) {
  Trpperc <- nrow(pwys_cov_Trp[pwys_cov_Trp$rxn.metacyc == rxni & pwys_cov_Trp$prediction == "FALSE"]) / nrow(pwys_cov_Trp[pwys_cov_Trp$rxn.metacyc == rxni])
  trpperc <- data.frame(Trpperc)
  colnames(trpperc) <- "perc"
  trpperc$rxn <- rxni
  trpperc$AA <- "Tryptophan"
  Trp[[k]] <- trpperc
  k <- k+1
}
trp <- rbindlist(Trp)


#Tyrosine
relGenomes <- Auxo_info[Tyr == 0, `Genomes`]
models <- fetch_model_collection("/mnt/nuuk/2021/HRGM/models/", IDs = relGenomes)
rxns <- fetch_model_collection("/mnt/nuuk/2021/HRGM/models", file.type = "reactions", IDs = relGenomes)
pwys <- fetch_model_collection("/mnt/nuuk/2021/HRGM/models", file.type = "pathways", IDs = relGenomes)
pwys_cov_Tyr <- get_pathway_coverage(models,rxns,pwys, pathways.of.interest = c("PWY-3461","PWY-6120","TYRSYN"))
relrxns <- unique(pwys_cov_Tyr$rxn.metacyc)
Tyr <- list()
k <- 1
for (rxni in relrxns) {
  Tyrperc <- nrow(pwys_cov_Tyr[pwys_cov_Tyr$rxn.metacyc == rxni & pwys_cov_Tyr$prediction == "FALSE"]) / nrow(pwys_cov_Tyr[pwys_cov_Tyr$rxn.metacyc == rxni])
  tyrperc <- data.frame(Tyrperc)
  colnames(tyrperc) <- "perc"
  tyrperc$rxn <- rxni
  tyrperc$AA <- "Tyrosine"
  Tyr[[k]] <- tyrperc
  k <- k+1
}
tyr <- rbindlist(Tyr)


#Valine
relGenomes <- Auxo_info[Val == 0, `Genomes`]
models <- fetch_model_collection("/mnt/nuuk/2021/HRGM/models/", IDs = relGenomes)
rxns <- fetch_model_collection("/mnt/nuuk/2021/HRGM/models", file.type = "reactions", IDs = relGenomes)
pwys <- fetch_model_collection("/mnt/nuuk/2021/HRGM/models", file.type = "pathways", IDs = relGenomes)
pwys_cov_Val <- get_pathway_coverage(models,rxns,pwys, pathways.of.interest = c("VALSYN-PWY"))
relrxns <- unique(pwys_cov_Val$rxn.metacyc)
Val <- list()
k <- 1
for (rxni in relrxns) {
  Valperc <- nrow(pwys_cov_Val[pwys_cov_Val$rxn.metacyc == rxni & pwys_cov_Val$prediction == "FALSE"]) / nrow(pwys_cov_Val[pwys_cov_Val$rxn.metacyc == rxni])
  valperc <- data.frame(Valperc)
  colnames(valperc) <- "perc"
  valperc$rxn <- rxni
  valperc$AA <- "Valine"
  Val[[k]] <- valperc
  k <- k+1
}
val <- rbindlist(Val)


#Bind all data together
pathway_completeness <- rbind(arg, asn, chor, cys, gln, gly, his, ile, leu, lys, met, phe, pro, ser, thr, trp, tyr, val)
pathway_completeness$perc <- pathway_completeness$perc * 100
pathway_completeness$perc <- round(pathway_completeness$perc)
pathway_completeness <- data.frame(pathway_completeness)
pathway_new <- merge(pathway_completeness, test, by.x = "rxn", by.y="rxn")
pathway_completeness <- pathway_new[order(pathway_new$AA),]
colnames(pathway_completeness) <- c("Reaction names", "Abundance[%]", "AA", "EC number")
install.packages("gt")
library(gt)

################################   visualization ###############################
t <- gt(pathway_completeness,
   rowname_col = "rowname",
   groupname_col = "AA",
   rownames_to_stub = FALSE,
   auto_align = TRUE,
   id = NULL)

t1 <- t %>%
  tab_header(title = "Completeness of the pathways",
             subtitle = "Abundance of missing enzymes in auxotrophic microbiota")
t1

t1 %>%
  gtsave("pathway.pdf", path = "/home/svenja/Documents")








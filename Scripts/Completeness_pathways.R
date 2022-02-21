################## get pathway coverage for all amino acids ####################
#Arginine
#get arginine auxotrophic genomes
relGenomes <- Auxo_info[Arg == 0, `Genomes`]
models <- fetch_model_collection("/mnt/nuuk/2021/HRGM/models/", IDs = relGenomes)
rxns <- fetch_model_collection("/mnt/nuuk/2021/HRGM/models", file.type = "reactions", IDs = relGenomes)
pwys <- fetch_model_collection("/mnt/nuuk/2021/HRGM/models", file.type = "pathways", IDs = relGenomes)
pwys_cov_Arg <- get_pathway_coverage(models,rxns,pwys, pathways.of.interest = c("ARGSYNBSUB-PWY","PWY-5154"))
View(pwys_cov_Arg)
relrxns <- unique(pwys_cov_Arg$rxn.metacyc)
#relpathway <- unique(pwys_cov_Arg$pathway.name)
Arg <- list()
k <- 1
for (rxni in relrxns) {
  #print(rxni)
  #for (pwayi in relpathway) {
  Argperc <- nrow(pwys_cov_Arg[pwys_cov_Arg$rxn.metacyc == rxni & pwys_cov_Arg$prediction == "FALSE" ]) / nrow(pwys_cov_Arg[pwys_cov_Arg$rxn.metacyc == rxni])
  argperc<- data.frame(Argperc)
  colnames(argperc) <- "perc"
  argperc$rxn <- rxni
  argperc$AA <- "Arginine"
  Arg[[k]] <- argperc
  k <- k+1    
  #}
 
}
arg <- rbindlist(Arg)
arg <- dplyr::mutate(arg, ID = row_number())

remove(arg)
#get the enzyme names for merging with the table later
head(rxns)
l <- rxns[[1]]
is.data.frame(l)
x <- l[,c(1,3)]
test <- distinct(x)
View(x)

#Asparagine
relGenomes <- Auxo_info[Asn == 0, `Genomes`]
models <- fetch_model_collection("/mnt/nuuk/2021/HRGM/models/", IDs = relGenomes)
rxns <- fetch_model_collection("/mnt/nuuk/2021/HRGM/models", file.type = "reactions", IDs = relGenomes)
pwys <- fetch_model_collection("/mnt/nuuk/2021/HRGM/models", file.type = "pathways", IDs = relGenomes)
pwys_cov_Asn <- get_pathway_coverage(models,rxns,pwys, pathways.of.interest = c("ASPARAGINE-BIOSYNTHESIS","ASPARAGINESYN-PWY"))
relrxns <- unique(pwys_cov_Asn$rxn.metacyc)
#relpathway <- unique(pwys_cov_Asn$pathway.name)
Asn <- list()
k <- 1
for (rxni in relrxns) {
  #print(rxni)
 # for (pwayi in relpathway) {
  Asnperc <- nrow(pwys_cov_Asn[pwys_cov_Asn$rxn.metacyc == rxni & pwys_cov_Asn$prediction == "FALSE"]) / nrow(pwys_cov_Asn[pwys_cov_Asn$rxn.metacyc == rxni])
  asnperc <- data.frame(Asnperc)
  colnames(asnperc) <- "perc"
  asnperc$rxn <- rxni
  asnperc$AA <- "Asparagine"
  Asn[[k]] <- asnperc
  k <- k+1
  #}
}
asn <- rbindlist(Asn)
asn <- dplyr::mutate(asn, ID = row_number())
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
chor <- dplyr::mutate(chor, ID = row_number())

#prepare the visualization part
chorperce <- pwys_cov_Chor[ pwys_cov_Chor$model == "HRGM_Genome_0003",]
chorperc <- chorperce[,c(3,5)]
distinct(chorperc)
chor <- merge(chorperc, chor, by.x="rxn.metacyc",by.y = "rxn")
chor_new <- distinct(chor)
chor <- data.frame(chor_new)
chor$perc <- chor$perc *100


##visualization completeness of the chor pathway with the right order
ch <- ggplot(chor, aes(x =`ec`, y =perc))+
  geom_bar(stat="identity") +
  ylab("Abundance of missing enzymes [%]") +
  xlab("Enzymes in the chorismate pathway") +
  theme(axis.line = element_line(size=0.2, colour = "black")) +
  theme(panel.background = element_rect(fill="white", colour= "white")) +
  theme(axis.title.y = element_text(colour = "black", size = 16, face = "bold", margin = margin(0,10,0,0)))+
  theme(axis.title.x = element_text(colour = "black", size = 16, face = "bold", margin = margin(10,0,0,0))) +
  theme(axis.text.x = element_text(size=14, colour = "black", hjust = 1,  angle = 45, margin = margin(10,0,0,0))) +
  theme(axis.text.y = element_text(size = 14, colour = "black")) +
  theme(plot.margin= margin(0.5,0.5,0.5,0.5, "cm")) +
  coord_cartesian(ylim=c(0,100)) +
  scale_x_discrete(limits = c("4.2.1.10","1.1.1.25","2.7.1.71","2.5.1.19","4.2.3.5/4.6.1.4"))

ggsave("output/plots/Completeness_Chor_pathway.pdf", plot = ch,
       width =7, height = 5.5)

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
cys <- filter(cys, perc !=0)
cys <- dplyr::mutate(cys, ID = row_number())
View(pwys_cov_Cys)


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
gln <- dplyr::mutate(gln, ID = row_number())

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
gly <- dplyr::mutate(gly, ID = row_number())


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
his <- dplyr::mutate(his, ID = row_number())

#prepare the visualization part
hisperce <- pwys_cov_His[ pwys_cov_His$model == "HRGM_Genome_0003",]
hisperc <- hisperce[,c(3,5)]
distinct(hisperc)
his <- merge(hisperc, his, by.x="rxn.metacyc",by.y = "rxn")
his_new <- distinct(his)
his <- data.frame(his_new)
his$perc <- his$perc *100

his <- his[-7,]

##visualization completeness of the chor pathway with the right order
hi <- ggplot(his, aes(x =`ec`, y =perc))+
  geom_bar(stat="identity") +
  ylab("Abundance of missing enzymes [%]") +
  xlab("Enzymes in the histidine pathway") +
  theme(axis.line = element_line(size=0.2, colour = "black")) +
  theme(panel.background = element_rect(fill="white", colour= "white")) +
  theme(axis.title.y = element_text(colour = "black", size = 16, face = "bold", margin = margin(0,10,0,0)))+
  theme(axis.title.x = element_text(colour = "black", size = 16, face = "bold", margin = margin(10,0,0,0))) +
  theme(axis.text.x = element_text(size=14, colour = "black", hjust = 1,  angle = 45, margin = margin(10,0,0,0))) +
  theme(axis.text.y = element_text(size = 14, colour = "black")) +
  theme(plot.margin= margin(0.5,0.5,0.5,0.5, "cm")) +
  coord_cartesian(ylim=c(0,100)) +
  scale_x_discrete(limits = c("2.4.2.17", "3.6.1.31", "3.5.4.19","5.3.1.16","4.3.2.10","4.2.1.19","2.6.1.9","3.1.3.15","1.1.1.23","1.1.1.23"))
hi
ggsave("output/plots/Completeness_His_pathway.pdf", plot = hi,
       width =7, height = 5.5)


#Isoleucine
relGenomes <- Auxo_info[Ile == 0, `Genomes`]
models <- fetch_model_collection("/mnt/nuuk/2021/HRGM/models/", IDs = relGenomes)
rxns <- fetch_model_collection("/mnt/nuuk/2021/HRGM/models", file.type = "reactions", IDs = relGenomes)
pwys <- fetch_model_collection("/mnt/nuuk/2021/HRGM/models", file.type = "pathways", IDs = relGenomes)
pwys_cov_Ile <- get_pathway_coverage(models,rxns,pwys, pathways.of.interest = c("ILEUSYN-PWY","PWY-5101","PWY-5103","PWY-5104","PWY-5108"))
View(pwys_cov_Ile)
#analyse only the isoleucine biosynthesis pathway from threonine
pwys_cov_Ile <- get_pathway_coverage(models,rxns,pwys, pathways.of.interest = c("ILEUSYN-PWY"))

relrxns <- unique(pwys_cov_Ile$rxn.metacyc)
Ile <- list()
k <- 1
for (rxni in relrxns) {
  Ileperc <- nrow(pwys_cov_Ile[pwys_cov_Ile$rxn.metacyc == rxni & pwys_cov_Ile$prediction == FALSE]) / nrow(pwys_cov_Ile[pwys_cov_Ile$rxn.metacyc == rxni])
  ileperc <- data.frame(Ileperc)
  colnames(ileperc) <- "perc"
  ileperc$rxn <- rxni
  ileperc$AA <- "Isoleucine"
  Ile[[k]] <- ileperc
  k <- k+1
}
ile <- rbindlist(Ile)
#delete 6th and 7th row because that reactions are spontaenous(zero)
ile <-ile[-c(6,7),]
ile <- dplyr::mutate(ile, ID = row_number())


#prepare the visualization part
ileperce <- pwys_cov_Ile[ pwys_cov_Ile$model == "HRGM_Genome_0167",]
ileperc <- ileperce[,c(3,5)]
distinct(ileperc)
ile1 <- merge(ileperc, ile, by.x="rxn.metacyc",by.y = "rxn")
ile_new1 <- distinct(ile1)
ile1 <- data.frame(ile_new1)
ile1$perc <- ile1$perc *100
ile1[5,2] <- "RXN-15121"
names(ile1) [names(ile1) == "ec"] <- "Enzymes"
##visualization completeness of the chor pathway with the right order
il1 <- ggplot(ile1, aes(x =Enzymes, y =perc, fill =factor(ifelse(Enzymes == "4.3.1.19/4.2.1.16", "not shared", "shared enzyme"))))+
  geom_bar(stat="identity") +
  ylab("Abundance of missing enzymes [%]") +
  xlab("Isoleucine biosynthesis I pathway") +
  theme(axis.line = element_line(size=0.2, colour = "black")) +
  theme(panel.background = element_rect(fill="white", colour= "white")) +
  theme(axis.title.y = element_text(colour = "black", size = 6, face = "bold", margin = margin(0,10,0,0)))+
  theme(axis.title.x = element_text(colour = "black", size = 6, face = "bold", margin = margin(10,0,0,0))) +
  theme(axis.text.x = element_text(size=6, colour = "black", hjust = 1,  angle = 45, margin = margin(10,0,0,0))) +
  theme(axis.text.y = element_text(size = 6, colour = "black")) +
  theme(plot.margin= margin(0.5,0.5,0.5,0.5, "cm")) +
  coord_cartesian(ylim=c(0,100)) +
  scale_x_discrete(limits = c("4.3.1.19/4.2.1.16","2.2.1.6/4.1.3.18","1.1.1.86/1.1.1.89","4.2.1.9","2.6.1.42")) +
  scale_fill_manual(name = "Enzymes", values=c("#fcbba1","#99000d")) +
  theme(legend.text = element_text(size=6))  +
  theme(legend.title = element_text(size=6))
il1
ggsave("output/plots/Completeness_Ile1_pathway.pdf", plot = il1,
       width =7, height = 5.5)


#PWY-5101
pwys_cov_Ile <- get_pathway_coverage(models,rxns,pwys, pathways.of.interest = c("PWY-5101"))

relrxns <- unique(pwys_cov_Ile$rxn.metacyc)
Ile <- list()
k <- 1
for (rxni in relrxns) {
  Ileperc <- nrow(pwys_cov_Ile[pwys_cov_Ile$rxn.metacyc == rxni & pwys_cov_Ile$prediction == FALSE]) / nrow(pwys_cov_Ile[pwys_cov_Ile$rxn.metacyc == rxni])
  ileperc <- data.frame(Ileperc)
  colnames(ileperc) <- "perc"
  ileperc$rxn <- rxni
  ileperc$AA <- "Isoleucine"
  Ile[[k]] <- ileperc
  k <- k+1
}
ile <- rbindlist(Ile)
#delete 6th and 7th row because that reactions are spontaenous(zero)
ile <-ile[-c(6,7),]
ile <- dplyr::mutate(ile, ID = row_number())


#prepare the visualization part
ileperce <- pwys_cov_Ile[pwys_cov_Ile$model == "HRGM_Genome_0167",]
ileperc <- ileperce[,c(3,5)]
distinct(ileperc)
ile2 <- merge(ileperc, ile, by.x="rxn.metacyc",by.y = "rxn")
ile_new2 <- distinct(ile2)
ile2 <- data.frame(ile_new2)
ile2$perc <- ile2$perc *100
names(ile2) [names(ile2) == "ec"] <- "Enzymes"
#give ec number in missing values
ile2[8,2] = "1.1.1.-"
ile2[7,2] = "RXN-7744"
il2 <- ggplot(ile2, aes(x =Enzymes, y =perc, fill =factor(ifelse(Enzymes == "2.2.1.6/4.1.3.18" | Enzymes == "4.2.1.9"| Enzymes =="2.6.1.42", "shared", "not shared"))))+
  geom_bar(stat="identity", position = "stack") +
  ylab("Abundance of missing enzymes [%]") +
  xlab("Isoleucine biosynthesis II pathway") +
  theme(axis.line = element_line(size=0.2, colour = "black")) +
  theme(panel.background = element_rect(fill="white", colour= "white")) +
  theme(axis.title.y = element_text(colour = "black", size = 6, face = "bold", margin = margin(0,10,0,0)))+
  theme(axis.title.x = element_text(colour = "black", size = 6, face = "bold", margin = margin(10,0,0,0))) +
  theme(axis.text.x = element_text(size=6, colour = "black", hjust = 1,  angle = 45, margin = margin(10,0,0,0))) +
  theme(axis.text.y = element_text(size = 6, colour = "black")) +
  theme(plot.margin= margin(0.5,0.5,0.5,0.5, "cm")) +
  coord_cartesian(ylim=c(0,100)) +
  scale_x_discrete(limits = c("2.3.1.182","4.2.1.35","RXN-7744","1.1.1.-","2.2.1.6/4.1.3.18", "1.1.1.383","4.2.1.9", "2.6.1.42")) +
  scale_fill_manual(name = "Enzymes", values=c("#fcbba1","#99000d")) +
  theme(legend.text = element_text(size=6))  +
  theme(legend.title = element_text(size=6))
il2
ggsave("output/plots/Completeness_Ile2_pathway.pdf", plot = il2,
       width =7, height = 5.5)

####PWY-5103
pwys_cov_Ile <- get_pathway_coverage(models,rxns,pwys, pathways.of.interest = c("PWY-5103"))

relrxns <- unique(pwys_cov_Ile$rxn.metacyc)
Ile <- list()
k <- 1
for (rxni in relrxns) {
  Ileperc <- nrow(pwys_cov_Ile[pwys_cov_Ile$rxn.metacyc == rxni & pwys_cov_Ile$prediction == FALSE]) / nrow(pwys_cov_Ile[pwys_cov_Ile$rxn.metacyc == rxni])
  ileperc <- data.frame(Ileperc)
  colnames(ileperc) <- "perc"
  ileperc$rxn <- rxni
  ileperc$AA <- "Isoleucine"
  Ile[[k]] <- ileperc
  k <- k+1
}
ile <- rbindlist(Ile)


#prepare the visualization part
ileperce <- pwys_cov_Ile[ pwys_cov_Ile$model == "HRGM_Genome_0167",]
ileperc <- ileperce[,c(3,5)]
distinct(ileperc)
ile3 <- merge(ileperc, ile, by.x="rxn.metacyc",by.y = "rxn")
ile_new3 <- distinct(ile3)
ile3 <- data.frame(ile_new3)
ile3$perc <- ile3$perc *100
View(pwys_cov_Ile)
names(ile3) [names(ile3) == "ec"] <- "Enzymes"
###fill in missing gaps with ec numbers
ile3[7,2] = "2.6.1.-"
ile3[6,2] = "RXN-7746"
##visualization completeness of the ile pathway with the right order
il3 <- ggplot(ile3, aes(x =Enzymes, y =perc, fill =factor(ifelse(Enzymes == "2.2.1.6/4.1.3.18" | Enzymes == "1.1.1.86/1.1.1.89"| Enzymes == "4.2.1.9" | Enzymes == "2.6.1.42", "shared", "not shared"))))+
  geom_bar(stat="identity") +
  ylab("Abundance of missing enzymes [%]") +
  xlab("Isoleucine biosynthesis III pathway") +
  theme(axis.line = element_line(size=0.2, colour = "black")) +
  theme(panel.background = element_rect(fill="white", colour= "white")) +
  theme(axis.title.y = element_text(colour = "black", size = 6, face = "bold", margin = margin(0,10,0,0)))+
  theme(axis.title.x = element_text(colour = "black", size = 6, face = "bold", margin = margin(10,0,0,0))) +
  theme(axis.text.x = element_text(size=6, colour = "black", hjust = 1,  angle = 45, margin = margin(10,0,0,0))) +
  theme(axis.text.y = element_text(size = 6, colour = "black")) +
  theme(plot.margin= margin(0.5,0.5,0.5,0.5, "cm")) +
  coord_cartesian(ylim=c(0,100)) +
  scale_x_discrete(limits = c("5.4.99.1","2.6.1.-","RXN-7746","2.2.1.6/4.1.3.18","1.1.1.86/1.1.1.89","4.2.1.9","2.6.1.42")) +
  scale_fill_manual(name = "Enzymes", values=c("#fcbba1","#99000d")) +
  theme(legend.text = element_text(size=6))  +
  theme(legend.title = element_text(size=6))
il3
ggsave("output/plots/Completeness_Ile3_pathway.pdf", plot = il3,
       width =7, height = 5.5)

####PWY-5104
pwys_cov_Ile <- get_pathway_coverage(models,rxns,pwys, pathways.of.interest = c("PWY-5104"))

relrxns <- unique(pwys_cov_Ile$rxn.metacyc)
Ile <- list()
k <- 1
for (rxni in relrxns) {
  Ileperc <- nrow(pwys_cov_Ile[pwys_cov_Ile$rxn.metacyc == rxni & pwys_cov_Ile$prediction == FALSE]) / nrow(pwys_cov_Ile[pwys_cov_Ile$rxn.metacyc == rxni])
  ileperc <- data.frame(Ileperc)
  colnames(ileperc) <- "perc"
  ileperc$rxn <- rxni
  ileperc$AA <- "Isoleucine"
  Ile[[k]] <- ileperc
  k <- k+1
}
ile <- rbindlist(Ile)


#prepare the visualization part
ileperce <- pwys_cov_Ile[ pwys_cov_Ile$model == "HRGM_Genome_0167",]
ileperc <- ileperce[,c(3,5)]
distinct(ileperc)
ile4 <- merge(ileperc, ile, by.x="rxn.metacyc",by.y = "rxn")
ile_new4 <- distinct(ile4)
ile4 <- data.frame(ile_new4)
ile4$perc <- ile4$perc *100
View(pwys_cov_Ile)
names(ile4) [names(ile4) == "ec"] <- "Enzymes"
###fill in missing gaps with ec numbers

##visualization completeness of the ile pathway with the right order
il4 <- ggplot(ile4, aes(x =Enzymes, y =perc, fill =factor(ifelse(Enzymes == "2.2.1.6/4.1.3.18"| Enzymes == "4.2.1.9" | Enzymes == "2.6.1.42", "shared", "not shared"))))+
  geom_bar(stat="identity") +
  ylab("Abundance of missing enzymes [%]") +
  xlab("Isoleucine biosynthesis IV pathway") +
  theme(axis.line = element_line(size=0.2, colour = "black")) +
  theme(panel.background = element_rect(fill="white", colour= "white")) +
  theme(axis.title.y = element_text(colour = "black", size = 6, face = "bold", margin = margin(0,10,0,0)))+
  theme(axis.title.x = element_text(colour = "black", size = 6, face = "bold", margin = margin(10,0,0,0))) +
  theme(axis.text.x = element_text(size=6, colour = "black", hjust = 1,  angle = 45, margin = margin(10,0,0,0))) +
  theme(axis.text.y = element_text(size = 6, colour = "black")) +
  theme(plot.margin= margin(0.5,0.5,0.5,0.5, "cm")) +
  coord_cartesian(ylim=c(0,100)) +
  scale_x_discrete(limits = c("6.2.1.17","1.2.7.7","2.2.1.6/4.1.3.18","1.1.1.383","4.2.1.9","2.6.1.42")) +
  scale_fill_manual(name = "Enzymes", values=c("#fcbba1","#99000d")) +
  theme(legend.text = element_text(size=6))  +
  theme(legend.title = element_text(size=6))
il4
ggsave("output/plots/Completeness_Ile_pathway.pdf", plot = il4,
       width =7, height = 5.5)

####PWY-5108
pwys_cov_Ile <- get_pathway_coverage(models,rxns,pwys, pathways.of.interest = c("PWY-5108"))

relrxns <- unique(pwys_cov_Ile$rxn.metacyc)
Ile <- list()
k <- 1
for (rxni in relrxns) {
  Ileperc <- nrow(pwys_cov_Ile[pwys_cov_Ile$rxn.metacyc == rxni & pwys_cov_Ile$prediction == FALSE]) / nrow(pwys_cov_Ile[pwys_cov_Ile$rxn.metacyc == rxni])
  ileperc <- data.frame(Ileperc)
  colnames(ileperc) <- "perc"
  ileperc$rxn <- rxni
  ileperc$AA <- "Isoleucine"
  Ile[[k]] <- ileperc
  k <- k+1
}
ile <- rbindlist(Ile)


#prepare the visualization part
ileperce <- pwys_cov_Ile[ pwys_cov_Ile$model == "HRGM_Genome_0167",]
ileperc <- ileperce[,c(3,5)]
distinct(ileperc)
ile5 <- merge(ileperc, ile, by.x="rxn.metacyc",by.y = "rxn")
ile_new5 <- distinct(ile5)
ile5 <- data.frame(ile_new5)
ile5$perc <- ile5$perc *100
names(ile5) [names(ile5) == "ec"] <- "Enzymes"


##visualization completeness of the ile pathway with the right order
il5 <- ggplot(ile5, aes(x =Enzymes, y =perc, fill = factor(ifelse(Enzymes == "2.6.1.42", "shared", "not shared"))))+
  geom_bar(stat="identity") +
  ylab("Abundance of missing enzymes [%]") +
  xlab("Isoleucine biosynthesis V pathway") +
  theme(axis.line = element_line(size=0.2, colour = "black")) +
  theme(panel.background = element_rect(fill="white", colour= "white")) +
  theme(axis.title.y = element_text(colour = "black", size = 6, face = "bold", margin = margin(0,10,0,0)))+
  theme(axis.title.x = element_text(colour = "black", size = 6, face = "bold", margin = margin(10,0,0,0))) +
  theme(axis.text.x = element_text(size=6, colour = "black", hjust = 1,  angle = 45, margin = margin(10,0,0,0))) +
  theme(axis.text.y = element_text(size = 6, colour = "black")) +
  theme(plot.margin= margin(0.5,0.5,0.5,0.5, "cm")) +
  coord_cartesian(ylim=c(0,100)) +
  scale_x_discrete(limits = c("6.2.1.1","1.2.7.7","2.6.1.42")) +
  scale_fill_manual(name = "Enzymes", values=c("#fcbba1","#99000d")) +
  theme(legend.text = element_text(size=6))  +
  theme(legend.title = element_text(size=6))
il5
ggsave("output/plots/Completeness_Ile_pathway.pdf", plot = il5,
       width =7, height = 5.5)


#Leucine
relGenomes <- Auxo_info[Leu == 0, `Genomes`]
models <- fetch_model_collection("/mnt/nuuk/2021/HRGM/models/", IDs = relGenomes)
rxns <- fetch_model_collection("/mnt/nuuk/2021/HRGM/models", file.type = "reactions", IDs = relGenomes)
pwys <- fetch_model_collection("/mnt/nuuk/2021/HRGM/models", file.type = "pathways", IDs = relGenomes)
pwys_cov_Leu <- get_pathway_coverage(models,rxns,pwys, pathways.of.interest = c("LEUSYN-PWY"))
View(pwys_cov_Leu)
relrxns <- unique(pwys_cov_Leu$rxn.metacyc)
Leu<- list()
k <- 1
for (rxni in relrxns) {
  Leuperc <- nrow(pwys_cov_Leu[pwys_cov_Leu$rxn.metacyc == rxni & pwys_cov_Leu$prediction == "FALSE" & pwys_cov_Leu$spontaneous == "FALSE"]) / nrow(pwys_cov_Leu[pwys_cov_Leu$rxn.metacyc == rxni])
  leuperc <- data.frame(Leuperc)
  colnames(leuperc) <- "perc"
  leuperc$rxn <- rxni
  leuperc$AA <- "Leucine"
  Leu[[k]] <- leuperc
  k <- k+1
}
leu <- rbindlist(Leu)
#delete RXN-7800 because its spontaneous and zero
leu <-leu[-5,]
leu <- dplyr::mutate(leu, ID = row_number())

#prepare the visualization part
leuperce <- pwys_cov_Leu[ pwys_cov_Leu$model == "HRGM_Genome_0167",]
leuperc <- leuperce[,c(3,5)]
distinct(leuperc)
leu <- merge(leuperc, leu, by.x="rxn.metacyc",by.y = "rxn")
leu_new <- distinct(leu)
leu <- data.frame(leu_new)
leu$perc <- leu$perc *100
names(leu) [names(leu) == "ec"] <- "Enzymes"

#delete one reaction meaning the same enzyme because the percentage is the same
leu <- leu[-5,]

##visualization completeness of the leu pathway with the right order
le <- ggplot(leu, aes(x =Enzymes, y =perc, fill =factor(ifelse(Enzymes == "2.6.1.6/2.6.1.42", "shared", "not shared"))))+
  geom_bar(stat="identity") +
  ylab("Abundance of missing enzymes [%]") +
  xlab("Leucine pathway") +
  theme(axis.line = element_line(size=0.2, colour = "black")) +
  theme(panel.background = element_rect(fill="white", colour= "white")) +
  theme(axis.title.y = element_text(colour = "black", size = 6, face = "bold", margin = margin(0,10,0,0)))+
  theme(axis.title.x = element_text(colour = "black", size = 6, face = "bold", margin = margin(10,0,0,0))) +
  theme(axis.text.x = element_text(size=6, colour = "black", hjust = 1,  angle = 45, margin = margin(10,0,0,0))) +
  theme(axis.text.y = element_text(size = 6, colour = "black")) +
  theme(plot.margin= margin(0.5,0.5,0.5,0.5, "cm")) +
  coord_cartesian(ylim=c(0,100)) +
  scale_x_discrete(limits = c("2.3.3.13/4.1.3.12","4.2.1.33","1.1.1.85","2.6.1.6/2.6.1.42")) +
  scale_fill_manual(name = "Enzymes", values=c("#fcbba1","#99000d")) +
  theme(legend.text = element_text(size=6))  +
  theme(legend.title = element_text(size=6))

le

ggsave("output/plots/Completeness_Leu_pathway.pdf", plot = le,
       width =7, height = 5.5)

#Lysine
relGenomes <- Auxo_info[Lys == 0, `Genomes`]
models <- fetch_model_collection("/mnt/nuuk/2021/HRGM/models/", IDs = relGenomes)
rxns <- fetch_model_collection("/mnt/nuuk/2021/HRGM/models", file.type = "reactions", IDs = relGenomes)
pwys <- fetch_model_collection("/mnt/nuuk/2021/HRGM/models", file.type = "pathways", IDs = relGenomes)
pwys_cov_Lys <- get_pathway_coverage(models,rxns,pwys, pathways.of.interest = c("PWY-2941", "PWY-2942", "PWY-3081","PWY-5097","DAPLYSINESYN-PWY"))
View(pwys_cov_Lys)
relrxns <- unique(pwys_cov_Lys$rxn.metacyc)
Lys <- list()
k <- 1
for (rxni in relrxns) {
  Lysperc <- nrow(pwys_cov_Lys[pwys_cov_Lys$rxn.metacyc == rxni & pwys_cov_Lys$prediction == "FALSE"& pwys_cov_Lys$spontaneous == "FALSE"]) / nrow(pwys_cov_Lys[pwys_cov_Lys$rxn.metacyc == rxni])
  lysperc <- data.frame(Lysperc)
  colnames(lysperc) <- "perc"
  lysperc$rxn <- rxni
  lysperc$AA <- "Lysine"
  Lys[[k]] <- lysperc
  k <- k+1
}
lys <- rbindlist(Lys)
lys <- filter(lys, perc != 0)
lys <- dplyr::mutate(lys, ID = row_number())


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
met <- filter(met, perc != 0)
met <- dplyr::mutate(met, ID = row_number())

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
phe <- dplyr::mutate(phe, ID = row_number())

#Proline
relGenomes <- Auxo_info[Pro == 0, `Genomes`]
models <- fetch_model_collection("/mnt/nuuk/2021/HRGM/models/", IDs = relGenomes)
rxns <- fetch_model_collection("/mnt/nuuk/2021/HRGM/models", file.type = "reactions", IDs = relGenomes)
pwys <- fetch_model_collection("/mnt/nuuk/2021/HRGM/models", file.type = "pathways", IDs = relGenomes)
pwys_cov_Pro <- get_pathway_coverage(models,rxns,pwys, pathways.of.interest = c("PROSYN-PWY","PWY-4981"))
View(pwys_cov_Pro)
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
pro <- filter(pro, perc != 0)
pro <- dplyr::mutate(pro, ID = row_number())

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
ser <- dplyr::mutate(ser, ID = row_number())


#prepare the visualization part
serperce <- pwys_cov_Ser[ pwys_cov_Ser$model == "HRGM_Genome_0002",]
serperc <- serperce[,c(3,5)]
distinct(serperc)
ser <- merge(serperc, ser, by.x="rxn.metacyc",by.y = "rxn")
ser_new <- distinct(ser)
ser <- data.frame(ser_new)
ser$perc <- ser$perc *100

##visualization completeness of the leu pathway with the right order
se <- ggplot(ser, aes(x =`ec`, y =perc))+
  geom_bar(stat="identity") +
  ylab("Abundance of missing enzymes [%]") +
  xlab("Enzymes in the serine pathway") +
  theme(axis.line = element_line(size=0.2, colour = "black")) +
  theme(panel.background = element_rect(fill="white", colour= "white")) +
  theme(axis.title.y = element_text(colour = "black", size = 16, face = "bold", margin = margin(0,10,0,0)))+
  theme(axis.title.x = element_text(colour = "black", size = 16, face = "bold", margin = margin(10,0,0,0))) +
  theme(axis.text.x = element_text(size=14, colour = "black", hjust = 1,  angle = 45, margin = margin(10,0,0,0))) +
  theme(axis.text.y = element_text(size = 14, colour = "black")) +
  theme(plot.margin= margin(0.5,0.5,0.5,0.5, "cm")) +
  coord_cartesian(ylim=c(0,100)) +
  scale_x_discrete(limits = c("1.1.1.95","2.6.1.52","3.1.3.3"))
se
ggsave("output/plots/Completeness_Ser_pathway.pdf", plot = se,
       width =7, height = 5.5)

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
thr <- filter(thr, perc != 0)
thr <- dplyr::mutate(thr, ID = row_number())

#Tryptophan
relGenomes <- Auxo_info[Trp == 0, `Genomes`]
models <- fetch_model_collection("/mnt/nuuk/2021/HRGM/models/", IDs = relGenomes)
rxns <- fetch_model_collection("/mnt/nuuk/2021/HRGM/models", file.type = "reactions", IDs = relGenomes)
pwys <- fetch_model_collection("/mnt/nuuk/2021/HRGM/models", file.type = "pathways", IDs = relGenomes)
pwys_cov_Trp <- get_pathway_coverage(models,rxns,pwys, pathways.of.interest = c("TRPSYN-PWY","TRPSYN-PWY2"))

# identify genomes WITH 4.2.1.20
tmpmods <- pwys_cov_Trp[ec == "4.2.1.20" & prediction == TRUE, model]
pwys_cov_Trp[model %in% tmpmods & ec %in% c("4.1.2.8","4.2.1.122"), prediction := TRUE]
rm(tmpmods)

relrxns <- unique(pwys_cov_Trp$rxn.metacyc)

Trp <- list()
k <- 1
for (rxni in relrxns) {
  Trpperc <- (nrow(pwys_cov_Trp[pwys_cov_Trp$rxn.metacyc == rxni & pwys_cov_Trp$prediction == "FALSE"]) / nrow(pwys_cov_Trp[pwys_cov_Trp$rxn.metacyc == rxni]))*100
  trpperc <- data.frame(Trpperc)
  colnames(trpperc) <- "perc"
  trpperc$rxn <- rxni
  trpperc$AA <- "Tryptophan"
  Trp[[k]] <- trpperc
  k <- k+1
}

trp <- rbindlist(Trp)
#prepare the visualization part
trpperce <- pwys_cov_Trp[ pwys_cov_Trp$model == "HRGM_Genome_0094",]
trpperc <- trpperce[,c(3,5)]
distinct(trpperc)
trp <- merge(trpperc, trp, by.x="rxn.metacyc",by.y = "rxn")
trp_new <- distinct(trp)
trp <- data.frame(trp_new)
trp <- dplyr::mutate(trp, ID = row_number())


trp_compl <- filter(pathway1, AA == "Tryptophan")

##visualization completeness of the trp pathway with the right order
tr <- ggplot(trp, aes(x =`ec`, y =perc))+
  geom_bar(stat="identity") +
  ylab("Abundance of missing enzymes [%]") +
  xlab("Enzymes in the tryptophan pathway") +
  theme(axis.line = element_line(size=0.2, colour = "black")) +
  theme(panel.background = element_rect(fill="white", colour= "white")) +
  theme(axis.title.y = element_text(colour = "black", size = 16, face = "bold", margin = margin(0,10,0,0)))+
  theme(axis.title.x = element_text(colour = "black", size = 16, face = "bold", margin = margin(10,0,0,0))) +
  theme(axis.text.x = element_text(size=14, colour = "black", hjust = 1,  angle = 45, margin = margin(10,0,0,0))) +
  theme(axis.text.y = element_text(size = 14, colour = "black")) +
  theme(plot.margin= margin(0.5,0.5,0.5,0.5, "cm")) +
  coord_cartesian(ylim=c(0,100)) +
  scale_x_discrete(limits = c("4.1.3.27","2.4.2.18","5.3.1.24","4.1.1.48","4.2.1.20","4.2.1.122","4.1.2.8"))
tr
ggsave("output/plots/Completeness_Trp_pathway.pdf", plot = tr,
       width =7, height = 5.5)


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
tyr <- dplyr::mutate(tyr, ID = row_number())

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
val <- dplyr::mutate(val, ID = row_number())

#prepare the visualization part
valperce <- pwys_cov_Val[ pwys_cov_Val$model == "HRGM_Genome_0097",]
valperc <- valperce[,c(3,5)]
distinct(valperc)
val <- merge(valperc, val, by.x="rxn.metacyc",by.y = "rxn")
val_new <- distinct(val)
val <- data.frame(val_new)
val$perc <- val$perc *100
names(val)[names(val) == "ec"] <- "Enzymes"

##visualization completeness of the leu pathway with the right order
va <- ggplot(val, aes(x =Enzymes, y =perc))+
  geom_bar(stat="identity", fill = "#99000d") +
  ylab("Abundance of missing enzymes [%]") +
  xlab("Valine pathway") +
  theme(axis.line = element_line(size=0.2, colour = "black")) +
  theme(panel.background = element_rect(fill="white", colour= "white")) +
  theme(axis.title.y = element_text(colour = "black", size = 6, face = "bold", margin = margin(0,10,0,0)))+
  theme(axis.title.x = element_text(colour = "black", size = 6, face = "bold", margin = margin(10,0,0,0))) +
  theme(axis.text.x = element_text(size=6, colour = "black", hjust = 1,  angle = 45, margin = margin(10,0,0,0))) +
  theme(axis.text.y = element_text(size = 6, colour = "black")) +
  theme(plot.margin= margin(0.5,0.5,0.5,0.5, "cm")) +
  coord_cartesian(ylim=c(0,100)) +
  scale_x_discrete(limits = c("2.2.1.6/4.1.3.18","1.1.1.86/1.1.1.89","4.2.1.9", "2.6.1.42")) +
  theme(legend.text = element_text(size=6))  +
  theme(legend.title = element_text(size=6))
va
ggsave("output/plots/Completeness_Val_pathway.pdf", plot = va,
       width =7, height = 5.5)

#Bind all data together
pathway_completeness <- rbind(arg, asn, chor, cys, gln, gly, his, ile, leu, lys, met, phe, pro, ser, thr, trp, tyr, val)
pathway_completeness$perc <- pathway_completeness$perc * 100
pathway_completeness$perc <- round(pathway_completeness$perc)
pathway_completeness <- data.frame(pathway_completeness)
pathway_new <- merge(pathway_completeness, test, by.x = "rxn", by.y="rxn")
pathway_completeness <- pathway_new[order(pathway_new$AA),]
colnames(pathway_completeness) <- c("Reaction names", "Abundance[%]", "AA","ID", "EC number")
pathway <- pathway_completeness[,c(3,4,1,2)]
pathway1 <- data.table(pathway_completeness)

remove(pathway_completeness)

 r <- ggplot(pathway1, aes(AA, ID, fill = `Abundance[%]`)) +
  geom_tile() +
 # geom_text(aes(label = `EC number`)) +
  scale_fill_gradientn(colours = met.brewer("OKeeffe2")) +
  theme(axis.text.x = element_text(angle =90, hjust =0.5, vjust = 0.5)) +
  theme(panel.background = element_blank()) +
  ylab("Missing enzymes") +
  xlab("Amino acid synthesis pathway") +
  theme(axis.text.y = element_blank()) +
  theme(axis.ticks.y = element_blank()) +
  theme(axis.ticks.x = element_blank()) +
  ggtitle("Missing Enzymes in the biosynthesis pathways of amino acidsin auxotrophic microbiota") +
  theme(plot.title = element_text(hjust = 0.5, vjust =0.5, size =12))


ggsave("output/plots/Barplot_Completeness.pdf", plot = r,
       width = 8, height = 6)


library(MetBrewer)
install.packages("patchwork")
library(patchwork)
library(ggpubr)
 a <- ggplot(pathway1[AA == "Arginine"], aes(AA, `EC number`, fill = `Abundance[%]`)) +
  geom_tile() +
  theme(panel.background = element_rect(fill="white", colour = "white")) +
  theme(panel.border = element_blank()) +
  theme(panel.grid = element_blank()) +
  ggtitle("Arginine") +
  theme(plot.title = element_text(hjust = 0.5, size = 8)) +
  theme(axis.title.x = element_blank()) +
  theme(axis.text.x = element_blank())  +
   theme(axis.ticks.x = element_blank()) +
   theme(axis.title.y = element_blank()) +
   theme(axis.text.y = element_text(size = 5)) +
   geom_text(aes(label = ID))

b <- ggplot(pathway1[AA == "Cysteine"], aes(AA, `EC number`, fill = `Abundance[%]`)) +
  geom_tile() +
  theme(panel.background = element_rect(fill="white", colour = "white")) +
  theme(panel.border = element_blank()) +
  theme(panel.grid = element_blank()) +
  ggtitle("Cysteine") +
  theme(plot.title = element_text(hjust = 0.5, size = 8)) +
  theme(axis.title.x = element_blank()) +
  theme(axis.text.x = element_blank()) +
  theme(axis.ticks.x = element_blank()) +
  theme(axis.title.y = element_blank()) +
  theme(axis.text.y = element_text(size = 5)) +
  theme()
b
c <- ggplot(pathway1[AA == "Tryptophan"], aes(AA, `EC number`, fill = `Abundance[%]`)) +
  geom_tile() +
  theme(panel.background = element_rect(fill="white", colour = "white")) +
  theme(panel.border = element_blank()) +
  theme(panel.grid = element_blank()) +
  ggtitle("Tryptophan") +
  theme(plot.title = element_text(hjust = 0.5, size = 8)) +
  theme(axis.title.x = element_blank()) +
  theme(axis.text.x = element_blank())  +
 theme(axis.ticks.x = element_blank()) +
  theme(axis.title.y = element_blank())+
  theme(axis.text.y = element_text(size = 5))

c <- ggplot(pathway1[AA == "Isoleucine"], aes(AA, `EC number`, fill = `Abundance[%]`)) +
  geom_tile() +
  theme(panel.background = element_rect(fill="white", colour = "white")) +
  theme(panel.border = element_blank()) +
  theme(panel.grid = element_blank()) +
  ggtitle("Isoleucine") +
  theme(plot.title = element_text(hjust = 0.5, size = 8)) +
  theme(axis.title.x = element_blank()) +
  theme(axis.text.x = element_blank())  +
  theme(axis.ticks.x = element_blank()) +
  theme(axis.title.y = element_blank()) +
  theme(axis.text.y = element_text(size = 5))

d <- ggplot(pathway1[AA == "Valine"], aes(AA, `EC number`, fill = `Abundance[%]`)) +
  geom_tile() +
  theme(panel.background = element_rect(fill="white", colour = "white")) +
  theme(panel.border = element_blank()) +
  theme(panel.grid = element_blank()) +
  ggtitle("Valine") +
  theme(plot.title = element_text(hjust = 0.5, size = 8)) +
  theme(axis.title.x = element_blank()) +
  theme(axis.text.x = element_blank()) +
  theme(axis.ticks.x = element_blank()) +
  theme(axis.title.y = element_blank()) +
  theme(axis.text.y = element_text(size = 5))

e <- ggplot(pathway1[AA == "Leucine"], aes(AA, `EC number`, fill = `Abundance[%]`)) +
  geom_tile() +
  theme(panel.background = element_rect(fill="white", colour = "white")) +
  theme(panel.border = element_blank()) +
  theme(panel.grid = element_blank()) +
  ggtitle("Leucine") +
  theme(plot.title = element_text(hjust = 0.5, size = 8)) +
  theme(axis.title.x = element_blank()) +
  theme(axis.text.x = element_blank()) +
  theme(axis.ticks.x = element_blank()) +
  theme(axis.title.y = element_blank()) +
  theme(axis.text.y = element_text(size = 5))


e <- ggplot(pathway1[AA == "Methionine"], aes(AA, `EC number`, fill = `Abundance[%]`)) +
  geom_tile() +
  theme(panel.background = element_rect(fill="white", colour = "white")) +
  theme(panel.border = element_blank()) +
  theme(panel.grid = element_blank()) +
  ggtitle("Methionine") +
  theme(plot.title = element_text(hjust = 0.5, size = 8)) +
  theme(axis.title.x = element_blank()) +
  theme(axis.text.x = element_blank()) +
  theme(axis.ticks.x = element_blank()) +
  theme(axis.ticks.x = element_blank()) +
  theme(axis.title.y = element_blank()) +
  theme(axis.text.y = element_text(size = 5))

a+b+c

ggarrange(a, b,c,d,e,
          ncol=5, nrow = 1,
          common.legend = TRUE,
          heights = c(0.8,0.1,1,0.1,0.1),
          legend = "bottom")



################################   visualization ###############################
library(gt)
t <- gt(pathway,
   rowname_col = "rowname",
   groupname_col = "AA",
   rownames_to_stub = FALSE,
   auto_align = TRUE,
   id = NULL)

t1 <- t %>%
  tab_header(title = md("**Completeness of the amino acid biosynthesis pathways**"),
             subtitle = "Abundance of missing enzymes in auxotrophic microbiota")
t2 <- t1 %>%
  tab_options(table.width = pct(70),
              data_row.padding = px(10),
              column_labels.font.weight = "bold",
              row_group.font.weight = "bold",
              table.font.size = px(12))
t2
t2 %>%
  gtsave("pathway.png", expand = 100,  path = "/home/svenja/Documents")

t2 %>%
  gtsave("pathway_12.02.2022.pdf", path = "/home/svenja/Documents")



############## combined plot for BCAA superpathway ##############################
BCAA <- ggarrange(le,il1,il2,il3,va,il4,il5,
                  labels = c("A","B","C", "D","E","F","G","H"),
                  ncol=4, nrow= 2, common.legend = TRUE, legend = "right")
BCAA


ggsave("output/plots/BCAA_completeness.pdf", plot = BCAA,
       width = 9, height = 6)



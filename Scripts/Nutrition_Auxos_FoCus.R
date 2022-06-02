##########  association of auxotrophic bacteria with nutritional data ##########

nutr_info <- fread("/mnt/nuuk/2021/HRGM/2020-031_FOC_nutrintake.csv")

#### analysis for amino acids 
nutr_info_AA <- nutr_info[,c(1:25,70)]
#deleting other non relevant columns
nutr_info_AA$EH <- NULL
nutr_info_AA$EP <- NULL
nutr_info_AA$EPF <- NULL
nutr_info_AA$ENA <- NULL
nutr_info_AA$EEA <- NULL

####   AA/day (E%)
library(dplyr)
nutr_info_AA[,3:20] <- nutr_info_AA[,3:20]* 17
nutr_info_AA <- nutr_info_AA %>% mutate(across(c(3:20), .fns = ~./GJ)*100)

#####  
data <- fread("/mnt/nuuk/2021/HRGM/FOCUS_HRGM_abundancies.csv.gz")
data <- data[focus.call == "BL"]

#correlation analysis
subj <- unique(data$subject)
relAA <- unique(Auxotrophy_2$Compound)
remove(e)
e <- list()
k <- 1

for (subi in subj) {
  print(subi) 
  for (AAi in relAA) {
    x <- data[subject == subi]
    y <- Auxotrophy_2[Compound == AAi]
    z <- merge(x,y, by.x = "model", by.y = "Genomes")
    t <- z[Prototrophy == 0]
    e[[k]] <- t
    k <- k +1
  }
}

nutr_info_auxo <- rbindlist(e) 
View(nutr_info_auxo)


sumfreq_nutr_auxo <- aggregate(nutr_info_auxo$freq, by=list(subject=nutr_info_auxo$subject, AA=nutr_info_auxo$Compound), FUN=sum)
describe(sumfreq_nutr_auxo$subject)
nutri_all_info <- merge(sumfreq_nutr_auxo,nutr_info_AA, by.x="subject", by.y="new_id")
FoCus_data <- fread("/mnt/nuuk/2021/HRGM/FOCUS_meta_info_v2.csv")
nutri_all_info <- merge(nutri_all_info,FoCus_data, by.x="subject", by.y="sample")
nutri_all_info <- data.table(nutri_all_info)

#exclude all people with diabetes
nutri_all_info <- nutri_all_info[diabetes!=1, ]
describe(nutri_all_info$subject)
View(nutri_all_info)


##### ALA
AAI
AA <- unique(nutri_all_info$AA)
rm(l)
ALA <- list()
k <- 1

for(AAi in AA) {
  print(AAi)
    n <- partial_Spearman(x|EALA ~ Gender + age + BMI, data = nutri_all_info[AA == AAi])
    nutri_all_info$Estimate_Part_Spear <- n$TS$TB$ts
    nutri_all_info$pvalue_Part_Spear <- n$TS$TB$pval
    table_corr <- data.table(Pvalue = n$TS$TB$pval,
                             AA = AAi,
                             day_int = "ALA",
                             Estimate = n$TS$TB$ts)
    ALA[[k]] <- table_corr
    k <- k+1
}

tmp_ALA <- rbindlist(ALA)
#### ARGININE
AA <- unique(nutri_all_info$AA)
rm(l)
ARG <- list()
k <- 1

for(AAi in AA) {
  print(AAi)
  n <- partial_Spearman(x|EARG ~ Gender + age + BMI, data = nutri_all_info[AA == AAi])
  nutri_all_info$Estimate_Part_Spear <- n$TS$TB$ts
  nutri_all_info$pvalue_Part_Spear <- n$TS$TB$pval
  table_corr <- data.table(Pvalue = n$TS$TB$pval,
                           AA = AAi,
                           day_int = "Arg",
                           Estimate = n$TS$TB$ts)
  ARG[[k]] <- table_corr
  k <- k+1
}
tmp_ARG <- rbindlist(ARG)
######Asparagine
AA <- unique(nutri_all_info$AA)
ASP <- list()
k <- 1

for(AAi in AA) {
  print(AAi)
  n <- partial_Spearman(x|EASP ~ Gender + age + BMI, data = nutri_all_info[AA == AAi])
  nutri_all_info$Estimate_Part_Spear <- n$TS$TB$ts
  nutri_all_info$pvalue_Part_Spear <- n$TS$TB$pval
  table_corr <- data.table(Pvalue = n$TS$TB$pval,
                           AA = AAi,
                           day_int = "Asp",
                           Estimate = n$TS$TB$ts)
  ASP[[k]] <- table_corr
  k <- k+1
}
tmp_ASP <- rbindlist(ASP)
############Cysteine
AA <- unique(nutri_all_info$AA)
CYS <- list()
k <- 1

for(AAi in AA) {
  print(AAi)
  n <- partial_Spearman(x|ECYS ~ Gender + age + BMI, data = nutri_all_info[AA == AAi])
  nutri_all_info$Estimate_Part_Spear <- n$TS$TB$ts
  nutri_all_info$pvalue_Part_Spear <- n$TS$TB$pval
  table_corr <- data.table(Pvalue = n$TS$TB$pval,
                           AA = AAi,
                           day_int = "Cys",
                           Estimate = n$TS$TB$ts)
  CYS[[k]] <- table_corr
  k <- k+1
}
tmp_CYS <- rbindlist(CYS)
###########Glutamate
AA <- unique(nutri_all_info$AA)
GLU <- list()
k <- 1

for(AAi in AA) {
  print(AAi)
  n <- partial_Spearman(x|EGLU ~ Gender + age + BMI, data = nutri_all_info[AA == AAi])
  nutri_all_info$Estimate_Part_Spear <- n$TS$TB$ts
  nutri_all_info$pvalue_Part_Spear <- n$TS$TB$pval
  table_corr <- data.table(Pvalue = n$TS$TB$pval,
                           AA = AAi,
                           day_int = "Glu",
                           Estimate = n$TS$TB$ts)
  GLU[[k]] <- table_corr
  k <- k+1
}
tmp_GLU <- rbindlist(GLU)
###########Glycin
AA <- unique(nutri_all_info$AA)
GLY<- list()
k <- 1

for(AAi in AA) {
  print(AAi)
  n <- partial_Spearman(x|EGLY ~ Gender + age + BMI, data = nutri_all_info[AA == AAi])
  nutri_all_info$Estimate_Part_Spear <- n$TS$TB$ts
  nutri_all_info$pvalue_Part_Spear <- n$TS$TB$pval
  table_corr <- data.table(Pvalue = n$TS$TB$pval,
                           AA = AAi,
                           day_int = "Gly",
                           Estimate = n$TS$TB$ts)
  GLY[[k]] <- table_corr
  k <- k+1
}
tmp_GLY <- rbindlist(GLY)
###########Histidine
AA <- unique(nutri_all_info$AA)
HIS <- list()
k <- 1

for(AAi in AA) {
  print(AAi)
  n <- partial_Spearman(x|EHIS ~ Gender + age + BMI, data = nutri_all_info[AA == AAi])
  nutri_all_info$Estimate_Part_Spear <- n$TS$TB$ts
  nutri_all_info$pvalue_Part_Spear <- n$TS$TB$pval
  table_corr <- data.table(Pvalue = n$TS$TB$pval,
                           AA = AAi,
                           day_int = "HIS",
                           Estimate = n$TS$TB$ts)
  HIS[[k]] <- table_corr
  k <- k+1
}
tmp_HIS <- rbindlist(HIS)
######Isoleucine
AA <- unique(nutri_all_info$AA)
ILE <- list()
k <- 1

for(AAi in AA) {
  print(AAi)
  n <- partial_Spearman(x|EILE ~ Gender + age + BMI, data = nutri_all_info[AA == AAi])
  nutri_all_info$Estimate_Part_Spear <- n$TS$TB$ts
  nutri_all_info$pvalue_Part_Spear <- n$TS$TB$pval
  table_corr <- data.table(Pvalue = n$TS$TB$pval,
                           AA = AAi,
                           day_int = "Ile",
                           Estimate = n$TS$TB$ts)
  ILE[[k]] <- table_corr
  k <- k+1
}
tmp_ILE <- rbindlist(ILE)
########Leucine
AA <- unique(nutri_all_info$AA)
LEU <- list()
k <- 1

for(AAi in AA) {
  print(AAi)
  n <- partial_Spearman(x|ELEU ~ Gender + age + BMI, data = nutri_all_info[AA == AAi])
  nutri_all_info$Estimate_Part_Spear <- n$TS$TB$ts
  nutri_all_info$pvalue_Part_Spear <- n$TS$TB$pval
  table_corr <- data.table(Pvalue = n$TS$TB$pval,
                           AA = AAi,
                           day_int = "Leu",
                           Estimate = n$TS$TB$ts)
  LEU[[k]] <- table_corr
  k <- k+1
}
tmp_LEU <- rbindlist(LEU)
##########lysine
AA <- unique(nutri_all_info$AA)
LYS <- list()
k <- 1

for(AAi in AA) {
  print(AAi)
  n <- partial_Spearman(x|ELYS ~ Gender + age + BMI, data = nutri_all_info[AA == AAi])
  nutri_all_info$Estimate_Part_Spear <- n$TS$TB$ts
  nutri_all_info$pvalue_Part_Spear <- n$TS$TB$pval
  table_corr <- data.table(Pvalue = n$TS$TB$pval,
                           AA = AAi,
                           day_int = "Lys",
                           Estimate = n$TS$TB$ts)
  LYS[[k]] <- table_corr
  k <- k+1
}
tmp_LYS <- rbindlist(LYS)
######Methionine
AA <- unique(nutri_all_info$AA)
MET <- list()
k <- 1

for(AAi in AA) {
  print(AAi)
  n <- partial_Spearman(x|EMET ~ Gender + age + BMI, data = nutri_all_info[AA == AAi])
  nutri_all_info$Estimate_Part_Spear <- n$TS$TB$ts
  nutri_all_info$pvalue_Part_Spear <- n$TS$TB$pval
  table_corr <- data.table(Pvalue = n$TS$TB$pval,
                           AA = AAi,
                           day_int = "Met",
                           Estimate = n$TS$TB$ts)
  MET[[k]] <- table_corr
  k <- k+1
}
tmp_MET <- rbindlist(MET)
###########Phenylalanine
AA <- unique(nutri_all_info$AA)
PHE <- list()
k <- 1

for(AAi in AA) {
  print(AAi)
  n <- partial_Spearman(x|EPHE ~ Gender + age + BMI, data = nutri_all_info[AA == AAi])
  nutri_all_info$Estimate_Part_Spear <- n$TS$TB$ts
  nutri_all_info$pvalue_Part_Spear <- n$TS$TB$pval
  table_corr <- data.table(Pvalue = n$TS$TB$pval,
                           AA = AAi,
                           day_int = "Phe",
                           Estimate = n$TS$TB$ts)
  PHE[[k]] <- table_corr
  k <- k+1
}
tmp_PHE <- rbindlist(PHE)
#######Proline
AA <- unique(nutri_all_info$AA)
PRO <- list()
k <- 1

for(AAi in AA) {
  print(AAi)
  n <- partial_Spearman(x|EPRO ~ Gender + age + BMI, data = nutri_all_info[AA == AAi])
  nutri_all_info$Estimate_Part_Spear <- n$TS$TB$ts
  nutri_all_info$pvalue_Part_Spear <- n$TS$TB$pval
  table_corr <- data.table(Pvalue = n$TS$TB$pval,
                           AA = AAi,
                           day_int = "Pro",
                           Estimate = n$TS$TB$ts)
  PRO[[k]] <- table_corr
  k <- k+1
}
tmp_PRO <- rbindlist(PRO)
#######Serine
AA <- unique(nutri_all_info$AA)
SER <- list()
k <- 1

for(AAi in AA) {
  print(AAi)
  n <- partial_Spearman(x|ESER ~ Gender + age + BMI, data = nutri_all_info[AA == AAi])
  nutri_all_info$Estimate_Part_Spear <- n$TS$TB$ts
  nutri_all_info$pvalue_Part_Spear <- n$TS$TB$pval
  table_corr <- data.table(Pvalue = n$TS$TB$pval,
                           AA = AAi,
                           day_int = "Ser",
                           Estimate = n$TS$TB$ts)
  SER[[k]] <- table_corr
  k <- k+1
}
tmp_SER <- rbindlist(SER)
#########threonine
AA <- unique(nutri_all_info$AA)
THR <- list()
k <- 1

for(AAi in AA) {
  print(AAi)
  n <- partial_Spearman(x|ETHR ~ Gender + age + BMI, data = nutri_all_info[AA == AAi])
  nutri_all_info$Estimate_Part_Spear <- n$TS$TB$ts
  nutri_all_info$pvalue_Part_Spear <- n$TS$TB$pval
  table_corr <- data.table(Pvalue = n$TS$TB$pval,
                           AA = AAi,
                           day_int = "Thr",
                           Estimate = n$TS$TB$ts)
  THR[[k]] <- table_corr
  k <- k+1
}
tmp_THR <- rbindlist(THR)
#######tryptophan
AA <- unique(nutri_all_info$AA)
TRP <- list()
k <- 1

for(AAi in AA) {
  print(AAi)
  n <- partial_Spearman(x|ETRP ~ Gender + age + BMI, data = nutri_all_info[AA == AAi])
  nutri_all_info$Estimate_Part_Spear <- n$TS$TB$ts
  nutri_all_info$pvalue_Part_Spear <- n$TS$TB$pval
  table_corr <- data.table(Pvalue = n$TS$TB$pval,
                           AA = AAi,
                           day_int = "Trp",
                           Estimate = n$TS$TB$ts)
  TRP[[k]] <- table_corr
  k <- k+1
}
tmp_TRP <- rbindlist(TRP)

##########tyrosine
AA <- unique(nutri_all_info$AA)
TYR <- list()
k <- 1

for(AAi in AA) {
  print(AAi)
  n <- partial_Spearman(x|ETYR ~ Gender + age + BMI, data = nutri_all_info[AA == AAi])
  nutri_all_info$Estimate_Part_Spear <- n$TS$TB$ts
  nutri_all_info$pvalue_Part_Spear <- n$TS$TB$pval
  table_corr <- data.table(Pvalue = n$TS$TB$pval,
                           AA = AAi,
                           day_int = "Tyr",
                           Estimate = n$TS$TB$ts)
  TYR[[k]] <- table_corr
  k <- k+1
}
tmp_TYR <- rbindlist(TYR)
########valine
AA <- unique(nutri_all_info$AA)
VAL <- list()
k <- 1

for(AAi in AA) {
  print(AAi)
  n <- partial_Spearman(x|EVAL ~ Gender + age + BMI, data = nutri_all_info[AA == AAi])
  nutri_all_info$Estimate_Part_Spear <- n$TS$TB$ts
  nutri_all_info$pvalue_Part_Spear <- n$TS$TB$pval
  table_corr <- data.table(Pvalue = n$TS$TB$pval,
                           AA = AAi,
                           day_int = "VAL",
                           Estimate = n$TS$TB$ts)
  VAL[[k]] <- table_corr
  k <- k+1
}
tmp_VAL <- rbindlist(VAL)

###combine all AAA tables in one table
Part_Spear_EAA <- rbind(tmp_ALA, tmp_ASP, tmp_CYS, tmp_GLU, tmp_VAL, tmp_HIS, tmp_ILE, tmp_LEU,
                        tmp_LYS, tmp_MET, tmp_PHE, tmp_PRO, tmp_SER, tmp_THR, tmp_TRP, tmp_TYR, 
                        tmp_ARG, tmp_GLY)

##adjust pvalue
Part_Spear_EAA$padjust = p.adjust(Part_Spear_EAA$Pvalue, method = "fdr")
Part_Spear_EAA[padjust < 0.05, sign.label1 := "P < 0.05"]
Part_Spear_EAA$padjust_2 = p.adjust(Part_Spear_EAA$Pvalue, method = "fdr")
Part_Spear_EAA[padjust_2 < 0.1, sign.label2 := "P < 0.1"]

View(Part_Spear_EAA)

heatmap_AA <- ggplot(Part_Spear_EAA, aes(AA,day_int, fill = Estimate)) +
  geom_tile(linetype = 1.5, colour ="grey", lwd=0.2) +
  geom_point(aes(shape=sign.label2), size=1) +
  scale_fill_gradient2(high = "#ca0020", mid="white", low = "#0571b0") +
  theme_minimal() +
  theme(axis.text.x = element_text(colour="black", size = 10)) +
  theme(axis.text.y = element_text(colour = "black", size=10)) +
  theme(axis.title.x = element_text(colour="black", size = 12, face = "bold")) +
  theme(axis.title.y = element_text(colour = "black", size =12, face = "bold")) +
  scale_shape_manual(values = 8, na.translate = FALSE)+
  xlab("Amino acid auxotrophies") +
  ylab("E(%)") +
  theme(panel.background = element_blank())
heatmap_AA


ggsave("output/plots/heatmap_nutr_AA.pdf", plot = heatmap_AA,
       width = 9, height = 6)


##################### B-VITAMINS ###############################
#### analysis for amino acids 
nutr_info_BVit <- nutr_info[,c("new_id","gramm","VB1","VB12","VB2","VB3","VB3A","VB5","VB6","VB7","VB9","VB9F","VB9G","GJ")]
####   AA/day (E%)
library(dplyr)
nutr_info_BVit <- nutr_info_BVit %>% mutate(across(c(3:13), .fns = ~./GJ)*100)

#####  
data <- fread("/mnt/nuuk/2021/HRGM/FOCUS_HRGM_abundancies.csv.gz")
data <- data[focus.call == "BL"]

#merge
nutri_all_info_BV <- merge(sumfreq_nutr_auxo,nutr_info_BVit, by.x="subject", by.y="new_id")


colnames <- colnames(nutri_all_info_BV)
nrow(nutri_all_info[nutri_all_info$AA=="Gln"])

newcol <- colnames[-c(1:4,16)]

AAI <- unique(newcol)
B_Vit <- unique(nutri_all_info_BV$AA)
rm(l)
l <- list()
k <- 1

for(AAi in AA) {
  print(AAi)
  for(AAIi in AAI){
    table <- nutri_all_info_BV[nutri_all_info_BV$AA == AAi, ]
    res <- cor.test(table$x, table[[AAIi]], method = "spearman", exact = FALSE)
    table_corr_B <- data.table(Pvalue = res$p.value,
                             AA = AAi,
                             day_int = AAIi,
                             Estimate = res$estimate)
    l[[k]] <- table_corr_B
    k <- k+1
  }
}


BVit_auxo <- rbindlist(l)
BVit_auxo[Pvalue< 0.05, sign.label1 := "not adj.P < 0.05"]

heatmap_BVit_auxo <- ggplot(BVit_auxo, aes(AA,day_int, fill = Estimate)) +
  geom_tile(linetype = 1.5, colour ="grey", lwd=0.2) +
  geom_point(aes(shape=sign.label1), size=1) +
  scale_fill_gradient2(high = "#ca0020", mid="white", low = "#0571b0") +
  theme_minimal() +
  theme(axis.text.x = element_text(colour="black", size = 10)) +
  theme(axis.text.y = element_text(colour = "black", size=10)) +
  theme(axis.title.x = element_text(colour="black", size = 12, face = "bold")) +
  theme(axis.title.y = element_text(colour = "black", size =12, face = "bold")) +
  scale_shape_manual(values = 8, na.translate = FALSE)+
  xlab("Amino acid auxotrophies") +
  ylab("E(%)") +
  theme(panel.background = element_blank())

heatmap_BVit_auxo
ggsave("output/plots/heatmap_nutr_BVit.pdf", plot = heatmap_BVit_auxo,
       width = 9, height = 6)













#analyzing the amount of amino acids
boxplot(nutr_info_AA$EALA)


boxplot(nutr_info_AA$EALA)

ggplot(nutr_info_AA[(EALA)], aes(x="ALA", y = EALA)) +
  geom_boxplot()

ARG <- ggplot(nutr_info_AA[nutr_info_AA$EARG], aes(x = "ARG", y = EARG)) +
  geom_boxplot()+
  labs(x="") +
  theme_bw()  +
  #ylim(0,25)+
  labs(y="")
ARG
ASP <- ggplot(nutr_info_AA[nutr_info_AA$EASP], aes(x = "ASP", y = EASP)) +
  geom_boxplot()+
  labs(x="") +
  theme_bw()  +
  #ylim(0,25)+
  labs(y="")

#ggplot(nutr_info[nutr_info$EEA], aes(x = "total_essAA", y = EEA)) +
#geom_boxplot()+
#labs(x="") +
#theme_bw()

CYS <- ggplot(nutr_info_AA[nutr_info_AA$ECYS], aes(x = "CYS", y = ECYS)) +
  geom_boxplot()+
  labs(x="") +
  theme_bw()  +
  #ylim(0,25)+
  labs(y="")

GLU <- ggplot(nutr_info_AA[nutr_info_AA$EGLU], aes(x = "GLU", y = EGLU)) +
  geom_boxplot()+
  labs(x="") +
  theme_bw()  +
  #ylim(0,25)+
  labs(y="")

GLY <- ggplot(nutr_info_AA[nutr_info_AA$EGLY], aes(x = "GLY", y = EGLY)) +
  geom_boxplot()+
  labs(x="") +
  theme_bw()  +
  #ylim(0,25)+
  labs(y="")

HIS <- ggplot(nutr_info_AA[nutr_info_AA$EHIS], aes(x = "HIS", y = EHIS)) +
  geom_boxplot()+
  labs(x="") +
  theme_bw()  +
  #ylim(0,25)+
  labs(y="")

ILE <- ggplot(nutr_info_AA[nutr_info_AA$EILE], aes(x = "ILE", y = EILE)) +
  geom_boxplot()+
  labs(x="") +
  theme_bw()  +
  # ylim(0,25)+
  labs(y="")

LEU <- ggplot(nutr_info_AA[nutr_info_AA$ELEU], aes(x = "LEU", y = ELEU)) +
  geom_boxplot()+
  labs(x="") +
  theme_bw()  +
  #ylim(0,25)+
  labs(y="")

LYS <- ggplot(nutr_info_AA[nutr_info_AA$ELYS], aes(x = "LYS", y = ELYS)) +
  geom_boxplot()+
  labs(x="") +
  theme_bw()  +
  #ylim(0,25)+
  labs(y="")

MET <- ggplot(nutr_info_AA[nutr_info_AA$EMET], aes(x = "MET", y = EMET)) +
  geom_boxplot()+
  labs(x="") +
  theme_bw()  +
  #ylim(0,25)+
  labs(y="")

#ggplot(nutr_info_AA[nutr_info$ENA], aes(x = "total_non_ess", y = ENA)) +
#geom_boxplot()+
#labs(x="") +
# theme_bw()

PHE <- ggplot(nutr_info_AA[nutr_info_AA$EPHE], aes(x = "PHE", y = EPHE)) +
  geom_boxplot()+
  labs(x="") +
  theme_bw()  +
  #ylim(0,25)+
  labs(y="")

boxplot(nutr_info_AA$EMET)


PRO <- ggplot(nutr_info_AA[nutr_info_AA$EPRO], aes(x = "PRO", y = EPRO)) +
  geom_boxplot()+
  labs(x="") +
  theme_bw()  +
  #ylim(0,25)+
  labs(y="")

SER <- ggplot(nutr_info_AA[nutr_info_AA$ESER], aes(x = "SER", y = ESER)) +
  geom_boxplot()+
  labs(x="") +
  theme_bw()  +
  #ylim(0,25)+
  labs(y="")

THR <- ggplot(nutr_info_AA[nutr_info_AA$ETHR], aes(x = "THR", y = ETHR)) +
  geom_boxplot()+
  labs(x="") +
  theme_bw()  +
  #ylim(0,25)+
  labs(y="")

TRP <- ggplot(nutr_info_AA[nutr_info_AA$ETRP], aes(x = "TRP", y = ETRP)) +
  geom_boxplot()+
  labs(x="") +
  theme_bw()
#ylim(0,25)

TYR <- ggplot(nutr_info_AA[nutr_info_AA$ETYR], aes(x = "TYR", y = ETYR)) +
  geom_boxplot()+
  labs(x="") +
  theme_bw()+
  labs(y="")

VAL <- ggplot(nutr_info_AA[nutr_info_AA$EVAL], aes(x = "VAL", y = EVAL)) +
  geom_boxplot()+
  labs(x="") +
  theme_bw()  +
  #ylim(0,25)+
  labs(y="")

part1_AA <- ggarrange(ALA,ARG,ASP,CYS,GLU,GLY,HIS,ILE,LEU,
                      ncol=4,nrow=4)

ggsave("output/plots/AA_uptake_day_1.pdf", plot = part1_AA,
       width = 9, height = 8)

part2_AA <- ggarrange(LYS,MET,PHE,PRO,SER,THR,TRP,TYR,VAL,
                      ncol=4,nrow=3)
ggsave("output/plots/AA_uptake_day_2.pdf", plot = part2_AA,
       width = 9, height = 8)








##########  association of auxotrophic bacteria with nutritional data ##########

nutr_info <- fread("/mnt/nuuk/2021/HRGM/2020-031_FOC_nutrintake.csv")

#### analysis for amino acids 
nutr_info_AA <- nutr_info[,c(1:25,70)]
#deleting other non relevant columns
nutr_info_AA$EH <- NULL
nutr_info_AA$EP <- NULL
nutr_info_AA$EPF <- NULL
remove(nutr_info_AA)
####   AA/day (E%)
library(dplyr)
nutr_info_AA[,3:22] <- nutr_info_AA[,3:22]* 17
nutr_info_AA <- nutr_info_AA %>% mutate(across(c(3:22), .fns = ~./GJ)*100)

#####  
data <- fread("/mnt/nuuk/2021/HRGM/FOCUS_HRGM_abundancies.csv.gz")
data <- data[focus.call == "BL"]

#correlation analysis
subj <- unique(data$subject)
relAA <- unique(Auxotrophy_2$Compound)
View(Auxotrophy_2)
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
rm(sumfreq_nutr_auxo)

View(nutr_info_auxo)
sumfreq_nutr_auxo <- aggregate(nutr_info_auxo$freq, by=list(subject=nutr_info_auxo$subject, AA=nutr_info_auxo$Compound), FUN=sum)
nutri_all_info <- merge(sumfreq_nutr_auxo,nutr_info_AA, by.x="subject", by.y="new_id")


colnames <- colnames(nutri_all_info)
newcol <- colnames[-c(1:4,25)]

AAI <- unique(newcol)
AA <- unique(nutri_all_info$AA)
rm(l)
l <- list()
k <- 1

for(AAi in AA) {
  print(AAi)
  for(AAIi in AAI){
    table <- nutri_all_info[nutri_all_info$AA == AAi, ]
    res <- cor.test(table$x, table[[AAIi]], method = "spearman", exact = FALSE)
    table_corr <- data.table(Pvalue = res$p.value,
                             AA = AAi,
                             day_int = AAIi,
                             Estimate = res$estimate)
    l[[k]] <- table_corr
    k <- k+1
  }
}


nutr_auxo <- rbindlist(l)
View(nutr_auxo)
rm(nutr_auxo)
nutr_auxo[Pvalue< 0.05, sign.label1 := "not adj.P < 0.05"]

nutr_auxo[, padj := p.adjust(Pvalue, method = "fdr")]
nutr_auxo[padj< 0.05, sign.label1 := "Padj < 0.05"]

heatmap_EAA_auxo <- ggplot(nutr_auxo, aes(AA,day_int, fill = Estimate)) +
  geom_tile() +
  geom_point(aes(shape=sign.label1), size=1) +
  scale_fill_gradient2(high = "#ca0020", mid="white", low = "#0571b0") +
  theme_minimal() +
  scale_shape_manual(values = 8, na.translate = FALSE)+
  xlab("Amino acid auxotrophies") +
  ylab("E(%)") 
  
Tryp <- nutri_all_info[nutri_all_info$AA == "Trp",]
cor.test(Tryp$x, Tryp$EALA, method = "spearman")














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








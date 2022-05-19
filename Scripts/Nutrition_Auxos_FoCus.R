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
nutri_all_info <- merge(sumfreq_nutr_auxo,nutr_info_AA, by.x="subject", by.y="new_id")
nutri_all_info <- merge(nutri_all_info,FoCus_data, by.x="subject", by.y="sample")
FoCus_data <- fread("/mnt/nuuk/2021/HRGM/FOCUS_meta_info_v2.csv")


colnames <- colnames(nutri_all_info2)
newcol <- colnames[-c(1:4,23)]

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

#nutr_auxo[Pvalue< 0.05, sign.label1 := "not adj.P < 0.05"]

nutr_auxo[, padj := p.adjust(Pvalue, method = "fdr")]
nutr_auxo[padj< 0.05, sign.label1 := "Padj < 0.05"]

heatmap_EAA_auxo <- ggplot(nutr_auxo, aes(AA,day_int, fill = Estimate)) +
  geom_tile(linetype = 1.5, colour ="grey", lwd=0.2) +
  geom_point(aes(shape=sign.label1), size=1, show.legend = FALSE) +
  scale_fill_gradient2(high = "#ca0020", mid="white", low = "#0571b0") +
  theme_minimal() +
  theme(axis.text.x = element_text(colour="black", size = 10)) +
  theme(axis.text.y = element_text(colour = "black", size=10)) +
  theme(axis.title.x = element_text(colour="black", size = 12, face = "bold")) +
  theme(axis.title.y = element_text(colour = "black", size =12, face = "bold")) +
  scale_shape_manual(values = 8, na.translate = FALSE)+
  xlab("Amino acid auxotrophies") +
  ylab("E(%)") +
  labs(shape = "") +
  theme(panel.background = element_blank())

heatmap_EAA_auxo

ggsave("output/plots/heatmap_nutr_AA_adj_pvalue.pdf", plot = heatmap_EAA_auxo,
       width = 9, height = 6)

Tryp <- nutri_all_info[nutri_all_info$AA == "Trp",]
cor.test(Tryp$x, Tryp$EALA, method = "spearman")


#linear model
colnames <- colnames(nutri_all_info2)
newcol <- colnames[-c(1:4,23)]

AAI <- unique(newcol)
AA <- unique(nutri_all_info2$AA)
rm(l)
l <- list()
k <- 1

for(AAi in AA) {
  print(AAi)
  m0 <- lm(formula = x ~ Gender + age + BMI + parodontitis, data = nutri_all_info2[AA == AAi])
  m0.sum <- summary(m0)
  lin_mod <- m0.sum$coefficients
  lin_mod <- data.table(lin_mod)
  lin_mod$AA <- AAi
  lin_mod$factor <- c("Intercept","sex", "age", "BMI","parodontitis")
  lin_nutr[[k]] <- lin_mod
  k <- k +1
}


linear_mod_parodontitis <- rbindlist(lin_parodontitis)
remove(lin_parodontitis,m0,m0.sum,lin_mod)
names(linear_mod_parodontitis)[names(linear_mod_parodontitis) == "Pr(>|t|)"] <- "pvalue"

linear_mod_parodontitis_adjust <- linear_mod_parodontitis[factor == "parodontitis"]
linear_mod_parodontitis_adjust$padjust = p.adjust(linear_mod_parodontitis_adjust$pvalue, method = "fdr")
linear_mod_parodontitis_adjust[padjust < 0.05, sign.label1 := "P < 0.05"]
linear_mod_parodontitis_adjust <- data.table(linear_mod_parodontitis_adjust)

parodontitis <- ggplot(linear_mod_parodontitis_adjust, aes(AA, factor, fill = `t value`))+
  geom_tile() +
  labs(x = "Auxotrophy", y = "", shape = "")+
  geom_point(aes(shape = sign.label1), size = 1) +
  scale_fill_gradient2(high = "#b2182b", mid = "white", low = "#2166ac") +
  scale_shape_manual(values = 8, na.translate = FALSE) +
  theme_minimal() +
  theme(legend.position = "bottom",
        legend.justification = 	1,
        axis.text.x = element_text(color = "black", angle = 45, hjust = 1, size = 10),
        axis.text.y = element_text(color = "black", size = 10)) +
  theme(axis.title.y = element_text(margin = margin(t =0, r = 20, b= 0, l = 0))) +
  theme(axis.title.x = element_text(face  = "bold"))+
  theme(panel.background = element_blank()) +
  theme(legend.title = element_text(size=9)) +
  labs(fill="") +
  theme(legend.position = "right",
        legend.justification = 	1) +
  theme(legend.text = element_text(size=9)) +
  theme(panel.grid.major = element_blank())
parodontitis

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








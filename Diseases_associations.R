###############     Auxotrophies and health parameters/diseases   ##############
data <- fread("/mnt/nuuk/2021/HRGM/FOCUS_HRGM_abundancies.csv.gz")
focus_info <- fread("/mnt/nuuk/2021/HRGM/FOCUS_meta_info.csv")
View(focus_info)
focus_all_info <- merge(data,focus_info, by.x="subject", by.y="subject")
#filter for BL
focus_all_info <- focus_all_info[focus.call == "BL"]
View(focus_all_info)
#filter
library(data.table)
Auxotrophy
Auxotrophy_2[,c(4:16)] <- NULL

############################# prepare data #####################################
sub <- unique(focus_all_info$subject)
relAA <- unique(Auxotrophy_2$Compound)
Auxotrophy_2
h <- list()
k <- 1

for (subi in sub) {
  print(subi) 
  for (AAi in relAA) {
    x <- focus_all_info[subject == subi]
    y <- Auxotrophy_2[Compound == AAi]
    z <- merge(x,y, by.x = "model", by.y = "Genomes")
    t <- z[Prototrophy == 0]
    h[[k]] <- t
    k <- k +1
  }
}
info_auxo <- rbindlist(h) 
View(info_auxo)
#########################          diabetes       ##############################
sumfreq_diabetes <- aggregate(info_auxo$freq, by=list(subject=info_auxo$subject, AA=info_auxo$Compound, diabetes=info_auxo$diabetes), FUN=sum)
View(sumfreq)
sumfreq_diabetes <- data.table(sumfreq_diabetes)
View(sumfreq_diabetes)
AA <- unique(sumfreq_diabetes$AA)
sumfreq_diabetes$diabetes.status <- ifelse(sumfreq_diabetes$diabetes ==1, "yes", "no")


l <- list()
k <- 1

for (AAi in AA) {
  print(AAi)
  wilcox_data <- sumfreq_diabetes[AA == AAi]
  wilcox <- wilcox_test(x~ diabetes, data = wilcox_data)
  wilcox.pvalue <- data.table(wilcox.p = wilcox$p)
  wilcox.pvalue$AA <- AAi
  new <- tapply(wilcox_data$x, wilcox_data$diabetes, mean)
  new5 <- tapply(wilcox_data$x, wilcox_data$diabetes, median)
  new1 <- t(new)
  new2 <- t(new5)
  mean1 <- new1[1,1]
  mean2 <- new1[1,2]
  wilcox.pvalue$mean1 <- mean1
  wilcox.pvalue$mean2 <- mean2
  median1 <- new2[1,1]
  median2 <- new2[1,2]
  wilcox.pvalue$median1 <- median1
  wilcox.pvalue$median2 <- median2
  
  l[[k]] <- wilcox.pvalue
  k <- k +1
}

diabetes_auxos <- rbindlist(l) 
diabetes_auxos

diabetes_auxos[wilcox.p < 0.05, sign.label := "P < 0.05"]
diabetes_auxos

diabetes <- ggplot(sumfreq_diabetes, aes(AA, x, fill = diabetes.status)) +
  geom_boxplot(outlier.shape = NA) +
  geom_signif(y_position = c(0.90,0.18,0.07), xmin = c(5.8,11.8,16.8), xmax = c(6.2,12.2,17.2), annotation = c("*", "*", "*"), tip_length = 0)+
  ylab("Frequence of auxotrophic bacteria[%]") +
  xlab("Amino acid auxotrophies") +
  theme(axis.line = element_line(size=0.2, colour = "black")) +
  theme(panel.background = element_rect(fill="white", colour= "white")) +
  theme(axis.title.y = element_text(colour = "black", size = 10, face = "bold", margin = margin(0,10,0,0)))+
  theme(axis.title.x = element_text(colour = "black", size = 10, face = "bold", margin = margin(10,0,0,0))) +
  theme(axis.text.x = element_text(size=8, colour = "black", hjust = 1,  angle = 45, margin = margin(10,0,0,0))) +
  theme(axis.text.y = element_text(size =8, colour = "black")) +
  theme(legend.title = element_text(size = 8)) +
  theme(legend.text = element_text(size=8)) +
  theme(plot.margin= margin(0.5,0.5,0.5,0.5, "cm")) +
  scale_fill_manual(values=c("#4575b4", "#d73027")) +
  labs(fill = "Diabetes status")
diabetes
ggsave("output/plots/diabetes.pdf", plot = diabetes,
       width = 9, height = 5)

#########################         IBD      ##############################
sumfreq_IBD <- aggregate(info_auxo$freq, by=list(subject=info_auxo$subject, AA=info_auxo$Compound, IBD=info_auxo$IBD), FUN=sum)
sumfreq_IBD <- data.table(sumfreq_IBD)
sumfreq_IBD$IBD.status <- ifelse(sumfreq_IBD$IBD ==1, "yes", "no")

AA <- unique(sumfreq_IBD$AA)
remove(l)
l <- list()
k <- 1

for (AAi in AA) {
  print(AAi)
  wilcox_data <- sumfreq_IBD[AA == AAi]
  wilcox <- wilcox_test(x~ IBD, data = wilcox_data)
  wilcox.pvalue <- data.table(wilcox.p = wilcox$p)
  wilcox.pvalue$AA <- AAi
  l[[k]] <- wilcox.pvalue
  k <- k +1
}

IBD_auxos <- rbindlist(l) 
IBD_auxos

IBD_auxos[wilcox.p < 0.05, sign.label := "P < 0.05"]

IBD <- ggplot(sumfreq_IBD, aes(AA, x, fill = IBD.status)) +
  geom_boxplot(outlier.shape = NA) +
  geom_signif(y_position = c(0.07), xmin = c(16.8), xmax = c(17.2), annotation = c("*"), tip_length = 0)+
  ylab("Frequence of auxotrophic bacteria[%]") +
  xlab("Amino acid auxotrophies") +
  theme(axis.line = element_line(size=0.2, colour = "black")) +
  theme(panel.background = element_rect(fill="white", colour= "white")) +
  theme(axis.title.y = element_text(colour = "black", size = 10, face = "bold", margin = margin(0,10,0,0)))+
  theme(axis.title.x = element_text(colour = "black", size = 10, face = "bold", margin = margin(10,0,0,0))) +
  theme(axis.text.x = element_text(size=8, colour = "black", hjust = 1,  angle = 45, margin = margin(10,0,0,0))) +
  theme(axis.text.y = element_text(size = 8, colour = "black")) +
  theme(plot.margin= margin(0.5,0.5,0.5,0.5, "cm")) +
  scale_fill_manual(values=c("#4575b4", "#d73027")) +
  theme(legend.title = element_text(size = 8)) +
  theme(legend.text = element_text(size=8))  +
  labs(fill = "IBD status")
  
  
IBD
ggsave("output/plots/IBD.pdf", plot = IBD,
       width = 9, height = 5)


#########################       chronic_diarrhea  ##############################
sumfreq_chronic_diarrhea <- aggregate(info_auxo$freq, by=list(subject=info_auxo$subject, AA=info_auxo$Compound, chronic_diarrhea=info_auxo$chronic_diarrhea), FUN=sum)
sumfreq_chronic_diarrhea <- data.table(sumfreq_chronic_diarrhea)
sumfreq_chronic_diarrhea$chronic_diarrhea.status <- ifelse(sumfreq_chronic_diarrhea$chronic_diarrhea ==1, "yes", "no")

AA <- unique(sumfreq_chronic_diarrhea$AA)
remove(l)
l <- list()
k <- 1

for (AAi in AA) {
  print(AAi)
  wilcox_data <- sumfreq_chronic_diarrhea[AA == AAi]
  wilcox <- wilcox_test(x~ chronic_diarrhea, data = wilcox_data)
  wilcox.pvalue <- data.table(wilcox.p = wilcox$p)
  wilcox.pvalue$AA <- AAi
  l[[k]] <- wilcox.pvalue
  k <- k +1
}

chronic_diarrhea_auxos <- rbindlist(l) 
chronic_diarrhea_auxos[wilcox.p < 0.05, sign.label := "P < 0.05"]
chronic_diarrhea_auxos

chron_diar <- ggplot(sumfreq_chronic_diarrhea, aes(AA, x, fill = chronic_diarrhea.status)) +
  geom_boxplot(outlier.shape = NA) +
  geom_signif(y_position = c(0.07), xmin = c(16.8), xmax = c(17.2), annotation = c("*"), tip_length = 0)+
  ylab("Frequence of auxotrophic bacteria[%]") +
  xlab("Amino acid auxotrophies") +
  theme(axis.line = element_line(size=0.2, colour = "black")) +
  theme(panel.background = element_rect(fill="white", colour= "white")) +
  theme(axis.title.y = element_text(colour = "black", size = 10, face = "bold", margin = margin(0,10,0,0)))+
  theme(axis.title.x = element_text(colour = "black", size = 10, face = "bold", margin = margin(10,0,0,0))) +
  theme(axis.text.x = element_text(size=8, colour = "black", hjust = 1,  angle = 45, margin = margin(10,0,0,0))) +
  theme(axis.text.y = element_text(size = 8, colour = "black")) +
  theme(plot.margin= margin(0.5,0.5,0.5,0.5, "cm")) +
  scale_fill_manual(values=c("#4575b4", "#d73027")) +
  theme(legend.title = element_text(size = 8)) +
  theme(legend.text = element_text(size=8)) +
  labs(fill = "Chronic diarrhea status")
chron_diar
ggsave("output/plots/chron_diar.pdf", plot = chron_diar,
       width = 9, height = 5)

#########################       IBS  ##############################
sumfreq_IBS <- aggregate(info_auxo$freq, by=list(subject=info_auxo$subject, AA=info_auxo$Compound, IBS=info_auxo$IBS), FUN=sum)
sumfreq_IBS <- data.table(sumfreq_IBS)
sumfreq_IBS$IBS.status <- ifelse(sumfreq_IBS$IBS == 1, "yes", "no")

AA <- unique(sumfreq_IBS$AA)
remove(l)
l <- list()
k <- 1

for (AAi in AA) {
  print(AAi)
  wilcox_data <- sumfreq_IBS[AA == AAi]
  wilcox <- wilcox_test(x~ IBS, data = wilcox_data)
  wilcox.pvalue <- data.table(wilcox.p = wilcox$p)
  wilcox.pvalue$AA <- AAi
  l[[k]] <- wilcox.pvalue
  k <- k +1
}

IBS_auxos <- rbindlist(l) 
IBS_auxos[wilcox.p < 0.05, sign.label := "P < 0.05"]
IBS_auxos

IBS <- ggplot(sumfreq_IBS, aes(AA, x, fill = IBS.status)) +
  geom_boxplot(outlier.shape = NA) +
  geom_signif(y_position = c(0.30,0.40), xmin = c(12.8,15.8), xmax = c(13.2,16.2), annotation = c("*","*"), tip_length = 0)+
  ylab("Frequence of auxotrophic bacteria[%]") +
  xlab("Amino acid auxotrophies") +
  theme(axis.line = element_line(size=0.2, colour = "black")) +
  theme(panel.background = element_rect(fill="white", colour= "white")) +
  theme(axis.title.y = element_text(colour = "black", size = 10, face = "bold", margin = margin(0,10,0,0)))+
  theme(axis.title.x = element_text(colour = "black", size = 10, face = "bold", margin = margin(10,0,0,0))) +
  theme(axis.text.x = element_text(size=8, colour = "black", hjust = 1,  angle = 45, margin = margin(10,0,0,0))) +
  theme(axis.text.y = element_text(size = 8, colour = "black")) +
  theme(plot.margin= margin(0.5,0.5,0.5,0.5, "cm"))  +
  scale_fill_manual(values=c("#4575b4", "#d73027")) +
  theme(legend.title = element_text(size = 8)) +
  theme(legend.text = element_text(size=8)) +
  labs(fill = "IBS status")

ggsave("output/plots/IBS.pdf", plot = IBS,
       width = 9, height = 5)


###############               correlation with BMI                ##############
sumfreq_BMI <- aggregate(info_auxo$freq, by=list(subject=info_auxo$subject, AA=info_auxo$Compound, BMI=info_auxo$BMI), FUN=sum)
sumfreq_BMI <- data.table(sumfreq_BMI)
ggplot(sumfreq_BMI[AA=="Cys"], aes(x, BMI)) +
  geom_point() +
  geom_smooth()


relAA <- unique(sumfreq_BMI$AA)
correlation <- list()
k <- 1
for (AAi in relAA) {
  print(AAi)
  test <- sumfreq_BMI[AA == AAi]
  corr <- cor.test(test$x, test$BMI, method = "spearman", exact = FALSE)
  corre <- data.table(pvalue = corr$p.value,
                     rho = corr$estimate,
                     Aminoacid = AAi,
                     factor = "BMI") 
  correlation[[k]] <- corre
  k <- k+1
}
correlation_BMI <- rbindlist(correlation)
correlation_BMI

ggplot(correlation_BMI, aes (Aminoacid, rho)) +
  geom_tile()




###############               correlation with weight      #####################
sumfreq_weight <- aggregate(info_auxo$freq, by=list(subject=info_auxo$subject, AA=info_auxo$Compound, Weight=info_auxo$weight), FUN=sum)
sumfreq_weight <- data.table(sumfreq_weight)
relAA <- unique(sumfreq_weight$AA)
correlation_weight <- list()
k <- 1
for (AAi in relAA) {
  print(AAi)
  test <- sumfreq_weight[AA == AAi]
  corr <- cor.test(test$x, test$Weight, method = "spearman", exact = FALSE)
  corre <- data.table(pvalue = corr$p.value,
                      rho = corr$estimate,
                      Aminoacid = AAi,
                      factor = "Weight") 
  correlation_weight[[k]] <- corre
  k <- k+1
}
correlation_weight <- rbindlist(correlation_weight)
correlation_weight

ggplot(sumfreq_weight[AA=="Cys"], aes(x, Weight)) +
  geom_point()+
  geom_smooth()



###############               correlation with age         #####################
sumfreq_age <- aggregate(info_auxo$freq, by=list(subject=info_auxo$subject, AA=info_auxo$Compound, Age=info_auxo$age), FUN=sum)
sumfreq_age <- data.table(sumfreq_age)
relAA <- unique(sumfreq_age$AA)
correlation_age <- list()
k <- 1
for (AAi in relAA) {
  print(AAi)
  test <- sumfreq_age[AA == AAi]
  corr <- cor.test(test$x, test$Age, method = "spearman", exact = FALSE)
  corre <- data.table(pvalue = corr$p.value,
                      rho = corr$estimate,
                      Aminoacid = AAi,
                      factor = "Age") 
  correlation_age[[k]] <- corre
  k <- k+1
}
correlation_age <- rbindlist(correlation_age)
correlation_age

###visualization
cys_age <- ggscatter(sumfreq_age[AA=="Cys"], x= "x", y = "Age", cor.method = "spearman", 
          cor.coef = TRUE, add = "reg.line")



######################## create big heatmap with all factors ###################
corr_all <- rbind(correlation_BMI, correlation_weight, correlation_age)
corr_all[pvalue < 0.05, sign.label1 := "Padj < 0.05"]
ggplot(corr_all, aes(Aminoacid, factor, fill = rho))+
  geom_tile() +
  labs(x = "Auxotrophy", y = "", shape = "")+
  geom_point(aes(shape = sign.label1), size = 0.5) +
  scale_fill_gradient2(high = "#ca0020", mid = "white", low = "#0571b0") +
  scale_shape_manual(values = 8, na.translate = FALSE) +
theme_classic() +
  theme(legend.position = "bottom",
        legend.justification = 	1,
        axis.text.x = element_text(color = "black", angle = 45, hjust = 1, size = 10),
        axis.text.y = element_text(color = "black", size = 10)) +
  theme(axis.title.y = element_text(margin = margin(t =0, r = 20, b= 0, l = 0))) +
  theme(panel.background = element_blank())








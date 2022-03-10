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
remove(wilcox_data, wilcox, wilcox.pvalue,new,new5,new1,new2, mean1,mean2,median1,median2,l)

diabetes_auxos$padjust = p.adjust(diabetes_auxos$wilcox.p, method = "fdr")
diabetes_auxos[padjust < 0.05, sign.label := "P < 0.05"]
pvalue_sumfreq_diabetes <- merge(sumfreq_diabetes, diabetes_auxos, by.x="AA", by.y="AA")
diabetes_auxos

diabetes <- ggplot(pvalue_sumfreq_diabetes, aes(AA, x, fill = diabetes.status)) +
  geom_boxplot(outlier.shape = NA) +
  ylab("Frequence of auxotrophic bacteria[%]") +
  xlab("Amino acid auxotrophies") +
  theme(axis.line = element_line(size=0.2, colour = "black")) +
  geom_point(aes(shape = sign.label), size = 0.5) +
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
remove(wilcox_data, wilcox, wilcox.pvalue,l)

IBD_auxos$padjust = p.adjust(IBD_auxos$wilcox.p, method = "fdr")
IBD_auxos[padjust < 0.05, sign.label := "P < 0.05"]
IBD_auxos
sumfreq_IBD
pvalue_sumfreq_IBD <- merge(sumfreq_IBD, IBD_auxos, by.x="AA", by.y="AA")

IBD <- ggplot(pvalue_sumfreq_IBD, aes(AA, x, fill = IBD.status)) +
  geom_boxplot(outlier.shape = NA) +
  geom_signif(y_position = c(0.07), xmin = c(16.8), xmax = c(17.2), annotation = c("*"), tip_length = 0)+
  ylab("Frequence of auxotrophic bacteria[%]") +
  geom_point(aes(shape = sign.label), size = 0.5) +
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
l <- list()
k <- 1

for (AAi in AA) {
  print(AAi)
  wilcox_data <- sumfreq_chronic_diarrhea[AA == AAi]
  wilcox <- wilcox_test(x~ chronic_diarrhea, data = wilcox_data)
  wilcox.pvalue <- data.table(wilcox.p = wilcox$p)
  wilcox.pvalue$AA <- AAi
  new <- tapply(wilcox_data$x, wilcox_data$chronic_diarrhea, mean)
  new5 <- tapply(wilcox_data$x, wilcox_data$chronic_diarrhea, median)
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

chronic_diarrhea_auxos <- rbindlist(l) 
remove(wilcox_data, wilcox, wilcox.pvalue,new,new5,new1,new2, mean1,mean2,median1,median2,l)
chronic_diarrhea_auxos$padjust = p.adjust(chronic_diarrhea_auxos$wilcox.p, method = "fdr")
chronic_diarrhea_auxos[padjust < 0.05, sign.label := "P < 0.05"]
chronic_diarrhea_auxos
pvalue_sumfreq_chronic_diarrhea <- merge(sumfreq_chronic_diarrhea, chronic_diarrhea_auxos, by.x="AA", by.y="AA")
View()
chron_diar <- ggplot(pvalue_sumfreq_chronic_diarrhea, aes(AA, x, fill = chronic_diarrhea.status)) +
  geom_boxplot(outlier.shape = NA) +
  geom_signif(y_position = c(0.20), xmin = c(1.8), xmax = c(2.2), annotation = c("*"), tip_length = 0)+
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
  labs(y = "Frequence of auxotrophic bacteria[%]", x = "Amino acid auxotrophies")+
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
remove(wilcox_data, wilcox, l)
IBS_auxos$padjust = p.adjust(IBS_auxos$wilcox.p, method = "fdr")
IBS_auxos[padjust < 0.05, sign.label := "P < 0.05"]
IBS_auxos
pvalue_sumfreq_IBS <- merge(sumfreq_IBS, IBS_auxos, by.x="AA", by.y="AA")

IBS <- ggplot(pvalue_sumfreq_IBS, aes(AA, x, fill = IBS.status)) +
  geom_boxplot(outlier.shape = NA) +
  ylab("Frequence of auxotrophic bacteria[%]") +
  xlab("Amino acid auxotrophies") +
  theme(axis.line = element_line(size=0.2, colour = "black")) +
  geom_point(aes(shape = sign.label), size = 0.2) +
  theme(panel.background = element_rect(fill="white", colour= "white")) +
  theme(axis.title.y = element_text(colour = "black", size = 10, face = "bold", margin = margin(0,10,0,0)))+
  theme(axis.title.x = element_text(colour = "black", size = 10, face = "bold", margin = margin(10,0,0,0))) +
  theme(axis.text.x = element_text(size=8, colour = "black", hjust = 1,  angle = 45, margin = margin(10,0,0,0))) +
  theme(axis.text.y = element_text(size = 8, colour = "black")) +
  theme(plot.margin= margin(0.5,0.5,0.5,0.5, "cm"))  +
  scale_shape_manual(values = 8, na.translate = FALSE)+
  scale_fill_manual(values=c("#4575b4", "#d73027")) +
  theme(legend.title = element_text(size = 8)) +
  theme(legend.text = element_text(size=8)) +
  labs(fill = "IBS status")
IBS
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
remove(correlation,test,corr,corre)
correlation_BMI$padjust = p.adjust(correlation_BMI$pvalue, method = "fdr")
correlation_BMI[padjust < 0.05, sign.label1 := "P < 0.05"]
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
remove(correlation,test,corr,corre)
correlation_weight$padjust = p.adjust(correlation_weight$pvalue, method = "fdr")
correlation_weight[padjust < 0.05, sign.label1 := "P < 0.05"]
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
remove(correlation,test,corr,corre)
correlation_age$padjust = p.adjust(correlation_age$pvalue, method = "fdr")
correlation_age[padjust < 0.05, sign.label1 := "P < 0.05"]

###visualization
cys_age <- ggscatter(sumfreq_age[AA=="Arg"], x= "x", y = "Age", cor.method = "spearman", 
                     cor.coef = TRUE, add = "reg.line")
cys_age


######################## create big heatmap with all factors ###################
corr_all <- rbind(correlation_BMI, correlation_weight, correlation_age)
corr_health <- ggplot(corr_all, aes(Aminoacid, factor, fill = rho))+
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
  labs(fill="Spearman correlation") +
  theme(legend.position = "top",
        legend.justification = 	1) +
  theme(legend.text = element_text(size=9)) +
  theme(panel.grid.major = element_blank())
corr_health
ggsave("output/plots/health_diseases_all_spearman.pdf", plot = corr_health,
       width = 6, height = 3)


##########################    linear modeling    ###############################
sumfreq_all <- aggregate(info_auxo$freq, by=list(subject=info_auxo$subject, AA=info_auxo$Compound,
                                                 BMI=info_auxo$BMI, sex=info_auxo$gender,age=info_auxo$age), FUN=sum)
sumfreq_all <- data.table(sumfreq_all)
test <- lm(formula = x ~ sex + age + BMI, data = sumfreq_all[AA == "Cys"])
summary(test)
relAA1 <- unique(sumfreq_all$AA)
lin <- list()
k <- 1
for (AAi in relAA1){
  print(AAi)
  m0 <- lm(formula = x ~ sex + age + BMI, data = sumfreq_all[AA == AAi])
  m0.sum <- summary(m0)
  lin_mod <- m0.sum$coefficients
  lin_mod <- data.table(lin_mod)
  lin_mod$AA <- AAi
  lin_mod$factor <- c("Intercept","sex", "age", "BMI")
  lin[[k]] <- lin_mod
  k <- k +1
  
}

linear_mod <- rbindlist(lin)
remove(lin,m0,m0.sum,lin_mod)
names(linear_mod)[names(linear_mod) == "Pr(>|t|)"] <- "pvalue"

linear_mod_sex <- linear_mod[factor == "sex"]
linear_mod_age <- linear_mod[factor == "age"]
linear_mod_BMI <- linear_mod[factor == "BMI"]

##sex
linear_mod_sex$padjust = p.adjust(linear_mod_sex$pvalue, method = "fdr")
linear_mod_sex[padjust < 0.05, sign.label1 := "P < 0.05"]

##age
linear_mod_age$padjust = p.adjust(linear_mod_age$pvalue, method = "fdr")
linear_mod_age[padjust < 0.05, sign.label1 := "P < 0.05"]

##BMI
linear_mod_BMI$padjust = p.adjust(linear_mod_BMI$pvalue, method = "fdr")
linear_mod_BMI[padjust < 0.05, sign.label1 := "P < 0.05"]

##bind all together
linear_all <- rbind(linear_mod_sex, linear_mod_age, linear_mod_BMI)

##visualization BMI and age
linear_wo_sex <- rbind(linear_mod_age, linear_mod_BMI)
linear_health_BMI_age <- ggplot(linear_wo_sex, aes(factor, AA, fill = Estimate))+
  geom_tile() +
  labs(y = "Auxotrophy", x = "Health parameter", shape = "")+
  geom_point(aes(shape = sign.label1), size = 1) +
  scale_fill_gradient2(high = "#b2182b", mid = "white", low = "#2166ac") +
  scale_shape_manual(values = 8, na.translate = FALSE) +
  theme_minimal() +
  theme(legend.position = "right",
        legend.justification = 	0.5,
        axis.text.x = element_text(color = "black", angle = 90, vjust = 0.2, size = 8),
        axis.text.y = element_text(color = "black", size = 8)) +
  theme(axis.title.y = element_text(margin = margin(t =0, r = 20, b= 0, l = 0))) +
  theme(axis.title.x = element_text(face  = "bold", size = 10, margin = margin(t=20, r = 0, b= 0, l = 0))) +
  theme(axis.title.y = element_text(face  = "bold", size = 10))+ 
  scale_x_discrete("Health parameter", labels = c("age" = "    Age", "BMI" = "    BMI")) +
  theme(panel.background = element_blank()) +
  theme(legend.title = element_text(size=7)) + 
  theme(legend.text = element_text(size=7)) +
  labs(fill="") +
  theme(legend.position = "top",
        legend.justification = 	1) +
  theme(legend.text = element_text(size=7)) +
  theme(panel.grid.major = element_blank())
linear_health_BMI_age + guides(shape = guide_legend(order = 1))
linear_health_BMI_age1 <- linear_health_BMI_age + guides(shape= "none")
linear_health_BMI_age1

ggsave("output/plots/health_diseases_lin_BMI_age.pdf", plot = linear_health_BMI_age,
       width = 6, height = 3)

###visualization for sex
linear_mod_sex <- ggplot(linear_mod_sex, aes(AA, factor, fill = Estimate))+
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
linear_mod_sex 
ggsave("output/plots/health_diseases_lin_sex.pdf", plot = linear_mod_sex,
       width = 6, height = 3)


#visualization for diabetes
sumfreq_diabetes_lin <- aggregate(info_auxo$freq, by=list(subject=info_auxo$subject, AA=info_auxo$Compound,
                                                 BMI=info_auxo$BMI, sex=info_auxo$gender,age=info_auxo$age, diabetes =info_auxo$diabetes), FUN=sum)
sumfreq_diabetes_lin <- data.table(sumfreq_diabetes_lin)
View(sumfreq_diabetes_lin)
m0 <- lm(formula = x ~ sex + age + BMI + diabetes, data = sumfreq_diabetes_lin[AA == "Val"])

relAA <- unique(sumfreq_diabetes_lin$AA)
lin_diab <- list()
k <- 1
for (AAi in relAA){
  print(AAi)
  m0 <- lm(formula = x ~ sex + age + BMI + diabetes, data = sumfreq_diabetes_lin[AA == AAi])
  m0.sum <- summary(m0)
  lin_mod <- m0.sum$coefficients
  lin_mod <- data.table(lin_mod)
  lin_mod$AA <- AAi
  lin_mod$factor <- c("Intercept","sex", "age", "BMI","diabetes")
  lin_diab[[k]] <- lin_mod
  k <- k +1
}
linear_mod_diab <- rbindlist(lin_diab)
remove(lin_diab,m0,m0.sum,lin_mod)
names(linear_mod_diab)[names(linear_mod_diab) == "Pr(>|t|)"] <- "pvalue"


linear_mod_diab_adjust <- linear_mod_diab[factor == "diabetes"]
linear_mod_diab_adjust$padjust = p.adjust(linear_mod_diab_adjust$pvalue, method = "fdr")
linear_mod_diab_adjust[padjust < 0.05, sign.label1 := "P < 0.05"]
linear_mod_diab_adjust <- data.table(linear_mod_diab_adjust)

d<- ggplot(linear_mod_diab_adjust[factor == "diabetes"], aes(AA, factor, fill = Estimate))+
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
d
ggsave("output/plots/health_diseases_lin_sex.pdf", plot = d,
       width = 6, height = 3)


#visualization for IBD
sumfreq_IBD_lin <- aggregate(info_auxo$freq, by=list(subject=info_auxo$subject, AA=info_auxo$Compound,
                                                          BMI=info_auxo$BMI, sex=info_auxo$gender,age=info_auxo$age, IBD =info_auxo$IBD), FUN=sum)
sumfreq_IBD_lin <- data.table(sumfreq_IBD_lin)
View(sumfreq_IBD_lin)

relAA <- unique(sumfreq_IBD_lin$AA)
lin_IBD <- list()
k <- 1
for (AAi in relAA){
  print(AAi)
  m0 <- lm(formula = x ~ sex + age + BMI + IBD, data = sumfreq_IBD_lin[AA == AAi])
  m0.sum <- summary(m0)
  lin_mod <- m0.sum$coefficients
  lin_mod <- data.table(lin_mod)
  lin_mod$AA <- AAi
  lin_mod$factor <- c("Intercept","sex", "age", "BMI","IBD")
  lin_IBD[[k]] <- lin_mod
  k <- k +1
}
linear_mod_IBD <- rbindlist(lin_IBD)
remove(lin_IBD,m0, m0.sum,lin_mod)
names(linear_mod_IBD)[names(linear_mod_IBD) == "Pr(>|t|)"] <- "pvalue"

linear_mod_IBD_adjust <- linear_mod_IBD[factor == "IBD"]
linear_mod_IBD_adjust$padjust = p.adjust(linear_mod_IBD_adjust$pvalue, method = "fdr")
linear_mod_IBD_adjust[padjust < 0.05, sign.label1 := "P < 0.05"]
linear_mod_IBD_adjust <- data.table(linear_mod_IBD_adjust)

i <- ggplot(linear_mod_IBD_adjust, aes(AA, factor, fill = Estimate))+
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
i
ggsave("output/plots/health_diseases_lin_IBD.pdf", plot = i,
       width = 6, height = 3)


#visualization for chronic diarrhea
sumfreq_chrond_lin <- aggregate(info_auxo$freq, by=list(subject=info_auxo$subject, AA=info_auxo$Compound,
                                                     BMI=info_auxo$BMI, sex=info_auxo$gender,age=info_auxo$age, chrond =info_auxo$chronic_diarrhea), FUN=sum)
sumfreq_chrond_lin <- data.table(sumfreq_chrond_lin)
View(sumfreq_chrond_lin)

relAA <- unique(sumfreq_chrond_lin$AA)
lin_chrond <- list()
k <- 1
for (AAi in relAA){
  print(AAi)
  m0 <- lm(formula = x ~ sex + age + BMI + chrond, data = sumfreq_chrond_lin[AA == AAi])
  m0.sum <- summary(m0)
  lin_mod <- m0.sum$coefficients
  lin_mod <- data.table(lin_mod)
  lin_mod$AA <- AAi
  lin_mod$factor <- c("Intercept","sex", "age", "BMI","chrond")
  lin_chrond[[k]] <- lin_mod
  k <- k +1
}
linear_mod_chrond <- rbindlist(lin_chrond)
remove(lin_chrond,m0,m0.sum,lin_mod)
names(linear_mod_chrond)[names(linear_mod_chrond) == "Pr(>|t|)"] <- "pvalue"

linear_mod_chrond_adjust <- linear_mod_chrond[factor == "chrond"]
linear_mod_chrond_adjust$padjust = p.adjust(linear_mod_chrond_adjust$pvalue, method = "fdr")
linear_mod_chrond_adjust[padjust < 0.05, sign.label1 := "P < 0.05"]
linear_mod_chrond_adjust <- data.table(linear_mod_chrond_adjust)

c <- ggplot(linear_mod_chrond_adjust, aes(AA, factor, fill = Estimate))+
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
c
ggsave("output/plots/health_diseases_lin_chronic_diarrhea.pdf", plot = c,
       width = 6, height = 3)

#visualization for IBS
sumfreq_IBS_lin <- aggregate(info_auxo$freq, by=list(subject=info_auxo$subject, AA=info_auxo$Compound,
                                                        BMI=info_auxo$BMI, sex=info_auxo$gender,age=info_auxo$age, IBS =info_auxo$IBS), FUN=sum)
sumfreq_IBS_lin <- data.table(sumfreq_IBS_lin)
View(sumfreq_chrond_lin)

relAA <- unique(sumfreq_IBS_lin$AA)
lin_IBS <- list()
k <- 1
for (AAi in relAA){
  print(AAi)
  m0 <- lm(formula = x ~ sex + age + BMI + IBS, data = sumfreq_IBS_lin[AA == AAi])
  m0.sum <- summary(m0)
  lin_mod <- m0.sum$coefficients
  lin_mod <- data.table(lin_mod)
  lin_mod$AA <- AAi
  lin_mod$factor <- c("Intercept","sex", "age", "BMI","IBS")
  lin_IBS[[k]] <- lin_mod
  k <- k +1
}
linear_mod_IBS <- rbindlist(lin_IBS)
remove(lin_IBS,m0,m0.sum,lin_mod)
names(linear_mod_IBS)[names(linear_mod_IBS) == "Pr(>|t|)"] <- "pvalue"

linear_mod_IBS_adjust <- linear_mod_IBS[factor == "IBS"]
linear_mod_IBS_adjust$padjust = p.adjust(linear_mod_IBS_adjust$pvalue, method = "fdr")
linear_mod_IBS_adjust[padjust < 0.05, sign.label1 := "P < 0.05"]
linear_mod_IBS_adjust <- data.table(linear_mod_IBS_adjust)

ibs <- ggplot(linear_mod_IBS_adjust, aes(AA, factor, fill = Estimate))+
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
ibs
ggsave("output/plots/health_diseases_IBS.pdf", plot = ibs,
       width = 6, height = 3)

################## all heatmap conbined of the linear modeling #################
all <- rbind(linear_mod_IBS_adjust, linear_mod_IBD_adjust, linear_mod_diab_adjust, linear_mod_chrond_adjust)
all1 <- ggplot(all, aes(factor, AA, fill = Estimate))+
  geom_tile() +
  labs(y = "", x = "Diseases", shape = "")+
  geom_point(aes(shape = sign.label1), size = 1) +
  scale_fill_gradient2(high = "#b2182b", mid = "white", low = "#2166ac") +
  scale_shape_manual(values = 8, na.translate = FALSE) +
  theme_minimal() +
  scale_x_discrete("Diseases", labels = c("IBS" = "IBS", "IBD" = "IBD", "diabetes" = "Diab.", "chrond" = "Chron.\nDiarr.")) +
  theme(legend.position = "top", legend.box = "horizontal",
        legend.justification = 	0.5,
        axis.text.x = element_text(color = "black", angle = 90, vjust = 0.2, size = 8),
        axis.text.y = element_text(color = "black", size = 8)) +
  theme(axis.title.y = element_text(margin = margin(t =0, r = 20, b= 0, l = 0))) +
  theme(axis.title.x = element_text(face  = "bold", size = 10, margin = margin(t=20, r = 0, b= 0, l = 0))) +
  theme(axis.title.y = element_text(face  = "bold", size = 10))+
  theme(panel.background = element_blank()) +
  theme(legend.title = element_text(size=9)) +
  labs(fill="", x = "") +
  theme(legend.text = element_text(size=7)) +
  theme(panel.grid.major = element_blank())
all1 + guides(shape = guide_legend(order = 1))
all2 <- all1 + guides(shape= "none")
all2

ggsave("output/plots/health_diseases_linear_mod_all.pdf", plot = all,
       width = 6, height = 3)












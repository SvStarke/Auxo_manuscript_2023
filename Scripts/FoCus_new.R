###############   new FoCus cohort data   ##################

FoCus_data <- fread("/mnt/nuuk/2021/HRGM/FOCUS_meta_info_v2.csv")
describe(FoCus_data$sample)
#checking distribution of data
hist(FoCus_data$IL6)
describe(FoCus_data$IL6)
FoCus_data$IL6
View(FoCus_data)

#create NA value for unrealistic value
FoCus_data$IL6[FoCus_data$IL6 == -999.0] <- NA

boxplot(FoCus_data$IL6)
#correcting for some values for IL6,CRP, HOMA
#IL6
IL6_1 <- (quantile(FoCus_data$IL6, prob = 0.75, na.rm = TRUE) - quantile(FoCus_data$IL6, prob = 0.25, na.rm = TRUE)) * 2
FoCus_data$IL6[FoCus_data$IL6 > IL6_1] <- NA

#CRP
CRP_1 <- (quantile(FoCus_data$CRP, prob = 0.75, na.rm = TRUE) - quantile(FoCus_data$CRP, prob = 0.25, na.rm = TRUE)) * 2
FoCus_data$CRP[FoCus_data$CRP > CRP_1] <- NA

#HOMA
HOMA_1 <- (quantile(FoCus_data$HOMA, prob = 0.75, na.rm = TRUE) - quantile(FoCus_data$HOMA, prob = 0.25, na.rm = TRUE)) * 2
FoCus_data$HOMA[FoCus_data$HOMA > HOMA_1] <- NA

#linear model adjusted for BMI,age,gender
data <- fread("/mnt/nuuk/2021/HRGM/FOCUS_HRGM_abundancies.csv.gz")
describe(data$subject)
describe(FoCus_data$sample)
describe(FoCus_data$hypertens)
#filter for BL
data <- data[focus.call == "BL"]

####merge data
FoCus_info <- merge(data, FoCus_data, by.x="subject", by.y="sample")
describe(FoCus_info$subject)

#include only relevant information
Auxotrophy_2[,c(4:16)] <- NULL


############################# prepare data #####################################
sub <- unique(FoCus_info$subject)
relAA <- unique(Auxotrophy_2$Compound)
Auxotrophy_2
h <- list()
k <- 1

for (subi in sub) {
  print(subi) 
  for (AAi in relAA) {
    x <- FoCus_info[subject == subi]
    y <- Auxotrophy_2[Compound == AAi]
    z <- merge(x,y, by.x = "model", by.y = "Genomes")
    t <- z[Prototrophy == 0]
    h[[k]] <- t
    k <- k +1
  }
}
info_auxo <- rbindlist(h) 
describe(info_auxo$subject)
View(info_auxo)

#############################    BMI,age,sex  ##################################
sumfreq_all <- aggregate(info_auxo$freq, by=list(subject=info_auxo$subject, AA=info_auxo$Compound,
                                                 BMI=info_auxo$BMI, sex=info_auxo$Gender,age=info_auxo$age), FUN=sum)
describe(sumfreq_all$subject)
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
linear_wo_sex <- rbind(linear_mod_age, linear_mod_BMI, linear_mod_sex)
linear_health_BMI_age <- ggplot(linear_wo_sex, aes(factor, AA, fill =`t value`))+
  geom_tile() +
  labs(y = "Frequency of auxotrophic bacteria [%]", x = "Health parameter", shape = "")+
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

###visualization for sex
# linear_mod_sex <- ggplot(linear_mod_sex, aes(AA, factor, fill = Estimate))+
#   geom_tile() +
#   labs(x = "Auxotrophy", y = "", shape = "")+
#   geom_point(aes(shape = sign.label1), size = 1) +
#   scale_fill_gradient2(high = "#b2182b", mid = "white", low = "#2166ac") +
#   scale_shape_manual(values = 8, na.translate = FALSE) +
#   theme_minimal() +
#   theme(legend.position = "bottom",
#         legend.justification = 	1,
#         axis.text.x = element_text(color = "black", angle = 45, hjust = 1, size = 10),
#         axis.text.y = element_text(color = "black", size = 10)) +
#   theme(axis.title.y = element_text(margin = margin(t =0, r = 20, b= 0, l = 0))) +
#   theme(axis.title.x = element_text(face  = "bold"))+
#   theme(panel.background = element_blank()) +
#   theme(legend.title = element_text(size=9)) +
#   labs(fill="") +
#   theme(legend.position = "right",
#         legend.justification = 	1) +
#   theme(legend.text = element_text(size=9)) +
#   theme(panel.grid.major = element_blank())
# linear_mod_sex 
# linear_mod_sex  + guides(shape = guide_legend(order = 1))
# linear_mod_sex  <- linear_mod_sex  + guides(shape= "none")
# linear_mod_sex 
# ggsave("output/plots/health_diseases_lin_sex.pdf", plot = linear_mod_sex,
#        width = 6, height = 3)


#visualization for diabetes
sumfreq_diabetes_lin <- aggregate(info_auxo$freq, by=list(subject=info_auxo$subject, AA=info_auxo$Compound,
                                                          BMI=info_auxo$BMI, sex=info_auxo$Gender,age=info_auxo$age, diabetes =info_auxo$diabetes), FUN=sum)
sumfreq_diabetes_lin <- data.table(sumfreq_diabetes_lin)
describe(sumfreq_diabetes_lin$subject)
View(info_auxo)
describe(sumfreq_diabetes_lin$diabetes)
View(sumfreq_diabetes_lin)
m0 <- lm(formula = x ~ sex + age + BMI + diabetes, data = sumfreq_diabetes_lin[AA == "Val"])
n <- partial_Spearman(x|diabetes~ sex + age + BMI, data = sumfreq_diabetes_lin[AA == "Trp"])
n

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
  n <- partial_Spearman(x|diabetes~ sex + age + BMI, data = sumfreq_diabetes_lin[AA == AAi])
  lin_mod$Estimate_Part_Spear <- n$TS$TB$ts
  lin_mod$pvalue_Part_Spear <- n$TS$TB$pval
  lin_diab[[k]] <- lin_mod
  k <- k +1
}

linear_mod_diab <- rbindlist(lin_diab)
remove(lin_diab,m0,m0.sum,lin_mod,n)
names(linear_mod_diab)[names(linear_mod_diab) == "Pr(>|t|)"] <- "pvalue"


linear_mod_diab_adjust <- linear_mod_diab[factor == "diabetes"]
linear_mod_diab_adjust$padjust = p.adjust(linear_mod_diab_adjust$pvalue, method = "fdr")
linear_mod_diab_adjust[padjust < 0.05, sign.label1 := "P < 0.05"]
linear_mod_diab_adjust$padjust_Spear = p.adjust(linear_mod_diab_adjust$pvalue_Part_Spear, method = "fdr")
linear_mod_diab_adjust[padjust_Spear < 0.05, sign.label2 := "P < 0.05"]
linear_mod_diab_adjust <- data.table(linear_mod_diab_adjust)

d<- ggplot(linear_mod_diab_adjust[factor == "diabetes"], aes(AA, factor, fill = `t value`))+
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
                                                     BMI=info_auxo$BMI, sex=info_auxo$Gender,age=info_auxo$age, IBD =info_auxo$IBD), FUN=sum)
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
  n <- partial_Spearman(x|IBD~ sex + age + BMI, data = sumfreq_IBD_lin[AA == AAi])
  lin_mod$Estimate_Part_Spear <- n$TS$TB$ts
  lin_mod$pvalue_Part_Spear <- n$TS$TB$pval
  lin_IBD[[k]] <- lin_mod
  k <- k +1
}
linear_mod_IBD <- rbindlist(lin_IBD)
remove(lin_IBD,m0, m0.sum,lin_mod,n)
names(linear_mod_IBD)[names(linear_mod_IBD) == "Pr(>|t|)"] <- "pvalue"

linear_mod_IBD_adjust <- linear_mod_IBD[factor == "IBD"]
linear_mod_IBD_adjust$padjust = p.adjust(linear_mod_IBD_adjust$pvalue, method = "fdr")
linear_mod_IBD_adjust[padjust < 0.05, sign.label1 := "P < 0.05"]
linear_mod_IBD_adjust$padjust_Spear = p.adjust(linear_mod_IBD_adjust$pvalue_Part_Spear, method = "fdr")
linear_mod_IBD_adjust[padjust_Spear < 0.05, sign.label2 := "P < 0.05"]
linear_mod_IBD_adjust <- data.table(linear_mod_IBD_adjust)

i <- ggplot(linear_mod_IBD_adjust, aes(AA, factor, fill = `t value`))+
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
                                                        BMI=info_auxo$BMI, sex=info_auxo$Gender,age=info_auxo$age, chrond =info_auxo$chr_diarrh), FUN=sum)
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
  n <- partial_Spearman(x|chrond~ sex + age + BMI, data = sumfreq_chrond_lin[AA == AAi])
  lin_mod$Estimate_Part_Spear <- n$TS$TB$ts
  lin_mod$pvalue_Part_Spear <- n$TS$TB$pval
  lin_chrond[[k]] <- lin_mod
  k <- k +1
}
linear_mod_chrond <- rbindlist(lin_chrond)
remove(lin_chrond,m0,m0.sum,lin_mod)
names(linear_mod_chrond)[names(linear_mod_chrond) == "Pr(>|t|)"] <- "pvalue"

linear_mod_chrond_adjust <- linear_mod_chrond[factor == "chrond"]
linear_mod_chrond_adjust$padjust = p.adjust(linear_mod_chrond_adjust$pvalue, method = "fdr")
linear_mod_chrond_adjust[padjust < 0.05, sign.label1 := "P < 0.05"]
linear_mod_chrond_adjust$padjust_Spear = p.adjust(linear_mod_chrond_adjust$pvalue_Part_Spear, method = "fdr")
linear_mod_chrond_adjust[padjust_Spear < 0.05, sign.label2 := "P < 0.05"]
linear_mod_chrond_adjust <- data.table(linear_mod_chrond_adjust)

c <- ggplot(linear_mod_chrond_adjust, aes(AA, factor, fill = `t value`))+
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
                                                     BMI=info_auxo$BMI, sex=info_auxo$Gender,age=info_auxo$age, IBS =info_auxo$IBS), FUN=sum)
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
  n <- partial_Spearman(x|IBS~ sex + age + BMI, data = sumfreq_IBS_lin[AA == AAi])
  lin_mod$Estimate_Part_Spear <- n$TS$TB$ts
  lin_mod$pvalue_Part_Spear <- n$TS$TB$pval
  lin_IBS[[k]] <- lin_mod
  k <- k +1
}
linear_mod_IBS <- rbindlist(lin_IBS)
remove(lin_IBS,m0,m0.sum,lin_mod)
names(linear_mod_IBS)[names(linear_mod_IBS) == "Pr(>|t|)"] <- "pvalue"

linear_mod_IBS_adjust <- linear_mod_IBS[factor == "IBS"]
linear_mod_IBS_adjust$padjust = p.adjust(linear_mod_IBS_adjust$pvalue, method = "fdr")
linear_mod_IBS_adjust[padjust < 0.05, sign.label1 := "P < 0.05"]
linear_mod_IBS_adjust$padjust_Spear = p.adjust(linear_mod_IBS_adjust$pvalue_Part_Spear, method = "fdr")
linear_mod_IBS_adjust[padjust_Spear < 0.05, sign.label2 := "P < 0.05"]
linear_mod_IBS_adjust <- data.table(linear_mod_IBS_adjust)

ibs <- ggplot(linear_mod_IBS_adjust, aes(AA, factor, fill = `t value`))+
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

#visualization for HOMA
sumfreq_HOMA_lin <- aggregate(info_auxo$freq, by=list(subject=info_auxo$subject, AA=info_auxo$Compound,
                                                     BMI=info_auxo$BMI, sex=info_auxo$Gender,age=info_auxo$age, HOMA =info_auxo$HOMA), FUN=sum)
sumfreq_HOMA_lin <- data.table(sumfreq_HOMA_lin)
View(sumfreq_chrond_lin)

relAA <- unique(sumfreq_HOMA_lin$AA)
lin_HOMA <- list()
k <- 1
for (AAi in relAA){
  print(AAi)
  m0 <- lm(formula = x ~ sex + age + BMI + HOMA, data = sumfreq_HOMA_lin[AA == AAi])
  m0.sum <- summary(m0)
  lin_mod <- m0.sum$coefficients
  lin_mod <- data.table(lin_mod)
  lin_mod$AA <- AAi
  lin_mod$factor <- c("Intercept","sex", "age", "BMI","HOMA")
  n <- partial_Spearman(x|HOMA ~ sex + age + BMI, data = sumfreq_HOMA_lin[AA == AAi])
  lin_mod$Estimate_Part_Spear <- n$TS$TB$ts
  lin_mod$pvalue_Part_Spear <- n$TS$TB$pval
  lin_HOMA[[k]] <- lin_mod
  k <- k +1
}
linear_mod_HOMA <- rbindlist(lin_HOMA)
remove(lin_HOMA,m0,m0.sum,lin_mod)
names(linear_mod_HOMA)[names(linear_mod_HOMA) == "Pr(>|t|)"] <- "pvalue"

linear_mod_HOMA_adjust <- linear_mod_HOMA[factor == "HOMA"]
linear_mod_HOMA_adjust$padjust = p.adjust(linear_mod_HOMA_adjust$pvalue, method = "fdr")
linear_mod_HOMA_adjust[padjust < 0.05, sign.label1 := "P < 0.05"]
linear_mod_HOMA_adjust$padjust_Spear = p.adjust(linear_mod_HOMA_adjust$pvalue_Part_Spear, method = "fdr")
linear_mod_HOMA_adjust[padjust_Spear < 0.05, sign.label2 := "P < 0.05"]
linear_mod_HOMA_adjust <- data.table(linear_mod_HOMA_adjust)

homa <- ggplot(linear_mod_HOMA_adjust, aes(AA, factor, fill = `t value`))+
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
homa


#visualization for IL6
sumfreq_IL6_lin <- aggregate(info_auxo$freq, by=list(subject=info_auxo$subject, AA=info_auxo$Compound,
                                                      BMI=info_auxo$BMI, sex=info_auxo$Gender,age=info_auxo$age, IL6 =info_auxo$IL6), FUN=sum)
sumfreq_IL6_lin <- data.table(sumfreq_IL6_lin)
View(sumfreq_chrond_lin)

relAA <- unique(sumfreq_IL6_lin$AA)
lin_IL6 <- list()
k <- 1
for (AAi in relAA){
  print(AAi)
  m0 <- lm(formula = x ~ sex + age + BMI + IL6, data = sumfreq_IL6_lin[AA == AAi])
  m0.sum <- summary(m0)
  lin_mod <- m0.sum$coefficients
  lin_mod <- data.table(lin_mod)
  lin_mod$AA <- AAi
  lin_mod$factor <- c("Intercept","sex", "age", "BMI","IL6")
  n <- partial_Spearman(x|IL6~ sex + age + BMI, data = sumfreq_IL6_lin[AA == AAi])
  lin_mod$Estimate_Part_Spear <- n$TS$TB$ts
  lin_mod$pvalue_Part_Spear <- n$TS$TB$pval
  lin_IL6[[k]] <- lin_mod
  k <- k +1
}
linear_mod_IL6 <- rbindlist(lin_IL6)
remove(lin_IL6,m0,m0.sum,lin_mod)
names(linear_mod_IL6)[names(linear_mod_IL6) == "Pr(>|t|)"] <- "pvalue"

linear_mod_IL6_adjust <- linear_mod_IL6[factor == "IL6"]
linear_mod_IL6_adjust$padjust = p.adjust(linear_mod_IL6_adjust$pvalue, method = "fdr")
linear_mod_IL6_adjust[padjust < 0.05, sign.label1 := "P < 0.05"]
linear_mod_IL6_adjust$padjust_Spear = p.adjust(linear_mod_IL6_adjust$pvalue_Part_Spear, method = "fdr")
linear_mod_IL6_adjust[padjust_Spear < 0.05, sign.label2 := "P < 0.05"]
linear_mod_IL6_adjust <- data.table(linear_mod_IL6_adjust)

IL6 <- ggplot(linear_mod_IL6_adjust, aes(AA, factor, fill = `t value`))+
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
IL6

#visualization for CRP
sumfreq_CRP_lin <- aggregate(info_auxo$freq, by=list(subject=info_auxo$subject, AA=info_auxo$Compound,
                                                      BMI=info_auxo$BMI, sex=info_auxo$Gender,age=info_auxo$age, CRP =info_auxo$CRP), FUN=sum)
sumfreq_CRP_lin <- data.table(sumfreq_CRP_lin)
View(sumfreq_chrond_lin)

relAA <- unique(sumfreq_CRP_lin$AA)
lin_CRP <- list()
k <- 1
for (AAi in relAA){
  print(AAi)
  m0 <- lm(formula = x ~ sex + age + BMI + CRP, data = sumfreq_CRP_lin[AA == AAi])
  m0.sum <- summary(m0)
  lin_mod <- m0.sum$coefficients
  lin_mod <- data.table(lin_mod)
  lin_mod$AA <- AAi
  lin_mod$factor <- c("Intercept","sex", "age", "BMI","CRP")
  n <- partial_Spearman(x|CRP~ sex + age + BMI, data = sumfreq_CRP_lin[AA == AAi])
  lin_mod$Estimate_Part_Spear <- n$TS$TB$ts
  lin_mod$pvalue_Part_Spear <- n$TS$TB$pval
  lin_CRP[[k]] <- lin_mod
  k <- k +1
}
linear_mod_CRP <- rbindlist(lin_CRP)
remove(lin_CRP,m0,m0.sum,lin_mod)
names(linear_mod_CRP)[names(linear_mod_CRP) == "Pr(>|t|)"] <- "pvalue"

linear_mod_CRP_adjust <- linear_mod_CRP[factor == "CRP"]
linear_mod_CRP_adjust$padjust = p.adjust(linear_mod_CRP_adjust$pvalue, method = "fdr")
linear_mod_CRP_adjust[padjust < 0.05, sign.label1 := "P < 0.05"]
linear_mod_CRP_adjust$padjust_Spear = p.adjust(linear_mod_CRP_adjust$pvalue_Part_Spear, method = "fdr")
linear_mod_CRP_adjust[padjust_Spear < 0.05, sign.label2 := "P < 0.05"]
linear_mod_CRP_adjust <- data.table(linear_mod_CRP_adjust)

CRP <- ggplot(linear_mod_CRP_adjust, aes(AA, factor, fill = `t value`))+
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
CRP
#visualization for triglyceride
sumfreq_TG_lin <- aggregate(info_auxo$freq, by=list(subject=info_auxo$subject, AA=info_auxo$Compound,
                                                          BMI=info_auxo$BMI, sex=info_auxo$Gender,age=info_auxo$age, TG =info_auxo$triglyc), FUN=sum)
sumfreq_TG_lin <- data.table(sumfreq_TG_lin)

relAA <- unique(sumfreq_TG_lin$AA)
lin_TG <- list()
k <- 1
for (AAi in relAA){
  print(AAi)
  m0 <- lm(formula = x ~ sex + age + BMI + TG, data = sumfreq_TG_lin[AA == AAi])
  m0.sum <- summary(m0)
  lin_mod <- m0.sum$coefficients
  lin_mod <- data.table(lin_mod)
  lin_mod$AA <- AAi
  lin_mod$factor <- c("Intercept","sex", "age", "BMI","TG")
  n <- partial_Spearman(x|TG~ sex + age + BMI, data = sumfreq_TG_lin[AA == AAi])
  lin_mod$Estimate_Part_Spear <- n$TS$TB$ts
  lin_mod$pvalue_Part_Spear <- n$TS$TB$pval
  lin_TG[[k]] <- lin_mod
  k <- k +1
}
linear_mod_TG <- rbindlist(lin_TG)
remove(lin_TG,m0,m0.sum,lin_mod)
names(linear_mod_TG)[names(linear_mod_TG) == "Pr(>|t|)"] <- "pvalue"


linear_mod_TG_adjust <- linear_mod_TG[factor == "TG"]
linear_mod_TG_adjust$padjust = p.adjust(linear_mod_TG_adjust$pvalue, method = "fdr")
linear_mod_TG_adjust[padjust < 0.05, sign.label1 := "P < 0.05"]
linear_mod_TG_adjust$padjust_Spear = p.adjust(linear_mod_TG_adjust$pvalue_Part_Spear, method = "fdr")
linear_mod_TG_adjust[padjust_Spear < 0.05, sign.label2 := "P < 0.05"]
linear_mod_TG_adjust <- data.table(linear_mod_TG_adjust)

d <- ggplot(linear_mod_TG_adjust[factor == "TG"], aes(AA, factor, fill = `t value`))+
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
################## all heatmap conbined of the linear modeling #################
all <- rbind(linear_mod_IBS_adjust, linear_mod_IBD_adjust, linear_mod_diab_adjust, linear_mod_chrond_adjust, linear_mod_HOMA_adjust, linear_mod_IL6_adjust, linear_mod_CRP_adjust, linear_mod_TG_adjust)
all1 <- ggplot(all, aes(factor, AA, fill = `t value`))+
  geom_tile() +
  labs(y = "", x = "Diseases", shape = "")+
  geom_point(aes(shape = sign.label1), size = 1) +
  scale_fill_gradient2(high = "#b2182b", mid = "white", low = "#2166ac") +
  scale_shape_manual(values = 8, na.translate = FALSE) +
  theme_minimal() +
  scale_x_discrete("Diseases and inflammatory markers", labels = c("IBS" = "IBS", "IBD" = "IBD", "diabetes" = "Diab.", "chrond" = "Chron.\nDiarr.")) +
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

ggsave("output/plots/health_diseases_linear_mod_all.pdf", plot = all2,
       width = 6, height = 3)






View(info_auxo)
#visualization for cancer
sumfreq_cancer_lin <- aggregate(info_auxo$freq, by=list(subject=info_auxo$subject, AA=info_auxo$Compound,
                                                     BMI=info_auxo$BMI, sex=info_auxo$Gender,age=info_auxo$age, Cancer =info_auxo$cancer), FUN=sum)
sumfreq_cancer_lin <- data.table(sumfreq_cancer_lin)
View(sumfreq_chrond_lin)

relAA <- unique(sumfreq_cancer_lin$AA)
lin_cancer <- list()
k <- 1
for (AAi in relAA){
  print(AAi)
  m0 <- lm(formula = x ~ sex + age + BMI + Cancer, data = sumfreq_cancer_lin[AA == AAi])
  m0.sum <- summary(m0)
  lin_mod <- m0.sum$coefficients
  lin_mod <- data.table(lin_mod)
  lin_mod$AA <- AAi
  lin_mod$factor <- c("Intercept","sex", "age", "BMI","Cancer")
  n <- partial_Spearman(x|cancer~ sex + age + BMI, data = sumfreq_cancer_lin[AA == AAi])
  lin_mod$Estimate_Part_Spear <- n$TS$TB$ts
  lin_mod$pvalue_Part_Spear <- n$TS$TB$pval
  lin_cancer[[k]] <- lin_mod
  k <- k +1
}
linear_mod_cancer <- rbindlist(lin_cancer)
remove(lin_cancer,m0,m0.sum,lin_mod)
names(linear_mod_cancer)[names(linear_mod_cancer) == "Pr(>|t|)"] <- "pvalue"

linear_mod_cancer_adjust <- linear_mod_cancer[factor == "Cancer"]
linear_mod_cancer_adjust$padjust = p.adjust(linear_mod_cancer_adjust$pvalue, method = "fdr")
linear_mod_cancer_adjust[padjust < 0.05, sign.label1 := "P < 0.05"]
linear_mod_cancer_adjust$padjust_Spear = p.adjust(linear_mod_cancer_adjust$pvalue_Part_Spear, method = "fdr")
linear_mod_cancer_adjust[padjust_Spear < 0.05, sign.label2 := "P < 0.05"]
linear_mod_cancer_adjust <- data.table(linear_mod_cancer_adjust)

cancer <- ggplot(linear_mod_cancer_adjust, aes(AA, factor, fill = `t value`))+
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
cancer


#visualization for liver disease
sumfreq_liverdis_lin <- aggregate(info_auxo$freq, by=list(subject=info_auxo$subject, AA=info_auxo$Compound,
                                                        BMI=info_auxo$BMI, sex=info_auxo$Gender,age=info_auxo$age, liverdis =info_auxo$liverdisease), FUN=sum)
sumfreq_liverdis_lin <- data.table(sumfreq_liverdis_lin)
View(sumfreq_chrond_lin)

relAA <- unique(sumfreq_liverdis_lin$AA)
lin_liverdis <- list()
k <- 1
for (AAi in relAA){
  print(AAi)
  m0 <- lm(formula = x ~ sex + age + BMI + liverdis, data = sumfreq_liverdis_lin[AA == AAi])
  m0.sum <- summary(m0)
  lin_mod <- m0.sum$coefficients
  lin_mod <- data.table(lin_mod)
  lin_mod$AA <- AAi
  lin_mod$factor <- c("Intercept","sex", "age", "BMI","liverdis")
  n <- partial_Spearman(x|liverdis~ sex + age + BMI, data = sumfreq_liverdis_lin[AA == AAi])
  lin_mod$Estimate_Part_Spear <- n$TS$TB$ts
  lin_mod$pvalue_Part_Spear <- n$TS$TB$pval
  lin_liverdis[[k]] <- lin_mod
  k <- k +1
}
linear_mod_liverdis <- rbindlist(lin_liverdis)
remove(lin_liverdis,m0,m0.sum,lin_mod)
names(linear_mod_liverdis)[names(linear_mod_liverdis) == "Pr(>|t|)"] <- "pvalue"

linear_mod_liverdis_adjust <- linear_mod_liverdis[factor == "liverdis"]
linear_mod_liverdis_adjust$padjust = p.adjust(linear_mod_liverdis_adjust$pvalue, method = "fdr")
linear_mod_liverdis_adjust[padjust < 0.05, sign.label1 := "P < 0.05"]
linear_mod_liverdis_adjust$padjust_Spear = p.adjust(linear_mod_liverdis_adjust$pvalue_Part_Spear, method = "fdr")
linear_mod_liverdis_adjust[padjust_Spear < 0.05, sign.label2 := "P < 0.05"]
linear_mod_liverdis_adjust <- data.table(linear_mod_liverdis_adjust)

liverdis <- ggplot(linear_mod_liverdis_adjust, aes(AA, factor, fill = `t value`))+
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
liverdis


#visualization for rheumatic disease
sumfreq_rheumato_lin <- aggregate(info_auxo$freq, by=list(subject=info_auxo$subject, AA=info_auxo$Compound,
                                                          BMI=info_auxo$BMI, sex=info_auxo$Gender,age=info_auxo$age, rheumato =info_auxo$rheumato), FUN=sum)
sumfreq_rheumato_lin <- data.table(sumfreq_rheumato_lin)
View(sumfreq_chrond_lin)

relAA <- unique(sumfreq_rheumato_lin$AA)
lin_rheumato <- list()
k <- 1
for (AAi in relAA){
  print(AAi)
  m0 <- lm(formula = x ~ sex + age + BMI + rheumato, data = sumfreq_rheumato_lin[AA == AAi])
  m0.sum <- summary(m0)
  lin_mod <- m0.sum$coefficients
  lin_mod <- data.table(lin_mod)
  lin_mod$AA <- AAi
  lin_mod$factor <- c("Intercept","sex", "age", "BMI","rheumato")
  n <- partial_Spearman(x|rheumato~ sex + age + BMI, data = sumfreq_rheumato_lin[AA == AAi])
  lin_mod$Estimate_Part_Spear <- n$TS$TB$ts
  lin_mod$pvalue_Part_Spear <- n$TS$TB$pval
  lin_rheumato[[k]] <- lin_mod
  k <- k +1
}
linear_mod_rheumato <- rbindlist(lin_rheumato)
remove(lin_rheumato,m0,m0.sum,lin_mod)
names(linear_mod_rheumato)[names(linear_mod_rheumato) == "Pr(>|t|)"] <- "pvalue"

linear_mod_rheumato_adjust <- linear_mod_rheumato[factor == "rheumato"]
linear_mod_rheumato_adjust$padjust = p.adjust(linear_mod_rheumato_adjust$pvalue, method = "fdr")
linear_mod_rheumato_adjust[padjust < 0.05, sign.label1 := "P < 0.05"]
linear_mod_rheumato_adjust$padjust_Spear = p.adjust(linear_mod_rheumato_adjust$pvalue_Part_Spear, method = "fdr")
linear_mod_rheumato_adjust[padjust_Spear < 0.05, sign.label2 := "P < 0.05"]

linear_mod_rheumato_adjust <- data.table(linear_mod_rheumato_adjust)

rheumato <- ggplot(linear_mod_rheumato_adjust, aes(AA, factor, fill = `t value`))+
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
rheumato


#visualization for hypertension
sumfreq_hypertens_lin <- aggregate(info_auxo$freq, by=list(subject=info_auxo$subject, AA=info_auxo$Compound,
                                                          BMI=info_auxo$BMI, sex=info_auxo$Gender,age=info_auxo$age, hypertens =info_auxo$hypertens), FUN=sum)
sumfreq_hypertens_lin <- data.table(sumfreq_hypertens_lin)
View(sumfreq_chrond_lin)
describe(info_auxo$hypertens)
relAA <- unique(sumfreq_hypertens_lin$AA)
lin_hypertens <- list()
k <- 1
for (AAi in relAA){
  print(AAi)
  m0 <- lm(formula = x ~ sex + age + BMI + hypertens, data = sumfreq_hypertens_lin[AA == AAi])
  m0.sum <- summary(m0)
  lin_mod <- m0.sum$coefficients
  lin_mod <- data.table(lin_mod)
  lin_mod$AA <- AAi
  lin_mod$factor <- c("Intercept","sex", "age", "BMI","hypertens")
  n <- partial_Spearman(x|hypertens~ sex + age + BMI, data = sumfreq_hypertens_lin[AA == AAi])
  lin_mod$Estimate_Part_Spear <- n$TS$TB$ts
  lin_mod$pvalue_Part_Spear <- n$TS$TB$pval
  lin_hypertens[[k]] <- lin_mod
  k <- k +1
}

linear_mod_hypertens <- rbindlist(lin_hypertens)
remove(lin_hypertens,m0,m0.sum,lin_mod)
names(linear_mod_hypertens)[names(linear_mod_hypertens) == "Pr(>|t|)"] <- "pvalue"

linear_mod_hypertens_adjust <- linear_mod_hypertens[factor == "hypertens"]
linear_mod_hypertens_adjust$padjust = p.adjust(linear_mod_hypertens_adjust$pvalue, method = "fdr")
linear_mod_hypertens_adjust[padjust < 0.05, sign.label1 := "P < 0.05"]
linear_mod_hypertens_adjust$padjust_Spear = p.adjust(linear_mod_hypertens_adjust$pvalue_Part_Spear, method = "fdr")
linear_mod_hypertens_adjust[padjust_Spear < 0.05, sign.label2 := "P < 0.05"]
linear_mod_hypertens_adjust <- data.table(linear_mod_hypertens_adjust)

hypertens <- ggplot(linear_mod_hypertens_adjust, aes(AA, factor, fill = `Estimate_Part_Spear`))+
  geom_tile() +
  labs(x = "Auxotrophy", y = "", shape = "")+
  geom_point(aes(shape = sign.label2), size = 1) +
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
hypertens


#visualization for parodontitis
sumfreq_parodontitis_lin <- aggregate(info_auxo$freq, by=list(subject=info_auxo$subject, AA=info_auxo$Compound,
                                                           BMI=info_auxo$BMI, sex=info_auxo$Gender,age=info_auxo$age, parodontitis =info_auxo$parodontitis), FUN=sum)
sumfreq_parodontitis_lin <- data.table(sumfreq_parodontitis_lin)
View(sumfreq_chrond_lin)
describe(info_auxo$parodontitis)
relAA <- unique(sumfreq_parodontitis_lin$AA)
lin_parodontitis <- list()
k <- 1
for (AAi in relAA){
  print(AAi)
  m0 <- lm(formula = x ~ sex + age + BMI + parodontitis, data = sumfreq_parodontitis_lin[AA == AAi])
  m0.sum <- summary(m0)
  lin_mod <- m0.sum$coefficients
  lin_mod <- data.table(lin_mod)
  lin_mod$AA <- AAi
  lin_mod$factor <- c("Intercept","sex", "age", "BMI","parodontitis")
  n <- partial_Spearman(x|parodontitis~ sex + age + BMI, data = sumfreq_parodontitis_lin[AA == AAi])
  lin_mod$Estimate_Part_Spear <- n$TS$TB$ts
  lin_mod$pvalue_Part_Spear <- n$TS$TB$pval
  lin_parodontitis[[k]] <- lin_mod
  k <- k +1
}
linear_mod_parodontitis <- rbindlist(lin_parodontitis)
remove(lin_parodontitis,m0,m0.sum,lin_mod)
names(linear_mod_parodontitis)[names(linear_mod_parodontitis) == "Pr(>|t|)"] <- "pvalue"

linear_mod_parodontitis_adjust <- linear_mod_parodontitis[factor == "parodontitis"]
linear_mod_parodontitis_adjust$padjust = p.adjust(linear_mod_parodontitis_adjust$pvalue, method = "fdr")
linear_mod_parodontitis_adjust[padjust < 0.05, sign.label1 := "P < 0.05"]
linear_mod_parodontitis_adjust$padjust_Spear = p.adjust(linear_mod_parodontitis_adjust$pvalue_Part_Spear, method = "fdr")
linear_mod_parodontitis_adjust[padjust_Spear < 0.05, sign.label2 := "P < 0.05"]
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

########   heatmap for different diseases and health parameters
#health markers
linear_health_marks <- rbind(linear_mod_age, linear_mod_BMI, linear_mod_sex, 
                       linear_mod_CRP_adjust,linear_mod_HOMA_adjust,linear_mod_hypertens_adjust, linear_mod_IL6_adjust, linear_mod_TG_adjust)
linear_health <- ggplot(linear_health_marks, aes(factor, AA, fill =`t value`))+
  geom_tile() +
  labs(y = "Frequency of auxotrophic bacteria [%]", x = "Health markers", shape = "")+
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
  scale_x_discrete("Health markers", labels = c("age" = "    Age", "BMI" = "    BMI", "hypertens" = "Hypertension", "sex" = "Sex", "TG" = "Triglycerids")) +
  theme(panel.background = element_blank()) +
  theme(legend.title = element_text(size=9)) +
  theme(legend.text = element_text(size=7)) +
  labs(fill="T-Value") +
  theme(legend.position = "top",
        legend.justification = 	1) +
  theme(legend.text = element_text(size=7)) +
  theme(panel.grid.major = element_blank())
linear_health + guides(shape = guide_legend(order = 1)) 
linear_health<- linear_health + guides(shape= "none")
linear_health


### diseases 
dis_health_FoCus <- rbind(linear_mod_IBS_adjust, linear_mod_IBD_adjust, linear_mod_diab_adjust, 
             linear_mod_chrond_adjust, linear_mod_cancer_adjust, linear_mod_liverdis_adjust, 
             linear_mod_rheumato_adjust,linear_mod_parodontitis_adjust)
dis_health <- ggplot(dis_health_FoCus, aes(factor, AA, fill = `t value`))+
  geom_tile() +
  labs(y = "", x = "Diseases", shape = "") +
  geom_point(aes(shape = sign.label1), size = 1) +
  scale_fill_gradient2(high = "#b2182b", mid = "white", low = "#2166ac") +
  scale_shape_manual(values = 8, na.translate = FALSE, name = "T-value") +
  theme_minimal() +
  scale_x_discrete("Diseases", labels = c("IBS" = "IBS", "IBD" = "IBD", "chrond" = "Chronic\nDiarrhea", "liverdis" = "Liver", "diabetes" = "Diabetes", "parodontitis" = "Parodontitis", "rheumato"="Rheumatism")) +
  theme(legend.position = "top", legend.box = "horizontal",
        legend.justification = 	0.5,
        axis.text.x = element_text(color = "black", angle = 90, vjust = 0.2, size = 8),
        axis.text.y = element_text(color = "black", size = 8)) +
  theme(axis.title.y = element_text(margin = margin(t =0, r = 20, b= 0, l = 0))) +
  theme(axis.title.x = element_text(face  = "bold", size = 9, margin = margin(t=20, r = 0, b= 0, l = 0))) +
  theme(axis.title.y = element_text(face  = "bold", size = 9))+
  theme(panel.background = element_blank()) +
  theme(legend.title = element_text(size=9)) +
  labs(fill="T-Value", x = "") +
  theme(legend.text = element_text(size=7)) +
  theme(panel.grid.major = element_blank())
dis_health + guides(shape = guide_legend(order = 1))
dis_health  <- dis_health + guides(shape= "none")
dis_health

ggsave("output/plots/health_diseases_linear_mod_diseases_healthparam.pdf", plot = dis_health_param,
       width = 7, height = 5)

#####    Partial Spearman correlation ######

dzhk_auxo <- merge(dzhk_info, dzhk_relabun, by.x= "sample", by.y="sample")
describe(dzhk_auxo$sample)
# hist(dzhk_auxo$Neutrophils)
# hist(dzhk_auxo$Leucocytes)
# hist(dzhk_auxo$Lymphocytes)
# hist(dzhk_auxo$`HDL-Cholesterin`)
# hist(hist(dzhk_auxo$`LDL-Cholesterin`))
# hist(hist(dzhk_auxo$`Triglyceride`))
# hist(hist(dzhk_auxo$Monocytes))
# hist(hist(dzhk_auxo$Eosinophils))
# hist(hist(dzhk_auxo$Basophils))


neutro <- (quantile(dzhk_auxo$Triglyceride, prob = 0.75, na.rm = TRUE) - quantile(dzhk_auxo$Triglyceride, prob = 0.25, na.rm = TRUE)) * 2
dzhk_auxo$Triglyceride[dzhk_auxo$Triglyceride > neutro] <- NA



sub <- unique(dzhk_auxo$sample)
relAA <- unique(Auxotrophy_2$Compound)
Auxotrophy_2
h <- list()
k <- 1

for (subi in sub) {
  print(subi) 
  for (AAi in relAA) {
    x <- dzhk_auxo[sample == subi]
    y <- Auxotrophy_2[Compound == AAi]
    z <- merge(x,y, by.x = "model", by.y = "Genomes")
    t <- z[Prototrophy == 0]
    h[[k]] <- t
    k <- k +1
  }
}
info_auxo_DZHK <- rbindlist(h) 
info_auxo_DZHK <- info_auxo_DZHK[Compound != "Gly"]
info_auxo_DZHK$Geschlecht <- ifelse(info_auxo_DZHK$Geschlecht == "m", 0, 1)
describe(info_auxo_DZHK$Geschlecht)
head(info_auxo_DZHK)

#### BMI and age
sumfreq_all <- aggregate(info_auxo_DZHK$prop, by=list(subject=info_auxo_DZHK$sample, AA=info_auxo_DZHK$Compound,
                                                 BMI=info_auxo_DZHK$BMI, sex=info_auxo_DZHK$Geschlecht,age=info_auxo_DZHK$Age), FUN=sum)

sumfreq_all <- data.table(sumfreq_all)
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
  n1 <- partial_Spearman(x|sex~ age + BMI, data = sumfreq_all[AA == AAi])
  n2 <- partial_Spearman(x|age~ sex + BMI, data = sumfreq_all[AA == AAi])
  n3 <- partial_Spearman(x|BMI~ age + sex, data = sumfreq_all[AA == AAi])
  lin_mod$Estimate_Part_Spear <- c("Intercept",n1$TS$TB$ts, n2$TS$TB$ts,n3$TS$TB$ts)
  lin_mod$pvalue_Part_Spear <- c("Intercept",n1$TS$TB$pval,n2$TS$TB$pval,n3$TS$TB$pval)
  lin[[k]] <- lin_mod
  k <- k +1
  
}

linear_mod <- rbindlist(lin)

names(linear_mod)[names(linear_mod) == "Pr(>|t|)"] <- "pvalue"


linear_mod_age <- linear_mod[factor == "age"]
linear_mod_BMI <- linear_mod[factor == "BMI"]


##BMI
linear_mod_BMI$padjust = p.adjust(linear_mod_BMI$pvalue, method = "fdr")
linear_mod_BMI[padjust < 0.05, sign.label1 := "P < 0.05"]
linear_mod_BMI$padjust_Spear = p.adjust(linear_mod_BMI$pvalue_Part_Spear, method = "fdr")
linear_mod_BMI[padjust_Spear < 0.05, sign.label2 := "P < 0.05"]

#age
linear_mod_age$padjust = p.adjust(linear_mod_age$pvalue, method = "fdr")
linear_mod_age[padjust < 0.05, sign.label1 := "P < 0.05"]
linear_mod_age$padjust_Spear = p.adjust(linear_mod_age$pvalue_Part_Spear, method = "fdr")
linear_mod_age[padjust_Spear < 0.05, sign.label2 := "P < 0.05"]


linear_wo_sex <- rbind(linear_mod_age, linear_mod_BMI)
linear_health_BMI_age <- ggplot(linear_wo_sex, aes(factor, AA, fill =`Estimate`))+
  geom_tile() +
  labs(y = "Frequency of auxotrophic bacteria [%]", x = "Health parameter", shape = "")+
  geom_point(aes(shape = sign.label2), size = 1) +
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

#visualization for LDL-Cholesterin
sumfreq_LDL <- aggregate(info_auxo_DZHK$prop, by=list(subject=info_auxo_DZHK$sample, AA=info_auxo_DZHK$Compound,
                                                      BMI=info_auxo_DZHK$BMI, sex=info_auxo_DZHK$Geschlecht,age=info_auxo_DZHK$Age, LDL=info_auxo_DZHK$`LDL-Cholesterin`), FUN=sum)

sumfreq_LDL <- data.table(sumfreq_LDL)
sumfreq_LDL <- sumfreq_LDL[,`LDL-Cholesterin` := LDL]
View(sumfreq_chrond_lin)

relAA <- unique(sumfreq_LDL$AA)
lin_LDL <- list()
k <- 1
for (AAi in relAA){
  print(AAi)
  m0 <- lm(formula = x ~ sex + age + BMI + LDL, data = sumfreq_LDL[AA == AAi])
  m0.sum <- summary(m0)
  lin_mod <- m0.sum$coefficients
  lin_mod <- data.table(lin_mod)
  lin_mod$AA <- AAi
  lin_mod$factor <- c("Intercept","sex", "age", "BMI","LDL")
  n <- partial_Spearman(x|LDL~ sex + age + BMI, data = sumfreq_LDL[AA == AAi])
  lin_mod$Estimate_Part_Spear <- n$TS$TB$ts
  lin_mod$pvalue_Part_Spear <- n$TS$TB$pval
  lin_LDL[[k]] <- lin_mod
  k <- k +1
}
corr_LDL <- rbindlist(lin_LDL)
remove(lin_HOMA,m0,m0.sum,lin_mod)


corr_LDL_adjust <- corr_LDL[factor == "LDL"]
names(corr_LDL_adjust)[names(corr_LDL_adjust) == "Pr(>|t|)"] <- "pvalue"
corr_LDL_adjust$padjust = p.adjust(corr_LDL_adjust$pvalue, method = "fdr")
corr_LDL_adjust[padjust < 0.05, sign.label1 := "P < 0.05"]
corr_LDL_adjust$padjust_Spear = p.adjust(corr_LDL_adjust$pvalue_Part_Spear, method = "fdr")
corr_LDL_adjust[padjust_Spear < 0.05, sign.label2 := "P < 0.05"]
corr_LDL_adjust <- data.table(corr_LDL_adjust)

LDL <- ggplot(corr_LDL_adjust, aes(AA, factor, fill = `Estimate`))+
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
LDL

#visualization for HDL-Cholesterin
sumfreq_HDL <- aggregate(info_auxo_DZHK$prop, by=list(subject=info_auxo_DZHK$sample, AA=info_auxo_DZHK$Compound,
                                                      BMI=info_auxo_DZHK$BMI, sex=info_auxo_DZHK$Geschlecht,age=info_auxo_DZHK$Age, HDL=info_auxo_DZHK$`HDL-Cholesterin`), FUN=sum)

sumfreq_HDL <- data.table(sumfreq_HDL)

View(sumfreq_chrond_lin)

relAA <- unique(sumfreq_HDL$AA)
lin_HDL <- list()
k <- 1
for (AAi in relAA){
  print(AAi)
  m0 <- lm(formula = x ~ sex + age + BMI + HDL, data = sumfreq_HDL[AA == AAi])
  m0.sum <- summary(m0)
  lin_mod <- m0.sum$coefficients
  lin_mod <- data.table(lin_mod)
  lin_mod$AA <- AAi
  lin_mod$factor <- c("Intercept","sex", "age", "BMI","HDL")
  n <- partial_Spearman(x|HDL~ sex + age + BMI, data = sumfreq_HDL[AA == AAi])
  lin_mod$Estimate_Part_Spear <- n$TS$TB$ts
  lin_mod$pvalue_Part_Spear <- n$TS$TB$pval
  lin_HDL[[k]] <- lin_mod
  k <- k +1
}
corr_HDL <- rbindlist(lin_HDL)
remove(lin_HDL,m0,m0.sum,lin_mod)


corr_HDL_adjust <- corr_HDL[factor == "HDL"]
names(corr_HDL_adjust)[names(corr_HDL_adjust) == "Pr(>|t|)"] <- "pvalue"
corr_HDL_adjust$padjust_Spear = p.adjust(corr_HDL_adjust$pvalue_Part_Spear, method = "fdr")
corr_HDL_adjust[padjust_Spear < 0.05, sign.label2 := "P < 0.05"]
corr_HDL_adjust <- data.table(corr_HDL_adjust)

HDL <- ggplot(corr_HDL_adjust, aes(AA, factor, fill = `Estimate`))+
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
HDL

#visualization for TG
sumfreq_TG <- aggregate(info_auxo_DZHK$prop, by=list(subject=info_auxo_DZHK$sample, AA=info_auxo_DZHK$Compound,
                                                      BMI=info_auxo_DZHK$BMI, sex=info_auxo_DZHK$Geschlecht,age=info_auxo_DZHK$Age, TG=info_auxo_DZHK$Triglyceride), FUN=sum)

sumfreq_TG <- data.table(sumfreq_TG)
sumfreq_TG <- sumfreq_TG[,`TG-Cholesterin` := TG]
#View(sumfreq_chrond_lin)

relAA <- unique(sumfreq_TG$AA)
lin_TG <- list()
k <- 1
for (AAi in relAA){
  print(AAi)
  m0 <- lm(formula = x ~ sex + age + BMI + TG, data = sumfreq_TG[AA == AAi])
  m0.sum <- summary(m0)
  lin_mod <- m0.sum$coefficients
  lin_mod <- data.table(lin_mod)
  lin_mod$AA <- AAi
  lin_mod$factor <- c("Intercept","sex", "age", "BMI","TG")
  n <- partial_Spearman(x|TG~ sex + age + BMI, data = sumfreq_TG[AA == AAi])
  lin_mod$Estimate_Part_Spear <- n$TS$TB$ts
  lin_mod$pvalue_Part_Spear <- n$TS$TB$pval
  lin_TG[[k]] <- lin_mod
  k <- k +1
}
corr_TG <- rbindlist(lin_TG)
remove(lin_TG,m0,m0.sum,lin_mod)


corr_TG_adjust <- corr_TG[factor == "TG"]
names(corr_TG_adjust)[names(corr_TG_adjust) == "Pr(>|t|)"] <- "pvalue"
corr_TG_adjust$padjust = p.adjust(corr_TG_adjust$pvalue, method = "fdr")
corr_TG_adjust[padjust < 0.05, sign.label1 := "P < 0.05"]
corr_TG_adjust$padjust_Spear = p.adjust(corr_TG_adjust$pvalue_Part_Spear, method = "fdr")
corr_TG_adjust[padjust_Spear < 0.05, sign.label2 := "P < 0.05"]
corr_TG_adjust <- data.table(corr_TG_adjust)

#leucocytes
sumfreq_LEU <- aggregate(info_auxo_DZHK$prop, by=list(subject=info_auxo_DZHK$sample, AA=info_auxo_DZHK$Compound,
                                                     BMI=info_auxo_DZHK$BMI, sex=info_auxo_DZHK$Geschlecht,age=info_auxo_DZHK$Age, LEU=info_auxo_DZHK$Leucocytes), FUN=sum)

sumfreq_LEU <- data.table(sumfreq_LEU)

relAA <- unique(sumfreq_LEU$AA)
lin_LEU <- list()
k <- 1
for (AAi in relAA){
  print(AAi)
  m0 <- lm(formula = x ~ sex + age + BMI + LEU, data = sumfreq_LEU[AA == AAi])
  m0.sum <- summary(m0)
  lin_mod <- m0.sum$coefficients
  lin_mod <- data.table(lin_mod)
  lin_mod$AA <- AAi
  lin_mod$factor <- c("Intercept","sex", "age", "BMI","LEU")
  n <- partial_Spearman(x|LEU~ sex + age + BMI, data = sumfreq_LEU[AA == AAi])
  lin_mod$Estimate_Part_Spear <- n$TS$TB$ts
  lin_mod$pvalue_Part_Spear <- n$TS$TB$pval
  lin_LEU[[k]] <- lin_mod
  k <- k +1
}
corr_LEU <- rbindlist(lin_LEU)
remove(lin_LEU,m0,m0.sum,lin_mod)


corr_LEU_adjust <- corr_LEU[factor == "LEU"]
names(corr_LEU_adjust)[names(corr_LEU_adjust) == "Pr(>|t|)"] <- "pvalue"
corr_LEU_adjust$padjust = p.adjust(corr_LEU_adjust$pvalue, method = "fdr")
corr_LEU_adjust[padjust < 0.05, sign.label1 := "P < 0.05"]
corr_LEU_adjust$padjust_Spear = p.adjust(corr_LEU_adjust$pvalue_Part_Spear, method = "fdr")
corr_LEU_adjust[padjust_Spear < 0.05, sign.label2 := "P < 0.05"]
corr_LEU_adjust <- data.table(corr_LEU_adjust)


#neutrophils
sumfreq_NEUT <- aggregate(info_auxo_DZHK$prop, by=list(subject=info_auxo_DZHK$sample, AA=info_auxo_DZHK$Compound,
                                                      BMI=info_auxo_DZHK$BMI, sex=info_auxo_DZHK$Geschlecht,age=info_auxo_DZHK$Age, NEUT=info_auxo_DZHK$Neutrophils), FUN=sum)

sumfreq_NEUT <- data.table(sumfreq_NEUT)

relAA <- unique(sumfreq_NEUT$AA)
lin_NEUT <- list()
k <- 1
for (AAi in relAA){
  print(AAi)
  m0 <- lm(formula = x ~ sex + age + BMI + NEUT, data = sumfreq_NEUT[AA == AAi])
  m0.sum <- summary(m0)
  lin_mod <- m0.sum$coefficients
  lin_mod <- data.table(lin_mod)
  lin_mod$AA <- AAi
  lin_mod$factor <- c("Intercept","sex", "age", "BMI","NEUT")
  n <- partial_Spearman(x|NEUT~ sex + age + BMI, data = sumfreq_NEUT[AA == AAi])
  lin_mod$Estimate_Part_Spear <- n$TS$TB$ts
  lin_mod$pvalue_Part_Spear <- n$TS$TB$pval
  lin_NEUT[[k]] <- lin_mod
  k <- k +1
}
corr_NEUT <- rbindlist(lin_NEUT)
remove(lin_NEUT,m0,m0.sum,lin_mod)


corr_NEUT_adjust <- corr_NEUT[factor == "NEUT"]
names(corr_NEUT_adjust)[names(corr_NEUT_adjust) == "Pr(>|t|)"] <- "pvalue"
corr_NEUT_adjust$padjust = p.adjust(corr_NEUT_adjust$pvalue, method = "fdr")
corr_NEUT_adjust[padjust < 0.05, sign.label1 := "P < 0.05"]
corr_NEUT_adjust$padjust_Spear = p.adjust(corr_NEUT_adjust$pvalue_Part_Spear, method = "fdr")
corr_NEUT_adjust[padjust_Spear < 0.05, sign.label2 := "P < 0.05"]
corr_NEUT_adjust <- data.table(corr_NEUT_adjust)


#lymphocytes
sumfreq_Lymph <- aggregate(info_auxo_DZHK$prop, by=list(subject=info_auxo_DZHK$sample, AA=info_auxo_DZHK$Compound,
                                                      BMI=info_auxo_DZHK$BMI, sex=info_auxo_DZHK$Geschlecht,age=info_auxo_DZHK$Age, Lymph=info_auxo_DZHK$Lymphocytes), FUN=sum)

sumfreq_Lymph <- data.table(sumfreq_Lymph)

relAA <- unique(sumfreq_Lymph$AA)
lin_Lymph <- list()
k <- 1
for (AAi in relAA){
  print(AAi)
  m0 <- lm(formula = x ~ sex + age + BMI + Lymph, data = sumfreq_Lymph[AA == AAi])
  m0.sum <- summary(m0)
  lin_mod <- m0.sum$coefficients
  lin_mod <- data.table(lin_mod)
  lin_mod$AA <- AAi
  lin_mod$factor <- c("Intercept","sex", "age", "BMI","Lymph")
  n <- partial_Spearman(x|Lymph~ sex + age + BMI, data = sumfreq_Lymph[AA == AAi])
  lin_mod$Estimate_Part_Spear <- n$TS$TB$ts
  lin_mod$pvalue_Part_Spear <- n$TS$TB$pval
  lin_Lymph[[k]] <- lin_mod
  k <- k +1
}
corr_Lymph <- rbindlist(lin_Lymph)
remove(lin_Lymph,m0,m0.sum,lin_mod)


corr_Lymph_adjust <- corr_Lymph[factor == "Lymph"]
names(corr_Lymph_adjust)[names(corr_Lymph_adjust) == "Pr(>|t|)"] <- "pvalue"
corr_Lymph_adjust$padjust = p.adjust(corr_Lymph_adjust$pvalue, method = "fdr")
corr_Lymph_adjust[padjust < 0.05, sign.label1 := "P < 0.05"]
corr_Lymph_adjust$padjust_Spear = p.adjust(corr_Lymph_adjust$pvalue_Part_Spear, method = "fdr")
corr_Lymph_adjust[padjust_Spear < 0.05, sign.label2 := "P < 0.05"]
corr_Lymph_adjust <- data.table(corr_Lymph_adjust)

#monocytes
sumfreq_Mono <- aggregate(info_auxo_DZHK$prop, by=list(subject=info_auxo_DZHK$sample, AA=info_auxo_DZHK$Compound,
                                                      BMI=info_auxo_DZHK$BMI, sex=info_auxo_DZHK$Geschlecht,age=info_auxo_DZHK$Age, Mono=info_auxo_DZHK$Eosinophils), FUN=sum)

sumfreq_Mono <- data.table(sumfreq_Mono)

relAA <- unique(sumfreq_Mono$AA)
lin_Mono <- list()
k <- 1
for (AAi in relAA){
  print(AAi)
  m0 <- lm(formula = x ~ sex + age + BMI + Mono, data = sumfreq_Mono[AA == AAi])
  m0.sum <- summary(m0)
  lin_mod <- m0.sum$coefficients
  lin_mod <- data.table(lin_mod)
  lin_mod$AA <- AAi
  lin_mod$factor <- c("Intercept","sex", "age", "BMI","Mono")
  n <- partial_Spearman(x|Mono~ sex + age + BMI, data = sumfreq_Mono[AA == AAi])
  lin_mod$Estimate_Part_Spear <- n$TS$TB$ts
  lin_mod$pvalue_Part_Spear <- n$TS$TB$pval
  lin_Mono[[k]] <- lin_mod
  k <- k +1
}
corr_Mono <- rbindlist(lin_Mono)
remove(lin_Mono,m0,m0.sum,lin_mod)


corr_Mono_adjust <- corr_Mono[factor == "Mono"]
names(corr_Mono_adjust)[names(corr_Mono_adjust) == "Pr(>|t|)"] <- "pvalue"
corr_Mono_adjust$padjust = p.adjust(corr_Mono_adjust$pvalue, method = "fdr")
corr_Mono_adjust[padjust < 0.05, sign.label1 := "P < 0.05"]
corr_Mono_adjust$padjust_Spear = p.adjust(corr_Mono_adjust$pvalue_Part_Spear, method = "fdr")
corr_Mono_adjust[padjust_Spear < 0.05, sign.label2 := "P < 0.05"]
corr_Mono_adjust <- data.table(corr_Mono_adjust)

####esoinophils
sumfreq_Eos <- aggregate(info_auxo_DZHK$prop, by=list(subject=info_auxo_DZHK$sample, AA=info_auxo_DZHK$Compound,
                                                       BMI=info_auxo_DZHK$BMI, sex=info_auxo_DZHK$Geschlecht,age=info_auxo_DZHK$Age, Eos=info_auxo_DZHK$Eosinophils), FUN=sum)

sumfreq_Eos <- data.table(sumfreq_Eos)

relAA <- unique(sumfreq_Eos$AA)
lin_Eos <- list()
k <- 1
for (AAi in relAA){
  print(AAi)
  m0 <- lm(formula = x ~ sex + age + BMI + Eos, data = sumfreq_Eos[AA == AAi])
  m0.sum <- summary(m0)
  lin_mod <- m0.sum$coefficients
  lin_mod <- data.table(lin_mod)
  lin_mod$AA <- AAi
  lin_mod$factor <- c("Intercept","sex", "age", "BMI","Eos")
  n <- partial_Spearman(x|Eos~ sex + age + BMI, data = sumfreq_Eos[AA == AAi])
  lin_mod$Estimate_Part_Spear <- n$TS$TB$ts
  lin_mod$pvalue_Part_Spear <- n$TS$TB$pval
  lin_Eos[[k]] <- lin_mod
  k <- k +1
}
corr_Eos <- rbindlist(lin_Eos)
remove(lin_Eos,m0,m0.sum,lin_mod)


corr_Eos_adjust <- corr_Eos[factor == "Eos"]
names(corr_Eos_adjust)[names(corr_Eos_adjust) == "Pr(>|t|)"] <- "pvalue"
corr_Eos_adjust$padjust = p.adjust(corr_Eos_adjust$pvalue, method = "fdr")
corr_Eos_adjust[padjust < 0.05, sign.label1 := "P < 0.05"]
corr_Eos_adjust$padjust_Spear = p.adjust(corr_Eos_adjust$pvalue_Part_Spear, method = "fdr")
corr_Eos_adjust[padjust_Spear < 0.05, sign.label2 := "P < 0.05"]
corr_Eos_adjust <- data.table(corr_Eos_adjust)

####Basophils
sumfreq_Baso <- aggregate(info_auxo_DZHK$prop, by=list(subject=info_auxo_DZHK$sample, AA=info_auxo_DZHK$Compound,
                                                       BMI=info_auxo_DZHK$BMI, sex=info_auxo_DZHK$Geschlecht,age=info_auxo_DZHK$Age, Baso=info_auxo_DZHK$Monocytes), FUN=sum)

sumfreq_Baso <- data.table(sumfreq_Baso)

relAA <- unique(sumfreq_Baso$AA)
lin_Baso <- list()
k <- 1
for (AAi in relAA){
  print(AAi)
  m0 <- lm(formula = x ~ sex + age + BMI + Baso, data = sumfreq_Baso[AA == AAi])
  m0.sum <- summary(m0)
  lin_mod <- m0.sum$coefficients
  lin_mod <- data.table(lin_mod)
  lin_mod$AA <- AAi
  lin_mod$factor <- c("Intercept","sex", "age", "BMI","Baso")
  n <- partial_Spearman(x|Baso~ sex + age + BMI, data = sumfreq_Baso[AA == AAi])
  lin_mod$Estimate_Part_Spear <- n$TS$TB$ts
  lin_mod$pvalue_Part_Spear <- n$TS$TB$pval
  lin_Baso[[k]] <- lin_mod
  k <- k +1
}
corr_Baso <- rbindlist(lin_Baso)
remove(lin_Baso,m0,m0.sum,lin_mod)


corr_Baso_adjust <- corr_Baso[factor == "Baso"]
names(corr_Baso_adjust)[names(corr_Baso_adjust) == "Pr(>|t|)"] <- "pvalue"
corr_Baso_adjust$padjust = p.adjust(corr_Baso_adjust$pvalue, method = "fdr")
corr_Baso_adjust[padjust < 0.05, sign.label1 := "P < 0.05"]
corr_Baso_adjust$padjust_Spear = p.adjust(corr_Baso_adjust$pvalue_Part_Spear, method = "fdr")
corr_Baso_adjust[padjust_Spear < 0.05, sign.label2 := "P < 0.05"]
corr_Baso_adjust <- data.table(corr_Baso_adjust)

####Thrombocytes
sumfreq_Thr <- aggregate(info_auxo_DZHK$prop, by=list(subject=info_auxo_DZHK$sample, AA=info_auxo_DZHK$Compound,
                                                       BMI=info_auxo_DZHK$BMI, sex=info_auxo_DZHK$Geschlecht,age=info_auxo_DZHK$Age, Thr=info_auxo_DZHK$Thrombocytes), FUN=sum)

sumfreq_Thr <- data.table(sumfreq_Thr)

relAA <- unique(sumfreq_Thr$AA)
lin_Thr <- list()
k <- 1
for (AAi in relAA){
  print(AAi)
  m0 <- lm(formula = x ~ sex + age + BMI + Thr, data = sumfreq_Thr[AA == AAi])
  m0.sum <- summary(m0)
  lin_mod <- m0.sum$coefficients
  lin_mod <- data.table(lin_mod)
  lin_mod$AA <- AAi
  lin_mod$factor <- c("Intercept","sex", "age", "BMI","Thr")
  n <- partial_Spearman(x|Thr~ sex + age + BMI, data = sumfreq_Thr[AA == AAi])
  lin_mod$Estimate_Part_Spear <- n$TS$TB$ts
  lin_mod$pvalue_Part_Spear <- n$TS$TB$pval
  lin_Thr[[k]] <- lin_mod
  k <- k +1
}
corr_Thr <- rbindlist(lin_Thr)
remove(lin_Thr,m0,m0.sum,lin_mod)


corr_Thr_adjust <- corr_Thr[factor == "Thr"]
names(corr_Thr_adjust)[names(corr_Thr_adjust) == "Pr(>|t|)"] <- "pvalue"
corr_Thr_adjust$padjust = p.adjust(corr_Thr_adjust$pvalue, method = "fdr")
corr_Thr_adjust[padjust < 0.05, sign.label1 := "P < 0.05"]
corr_Thr_adjust$padjust_Spear = p.adjust(corr_Thr_adjust$pvalue_Part_Spear, method = "fdr")
corr_Thr_adjust[padjust_Spear < 0.05, sign.label2 := "P < 0.05"]
corr_Thr_adjust <- data.table(corr_Thr_adjust)

info_auxo_DZHK

####Erythrocytes
sumfreq_Ery <- aggregate(info_auxo_DZHK$prop, by=list(subject=info_auxo_DZHK$sample, AA=info_auxo_DZHK$Compound,
                                                       BMI=info_auxo_DZHK$BMI, sex=info_auxo_DZHK$Geschlecht,age=info_auxo_DZHK$Age, Ery=info_auxo_DZHK$Erythrocytes), FUN=sum)

sumfreq_Ery <- data.table(sumfreq_Ery)

relAA <- unique(sumfreq_Ery$AA)
lin_Ery <- list()
k <- 1
for (AAi in relAA){
  print(AAi)
  m0 <- lm(formula = x ~ sex + age + BMI + Ery, data = sumfreq_Ery[AA == AAi])
  m0.sum <- summary(m0)
  lin_mod <- m0.sum$coefficients
  lin_mod <- data.table(lin_mod)
  lin_mod$AA <- AAi
  lin_mod$factor <- c("Intercept","sex", "age", "BMI","Ery")
  n <- partial_Spearman(x|Ery~ sex + age + BMI, data = sumfreq_Ery[AA == AAi])
  lin_mod$Estimate_Part_Spear <- n$TS$TB$ts
  lin_mod$pvalue_Part_Spear <- n$TS$TB$pval
  lin_Ery[[k]] <- lin_mod
  k <- k +1
}
corr_Ery <- rbindlist(lin_Ery)
remove(lin_Ery,m0,m0.sum,lin_mod)


corr_Ery_adjust <- corr_Ery[factor == "Ery"]
names(corr_Ery_adjust)[names(corr_Ery_adjust) == "Pr(>|t|)"] <- "pvalue"
corr_Ery_adjust$padjust = p.adjust(corr_Ery_adjust$pvalue, method = "fdr")
corr_Ery_adjust[padjust < 0.05, sign.label1 := "P < 0.05"]
corr_Ery_adjust$padjust_Spear = p.adjust(corr_Ery_adjust$pvalue_Part_Spear, method = "fdr")
corr_Ery_adjust[padjust_Spear < 0.05, sign.label2 := "P < 0.05"]
corr_Ery_adjust <- data.table(corr_Ery_adjust)


####Erythrocytes
sumfreq_Hae <- aggregate(info_auxo_DZHK$prop, by=list(subject=info_auxo_DZHK$sample, AA=info_auxo_DZHK$Compound,
                                                      BMI=info_auxo_DZHK$BMI, sex=info_auxo_DZHK$Geschlecht,age=info_auxo_DZHK$Age, Hae=info_auxo_DZHK$Haematocrit), FUN=sum)

sumfreq_Hae <- data.table(sumfreq_Hae)

relAA <- unique(sumfreq_Hae$AA)
lin_Hae <- list()
k <- 1
for (AAi in relAA){
  print(AAi)
  m0 <- lm(formula = x ~ sex + age + BMI + Hae, data = sumfreq_Hae[AA == AAi])
  m0.sum <- summary(m0)
  lin_mod <- m0.sum$coefficients
  lin_mod <- data.table(lin_mod)
  lin_mod$AA <- AAi
  lin_mod$factor <- c("Intercept","sex", "age", "BMI","Hae")
  n <- partial_Spearman(x|Hae~ sex + age + BMI, data = sumfreq_Hae[AA == AAi])
  lin_mod$Estimate_Part_Spear <- n$TS$TB$ts
  lin_mod$pvalue_Part_Spear <- n$TS$TB$pval
  lin_Hae[[k]] <- lin_mod
  k <- k +1
}
corr_Hae <- rbindlist(lin_Hae)
remove(lin_Hae,m0,m0.sum,lin_mod)


corr_Hae_adjust <- corr_Hae[factor == "Hae"]
names(corr_Hae_adjust)[names(corr_Hae_adjust) == "Pr(>|t|)"] <- "pvalue"
corr_Hae_adjust$padjust = p.adjust(corr_Hae_adjust$pvalue, method = "fdr")
corr_Hae_adjust[padjust < 0.05, sign.label1 := "P < 0.05"]
corr_Hae_adjust$padjust_Spear = p.adjust(corr_Hae_adjust$pvalue_Part_Spear, method = "fdr")
corr_Hae_adjust[padjust_Spear < 0.05, sign.label2 := "P < 0.05"]
corr_Hae_adjust <- data.table(corr_Hae_adjust)

###summarize all results in one heatmap

corr_health <- rbind(linear_mod_age, linear_mod_BMI, #linear_mod_sex, 
                             corr_LDL_adjust,corr_LDL_adjust,corr_TG_adjust,corr_LEU_adjust,
                      corr_NEUT_adjust, corr_Lymph_adjust, corr_Mono_adjust,corr_Eos_adjust,corr_Baso_adjust,corr_Thr_adjust,
                      corr_Ery_adjust, corr_Hae_adjust)


corr_health$Estimate_Part_Spear = as.numeric(corr_health$Estimate_Part_Spear)
corr_health[corr_health == "P < 0.05"] <- "Padj <0.05"

corr_health_plot <- ggplot(corr_health, aes(x = factor,y= AA, fill =`Estimate_Part_Spear`))+
  geom_tile() +
  labs(y = "Auxotrophy", x = "Health markers", shape = "")+
  geom_point(aes(shape = sign.label2), size = 1) +
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
  scale_x_discrete("Health markers", labels = c("age" = "Age", "BMI" = "BMI", "hypertens" = "Hypertension", "sex" = "Sex", "TG" = "Triglycerids",
                                                "Baso" = "Basophils", "Eos" = "Eosinophils", "Ery" = "Erythrocytes", "Hae" = "Haemocrit",
                                                "LEU" = "Leucocytes", "Lymph" = "Lymphocytes", "Mono" = "Monocytes", "NEUT" = "Neutrophils", 
                                                "Thr"="Thrombocytes")) +
  theme(panel.background = element_blank()) +
  theme(legend.title = element_text(size=9)) +
  theme(legend.text = element_text(size=7)) +
  labs(fill="Estimate") +
  theme(legend.position = "top",
        legend.justification = 	1) +
  theme(legend.text = element_text(size=7)) +
  theme(panel.grid.major = element_blank())

corr_health_plot + guides(shape = guide_legend(order = 1)) 
corr_health_plot <- corr_health_plot + guides(shape= "none")
corr_health_plot

ggsave("output/plots/dis_health_DHZK.pdf", plot = corr_health_plot,
       width = 7, height = 5)


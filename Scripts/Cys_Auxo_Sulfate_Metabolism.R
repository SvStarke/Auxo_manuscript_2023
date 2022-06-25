###############     Auxotrophies and health parameters/diseases   ##############
data <- fread("/Users/Svenja/Downloads/FOCUS_HRGM_abundancies_2")
View(data)
###ONLY the FoCUs cohort with obese individuals!!!
focus_info <- fread("/Users/Svenja/Downloads/FOCUS_meta_info.csv")
summary
nrow(focus_info[focus_info$BMI>30])
View(focus_all_info)
focus_all_info <- merge(data,focus_info, by.x="subject", by.y="subject")
summarise_all(focus_all_info)
focus_info <- data.table(focus_info)
summarise_all(focus_info["BMI"], funs(nlevels(.), nmiss=sum(is.na(.))))


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

sumfreq_sulf <- aggregate(info_auxo$freq, by=list(subject=info_auxo$subject, AA=info_auxo$Compound, BMI=info_auxo$BMI), FUN=sum)
sumfreq_sulf <- data.table(sumfreq_sulf)
sumfreq_cys <- sumfreq_cys[AA == "Cys"]
sumfreq_cys[, tmp_BMI := BMI]
sumfreq_cys$tmp_BMI <- ifelse(sumfreq_cys$tmp_BMI <30, "non-obese", "obese")
nrow(sumfreq_cys[sumfreq_cys$BMI<30])
nrow(sumfreq_cys[sumfreq_cys$BMI>30])

sumfreq_cys

wilcox <- wilcox.test(sumfreq_cys$x ~ sumfreq_cys$tmp_BMI)

cbPalette <- c("#999999", "#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7")

ggplot(sumfreq_cys, aes(x=tmp_BMI, y=x, fill = tmp_BMI)) +
  geom_boxplot() +
  ylab("Frequence of cys auxotrophic bacteria") +
  xlab("") +
  theme_minimal() +
  theme(legend.position = "none") +
  scale_fill_manual(values = cbPalette) +
  theme(axis.text.x = element_text(size=12, color = "black"),
        axis.text.y = element_text(size=12, , color ="black"),
        axis.title.y = element_text(size=14, face="bold", color ="black", margin = margin(t=0, r = 20, b= 0, l = 0)))  +
  annotate("text",x=1.5, y=0.75, label="Wilcoxon-Test p <0.05") +
  scale_x_discrete(labels=c("non-obese" = "non-obese(n=163)", "obese"="obese(n=246)"))

cbPalette <- c( "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7")

##methionine auxotrophic bacteria
sumfreq_met <- sumfreq_sulf[AA == "Met"]
sumfreq_met[, tmp_BMI := BMI]
sumfreq_met$tmp_BMI <- ifelse(sumfreq_met$tmp_BMI <30, "non-obese", "obese")
nrow(sumfreq_met[sumfreq_cys$BMI<30])
nrow(sumfreq_met[sumfreq_cys$BMI>30])


wilcox <- wilcox.test(sumfreq_met$x ~ sumfreq_met$tmp_BMI)


ggplot(sumfreq_met, aes(x=tmp_BMI, y=x, fill = tmp_BMI)) +
  geom_boxplot() +
  ylab("Frequence of cys auxotrophic bacteria") +
  xlab("") +
  theme_minimal() +
  theme(legend.position = "none") +
  scale_fill_manual(values = cbPalette) +
  theme(axis.text.x = element_text(size=12, color = "black"),
        axis.text.y = element_text(size=12, , color ="black"),
        axis.title.y = element_text(size=14, face="bold", color ="black", margin = margin(t=0, r = 20, b= 0, l = 0)))  +
  annotate("text",x=1.5, y=0.75, label="Wilcoxon-Test p <0.05") +
  scale_x_discrete(labels=c("non-obese" = "non-obese(n=163)", "obese"="obese(n=246)"))

cbPalette <- c( "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7")


###H2S production by cys auxotrophic bacteria

cutoff_prodrate <- 1 # at which mmol/gDW the rate is considered as 'real'production

exchange <- get_exchanges(models)
names <- unique(exchange$name)
fwrite(exchange, file = "exchange.csv")

relCompounds <- c("H2S")

H2S <- exchange[name %in% relCompounds]
View(H2S)
m_gr <- lapply(models, FUN = get_growth)
head(m_gr)
m_growth <- data.table(Genome = names(m_gr),
                       Growth = unlist(m_gr))


#merge the files
relAuxos <- unique(Auxotrophy_2$Compound)

H2S_prod1 <- merge(H2S, m_growth, by.x = "model",
                    by.y = "Genome")
H2S_prod1[, prod_rate := flux / Growth]

H2S_all <- merge(H2S_prod1,info_auxo, by.x="model", by.y="model")
H2S_Prod <- H2S_all[prod_rate <0, prod_rate:=0]
H2S_Prod$realH2S_prod <- H2S_Prod$prod_rate*H2S_Prod$freq
View(H2S_Prod)
sumfreq_H2S <- aggregate(H2S_Prod$realH2S_prod, by=list(subject=H2S_Prod$subject, AA=H2S_Prod$Compound, BMI=H2S_Prod$BMI), FUN=sum)
sumfreq_H2S <- data.table(sumfreq_H2S)
sumfreq_H2S <- sumfreq_H2S[AA == "Cys"]
sumfreq_H2S[, tmp_BMI := BMI]
sumfreq_H2S$tmp_BMI <- ifelse(sumfreq_H2S$tmp_BMI <30, "non-obese", "obese")
View(sumfreq_H2S)
nrow(sumfreq_H2S[sumfreq_H2S$BMI<30])
nrow(sumfreq_H2S[sumfreq_H2S$BMI>30])



wilcox_h2S <- wilcox.test(sumfreq_H2S$x ~ sumfreq_H2S$tmp_BMI)

cpPalette <- c( "#F0E442", "#0072B2", "#D55E00", "#CC79A7")

ggplot(sumfreq_H2S, aes(x=tmp_BMI, y=x, fill = tmp_BMI)) +
  geom_boxplot() +
  ylab("Production of H2S by cys auxotrophic bacteria") +
  xlab("") +
  theme_minimal() +
  theme(legend.position = "none") +
  scale_fill_manual(values = cpPalette) +
  theme(axis.text.x = element_text(size=12, color = "black"),
        axis.text.y = element_text(size=12, , color ="black"),
        axis.title.y = element_text(size=14, face="bold", color ="black", margin = margin(t=0, r = 20, b= 0, l = 0)))  +
  annotate("text",x=1.5, y=0.75, label="Wilcoxon-Test p <0.05") +
  scale_x_discrete(labels=c("non-obese" = "non-obese(n=163)", "obese"="obese(n=246)"))




###############     Auxotrophies and health parameters/diseases   ##############
data <- fread("/mnt/nuuk/2021/HRGM/FOCUS_HRGM_abundancies.csv.gz")

View(data)
###ONLY the FoCUs cohort with obese individuals!!!
focus_info <- fread("/mnt/nuuk/2021/HRGM/FOCUS_meta_info.csv")
summary
nrow(focus_info[focus_info$BMI<30])
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
sumfreq_cys <- sumfreq_sulf[AA == "Cys"]
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
        axis.text.y = element_text(size=12, color ="black"),
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
        axis.text.y = element_text(size=12, color ="black"),
        axis.title.y = element_text(size=14, face="bold", color ="black", margin = margin(t=0, r = 20, b= 0, l = 0)))  +
  annotate("text",x=1.5, y=0.75, label="Wilcoxon-Test p <0.05") +
  scale_x_discrete(labels=c("non-obese" = "non-obese(n=163)", "obese"="obese(n=246)"))

cbPalette <- c( "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7")


###H2S production by cys auxotrophic bacteria
production <- get_produced_metabolites(models)
prod <- data.table(Genome = names(production),
                       ex = unlist(production["ex"]))
test <- rbindlist(production, idcol = TRUE)
View(test)
H2Sprod <- test[ex == "EX_cpd00239_e0"]
#get growth rates
m_gr <- get_growth(models)
head(m_gr)
m_growth <- data.table(Genome = names(m_gr),
                       Growth = unlist(m_gr))

#merge all together with info about auxotrophies
H2Stmp <- merge(m_growth, H2Sprod, by.x="Genome", by.y=".id")
describe(info_auxo$subject)
H2S_all <- merge(H2Stmp,info_auxo, by.x="Genome", by.y="model")
describe(H2S_all$subject)
H2S_all$prodrate <- H2S_all$mtf.flux/H2S_all$Growth
H2S_all$realH2S_prod <- H2S_all$prodrate*H2S_all$freq
View(H2S_all)
sumfreq_H2S <- aggregate(H2S_all$realH2S_prod, by=list(subject=H2S_all$subject, AA=H2S_all$Compound, BMI=H2S_all$BMI), FUN=sum)
sumfreq_H2S <- data.table(sumfreq_H2S)
sumfreq_H2S <- sumfreq_H2S[AA == "Cys"]
sumfreq_H2S[, tmp_BMI := BMI]
sumfreq_H2S$tmp_BMI <- ifelse(sumfreq_H2S$tmp_BMI <30, "non-obese", "obese")
View(sumfreq_H2S)
nrow(sumfreq_H2S[sumfreq_H2S$BMI<30])
nrow(sumfreq_H2S[sumfreq_H2S$BMI>30])
mean(sumfreq_H2S[sumfreq_H2S$tmp_BMI == "non-obese", x])
mean(sumfreq_H2S[sumfreq_H2S$tmp_BMI == "obese", x])

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
        axis.text.y = element_text(size=12,color ="black"),
        axis.title.y = element_text(size=14, face="bold", color ="black", margin = margin(t=0, r = 20, b= 0, l = 0)))  +
  annotate("text",x=1.5, y=0.75, label="Wilcoxon-Test p <0.05") +
  scale_x_discrete(labels=c("non-obese" = "non-obese(n=163)", "obese"="obese(n=246)"))

##### pathways for generation of cys auxotrophic bacteria

testx <- info_auxo[Compound == "Cys"]
describe(info_auxo$subject)
View(testx)
relGenomes <- unique(testx$model)
models <- fetch_model_collection("/mnt/nuuk/2021/HRGM/models/", IDs = relGenomes)
rxns <- fetch_model_collection("/mnt/nuuk/2021/HRGM/models", file.type = "reactions", IDs = relGenomes)
pwys <- fetch_model_collection("/mnt/nuuk/2021/HRGM/models", file.type = "pathways", IDs = relGenomes)
pwys_cov_sulf <- get_pathway_coverage(models,rxns,pwys, pathways.of.interest = c("DISSULFRED-PWY","PWY-5329", 
                                                                                 "PWY-6683", "SO4ASSIM-PWY","SULFMETII-PWY"))

pwys_cov_sulf$tmp_BMI <- ifelse(pwys_cov_sulf$tmp_BMI <30, "non-obese", "obese")



ggplot(pwys_cov_sulf[prediction =="TRUE"], aes(x=rxn.name, fill= pathway))+
         geom_bar()+
  theme(axis.text.x= element_text(angle=90))



###statistical analysis
relEC <- c("2.6.1.3","2.8.1.2")
relEC <- "2.6.1.3"
relPathways <- "PWY-5329"


cys_auxo_x <- list()
k <- 1

for(pi in relPathways) {
  print(pi)
  for (eci in relEC){
    tmp_cys_x <- pwys_cov_sulf[pathway == pi]
    tmp_cys <- tmp_cys_x[ec == eci]
    
    # Fisher Test
    cont_tab <- table(tmp_cys$Cys,tmp_cys$prediction)
      col_order <- c("TRUE", "FALSE")
      cont_tab <- cont_tab[, col_order]
      test_fish <- fisher.test(cont_tab)
      dttmp <- data.table(EC = eci,
                          fisher.p = test_fish$p.value,
                          fisher.or = test_fish$estimate,
                          Pathway = pi,
                          AA = "Cys")
      
      cys_auxo_x[[k]] <- dttmp
      k <- k + 1
  }
}


cys_auxos <- rbindlist(cys_auxo_x)
cys_auxos[, fisher.padj := p.adjust(fisher.p, method = "fdr")]
cys_auxos[, fisher.or.log2 := -log2(fisher.or)]
cys_auxos[fisher.padj < 0.05, sign.label1 := "Padj < 0.05"]


p <- ggplot(cys_auxos, aes(AA, EC,fill = -log2(fisher.or))) +
  geom_tile() +
  geom_point(aes(shape = sign.label1), size = 1) +
  scale_fill_gradient2(high = "#ca0020", mid = "white", low = "#0571b0") +
  scale_shape_manual(values = 8, na.translate = FALSE) +
  facet_grid(Pathway ~., scales = "free_y", space= "free_y") +
  labs(x = "", y = "Enzymatic reaction", shape = "",
       fill = expression(log[2]~'(odds ratio)')) +
  theme_bw() +
  theme(legend.position = "right",
        legend.justification = 	1,
        axis.text.x = element_blank(),
        axis.text.y = element_text(color = "black", size = 8)) +
  theme(axis.title.y=  element_text(size=9, margin = margin(t = 0, r = 20, b= 0, l = 0), face = "bold")) +
  theme(panel.background = element_blank()) +
  theme(strip.text.y = element_text(size=6))

p

ggsave("output/plots/TRP_auxo_Trp_degradation.pdf", plot = p,
       width = , height = 5)

met <- "EX_cpd00239_e0"
#### get reactions with metabolite

get_reactions_with_metabolite(models, met = "EX_cpd00239[e0]")


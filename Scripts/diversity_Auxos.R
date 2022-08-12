######   
sub <- unique(dzhk_relabun$sample)
relAA <- unique(Auxotrophy_2$Compound)
Auxotrophy_2

p <- list()
k <- 1


for (subi in sub) {
  print(subi) 
  for (AAi in relAA) {
    x <- dzhk_relabun[sample == subi]
    y <- Auxotrophy_2[Compound == AAi]
    z <- merge(x,y, by.x = "model", by.y = "Genomes")
    t <- z[Prototrophy == 0]
    p[[k]] <- t
    k <- k +1
  }
}

u <- rbindlist(p) 

sumfreq <- aggregate(u$prop, by=list(sample=u$sample, AA=u$Compound), FUN=sum)

divers_auxos <- merge(dzhk_diversity, sumfreq, by.x="sample", by.y="sample")
divers_auxos <- data.table(divers_auxos)

##### Shannon
d <- list()
k <- 1
AA <- unique(divers_auxos$AA)
D <- c("D.Shannon", "D.Simpson","D.invSimpson", "D.richness")
for(AAi in AA){
  print(AAi) 
  for(di in D) {
  divers <- divers_auxos[AA == AAi]
  cor1 <- cor.test(divers$x, divers[[di]], method = "spearman", exact = FALSE)
  div <- data.table(AA= AAi,
                    Estimate = cor1$estimate,
                    pvalue = cor1$p.value,
                    index = di)
  d[[k]] <- div
  k <- k+1
  }
}
d1 <- rbindlist(d)
d1$padjust = p.adjust(d1$pvalue, method = "fdr")
d1[padjust < 0.05, sign.label1 := "Padj <0.05"]
d1

d1 <- d1[AA != "Gly"]

diversity_auxo <- ggplot(d1, aes(index, AA, fill = Estimate))+
  geom_tile() +
  labs(y = "", x = "Diseases", shape = "") +
  geom_point(aes(shape = sign.label1), size = 1) +
  scale_fill_gradient2(high = "#b2182b", mid = "white", low = "#2166ac") +
  scale_shape_manual(values = 8, na.translate = FALSE) +
  theme_minimal() +
  theme(legend.position = "right",
        legend.justification = 	0.5,
        axis.text.x = element_text(color = "black", angle =0, vjust = 0.2, size = 8),
        axis.text.y = element_text(color = "black", size = 8)) +
  theme(axis.title.y = element_text(margin = margin(t =0, r = 20, b= 0, l = 0))) +
  theme(axis.title.x = element_text(face  = "bold", size = 10, margin = margin(t=20, r = 0, b= 0, l = 0))) +
  theme(axis.title.y = element_text(face  = "bold", size = 10))+
  theme(panel.background = element_blank()) +
  theme(legend.title = element_text(size=9)) +
  labs(fill="Estimate", x = "Diversity indices", y="Auxotrophy") +
  theme(legend.text = element_text(size=7)) +
  theme(panel.grid.major = element_blank()) +
  theme(legend.position = "top",
        legend.justification = 	1)
diversity_auxo + guides(shape = guide_legend(order = 1))
diversity_auxo <- diversity_auxo + guides(shape= "none")
diversity_auxo

ggsave("output/plots/diversity_auxos_DHZK.pdf", plot = diversity_auxo,
       width = 4, height = 6)


###only shannon diversity
diversity_auxo <- ggplot(d1[index == "D.Shannon"], aes(index, AA, fill = Estimate))+
  geom_tile() +
  labs(y = "", x = "Diseases", shape = "") +
  geom_point(aes(shape = sign.label1), size = 1) +
  scale_fill_gradient2(high = "#b2182b", mid = "white", low = "#2166ac") +
  scale_shape_manual(values = 8, na.translate = FALSE) +
  theme_minimal() +
  theme(legend.position = "right",
        legend.justification = 	0.5,
        axis.text.x = element_text(color = "black", angle =0, vjust = 0.2, size = 13),
        axis.text.y = element_text(color = "black", size = 13)) +
  theme(axis.title.y = element_text(margin = margin(t =0, r = 20, b= 0, l = 0))) +
  theme(axis.title.x = element_text(size = 14, margin = margin(t=20, r = 0, b= 0, l = 0))) +
  theme(axis.title.y = element_text(size = 14))+
  theme(panel.background = element_blank()) +
  theme(legend.title = element_text(size=12)) +
  labs(fill="Estimate", x = "Diversity index", y="Auxotrophy") +
  theme(legend.text = element_text(size=10)) +
  theme(panel.grid.major = element_blank()) +
  theme(legend.position = "top",
        legend.justification = 	1)
diversity_auxo + guides(shape = guide_legend(order = 1))
diversity_auxo <- diversity_auxo + guides(shape= "none")
diversity_auxo

ggsave("output/plots/diversity_auxos_DHZK_only_Shannon.pdf", plot = diversity_auxo,
       width = 3, height = 7)



###add diversity results to big figure of healthmarkers and freq of auxotrophies
##preparing data
corr_health_div <- corr_health[, c(5,6,7,8,11,12)]
head(corr_health_div)
d1 <- setcolorder(d1, c("AA", "index","Estimate","pvalue","padjust", "sign.label1"))
d1 <- d1[d1$index == "D.Shannon",]
colnames(corr_health_div) <- c("AA", "index","Estimate", "pvalue", "padjust","sign.label1")

corr_health_div_all <- rbind(corr_health_div, d1)
corr_health_div_all$parameter <- corr_health_div_all$index
corr_health_div_all$parameter[corr_health_div_all$parameter == "age"] <- "Health"
corr_health_div_all$parameter[corr_health_div_all$parameter == "LDL"] <- "Lipids"
corr_health_div_all$parameter[corr_health_div_all$parameter == "BMI"] <- "Health"
corr_health_div_all$parameter[corr_health_div_all$parameter == "D.Shannon"] <- "Diversity"
corr_health_div_all$parameter[corr_health_div_all$parameter == "Eos"] <- "Hematology"
corr_health_div_all$parameter[corr_health_div_all$parameter == "Ery"] <- "Hematology"
corr_health_div_all$parameter[corr_health_div_all$parameter == "Hae"] <- "Hematology"
corr_health_div_all$parameter[corr_health_div_all$parameter == "TG"] <- "Lipids"
corr_health_div_all$parameter[corr_health_div_all$parameter == "LEU"] <- "Hematology"
corr_health_div_all$parameter[corr_health_div_all$parameter == "NEUT"] <- "Hematology"
corr_health_div_all$parameter[corr_health_div_all$parameter == "Lymph"] <- "Hematology"
corr_health_div_all$parameter[corr_health_div_all$parameter == "Mono"] <- "Hematology"
corr_health_div_all$parameter[corr_health_div_all$parameter == "Baso"] <- "Hematology"
corr_health_div_all$parameter[corr_health_div_all$parameter == "Thr"] <- "Hematology"


unique(corr_health_div_all$index)

##visualization
corr_health_div_plot <- ggplot(corr_health_div_all[AA != "Gly"], aes(x = index,y= AA, fill =`Estimate`))+
  geom_tile() +
  labs(y = "Auxotrophy", x = "", shape = "")+
  geom_point(aes(shape = sign.label1), size = 1) +
  scale_fill_gradient2(high = "#b2182b", mid = "white", low = "#2166ac") +
  scale_shape_manual(values = 8, na.translate = FALSE) +
  theme_minimal() +
  theme(legend.position = "right",
        legend.justification = 	0.5,
        axis.text.x = element_text(color = "black", angle = 90, hjust = 1, size = 8),
        axis.text.y = element_text(color = "black", size = 8)) +
  theme(axis.title.y = element_text(size = 10,margin = margin(t =0, r = 10, b= 0, l = 0))) +
  theme(axis.title.x = element_blank()) +
 scale_x_discrete("Health markers", labels = c("age" = "Age", "BMI" = "BMI", " D.Shannon" = "Diversity","hypertens" = "Hypertension", "sex" = "Sex", "TG" = "Triglycerids",
                                                "Baso" = "Basophils", "Eos" = "Eosinophils", "Ery" = "Erythrocytes", "Hae" = "Haematocrit",
                                                "LEU" = "Leucocytes", "Lymph" = "Lymphocytes", "Mono" = "Monocytes", "NEUT" = "Neutrophils", 
                                                "Thr"="Thrombocytes")) +
  theme(panel.background = element_blank()) +
  theme(legend.title = element_text(size=9)) +
  theme(legend.text = element_text(size=7)) +
  facet_grid(. ~ parameter, scales = "free_x", space = "free_x") +
  theme(strip.background.x = element_rect(fill = "grey"))+
  theme(strip.text.x = element_text(size=9)) +
  labs(fill="Estimate") +
  theme(legend.position = "top",
        legend.justification = 	1) +
  theme(legend.text = element_text(size=7)) +
  theme(panel.grid.major = element_blank())

corr_health_div_plot + guides(shape = guide_legend(order = 1)) 
corr_health_div_plot <- corr_health_div_plot + guides(shape= "none")
corr_health_div_plot

ggsave("output/plots/dis_health_DHZK.pdf", plot = corr_health_div_plot,
       width = 8, height = 5)

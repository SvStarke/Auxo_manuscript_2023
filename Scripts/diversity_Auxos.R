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
d1
cortet <- cor.test(divers_auxos_Trp$x, divers_auxos_Trp$D.Shannon, method = "spearman", exact = FALSE)

d1
d1$padjust = p.adjust(d1$pvalue, method = "fdr")
d1[padjust < 0.05, sign.label1 := "P < 0.05"]
d1
divers_auxos_Trp <- divers_auxos[AA == "Trp"]
ggplot(divers_auxos[AA == "Asn"], aes(D.Shannon,x)) +
  geom_point() +
  geom_smooth()

diversity_auxo <- ggplot(d1, aes(index, AA, fill = Estimate))+
  geom_tile() +
  labs(y = "", x = "Diseases", shape = "") +
  geom_point(aes(shape = sign.label1), size = 1) +
  scale_fill_gradient2(high = "#b2182b", mid = "white", low = "#2166ac") +
  scale_shape_manual(values = 8, na.translate = FALSE) +
  theme_minimal() +
 # scale_x_discrete("Diseases", labels = c("IBS" = "IBS", "IBD" = "IBD", "chrond" = "Chr.Diarrhea", "liverdis" = "Liverdisease", "diabetes" = "Diabetes", "parodontitis" = "Periodontitis", "rheumato"="Rheumatism", "hypertens" ="Hypertension")) +
  theme(legend.position = "right",
        legend.justification = 	0.5,
        axis.text.x = element_text(color = "black", angle = 90, vjust = 0.2, size = 8),
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

citation("vegan")

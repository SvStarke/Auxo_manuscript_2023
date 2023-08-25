########### Abundancies of auxotrophic bacteria in the gut ###########


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
sumfreq <- as.data.table(sumfreq)

sumfreq[AA %in% c("Val","Met","Leu","Ile","Trp","Phe","Lys","His","Thr"), is.essential := "essential (human)"]
sumfreq[is.na(is.essential), is.essential := "not essential (human)"]
numb_samples <- dzhk_info$sample
sumfreq <- sumfreq[sumfreq$sample %in% numb_samples]
level_order <- c("His", "Ile", "Leu", "Lys","Met", "Phe", "Thr", "Trp", "Val", "Arg", "Asn", "Chor", "Cys", "Gln", "Pro", "Ser", "Tyr")
#boxplot
ü <- ggplot(sumfreq[AA != "Gly"], aes(x = factor(AA, level = level_order), x*100, fill = is.essential)) +
  geom_boxplot(outlier.shape = NA) +
  #geom_jitter(alpha = 0.05, width = 0.2, color = "black") +
  ylab("Relative abundance\n of auxotrophies [%]")+
  xlab("Amino acids") +
  scale_y_continuous(breaks = c(0,20,40,60,80),limits = c(0,90))+
  theme_bw()+
  theme(axis.line = element_line(size=0.2, colour = "black")) +
  #theme(panel.background = element_rect(fill="white", colour= "white")) +
  theme(axis.title.y = element_text(colour = "black", margin = margin(0,10,0,0)))+
  theme(axis.title.x = element_text(colour = "black", margin = margin(10,0,0,0))) +
  theme(axis.text.x = element_text( colour = "black", hjust = 1,  angle = 45, margin = margin(10,0,0,0))) +
  theme(axis.text.y = element_text( colour = "black")) +
  theme(plot.margin= margin(0.5,0.5,0.5,0.5, "cm")) +
  guides(fill = guide_legend(title = "Essentiality")) +
  scale_fill_manual(values = c("#fdb863", "white")) +
  theme(legend.position = "none",
        panel.grid.minor.y = element_blank(),
        panel.grid.major.x = element_blank()) +
  #theme(legend.text = element_text(size=6)) +
  #theme(legend.title = element_text(size =8, face = "bold")) +
  facet_grid(.~ is.essential, scales = "free_x", space = "free_x") 
  #theme(strip.text.x = element_text(size=8)) 

ü


ggsave("output/plots/boxplot_frequency_gut_DZHK_85.pdf", plot = ü,
       width = 6, height = 4)

##get plot for different compelteness cutoff values 
if(completeness_cuttoff == 80) {
ü_80 <- ü + ggh4x::facet_nested( .~ "80% completeness" +is.essential, scales = "free_x", space = "free_x")
}

ggsave("output/plots/boxplot_frequency_gut_DZHK_80.pdf", plot = ü_80,
       width = 6, height = 4)


if(completeness_cuttoff == 85) {
  ü_85 <- ü + ggh4x::facet_nested( .~ "85% completeness" +is.essential, scales = "free_x", space = "free_x")
}

ggsave("output/plots/boxplot_frequency_gut_DZHK_85.pdf", plot = ü_85,
       width = 6, height = 4)


if(completeness_cuttoff == 90) {
  ü_90 <- ü + ggh4x::facet_nested( .~ "90% completeness" +is.essential, scales = "free_x", space = "free_x")
}

ggsave("output/plots/boxplot_frequency_gut_DZHK_90.pdf", plot = ü_90,
       width = 6, height = 4)


if(completeness_cuttoff == 95) {
  ü_95 <- ü + ggh4x::facet_nested( .~ "95% completeness" +is.essential, scales = "free_x", space = "free_x")
}

ggsave("output/plots/boxplot_frequency_gut_DZHK_95.pdf", plot = ü_95,
       width = 6, height = 4)

fi_compl <- ggarrange(ü_80, ü_85, ü_90, ü_95, 
                 ncol=2,
                 nrow=2, 
                 labels = c("A","B", "C", "D"))
fi_compl
# test2
ggsave("output/plots/Auxo_freq_Completeness_levels.pdf", plot = fi_compl,
       width = 9.5, height = 8)


##get median
sumfreqAA <- aggregate(sumfreq$x, by=list(Aminoacid = sumfreq$AA), FUN=median)



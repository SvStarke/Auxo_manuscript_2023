

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


#boxplot
ü <- ggplot(sumfreq[AA != "Gly"], aes(AA, x*100, fill = is.essential)) +
  geom_boxplot(outlier.shape = NA) +
  #geom_jitter(alpha = 0.05, width = 0.2, color = "black") +
  ylab("Relative abundance of auxotrophies [%]")+
  xlab("Amino acids") +
  theme(axis.line = element_line(size=0.2, colour = "black")) +
  theme(panel.background = element_rect(fill="white", colour= "white")) +
  theme(axis.title.y = element_text(colour = "black", size = 10, margin = margin(0,10,0,0)))+
  theme(axis.title.x = element_text(colour = "black", size = 10, margin = margin(10,0,0,0))) +
  theme(axis.text.x = element_text(size=8, colour = "black", hjust = 1,  angle = 45, margin = margin(10,0,0,0))) +
  theme(axis.text.y = element_text(size = 8, colour = "black")) +
  theme(plot.margin= margin(0.5,0.5,0.5,0.5, "cm")) +
  guides(fill = guide_legend(title = "Essentiality")) +
  scale_fill_manual(values = c("#fdb863", "white")) +
  theme(legend.position = "none") +
  theme(legend.text = element_text(size=8)) +
  theme(legend.title = element_text(size =10, face = "bold")) +
  facet_grid(.~ is.essential, scales = "free_x", space = "free_x") +
  theme(strip.text.x = element_text(size=10))

ü


ggsave("output/plots/boxplot_frequency_gut_DZHK.pdf", plot = ü,
       width = 6, height = 4)

##get median 
sumfreqAA <- aggregate(sumfreq$x, by=list(Aminoacid = sumfreq$AA), FUN=median)



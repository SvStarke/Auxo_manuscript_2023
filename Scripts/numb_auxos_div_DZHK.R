##### number of auxotrophies and diversity
Auxotrophy$count <- rowSums(Auxotrophy == 0)

numb_auxos_DZHK <- merge(Auxotrophy, dzhk_relabun, by.x= "Genomes", by.y="model")
numb_auxos_DZHK <- data.table(numb_auxos_DZHK)
numb_auxos_DZHK[is.na(count), count:= 0]

new <- numb_auxos_DZHK[ ,sum(count*prop), by = sample]
dzhk_div_numb_auxos <- merge(new, dzhk_diversity, by.x="sample", by.y="sample")
numb_samples <- dzhk_info$sample
dzhk_div_numb_auxos <- dzhk_div_numb_auxos[dzhk_div_numb_auxos$sample %in% numb_samples]
describe(dzhk_div_numb_auxos$sample)

##spearman correlation
cor.test(dzhk_div_numb_auxos$D.Shannon, dzhk_div_numb_auxos$V1, method = "spearman", exact = FALSE)

###visualization
div_auxos <- ggplot(dzhk_div_numb_auxos, aes(D.Shannon, V1)) +
  geom_point() +
  geom_smooth(method=lm) +
  theme_bw() +
  xlab("Shannon diversity") +
  ylab("Abundance-weighted average of auxotrophies per MAG") +
  theme(axis.text.x = element_text(colour="black")) +
  theme(axis.text.y = element_text(colour= "black")) +
  theme(axis.title.y = element_text(size = 8, margin = margin(r = 10))) +
  theme(axis.title.x = element_text(size = 8, margin = margin(t = 10))) 
div_auxos
ggsave("output/plots/Numb_Auxos_Div.pdf", plot = div_auxos,
       width = 6, height = 5)


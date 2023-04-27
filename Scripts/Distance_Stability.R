#### Hamming distance and stability

Auxotroph_tmp <- copy(Auxotroph)
Auxotroph_tmp[which(is.na(Auxotroph), arr.ind = T)] <- 1

popgen_relabun <- popgen_relabun[popgen_relabun$model %in% popgen_rel_mags]
#numb_samples <- dzhk_info$sample
#dzhk_relabun <- dzhk_relabun[dzhk_relabun$sample %in% numb_samples]


HammingDT <- list()
for(spli in unique(popgen_relabun$sample)) {
  relmags <- popgen_relabun[sample == spli & prop > 0, model]
  
  Hamming_dist <- data.table(as.matrix(dist(Auxotroph_tmp[relmags,], method = "manhattan")))
  Hamming_dist$MAG1 <- colnames(Hamming_dist)
  Hamming_dist <- melt(Hamming_dist, id.vars = "MAG1", value.name = "Hamming", variable.name = "MAG2")
  
  Hamming_dist <- merge(Hamming_dist,
                        popgen_relabun[sample == spli, .(MAG1 = model, prop1 = prop)],
                        by = "MAG1")
  Hamming_dist <- merge(Hamming_dist,
                        popgen_relabun[sample == spli, .(MAG2 = model, prop2 = prop)],
                        by = "MAG2")
  
  HammingDT[[spli]] <- data.table(sample = spli,
                                  avgHamming = Hamming_dist[, sum(prop1 * prop2 * Hamming)])
  
}

HammingDT <- rbindlist(HammingDT)

###Popgen_stability Scripts

#source("Scripts/Popgen_stability.R")

Stability_hamming <- merge(HammingDT, popgen_bray_numb_auxos2, by.x= "sample", by.y="sample")

cor.test(Stability_hamming$Bray_distance, Stability_hamming$avgHamming, method = "spearman", exact = FALSE)

##visualization
stability_Hamming <- ggplot(Stability_hamming, aes(avgHamming, Bray_distance)) +
  geom_point() +
  geom_smooth(method=lm, se = FALSE) +
  theme_bw() +
  xlab("Average Hamming distance") +
  ylab("Bray Curtis distance") +
  theme(axis.text.x = element_text(colour="black")) +
  theme(axis.text.y = element_text(colour= "black")) +
  theme(axis.title.y = element_text(size = 12, margin = margin(r = 10))) +
  theme(axis.title.x = element_text(size = 12, margin = margin(t = 10))) 
stability_Hamming

ggsave("output/plots/stability_hamming.pdf", plot = stability_Hamming,
       width = 4, height = 4)

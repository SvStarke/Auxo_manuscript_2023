#–––––––––––––––––––––––––––––––––––––––#
# Weighted Hamming-Distances per sample #
#–––––––––––––––––––––––––––––––––––––––#
Auxotroph_tmp <- copy(Auxotroph)
Auxotroph_tmp[which(is.na(Auxotroph), arr.ind = T)] <- 1


HammingDT <- list()
for(spli in unique(dzhk_relabun$sample)) {
  relmags <- dzhk_relabun[sample == spli & prop > 0, model]
  
  Hamming_dist <- data.table(as.matrix(dist(Auxotroph_tmp[relmags,], method = "manhattan")))
  Hamming_dist$MAG1 <- colnames(Hamming_dist)
  Hamming_dist <- melt(Hamming_dist, id.vars = "MAG1", value.name = "Hamming", variable.name = "MAG2")
  
  Hamming_dist <- merge(Hamming_dist,
                        dzhk_relabun[sample == spli, .(MAG1 = model, prop1 = prop)],
                        by = "MAG1")
  Hamming_dist <- merge(Hamming_dist,
                        dzhk_relabun[sample == spli, .(MAG2 = model, prop2 = prop)],
                        by = "MAG2")
  
  HammingDT[[spli]] <- data.table(sample = spli,
                                  avgHamming = Hamming_dist[, sum(prop1 * prop2 * Hamming)])
  
}

HammingDT <- rbindlist(HammingDT)

Div_hamming <- merge(HammingDT, dzhk_diversity, by.x= "sample", by.y="sample")


###### visualization ###########
Hamming_shannon <- ggplot(Div_hamming, aes(D.Shannon, avgHamming)) +
  geom_point() +
  geom_smooth(method=lm) +
  theme_bw() +
  xlab("Shannon diversity") +
  ylab("Average Hamming distance") +
  theme(axis.text.x = element_text(colour="black")) +
  theme(axis.text.y = element_text(colour= "black")) +
  theme(axis.title.y = element_text(size = 10, margin = margin(r = 10))) +
  theme(axis.title.x = element_text(size = 10, margin = margin(t = 10))) 

ggsave("output/plots/hamming_Shannon.pdf", plot = Hamming_shannon,
       width = 6, height = 5)

##correlation
cor.test(Div_hamming$D.Shannon, Div_hamming$avgHamming, method = "spearman", exact = FALSE)


##### Correlation between hamming distance and abunance-weight average of auxotrophies

Auxotrophy$count <- rowSums(Auxotrophy == 0)

numb_auxos_DZHK <- merge(Auxotrophy, dzhk_relabun, by.x= "Genomes", by.y="model")
numb_auxos_DZHK <- data.table(numb_auxos_DZHK)

numb_auxos_DZHK[is.na(count), count:= 0]

new <- numb_auxos_DZHK[ ,sum(count*prop), by = sample]
Hamming_numb_auxos <- merge(new, HammingDT, by.x="sample", by.y="sample")


###### visualization ###########
Hamming_numb <- ggplot(Hamming_numb_auxos, aes(V1, avgHamming)) +
  geom_point() +
  geom_smooth(method=lm) +
  theme_bw() +
  xlab("Abundance-weighted average of auxotrophies per MAG") +
  ylab("Average Hamming distance") +
  theme(axis.text.x = element_text(colour="black")) +
  theme(axis.text.y = element_text(colour= "black")) +
  theme(axis.title.y = element_text(size = 10, margin = margin(r = 10))) +
  theme(axis.title.x = element_text(size = 10, margin = margin(t = 10))) 

ggsave("output/plots/hamming_numb.pdf", plot = Hamming_numb,
       width = 6, height = 5)

##correlation
cor.test(Div_hamming$D.Shannon, Div_hamming$avgHamming, method = "spearman", exact = FALSE)

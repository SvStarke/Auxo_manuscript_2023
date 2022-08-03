###### hamming distance for studying occurence of auxotrohpies together ###
library(cultevo)

Auxo_dzhk <- merge(Auxotrophy, dzhk_relabun, by.x = "Genomes", by.y="model")


Auxo_dzhk <- data.table(Auxo_dzhk)
Auxo_dzhk
rm(hd_auxo)
sample <- unique(Auxo_dzhk$sample)
#genome <- unique(Auxo_dzhk$Genomes)
#genome2 <- unique(Auxo_dzhk$Genomes)
hd_auxo <- list()
k <- 1

for (si in sample) {
  print(si)
  sample_dzhk_tmp <- Auxo_dzhk[sample == si]
  sample_dzhk_tmp <- data.table(sample_dzhk_tmp)
  sample_dzhk <- sample_dzhk_tmp[prop !=0, ]
  sample_dzhk <- sample_dzhk[, c(1,23,24):= NULL]
  t <- hammingdists(sample_dzhk)
  r <- data.table(t)
  mean_hd <- mean(r$t, na.rm = TRUE)
  hd <- data.table(Hamming = mean_hd,
                   sample = si)
  hd_auxo[[k]] <- hd
  k <- k+1
}

p <- rbindlist(hd_auxo)
p
t <- data.table(t)
t

###diversity
tz <- merge(p,dzhk_diversity, by.x="sample", by.y="sample")
tz


merge(divers_auxos, )

### visualization
Hamming_shannon <- ggplot(tz, aes(Hamming, D.Shannon)) +
                   geom_point() +
                    geom_smooth(method=lm) +
                    theme_bw() +
                    ylab("Shannon diversity") +
                    xlab("Average Hamming distance") +
                    theme(axis.text.x = element_text(colour="black")) +
                    theme(axis.text.y = element_text(colour= "black")) +
                    theme(axis.title.y = element_text(size = 10, margin = margin(r = 10))) +
                    theme(axis.title.x = element_text(size = 10, margin = margin(t = 10))) 

ggsave("output/plots/hamming_Shannon.pdf", plot = Hamming_shannon,
       width = 6, height = 5)

##correlation
cor.test(tz$D.Shannon, tz$Hamming, method = "spearman", exact = FALSE)














# for (si in sample) {
#   print(si)
#   sample_dzhk_tmp <- Auxo_dzhk[sample == si]
#   sample_dzhk_tmp <- data.table(sample_dzhk_tmp)
#   sample_dzhk <- sample_dzhk_tmp[prop !=0, ]
#   genome <- unique(sample_dzhk$Genomes)
#   genome2 <- unique(sample_dzhk$Genomes)
#         for(gi in genome) {
#           for(gi2 in genome2){
#             tmp_hd1 <- sample_dzhk[Genomes == gi]
#             tmp_hd2 <- sample_dzhk[Genomes == gi2]
#             tmp_hd1 <- tmp_hd1[1,]
#             tmp_hd2 <- tmp_hd2[1,]
#             tmp_hd1 <- tmp_hd1[, c(1,23,24):= NULL]
#             tmp_hd2 <- tmp_hd2[, c(1,23,24):= NULL]
#             if(gi != gi2) {
#             ##calculating hamming distance
#               tmp_hd <- sum(tmp_hd1 != tmp_hd2)
#               hd <- data.table(Hamming = tmp_hd,
#                                Genome1 = gi,
#                                Genome2 = gi2)
#               
#               hd_auxo[[k]] <- hd
#               k <- k+1
#         }      
#     }
#   }
# }
# 
# p <- rbindlist(hd_auxo)
# 
# for (si in sample) {
#   print(si)
#   sample_dzhk_tmp <- Auxo_dzhk[sample == si]
#   sample_dzhk_tmp <- data.table(sample_dzhk_tmp)
#   sample_dzhk <- sample_dzhk_tmp[prop !=0, ]
#   sample_dzhk <- sample_dzhk[, c(1,23,24):= NULL]
#   t <- hammingdists(Auxo_dzhk_t)
#   r <- data.table(t)
#   mean_hd <- mean(r$t, na.rm = TRUE)
#   hd <- data.table(Hamming = mean_hd,
#                    sample = si)
#   hd_auxo[[k]] <- hd
#   k <- k+1
# }
# 
# 
# 
# 
# 
# 
# 
# 
# 
# t <- hammingdists(Auxo_dzhk_t)
# View(t)
# fwrite(p, file = "hamming_distance_all_samples.csv")
# View(p)
# ###calculate the mean
# mean <- aggregate(p$Hamming, by=list(sample = p$sample), FUN = mean, na.rm = TRUE)
# mean
# 
# ggplot()
# 
# 
# pt <- p[sample== "I12311"]
# mean(pt$Hamming, na.rm = TRUE)
# 
# Auxo_dzhk_test <- data.table(Auxo_dzhk)
# Auxo_dzhk_te <- Auxo_dzhk_test[sample == "I12311"]
# Auxo_dzhk_te <- Auxo_dzhk_te[prop != 0]
# Auxo_dzhk_t <- Auxo_dzhk_te[, c(1,23,24):= NULL]
# 
# t <- expand_grid(Auxo_dzhk_t)
# r <- data.table(t)
# mean(r$t, na.rm = TRUE)
# mean(p$Hamming, na.rm= TRUE)
# new <- if (p$Genome1 == p$Genome2){
#   
# } 
# 
# 
# 
# 
# 
# 
# 
# 
# 
# 
# 
# 
# rm(tmp_hd1, tmp_hd2, tmp_hd, hd, hd_auxo)
# ###calculate average hamming distances per individual
# si <- "I12311"
# for(si in subj) {
#   print(si) 
#   for(gi in genome) {
#     for(gi2 in genome2) {
#       #p1 <- p[Genome1 == gi & Genome2 == gi2]
#       u <- dzhk_relabun[sample == si]
#       u1 <- merge(p, u, by.x="Genome1", by.y="model")
#       u2 <- merge(p, u, by.x="Genome2", by.y="model")
#     }
#   }
# }
# 
# 
# si <- "I12311"
# genome <- unique(Auxo_dzhk$Genomes)
# genome2 <- unique(Auxo_dzhk$Genomes)
# for(si in subj) {
#   print(si) 
#   for(gi in genome) {
#     for(gi2 in genome2) {
#       p1 <- p[Genome1 == gi & Genome2 == gi2]
#       u <- dzhk_relabun[sample == si]
#       u1 <- merge(p1, u, by.x="Genome1", by.y="model")
#       u2 <- merge(p1, u1, by.x="Genome2", by.y="model")
#     }
#   }
# }
# new <- merge(p, dzhk_relabun, by.x = )
# 
# 
# 
# 
# Auxo_dzhk <- data.table(Auxo_dzhk)
# Auxo_dzhk <- Auxo_dzhk[prop !=0]
# subj <- unique(Auxo_dzhk$sample)
# subj <- c("I12373", "I12381")
# genome <- unique(Auxo_dzhk$Genomes)
# genome2 <- unique(Auxo_dzhk$Genomes)
# l <- list()
# k <- 1
# for(si in subj) {
#   print(si) 
#   for(gi in genome) {
#     for(gi2 in genome2){
#     z <- Auxo_dzhk[sample == si]
#     z1 <- z[prop !=0]
#     if(gi != gi2) {
#       
#     tmp_hd1 <- z1[Genomes == gi]
#     tmp_hd2 <- z1[Genomes == gi2]
#     tmp_hd1 <- tmp_hd1[, c(1,23,24):= NULL]
#     tmp_hd2 <- tmp_hd2[, c(1,23,24):= NULL]
#     
#       ##calculating hamming distance
#     if(length(tmp_hd1) != length(tmp_hd2)) {
#       
#      
#       tmp_hd <- sum(tmp_hd1 != tmp_hd2)
#       hd <- data.table(Hamming = tmp_hd,
#                        Genome1 = gi,
#                        Genome2 = gi2,
#                        sample = si) 
#       l[[k]] <- hd
#       k <- k+1
#       }
#      }
#     }
#    
#   }
# }
# p <- rbindlist(l)
# p
# l
# tmp_hd1
# tmp_hd2
# z1[Genomes ]

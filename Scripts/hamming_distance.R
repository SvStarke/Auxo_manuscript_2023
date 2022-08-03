###### hamming distance for studying occurence of auxotrohpies together ###
library(cultevo)

Auxo_dzhk <- merge(Auxotrophy, dzhk_relabun, by.x = "Genomes", by.y="model")

Auxo_dzhk <- data.table(Auxo_dzhk)

sample <- unique(Auxo_dzhk$sample)
sample <- c("I12311","I12312")
#genome <- unique(Auxo_dzhk$Genomes)
#genome2 <- unique(Auxo_dzhk$Genomes)
hd_auxo <- list()
k <- 1

for(si in sample) {
  print(si)
  sample_dzhk_tmp <- Auxo_dzhk[sample == si]
  sample_dzhk <- sample_dzhk_tmp[prop !=0, ]
  genome <- unique(sample_dzhk$Genomes)
  genome2 <- unique(sample_dzhk$Genomes)
        for(gi in genome) {
          for(gi2 in genome2){
            if (gi !=gi2) {
            tmp_hd1 <- sample_dzhk[Genomes == gi]
            tmp_hd2 <- sample_dzhk[Genomes == gi2]
            tmp_hd1 <- tmp_hd1[1,]
            tmp_hd2 <- tmp_hd2[1,]
            tmp_hd1 <- tmp_hd1[, c(1,23,24):= NULL]
            tmp_hd2 <- tmp_hd2[, c(1,23,24):= NULL]
            
            ##calculating hamming distance
              tmp_hd <- sum(tmp_hd1 != tmp_hd2)
              hd <- data.table(Hamming = tmp_hd,
                               Genome1 = gi,
                               Genome2 = gi2,
                               sample = si) 
              hd_auxo[[k]] <- hd
              k <- k+1
        }      
    }
  }
}

p <- rbindlist(hd_auxo)


View(p)
###calculate the mean
mean <- aggregate(p$Hamming, by=list(sample = p$sample), FUN = mean)
mean

pt <- p[sample== "I12311"]
mean(pt$Hamming)

















#rm(tmp_hd1, tmp_hd2, tmp_hd, hd, hd_auxo)
###calculate average hamming distances per individual
si <- "I12311"
for(si in subj) {
  print(si) 
  for(gi in genome) {
    for(gi2 in genome2) {
      #p1 <- p[Genome1 == gi & Genome2 == gi2]
      u <- dzhk_relabun[sample == si]
      u1 <- merge(p, u, by.x="Genome1", by.y="model")
      u2 <- merge(p, u, by.x="Genome2", by.y="model")
    }
  }
}


si <- "I12311"
genome <- unique(Auxo_dzhk$Genomes)
genome2 <- unique(Auxo_dzhk$Genomes)
for(si in subj) {
  print(si) 
  for(gi in genome) {
    for(gi2 in genome2) {
      p1 <- p[Genome1 == gi & Genome2 == gi2]
      u <- dzhk_relabun[sample == si]
      u1 <- merge(p1, u, by.x="Genome1", by.y="model")
      u2 <- merge(p1, u1, by.x="Genome2", by.y="model")
    }
  }
}
new <- merge(p, dzhk_relabun, by.x = )




Auxo_dzhk <- data.table(Auxo_dzhk)
Auxo_dzhk <- Auxo_dzhk[prop !=0]
subj <- unique(Auxo_dzhk$sample)
subj <- c("I12373", "I12381")
genome <- unique(Auxo_dzhk$Genomes)
genome2 <- unique(Auxo_dzhk$Genomes)
l <- list()
k <- 1
for(si in subj) {
  print(si) 
  for(gi in genome) {
    for(gi2 in genome2){
    z <- Auxo_dzhk[sample == si]
    z1 <- z[prop !=0]
    if(gi != gi2) {
      
    tmp_hd1 <- z1[Genomes == gi]
    tmp_hd2 <- z1[Genomes == gi2]
    tmp_hd1 <- tmp_hd1[, c(1,23,24):= NULL]
    tmp_hd2 <- tmp_hd2[, c(1,23,24):= NULL]
    
      ##calculating hamming distance
    if(length(tmp_hd1) != length(tmp_hd2)) {
      
     
      tmp_hd <- sum(tmp_hd1 != tmp_hd2)
      hd <- data.table(Hamming = tmp_hd,
                       Genome1 = gi,
                       Genome2 = gi2,
                       sample = si) 
      l[[k]] <- hd
      k <- k+1
      }
     }
    }
   
  }
}
p <- rbindlist(l)
p
l
tmp_hd1
tmp_hd2
z1[Genomes ]

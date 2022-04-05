library(RaschSampler)

source("Scripts/init_models_filtered.R")

###### Predict auxotrophies

source("Scripts/predict_auxos.R")

View(Auxotroph)
mat <- Auxotroph
View(mat)
auxomat <- na.omit(mat)
nrow(auxomat)
ctrl <- rsctrl(burn_in = 100, n_eff = 1000, step = 16, seed = 0, tfixed = FALSE)
res <- rsampler(auxomat, ctrl)

rsextrmat(res, 2)
permut_mat <- list()
for(i in 2:1001) {
  permut_mat[[i-1]] <- rsextrmat(res, i)
}

remove(permut_mat)
###take very long (n=1000)
a <- list() # for distribution
b <- list() # for median
AAx1 <- colnames(Auxotroph)
AAx <- colnames(Auxotroph)
remove <- c("Ala","Asp","Glu")
AAx <- setdiff(AAx, remove)

for(aa1_ind in 1:(length(AAx)-1)) {
  aa1 <- AAx[aa1_ind]
  print(aa1)
  a[[aa1]] <- list()
  b[[aa1]] <- list()
  for(aa2_ind in (aa1_ind+1):length(AAx)) {
    aa2 <- AAx[aa2_ind]
    print(paste("  ",aa2))
    a[[aa1]][[aa2]] <- c()
    for(i in 1:1000) {
      #print(i)
      tmpmat <- permut_mat[[i]]
      colnames(tmpmat) <- AAx1
      coexist <- apply(tmpmat[,c(aa1,aa2)],
                       1, function(x) x[1] == 0 & x[2] == 0)
      
      a[[aa1]][[aa2]] <- c(a[[aa1]][[aa2]],
             sum(coexist) / length(coexist))
    }
    b[[aa1]][[aa2]] <- median(a[[aa1]][[aa2]])
  }
}

a
library(plyr)
test2 <- rbind.fill(b)
View(test2)
test2 <- rbindlist(b, fill = TRUE)
AAx2 <- AAx
remove1 <- c("Chor")
AAx2 <- setdiff(AAx2, remove1)
test3 <- data.frame(test2)
write.csv(test3,"result_rasch.csv")
test5 <- data.frame(x )
rownames(test3) <- AAx2
rasch_freq <- data.frame(x = unlist(test3))
rasch_freq$AA1 <-  AAx2
rasch_freq$AA2 <- rownames(rasch_freq)
head(rasch_freq)
colnames(rasch_freq) <- c("freq", "AA2", "AA1")
#reorder columns
rasch_freq[,c(3,2,1)]
rasch_freq$AA1 <- gsub("[0-9]+", "", rasch_freq$AA1) 
nrow(rasch_freq)
#test for the number of rows with actual numbers (deleted NAs)
test20 <- rasch_freq[complete.cases(rasch_freq), ]
nrow(test20)

# merge the results of the rasch sampler and the occurence of auxotrophies together
occurence3 <- occurence2[,c(4,5,7)]
rasch_freq <- rasch_freq[,c(3,2,1)]
all_freq <- merge(x = occurence3, y= rasch_freq, by.x= c("A1","A2"), by.y=c("AA1","AA2"))
all_freq <- all_freq[complete.cases(all_freq),]
all_freq[,log2FC := log2(Freq_all/freq)]
all_freq

#####################      visualization    ####################################
t <- ggplot(all_freq, aes(A1,A2, fill = log2FC)) +
  geom_tile(color ="white", lwd= 0.5, linetype = 1.5) +
  scale_fill_gradient2(high = "#ca0020", mid = "white", low = "#0571b0") +
  guides(fill = guide_colourbar(barwidth = 15, barheight = 1, title ="log2FoldChange",
                                label = TRUE, ticks = FALSE)) +
  theme(legend.position = "bottom") +
  theme(axis.text.x = element_text(colour = "Black", size = 10, angle = 45, vjust = 0.5, hjust=0.5)) +
  theme(axis.text.y = element_text(colour = "black", size = 10, vjust = 0.5, hjust=0.5)) +
  theme(axis.title.x = element_text(colour = "Black", face = "bold", size = 10, margin = margin(5,0,0,0))) +
  theme(axis.title.y = element_text(colour = "Black", face = "bold", size = 10, margin = margin(0,5,0,0))) +
  theme(legend.title = element_text(colour = "black", size = 10, face = "bold", vjust = 1, margin = margin(0,10,1,0))) +
  theme(legend.text = element_text(size=6)) +
  xlab("Amino acid auxotrophy 1") +
  ylab("Amino acid auxotrophy 2") + 
  theme(panel.background = element_blank()) +
  theme(plot.margin= margin(0.5,0.5,0.5,0.5, "cm"))
t

ggsave("output/plots/Rasch_Sampler.pdf", plot = t,
       width = 5, height = 5)


##### statistical analysis
all_results <- rbindlist(a, fill = TRUE)
remove(tmp_wilcox_rasch)
tmp_wilcox_rasch <- list()
k <- 1
#get third element of the first list 
for (i1 in  1:(length(AAx)-1)) {
  AA1 <- AAx[i1]
  print(AA1) 
  for(i2 in (i1+1):length(AAx)) {
    AA2 <- AAx[i2]
    tmp_occu <- a[[c(AA1,AA2)]] 
    #tmp_occu1 <- data.table(tmp_occu)
    #colnames(tmp_occu1) <- "perc"
    mu_occu <- occurence3[occurence3$A1 == AA1 & occurence3$A2 == AA2, Freq_all]
    res <- wilcox.test(tmp_occu, mu = mu_occu)
    tmp_wilcox <- data.table(A1 = AA1, A2 = AA2,
                             p.value = res$p.value, 
                             obs_freq = mu_occu, 
                             exp_freq_median = median(tmp_occu))
    tmp_wilcox_rasch[[k]] <- tmp_wilcox
    k <- k+1
  }
}
new_table <- rbindlist(tmp_wilcox_rasch)
new_table[, padj := p.adjust(p.value, method = "fdr")]
new_table[padj < 0.05, sign.label1 := "Padj < 0.05"]
new_table

# View(new_table)
# View(tmp_occu1)
# i1 <- "Ser"
# i2 <- "Chor"
# 
# a[[c("Ser","Chor")]]
# a[[c(i1,i2)]]
# nrow(a[["Met"]])
# tet <- a$Met
# nrow(a[a$Met])
# res <- wilcox.test(all)
# View(all_results)
# occurence4 <- occurence3
# occurence
# occurence5$A2 <- NULL
# occurence5$Freq_all <- NULL
# test11 <- c(i1,i2)
# colnames(occurence5) <- test11


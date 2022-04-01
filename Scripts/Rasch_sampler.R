library(RaschSampler)

source("Scripts/init_models_filtered.R")

###### Predict auxotrophies

source("Scripts/predict_auxos.R")

View(Auxotroph)
mat <- Auxotroph
View(mat)
auxomat <- na.omit(mat)
nrow(auxomat)
ctrl <- rsctrl(burn_in = 100, n_eff = 100, step = 16, seed = 0, tfixed = FALSE)
res <- rsampler(auxomat, ctrl)

rsextrmat(res, 2)

a <- list() # for distribution
b <- list() # for median
AAx1 <- colnames(Auxotroph)
remove <- c("Ala","Asp","Glu")
AAx <- setdiff(AAx, remove)
for(aa1 in AAx) {
  print(aa1)
  a[[aa1]] <- list()
  b[[aa1]] <- list()
  for(aa2 in AAx) {
    print(paste("  ",aa2))
    a[[aa1]][[aa2]] <- c()
    for(i in 2:11) {
      tmpmat <- rsextrmat(res, i)
      colnames(tmpmat) <- AAx1
      
      
      coexist <- apply(tmpmat[,c(aa1,aa2)],
                       1, function(x) x[1] == 0 & x[2] == 0)
      
      a[[aa1]][[aa2]] <- c(a[[aa1]][[aa2]],
             sum(coexist) / length(coexist))
    }
    b[[aa1]][[aa2]] <- median(a[[aa1]][[aa2]])
  }
}
remove(a)
test2 <- rbindlist(b)
rasch_freq <- data.frame(x=unlist(test2))
rasch_freq$AA2 <- AAx
rasch_freq$AA1 <- rownames(rasch_freq)
head(rasch_freq)
colnames(rasch_freq) <- c("freq", "AA2", "AA1")
#reorder columns
rasch_freq[,c(3,2,1)]
colnames(new_dataframe) <- c("AA1")
rasch_freq$AA1 <- gsub("[0-9]+", "", rasch_freq$AA1) 

ggplot(rasch_freq, aes(AA1,AA2, fill = freq)) +
  geom_tile()


rasch_freq1 <- as.data.table(rasch_freq)
occurence3 <- occurence2[,c(4,5,7)]
rasch_freq1 <- rasch_freq1[,c(3,2,1)]
all_freq <- merge(x = occurence3, y= rasch_freq1, by.x= c("A1","A2"), by.y=c("AA1","AA2"))
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

o



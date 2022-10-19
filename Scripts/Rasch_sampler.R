library(RaschSampler)

source("Scripts/init_models_filtered.R")

###### Predict auxotrophies

source("Scripts/predict_auxos.R")


mat <- Auxotroph
auxomat <- na.omit(mat)
nrow(auxomat)
ctrl <- rsctrl(burn_in = 100, n_eff = 1000, step = 16, seed = 0, tfixed = FALSE)
res <- rsampler(auxomat, ctrl)

rsextrmat(res, 2)
permut_mat <- list()
for(i in 2:1001) {
  permut_mat[[i-1]] <- rsextrmat(res, i)
}


###take very long (n=1000)
a <- list() # for distribution
b <- list() # for median
AAx1 <- colnames(Auxotroph)
AAx <- colnames(Auxotroph)
AAx <- sort(AAx)
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

occurence3 <- occurence2[,c(4,5,7)]

#########  new evaluation of the pvlaue  ###########

nsim <- 1000



tmp_pvalue_rasch <- list()
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
    pvalue <- (sum(tmp_occu >= abs(mu_occu)) + sum(tmp_occu <= -abs(mu_occu)))/nsim
    tmp_p <- data.table(A1 = AA1, A2 = AA2,
                             p.value = pvalue, 
                             obs_freq = mu_occu, 
                             exp_freq_median = median(tmp_occu))
    tmp_pvalue_rasch[[k]] <- tmp_p
    k <- k+1
  }
}

new_table_p <- rbindlist(tmp_pvalue_rasch)
new_table_p[, padj := p.adjust(p.value, method = "fdr")]
new_table_p[padj < 0.05, sign.label1 := "Padj < 0.05"]
new_table_p[,log2FC := log2(obs_freq/exp_freq_median)]
new_table_p

###exclude glycine
new_table_p <- new_table_p[A1 !="Gly"]
new_table_p <- new_table_p[A2 !="Gly"]

fwrite(new_table_p, file = "/home/svenja/Documents/Rasch_Sampler.csv")

#####################      visualization    ####################################
t <- ggplot(new_table_p, aes(A1,A2, fill = log2FC)) +
  geom_tile(color ="white", lwd= 0.5, linetype = 1.5) +
  scale_fill_gradient2(high = "#ca0020", mid = "white", low = "#0571b0") +
  guides(fill = guide_colourbar(barwidth = 12, barheight = 0.8, title ="log2FoldChange",
                                label = TRUE, ticks = FALSE)) +
  geom_point(aes(shape = sign.label1), size = 1, show.legend = FALSE) +
  theme(legend.position = "bottom") +
  scale_shape_manual(values = 8, na.translate = FALSE) +
  theme(axis.text.x = element_text(colour = "Black", size = 10, angle = 45, vjust = 0.5, hjust=0.5)) +
  theme(axis.text.y = element_text(colour = "black", size = 10, vjust = 0.5, hjust=0.5)) +
  theme(axis.title.x = element_text(colour = "Black", size = 10, margin = margin(5,0,0,0))) +
  theme(axis.title.y = element_text(colour = "Black", size = 10, margin = margin(0,5,0,0))) +
  theme(legend.title = element_text(colour = "black", size = 8, vjust = 1, margin = margin(0,10,1,0))) +
  theme(legend.text = element_text(size=8)) +
  xlab("Amino acid auxotrophy 1") +
  ylab("Amino acid auxotrophy 2") + 
  labs(shape = "") +
  scale_x_discrete(position = "top") +
  theme(panel.background = element_blank()) +
  theme(plot.margin= margin(0.5,0.5,0.5,0.5, "cm"))
t

ggsave("output/plots/Rasch_Sampler.pdf", plot = t,
       width = 5.5, height = 5)




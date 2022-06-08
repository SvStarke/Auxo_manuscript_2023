######  change in the abundance of amino acid auxotrophies  ###########

data <- fread("/mnt/nuuk/2021/HRGM/FOCUS_HRGM_abundancies.csv.gz")
View(data)
data[,.(subject_samples = .N), by = .(subject,model,freq)]

data_BL <- data[focus.call == "BL", unique(subject)]
data_F1 <- data[focus.call == "F1", unique(subject)]
data_F2 <- data[focus.call == "F2", unique(subject)]

x <- intersect(data_BL, intersect(data_F1,data_F2))

data_new <- data[subject %in% x]


##add information about auxotrophies
Auxotrophy_2[,c(4:17)] <- NULL

samp <- unique(data_new$sample)
relAA <- unique(Auxotrophy_2$Compound)
Auxotrophy_2

p <- list()
k <- 1

for (sampi in samp) {
  print(sampi) 
  for (AAi in relAA) {
    x <- data_new[sample == sampi]
    y <- Auxotrophy_2[Compound == AAi]
    z <- merge(x,y, by.x = "model", by.y = "Genomes")
    t <- z[Prototrophy == 0]
    p[[k]] <- t
    k <- k +1
  }
}

u <- rbindlist(p) 
View(u)

####sum per sample
sumfreq_type_Focus <- aggregate(u$freq, by=list(sample=u$sample, AA=u$Compound, FoCus_call=u$focus.call, subject=u$subject), FUN=sum)
sumfreq_type_Focus <- data.table(sumfreq_type_Focus)
View(sumfreq_type_Focus)


###visualization
auxos_time <- ggplot(sumfreq_type_Focus[AA == "Leu"], aes(x = subject, y=x, group = FoCus_call)) +
  geom_line(aes(color=FoCus_call)) +
  geom_point(aes(color=FoCus_call)) +
  theme(legend.position = "right",
        legend.justification = 	0.5,
        axis.text.x = element_text(color = "black", angle = 90, vjust = 0.2, size = 8),
        axis.text.y = element_text(color = "black", size = 8)) +
  #theme(axis.text.x = element_blank()) +
  #theme(axis.ticks.x = element_blank() ) +
  theme(panel.background = element_blank()) 
auxos_time + annotate("text", x=60, y=0.9, label = "Leucine, padj <0.05")

ggboxplot(sumfreq_type_Focus[AA=="Met"], x="FoCus_call", y="x", add="jitter")

View(sumfreq_type_Focus)

####statistical analysis
library(psych)
install.packages("datarium")
data("selfesteem", package="datarium")
AA <- unique(sumfreq_type_Focus$AA)
f<- list()
k <- 1

for(AAi in AA) {
  print(AAi)
    tmp_Focus <- sumfreq_type_Focus[AA == AAi]
    tmp_Focus <- data.table(tmp_Focus)
    fried <- friedman.test(y=tmp_Focus$x, groups=tmp_Focus$FoCus_call, blocks=tmp_Focus$subject)
    tmp_fried <- data.table(AA = AAi,
                             chisquare = fried$statistic,
                             pvalue = fried$p.value)
    f[[k]] <- tmp_fried
    k <- k +1
}

Fried_Focus <- rbindlist(f)

### pvalue adjustment
Fried_Focus$padjust = p.adjust(Fried_Focus$pvalue, method = "fdr")
Fried_Focus[padjust < 0.05, sign.label1 := "P < 0.05"]


### visualization

auxos_time_Val <- ggplot(sumfreq_type_Focus[AA == "Val"], aes(x = subject, y=x, group = FoCus_call)) +
  geom_line(aes(color=FoCus_call)) +
  geom_point(aes(color=FoCus_call)) +
  theme(legend.position = "right",
        legend.justification = 	0.5,
        axis.text.x = element_text(color = "black", angle = 90, vjust = 0.2, size = 8),
        axis.text.y = element_text(color = "black", size = 8)) +
  theme(axis.text.x = element_blank()) +
  theme(axis.ticks.x = element_blank() ) +
  theme(panel.background = element_blank()) 
auxos_time_Val + annotate("text", x=60, y=0.9, label = "Leucine, padj <0.05")




### pairwise wilcox test
c <- list()
k <- 1
relAA <- unique(sumfreq_type_Focus$AA)

for (AAi in relAA) {
  print(AAi)
  wilcox_auxo <- sumfreq_type_Focus[AA == AAi]
  pair_wilcox_all <- pairwise.wilcox.test(wilcox_auxo$x, wilcox_auxo$FoCus_call, p.adjust = "fdr", exact = FALSE)
  pairwise_wilcox <- pair_wilcox_all$p.value
  test <- data.frame(pairwise_wilcox)
  test$AA <- AAi
  test$BL <- ifelse(test$BL < 0.05, "*", "N.S")
  test$F1 <- ifelse(test$F1 < 0.05, "*", "N.S")
  test$F2 <- ifelse(test$F2 < 0.05, "*", "N.S")
  c[[k]] <- test
  k <- k+1
}

c <- list()
k <- 1
relAA <- unique(sumfreq_type_Focus$AA)

for(AAi in relAA) {
  print(AAi)
  wilcox_auxo <- sumfreq_type_Focus[AA == AAi]
  pwc <- wilcox_auxo %>% 
    wilcox_test(x ~ FoCus_call, paired = TRUE, p.adjust.method = "fdr")
  pwc$AA <- AAi
  pwc_new <- data.frame(pwc)
  c[[k]] <- pwc_new
  k <- k+1
}

Wil_Focus <- rbindlist(c)


ggplot(Wil_Focus, aes(group1,group2))
###alpha und beta-diversity



#changing healthparameters




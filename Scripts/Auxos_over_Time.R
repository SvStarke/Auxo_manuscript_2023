######  change in the abundance of amino acid auxotrophies  ###########

data <- fread("/mnt/nuuk/2021/HRGM/FOCUS_HRGM_abundancies.csv.gz")
View(data)
data[,.(subject_samples = .N), by = .(subject,model,freq)]

data_BL <- data[focus.call == "BL", unique(subject)]
data_F1 <- data[focus.call == "F1", unique(subject)]
data_F2 <- data[focus.call == "F2", unique(subject)]

x <- intersect(data_BL, intersect(data_F1,data_F2))
View(x)
data_new <- data[subject %in% x]
View(data_new)

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

sumfreq_type_Focus[order(sumfreq_type_Focus$x, sumfreq_type_Focus$FoCus_call, decreasing = TRUE), ]
reorder(sumfreq_type_Focus, x)
unique(sumfreq_type_Focus$FoCus_call)

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
##val
auxos_time_Val <- ggplot(sumfreq_type_Focus[AA == "Val"], aes(x = subject, y=x, group = FoCus_call)) +
  geom_line(aes(color=FoCus_call), size =0.2) +
  geom_point(aes(color=FoCus_call), size =1) +
  theme(legend.position = "right",
        legend.justification = 	0.5,
        axis.text.x = element_text(color = "black", angle = 90, vjust = 0.2, size = 8),
        axis.text.y = element_text(color = "black", size = 8)) +
  theme(axis.text.x = element_blank()) +
  theme(axis.ticks.x = element_blank() ) +
  theme(panel.background = element_blank()) +
  theme(axis.title.y = element_text(size=6)) +
  theme(axis.title.x = element_text(size=6)) +
  labs(y="Frequence of\n auxotrophies",x = "") +
  theme(axis.line.x = element_line(color="black", size = 0.4),
        axis.line.y = element_line(color = "black", size =0.4)) +
  annotate("text", x=60, y=0.9, size =2, label = "Valine, padj <0.05") +
  scale_color_manual(values = c("BL" = "#E69F00", "F1" = "#56B4E9", "F2" = "#009E73"))

##leu
auxos_time_Leu <- ggplot(sumfreq_type_Focus[AA == "Leu"], aes(x = subject, y=x, group = FoCus_call)) +
  geom_line(aes(color=FoCus_call), size =0.2) +
  geom_point(aes(color=FoCus_call), size =1) +
  theme(legend.position = "right",
        legend.justification = 	0.5,
        axis.text.x = element_text(color = "black", angle = 90, vjust = 0.2, size = 8),
        axis.text.y = element_text(color = "black", size = 8)) +
  theme(axis.text.x = element_blank()) +
  theme(axis.ticks.x = element_blank() ) +
  theme(axis.title.y = element_text(size=6)) +
  theme(axis.title.x = element_text(size=6)) +
  theme(panel.background = element_blank()) +
  labs(y="Frequence of\n auxotrophies",x = "") +
  theme(axis.line.x = element_line(color="black", size = 0.4),
        axis.line.y = element_line(color = "black", size =0.4)) +
  annotate("text", x=60, y=0.9, size =2, label = "Leucine, padj <0.05")+
  scale_color_manual(values = c("BL" = "#E69F00", "F1" = "#56B4E9", "F2" = "#009E73"))


#phe
auxos_time_Phe <- ggplot(sumfreq_type_Focus[AA == "Phe"], aes(x = subject, y=x, group = FoCus_call)) +
  geom_line(aes(color=FoCus_call), size =0.2) +
  geom_point(aes(color=FoCus_call), size =1) +
  theme(legend.position = "right",
        legend.justification = 	0.5,
        axis.text.x = element_text(color = "black", angle = 90, vjust = 0.2, size = 8),
        axis.text.y = element_text(color = "black", size = 8)) +
  theme(axis.text.x = element_blank()) +
  theme(axis.ticks.x = element_blank() ) +
  theme(axis.title.y = element_text(size=6)) +
  theme(axis.title.x = element_text(size=6)) +
  theme(panel.background = element_blank()) +
  labs(y="Frequence of\n auxotrophies",x = "") +
  theme(axis.line.x = element_line(color="black", size = 0.4),
        axis.line.y = element_line(color = "black", size =0.4)) +
  annotate("text", x=60, y=0.9, size =2, label = "Phenylalanine, padj <0.05") +
  scale_color_manual(values = c("BL" = "#E69F00", "F1" = "#56B4E9", "F2" = "#009E73"))

#Tyr
auxos_time_Tyr <- ggplot(sumfreq_type_Focus[AA == "Tyr"], aes(x = subject, y=x, group = FoCus_call)) +
  geom_line(aes(color=FoCus_call), size =0.2) +
  geom_point(aes(color=FoCus_call), size =1) +
  theme(legend.position = "right",
        legend.justification = 	0.5,
        axis.text.x = element_text(color = "black", angle = 90, vjust = 0.2, size = 8),
        axis.text.y = element_text(color = "black", size = 8)) +
  theme(axis.text.x = element_blank()) +
  theme(axis.ticks.x = element_blank() ) +
  theme(axis.title.y = element_text(size=6)) +
  theme(axis.title.x = element_text(size=6)) +
  theme(panel.background = element_blank()) +
  labs(y="Frequence of\n auxotrophies",x = "") +
  theme(axis.line.x = element_line(color="black", size = 0.4),
        axis.line.y = element_line(color = "black", size =0.4)) +
annotate("text", x=60, y=0.9, size =2, label = "Tyrosine, padj <0.05") +
  scale_color_manual(values = c("BL" = "#E69F00", "F1" = "#56B4E9", "F2" = "#009E73"))

##Gln
auxos_time_Gln <- ggplot(sumfreq_type_Focus[AA == "Gln"], aes(x = subject, y=x, group = FoCus_call)) +
  geom_line(aes(color=FoCus_call), size =0.2) +
  geom_point(aes(color=FoCus_call), size =1) +
  theme(legend.position = "right",
        legend.justification = 	0.5,
        axis.text.x = element_text(color = "black", angle = 90, vjust = 0.2, size = 8),
        axis.text.y = element_text(color = "black", size = 8)) +
  theme(axis.text.x = element_blank()) +
  theme(axis.ticks.x = element_blank() ) +
  theme(axis.title.y = element_text(size=6)) +
  theme(axis.title.x = element_text(size=6)) +
  theme(panel.background = element_blank()) +
  labs(y="Frequence of\n auxotrophies",x = "") +
  theme(axis.line.x = element_line(color="black", size = 0.4),
        axis.line.y = element_line(color = "black", size =0.4)) +
  annotate("text", x=60, y=0.9, size =2, label = "Glutamine, padj <0.05") +
  scale_color_manual(values = c("BL" = "#E69F00", "F1" = "#56B4E9", "F2" = "#009E73"))

##Ser
auxos_time_Ser <- ggplot(sumfreq_type_Focus[AA == "Ser"], aes(x = subject, y=x, group = FoCus_call)) +
  geom_line(aes(color=FoCus_call), size =0.2) +
  geom_point(aes(color=FoCus_call), size =1) +
  theme(legend.position = "right",
        legend.justification = 	0.5,
        axis.text.x = element_text(color = "black", angle = 90, vjust = 0.2, size = 8),
        axis.text.y = element_text(color = "black", size = 8)) +
  theme(axis.text.x = element_blank()) +
  theme(axis.ticks.x = element_blank() ) +
  theme(axis.title.y = element_text(size=6)) +
  theme(axis.title.x = element_text(size=6)) +
  theme(panel.background = element_blank()) +
  labs(y="Frequence of\n auxotrophies",x = "") +
  theme(axis.line.x = element_line(color="black", size = 0.4),
        axis.line.y = element_line(color = "black", size =0.4)) +
  annotate("text", x=60, y=0.9, size =2, label = "Serine, padj <0.05") +
  scale_color_manual(values = c("BL" = "#E69F00", "F1" = "#56B4E9", "F2" = "#009E73"))

##Cys
auxos_time_Cys <- ggplot(sumfreq_type_Focus[AA == "Cys"], aes(x = subject, y=x, group = FoCus_call)) +
  geom_line(aes(color=FoCus_call), size =0.2) +
  geom_point(aes(color=FoCus_call), size =1) +
  theme(legend.position = "right",
        legend.justification = 	0.5,
        axis.text.x = element_text(color = "black", angle = 90, vjust = 0.2, size = 8),
        axis.text.y = element_text(color = "black", size = 8)) +
  theme(axis.text.x = element_blank()) +
  theme(axis.ticks.x = element_blank() ) +
  theme(axis.title.y = element_text(size=6)) +
  theme(axis.title.x = element_text(size=6)) +
  theme(panel.background = element_blank()) +
  labs(y="Frequence of\n auxotrophies",x = "") +
  theme(axis.line.x = element_line(color="black", size = 0.4),
        axis.line.y = element_line(color = "black", size =0.4)) +
  annotate("text", x=60, y=0.9, size =2, label = "Cysteine, padj <0.05") +
  scale_color_manual(values = c("BL" = "#E69F00", "F1" = "#56B4E9", "F2" = "#009E73"))

##Asn
auxos_time_Asn <- ggplot(sumfreq_type_Focus[AA == "Asn"], aes(x = subject, y=x, group = FoCus_call)) +
  geom_line(aes(color=FoCus_call), size =0.2) +
  geom_point(aes(color=FoCus_call), size =1) +
  theme(legend.position = "right",
        legend.justification = 	0.5,
        axis.text.x = element_text(color = "black", angle = 90, vjust = 0.2, size = 8),
        axis.text.y = element_text(color = "black", size = 8)) +
  theme(axis.text.x = element_blank()) +
  theme(axis.ticks.x = element_blank() ) +
  theme(axis.title.y = element_text(size=6)) +
  theme(axis.title.x = element_text(size=6)) +
  theme(panel.background = element_blank()) +
  labs(y="Frequence of\n auxotrophies",x = "") +
  theme(axis.line.x = element_line(color="black", size = 0.4),
        axis.line.y = element_line(color = "black", size =0.4)) +
  annotate("text", x=60, y=0.9, size =2, label = "Asparagine, padj <0.05") +
  scale_color_manual(values = c("BL" = "#E69F00", "F1" = "#56B4E9", "F2" = "#009E73"))


###ggplot
auxos_over_time <- ggarrange(auxos_time_Val, auxos_time_Leu, auxos_time_Phe,auxos_time_Tyr,auxos_time_Gln,auxos_time_Ser,auxos_time_Cys,auxos_time_Asn,
                 ncol=2,
                 nrow=4, heights = c(1,1,1,1), widths= c(1,1),
                 labels = "AUTO", common.legend = TRUE, legend = c("bottom"))
auxos_over_time
ggsave("output/plots/Auxos_Over_time_only_signf_V2.pdf", plot = auxos_over_time,
       width = 10, height = 6)


### pairwise wilcox test
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

###alpha und beta-diversity



#changing healthparameters




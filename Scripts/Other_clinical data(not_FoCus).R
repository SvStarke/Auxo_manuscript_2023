############      data about patients with chronic inflammations    ############

frequencies <- fread("/mnt/nuuk/2021/HRGM/TrypCID_16S_abundancies.csv")
metainfo <- fread("/mnt/nuuk/2021/HRGM/TrypCID_16S_metaInfo.csv")

gut_dis <- merge(frequencies, metainfo, by.x= "sample", by.y="sample")
View(gut_dis)
Auxotrophy_2[,c(4:16)] <- NULL


sap <- unique(gut_dis$sample)
relAA <- unique(Auxotrophy_2$Compound)
Auxotrophy_2
h <- list()
k <- 1

for (sapi in sap) {
  print(sapi) 
  for (AAi in relAA) {
    x <- gut_dis[sample == sapi]
    y <- Auxotrophy_2[Compound == AAi]
    z <- merge(x,y, by.x = "model", by.y = "Genomes")
    t <- z[Prototrophy == 0]
    h[[k]] <- t
    k <- k +1
  }
}
gut_dis_auxos <- rbindlist(h) 
View(gut_dis_auxos)
remove(h,x,y,z,t)

sumfreq_all_diseases <- aggregate(gut_dis_auxos$freq, by=list(sample=gut_dis_auxos$sample, AA=gut_dis_auxos$Compound, diseases = gut_dis_auxos$short_desc), FUN=sum)
View(sumfreq_all_diseases)
sumfreq_all_diseases <- data.table(sumfreq_all_diseases)



#####################        visualization        ##############################
#####################           Psoriasis         ##############################
sumfreq_psoriasis <- sumfreq_all_diseases[diseases == "Psoriasis"]
median_psoriasis <- aggregate(sumfreq_psoriasis$x, by = list(Aminoacid =sumfreq_psoriasis$AA), FUN = median)
a <- ggplot(median_psoriasis, aes(Aminoacid, x)) +
  geom_bar(stat="identity", fill = "#7E0018DD",colour = "black", width= 1) +
  ylab("Abundance") +
  xlab("Amino acid auxotrophy") +
  theme(axis.line = element_line(size=0.2, colour = "black")) +
  theme(panel.background = element_rect(fill="white", colour= "white")) +
  theme(axis.title.y = element_text(colour = "black", size = 10, face = "bold", margin = margin(0,10,0,0)))+
  theme(axis.title.x = element_text(colour = "black", size = 10, face = "bold", margin = margin(10,0,0,0))) +
  theme(axis.text.x = element_text(size=8, colour = "black", hjust = 1,  angle = 45, margin = margin(10,0,0,0))) +
  theme(axis.text.y = element_text(size = 8, colour = "black")) +
  theme(plot.margin= margin(0.5,0.5,0.5,0.5, "cm")) +
  coord_cartesian(ylim = c(0,0.4)) +
  ggtitle("Psoriasis") +
  theme(title = element_text(size = 10))
  

a


ggsave("output/plots/barplot_frequence_gut_psoriasris.pdf", plot = a,
       width = 6, height = 4)

############################    Crohns disease   ###############################
sumfreq_crohn <- sumfreq_all_diseases[diseases == "Crohn's disease [regional enteritis]"]
median_crohn <- aggregate(sumfreq_crohn$x, by = list(Aminoacid =sumfreq_crohn$AA), FUN = median)
b <- ggplot(median_crohn, aes(Aminoacid, x)) +
  geom_bar(stat="identity", fill = "#7E0018DD",colour = "black", width= 1) +
  ylab("Abundance") +
  xlab("Amino acid auxotrophy") +
  theme(axis.line = element_line(size=0.2, colour = "black")) +
  theme(panel.background = element_rect(fill="white", colour= "white")) +
  theme(axis.title.y = element_text(colour = "black", size = 10, face = "bold", margin = margin(0,10,0,0)))+
  theme(axis.title.x = element_text(colour = "black", size = 10, face = "bold", margin = margin(10,0,0,0))) +
  theme(axis.text.x = element_text(size=8, colour = "black", hjust = 1,  angle = 45, margin = margin(10,0,0,0))) +
  theme(axis.text.y = element_text(size = 8, colour = "black")) +
  theme(plot.margin= margin(0.5,0.5,0.5,0.5, "cm")) +
  coord_cartesian(ylim = c(0,0.4)) +
  ggtitle("Crohns disease") +
  theme(title = element_text(size = 10))


b

ggsave("output/plots/barplot_frequence_gut_crohns.pdf", plot = b,
       width = 6, height = 4)

################# Rheumatoid arthritis with rheumatoid factor ##################
sumfreq_rheuma <- sumfreq_all_diseases[diseases == "Rheumatoid arthritis with rheumatoid factor"]
median_rheuma <- aggregate(sumfreq_rheuma$x, by = list(Aminoacid =sumfreq_rheuma$AA), FUN = median)
c <- ggplot(median_rheuma, aes(Aminoacid, x)) +
  geom_bar(stat="identity", fill = "#7E0018DD",colour = "black", width= 1) +
  ylab("Abundance") +
  xlab("Amino acid auxotrophy") +
  theme(axis.line = element_line(size=0.2, colour = "black")) +
  theme(panel.background = element_rect(fill="white", colour= "white")) +
  theme(axis.title.y = element_text(colour = "black", size = 10, face = "bold", margin = margin(0,10,0,0)))+
  theme(axis.title.x = element_text(colour = "black", size = 10, face = "bold", margin = margin(10,0,0,0))) +
  theme(axis.text.x = element_text(size=8, colour = "black", hjust = 1,  angle = 45, margin = margin(10,0,0,0))) +
  theme(axis.text.y = element_text(size = 8, colour = "black")) +
  theme(plot.margin= margin(0.5,0.5,0.5,0.5, "cm")) +
  coord_cartesian(ylim = c(0,0.4))   +
  ggtitle("Rheumatoid arthritis with\n rheumatoid factor") +
  theme(title = element_text(size = 10))


c

ggsave("output/plots/barplot_frequence_gut_rheuma.pdf", plot = c,
       width = 6, height = 4)


#################           Ulcerative colitis                ##################
sumfreq_uc <- sumfreq_all_diseases[diseases == "Ulcerative colitis"]
median_uc <- aggregate(sumfreq_uc$x, by = list(Aminoacid =sumfreq_uc$AA), FUN = median)
d <- ggplot(median_uc, aes(Aminoacid, x)) +
  geom_bar(stat="identity", fill = "#7E0018DD",colour = "black", width= 1) +
  ylab("Abundance") +
  xlab("Amino acid auxotrophy") +
  theme(axis.line = element_line(size=0.2, colour = "black")) +
  theme(panel.background = element_rect(fill="white", colour= "white")) +
  theme(axis.title.y = element_text(colour = "black", size = 10, face = "bold", margin = margin(0,10,0,0)))+
  theme(axis.title.x = element_text(colour = "black", size = 10, face = "bold", margin = margin(10,0,0,0))) +
  theme(axis.text.x = element_text(size=8, colour = "black", hjust = 1,  angle = 45, margin = margin(10,0,0,0))) +
  theme(axis.text.y = element_text(size = 8, colour = "black")) +
  theme(plot.margin= margin(0.5,0.5,0.5,0.5, "cm")) +
  coord_cartesian(ylim = c(0,0.4)) +
  ggtitle("Ulcerative colitis") +
  theme(title = element_text(size = 10))


d

ggsave("output/plots/barplot_frequence_gut_uc.pdf", plot = d,
       width = 6, height = 4)


#################           Ankylosing spondylitis              ################
sumfreq_as <- sumfreq_all_diseases[diseases == "Ankylosing spondylitis"]
median_as <- aggregate(sumfreq_as$x, by = list(Aminoacid =sumfreq_as$AA), FUN = median)
e <- ggplot(median_as, aes(Aminoacid, x)) +
  geom_bar(stat="identity", fill = "#7E0018DD",colour = "black", width= 1) +
  ylab("Abundance") +
  xlab("Amino acid auxotrophy") +
  theme(axis.line = element_line(size=0.2, colour = "black")) +
  theme(panel.background = element_rect(fill="white", colour= "white")) +
  theme(axis.title.y = element_text(colour = "black", size = 10, face = "bold", margin = margin(0,10,0,0)))+
  theme(axis.title.x = element_text(colour = "black", size = 10, face = "bold", margin = margin(10,0,0,0))) +
  theme(axis.text.x = element_text(size=8, colour = "black", hjust = 1,  angle = 45, margin = margin(10,0,0,0))) +
  theme(axis.text.y = element_text(size = 8, colour = "black")) +
  theme(plot.margin= margin(0.5,0.5,0.5,0.5, "cm")) +
  coord_cartesian(ylim = c(0,0.4)) +
  ggtitle("Ankylosing spondylitis") +
  theme(title = element_text(size = 10))


e

ggsave("output/plots/barplot_frequence_gut_ankylosing_spond.pdf", plot = e,
       width = 6, height = 4)

#################       Systemic lupus erythematosus (SLE)      ################
sumfreq_sle <- sumfreq_all_diseases[diseases == "Systemic lupus erythematosus (SLE)"]
median_sle <- aggregate(sumfreq_sle$x, by = list(Aminoacid =sumfreq_sle$AA), FUN = median)
f <- ggplot(median_sle, aes(Aminoacid, x)) +
  geom_bar(stat="identity", fill = "#7E0018DD",colour = "black", width= 1) +
  ylab("Abundance") +
  xlab("Amino acid auxotrophy") +
  theme(axis.line = element_line(size=0.2, colour = "black")) +
  theme(panel.background = element_rect(fill="white", colour= "white")) +
  theme(axis.title.y = element_text(colour = "black", size = 10, face = "bold", margin = margin(0,10,0,0)))+
  theme(axis.title.x = element_text(colour = "black", size = 10, face = "bold", margin = margin(10,0,0,0))) +
  theme(axis.text.x = element_text(size=8, colour = "black", hjust = 1,  angle = 45, margin = margin(10,0,0,0))) +
  theme(axis.text.y = element_text(size = 8, colour = "black")) +
  theme(plot.margin= margin(0.5,0.5,0.5,0.5, "cm")) +
  coord_cartesian(ylim = c(0,0.4)) +
  ggtitle("Systemic lupus erythematosus (SLE)") +
  theme(title = element_text(size = 10))


f

ggsave("output/plots/barplot_frequence_gut_SLE.pdf", plot = f,
       width = 6, height = 4)


########################             IBD             ###########################
sumfreq_IBD <- sumfreq_all_diseases[diseases == "Crohn's disease [regional enteritis]" | diseases == "Ulcerative colitis"
                                    | diseases == "Other and unsp noninfective gastroenteritis and colitis"]
median_IBD <- aggregate(sumfreq_IBD$x, by = list(Aminoacid =sumfreq_IBD$AA), FUN = median)
g <- ggplot(median_IBD, aes(Aminoacid, x)) +
  geom_bar(stat="identity", fill = "#7E0018DD",colour = "black", width= 1) +
  ylab("Abundance") +
  xlab("Amino acid auxotrophy") +
  theme(axis.line = element_line(size=0.2, colour = "black")) +
  theme(panel.background = element_rect(fill="white", colour= "white")) +
  theme(axis.title.y = element_text(colour = "black", size = 10, face = "bold", margin = margin(0,10,0,0)))+
  theme(axis.title.x = element_text(colour = "black", size = 10, face = "bold", margin = margin(10,0,0,0))) +
  theme(axis.text.x = element_text(size=8, colour = "black", hjust = 1,  angle = 45, margin = margin(10,0,0,0))) +
  theme(axis.text.y = element_text(size = 8, colour = "black")) +
  theme(plot.margin= margin(0.5,0.5,0.5,0.5, "cm")) +
  coord_cartesian(ylim = c(0,0.4)) +
  ggtitle("IBD") +
  theme(title = element_text(size = 10))


g

ggsave("output/plots/barplot_frequence_gut_IBD.pdf", plot = g,
       width = 6, height = 4)

autoimmun <- ggarrange(a,c,e,f,
          ncol = 2, nrow =2)
autoimmun
ggsave("output/plots/barplot_frequence_gut_autoimmun_diseases.pdf", plot = autoimmun,
       width = 6, height = 5)

IBD <- ggarrange(b,d,g,
          ncol = 2, nrow =2)
IBD
ggsave("output/plots/barplot_frequence_gut_IBD_combined.pdf", plot = IBD,
       width = 6, height = 5)

######################   Differences in the two IBD groups   ###################
###wilcoxon test

sumfreq_IBD_wilcox <- sumfreq_all_diseases[diseases == "Crohn's disease [regional enteritis]" | diseases == "Ulcerative colitis"]
ggplot(sumfreq)
AA <- unique(sumfreq_IBD_wilcox$AA)
l <- list()
k <- 1

for (AAi in AA) {
  print(AAi)
  wilcox_data <- sumfreq_IBD_wilcox[AA == AAi]
  wilcox <- wilcox_test(x~ diseases, data = wilcox_data)
  wilcox.pvalue <- data.table(wilcox.p = wilcox$p)
  wilcox.pvalue$AA <- AAi
  l[[k]] <- wilcox.pvalue
  k <- k +1
}

comp_IBD_groups <- rbindlist(l) 
remove(wilcox_data, wilcox, wilcox.pvalue,l)

comp_IBD_groups$padjust = p.adjust(comp_IBD_groups$wilcox.p, method = "fdr")
comp_IBD_groups[padjust < 0.05, sign.label := "P < 0.05"]
comp_IBD_groups
sumfreq_IBD
pvalue_sumfreq_IBD <- merge(sumfreq_IBD, IBD_auxos, by.x="AA", by.y="AA")

h <- ggplot(sumfreq_IBD_wilcox, aes(AA, x, fill = diseases)) +
  geom_boxplot(outlier.shape = NA) +
  ylab("Frequence of auxotrophic bacteria[%]") +
  xlab("Amino acid auxotrophies") +
  theme(axis.line = element_line(size=0.2, colour = "black")) +
  theme(panel.background = element_rect(fill="white", colour= "white")) +
  theme(axis.title.y = element_text(colour = "black", size = 10, face = "bold", margin = margin(0,10,0,0)))+
  theme(axis.title.x = element_text(colour = "black", size = 10, face = "bold", margin = margin(10,0,0,0))) +
  theme(axis.text.x = element_text(size=8, colour = "black", hjust = 1,  angle = 45, margin = margin(10,0,0,0))) +
  theme(axis.text.y = element_text(size =8, colour = "black")) +
  theme(legend.title = element_text(size = 8)) +
  theme(legend.text = element_text(size=8)) +
  theme(plot.margin= margin(0.5,0.5,0.5,0.5, "cm")) +
  scale_fill_manual(values=c("#4575b4", "#d73027")) +
  labs(fill = "Disease") +
  theme(legend.position = "bottom")
h
tbl <- tableGrob(comp_IBD_groups, rows = NULL, theme=tt)
grid.arrange(h,tbl,ncol =  2, heights = c(14,1,0.2))
#######################################################################
autoimmun1_wilcox <- sumfreq_all_diseases[diseases == "Psoriasis" | diseases == "Ankylosing spondylitis"]

AA <- unique(autoimmun1_wilcox$AA)
l <- list()
k <- 1

for (AAi in AA) {
  print(AAi)
  wilcox_data <- autoimmun1_wilcox[AA == AAi]
  wilcox <- wilcox_test(x~ diseases, data = wilcox_data)
  wilcox.pvalue <- data.table(wilcox.p = wilcox$p)
  wilcox.pvalue$AA <- AAi
  l[[k]] <- wilcox.pvalue
  k <- k +1
}

comp_autoimmun1_groups <- rbindlist(l) 
remove(wilcox_data, wilcox, wilcox.pvalue,l)

comp_autoimmun1_groups$padjust = p.adjust(comp_autoimmun1_groups$wilcox.p, method = "fdr")
comp_autoimmun1_groups[padjust < 0.05, sign.label := "P < 0.05"]
comp_autoimmun1_groups
sumfreq_IBD

j <- ggplot(autoimmun1_wilcox, aes(AA, x, fill = diseases)) +
  geom_boxplot(outlier.shape = NA) +
  ylab("Frequence of auxotrophic bacteria[%]") +
  xlab("Amino acid auxotrophies") +
  theme(axis.line = element_line(size=0.2, colour = "black")) +
  theme(panel.background = element_rect(fill="white", colour= "white")) +
  theme(axis.title.y = element_text(colour = "black", size = 10, face = "bold", margin = margin(0,10,0,0)))+
  theme(axis.title.x = element_text(colour = "black", size = 10, face = "bold", margin = margin(10,0,0,0))) +
  theme(axis.text.x = element_text(size=8, colour = "black", hjust = 1,  angle = 45, margin = margin(10,0,0,0))) +
  theme(axis.text.y = element_text(size =8, colour = "black")) +
  theme(legend.title = element_text(size = 8)) +
  theme(legend.text = element_text(size=8)) +
  theme(plot.margin= margin(0.5,0.5,0.5,0.5, "cm")) +
  scale_fill_manual(values=c("#4575b4", "#d73027")) +
  labs(fill = "Disease") +
  theme(legend.position = "bottom")
j

tbl1 <- tableGrob(comp_autoimmun1_groups, rows = NULL, theme=tt)
grid.arrange(j,tbl1,ncol =  2, heights = c(14,1,0.2))
########################################################################
autoimmun2_wilcox <- sumfreq_all_diseases[diseases == "Psoriasis" | diseases == "Rheumatoid arthritis with rheumatoid factor"]
View(autoimmun2_wilcox)
AA <- unique(autoimmun2_wilcox$AA)
l <- list()
k <- 1

for (AAi in AA) {
  print(AAi)
  wilcox_data <- autoimmun2_wilcox[AA == AAi]
  wilcox <- wilcox_test(x~ diseases, data = wilcox_data)
  wilcox.pvalue <- data.table(wilcox.p = wilcox$p)
  wilcox.pvalue$AA <- AAi
  l[[k]] <- wilcox.pvalue
  k <- k +1
}

comp_autoimmun2_groups <- rbindlist(l) 
remove(wilcox_data, wilcox, wilcox.pvalue,l)

comp_autoimmun2_groups$padjust = p.adjust(comp_autoimmun2_groups$wilcox.p, method = "fdr")
comp_autoimmun2_groups[padjust < 0.05, sign.label := "P < 0.05"]
comp_autoimmun2_groups
sumfreq_IBD

k <- ggplot(autoimmun2_wilcox, aes(AA, x, fill = diseases)) +
  geom_boxplot(outlier.shape = NA) +
  ylab("Frequence of auxotrophic bacteria[%]") +
  xlab("Amino acid auxotrophies") +
  theme(axis.line = element_line(size=0.2, colour = "black")) +
  theme(panel.background = element_rect(fill="white", colour= "white")) +
  theme(axis.title.y = element_text(colour = "black", size = 10, face = "bold", margin = margin(0,10,0,0)))+
  theme(axis.title.x = element_text(colour = "black", size = 10, face = "bold", margin = margin(10,0,0,0))) +
  theme(axis.text.x = element_text(size=8, colour = "black", hjust = 1,  angle = 45, margin = margin(10,0,0,0))) +
  theme(axis.text.y = element_text(size =8, colour = "black")) +
  theme(legend.title = element_text(size = 8)) +
  theme(legend.text = element_text(size=8)) +
  theme(plot.margin= margin(0.5,0.5,0.5,0.5, "cm")) +
  scale_fill_manual(values=c("#4575b4", "#d73027")) +
  labs(fill = "Disease")  +
  theme(legend.position = "bottom")
k

tbl2 <- tableGrob(comp_autoimmun2_groups, rows = NULL, theme=tt)
grid.arrange(k,tbl2,ncol =  2, heights = c(14,1,0.2))
#################################################################################
autoimmun3_wilcox <- sumfreq_all_diseases[diseases == "Psoriasis" | diseases == "Systemic lupus erythematosus (SLE)"]

AA <- unique(autoimmun3_wilcox$AA)
l <- list()
k <- 1

for (AAi in AA) {
  print(AAi)
  wilcox_data <- autoimmun3_wilcox[AA == AAi]
  wilcox <- wilcox_test(x~ diseases, data = wilcox_data)
  wilcox.pvalue <- data.table(wilcox.p = wilcox$p)
  wilcox.pvalue$AA <- AAi
  l[[k]] <- wilcox.pvalue
  k <- k +1
}

comp_autoimmun3_groups <- rbindlist(l) 
remove(wilcox_data, wilcox, wilcox.pvalue,l)

comp_autoimmun3_groups$padjust = p.adjust(comp_autoimmun3_groups$wilcox.p, method = "fdr")
comp_autoimmun3_groups[padjust < 0.05, sign.label := "P < 0.05"]
comp_autoimmun3_groups


l <- ggplot(autoimmun3_wilcox, aes(AA, x, fill = diseases)) +
  geom_boxplot(outlier.shape = NA) +
  ylab("Frequence of auxotrophic bacteria[%]") +
  xlab("Amino acid auxotrophies") +
  theme(axis.line = element_line(size=0.2, colour = "black")) +
  theme(panel.background = element_rect(fill="white", colour= "white")) +
  theme(axis.title.y = element_text(colour = "black", size = 10, face = "bold", margin = margin(0,10,0,0)))+
  theme(axis.title.x = element_text(colour = "black", size = 10, face = "bold", margin = margin(10,0,0,0))) +
  theme(axis.text.x = element_text(size=8, colour = "black", hjust = 1,  angle = 45, margin = margin(10,0,0,0))) +
  theme(axis.text.y = element_text(size =8, colour = "black")) +
  theme(legend.title = element_text(size = 8)) +
  theme(legend.text = element_text(size=8)) +
  theme(plot.margin= margin(0.5,0.5,0.5,0.5, "cm")) +
  scale_fill_manual(values=c("#4575b4", "#d73027")) +
  labs(fill = "Disease")  +
  theme(legend.position = "bottom")
l

tbl3 <- tableGrob(comp_autoimmun3_groups, rows = NULL, theme=tt)
grid.arrange(l,tbl3,ncol =  2, heights = c(14,1,0.2))
################################################################################
autoimmun4_wilcox <- sumfreq_all_diseases[diseases == "Ankylosing spondylitis" | diseases == "Rheumatoid arthritis with rheumatoid factor"]

AA <- unique(autoimmun4_wilcox$AA)
l <- list()
k <- 1

for (AAi in AA) {
  print(AAi)
  wilcox_data <- autoimmun4_wilcox[AA == AAi]
  wilcox <- wilcox_test(x~ diseases, data = wilcox_data)
  wilcox.pvalue <- data.table(wilcox.p = wilcox$p)
  wilcox.pvalue$AA <- AAi
  l[[k]] <- wilcox.pvalue
  k <- k +1
}

comp_autoimmun4_groups <- rbindlist(l) 
remove(wilcox_data, wilcox, wilcox.pvalue,l)

comp_autoimmun4_groups$padjust = p.adjust(comp_autoimmun4_groups$wilcox.p, method = "fdr")
comp_autoimmun4_groups[padjust < 0.05, sign.label := "P < 0.05"]
comp_autoimmun4_groups


m <- ggplot(autoimmun4_wilcox, aes(AA, x, fill = diseases)) +
  geom_boxplot(outlier.shape = NA) +
  ylab("Frequence of auxotrophic bacteria[%]") +
  xlab("Amino acid auxotrophies") +
  theme(axis.line = element_line(size=0.2, colour = "black")) +
  theme(panel.background = element_rect(fill="white", colour= "white")) +
  theme(axis.title.y = element_text(colour = "black", size = 10, face = "bold", margin = margin(0,10,0,0)))+
  theme(axis.title.x = element_text(colour = "black", size = 10, face = "bold", margin = margin(10,0,0,0))) +
  theme(axis.text.x = element_text(size=8, colour = "black", hjust = 1,  angle = 45, margin = margin(10,0,0,0))) +
  theme(axis.text.y = element_text(size =8, colour = "black")) +
  theme(legend.title = element_text(size = 8)) +
  theme(legend.text = element_text(size=8)) +
  theme(plot.margin= margin(0.5,0.5,0.5,0.5, "cm")) +
  scale_fill_manual(values=c("#4575b4", "#d73027")) +
  labs(fill = "Disease") +
  theme(legend.position = "bottom")
m
tbl4 <- tableGrob(comp_autoimmun4_groups, rows = NULL, theme=tt)
grid.arrange(m,tbl4,ncol =  2, heights = c(14,1,0.2))
##############################################################################
autoimmun5_wilcox <- sumfreq_all_diseases[diseases == "Ankylosing spondylitis" | diseases == "Systemic lupus erythematosus (SLE)"]

AA <- unique(autoimmun5_wilcox$AA)
l <- list()
k <- 1

for (AAi in AA) {
  print(AAi)
  wilcox_data <- autoimmun5_wilcox[AA == AAi]
  wilcox <- wilcox_test(x~ diseases, data = wilcox_data)
  wilcox.pvalue <- data.table(wilcox.p = wilcox$p)
  wilcox.pvalue$AA <- AAi
  l[[k]] <- wilcox.pvalue
  k <- k +1
}

comp_autoimmun5_groups <- rbindlist(l) 
remove(wilcox_data, wilcox, wilcox.pvalue,l)

comp_autoimmun5_groups$padjust = p.adjust(comp_autoimmun5_groups$wilcox.p, method = "fdr")
comp_autoimmun5_groups[padjust < 0.05, sign.label := "P < 0.05"]
comp_autoimmun5_groups



n <- ggplot(autoimmun5_wilcox, aes(AA, x, fill = diseases)) +
  geom_boxplot(outlier.shape = NA) +
  ylab("Frequence of auxotrophic bacteria[%]") +
  xlab("Amino acid auxotrophies") +
  theme(axis.line = element_line(size=0.2, colour = "black")) +
  theme(panel.background = element_rect(fill="white", colour= "white")) +
  theme(axis.title.y = element_text(colour = "black", size = 10, face = "bold", margin = margin(0,10,0,0)))+
  theme(axis.title.x = element_text(colour = "black", size = 10, face = "bold", margin = margin(10,0,0,0))) +
  theme(axis.text.x = element_text(size=8, colour = "black", hjust = 1,  angle = 45, margin = margin(10,0,0,0))) +
  theme(axis.text.y = element_text(size =8, colour = "black")) +
  theme(legend.title = element_text(size = 8)) +
  theme(legend.text = element_text(size=8)) +
  theme(plot.margin= margin(0.5,0.5,0.5,0.5, "cm")) +
  scale_fill_manual(values=c("#4575b4", "#d73027")) +
  labs(fill = "Disease") +
  theme(legend.position = "bottom")
n
tbl5 <- tableGrob(comp_autoimmun5_groups, rows = NULL, theme=tt)
grid.arrange(n,tbl5,ncol =  2, heights = c(14,1,0.2))
################################################################################
autoimmun6_wilcox <- sumfreq_all_diseases[diseases == "Rheumatoid arthritis with rheumatoid factor" | diseases == "Systemic lupus erythematosus (SLE)"]

AA <- unique(autoimmun6_wilcox$AA)
l <- list()
k <- 1

for (AAi in AA) {
  print(AAi)
  wilcox_data <- autoimmun6_wilcox[AA == AAi]
  wilcox <- wilcox_test(x~ diseases, data = wilcox_data)
  wilcox.pvalue <- data.table(wilcox.p = wilcox$p)
  wilcox.pvalue$AA <- AAi
  l[[k]] <- wilcox.pvalue
  k <- k +1
}

comp_autoimmun6_groups <- rbindlist(l) 
remove(wilcox_data, wilcox, wilcox.pvalue,l)

comp_autoimmun6_groups$padjust = p.adjust(comp_autoimmun6_groups$wilcox.p, method = "fdr")
comp_autoimmun6_groups[padjust < 0.05, sign.label := "P < 0.05"]
comp_autoimmun6_groups
library(grid)
library(gridExtra)
library(scales)

o <- ggplot(autoimmun6_wilcox, aes(AA, x, fill = diseases)) +
  geom_boxplot(outlier.shape = NA) +
  ylab("Frequence of auxotrophic bacteria[%]") +
  xlab("Amino acid auxotrophies") +
  theme(axis.line = element_line(size=0.2, colour = "black")) +
  theme(panel.background = element_rect(fill="white", colour= "white")) +
  theme(axis.title.y = element_text(colour = "black", size = 10, face = "bold", margin = margin(0,10,0,0)))+
  theme(axis.title.x = element_text(colour = "black", size = 10, face = "bold", margin = margin(10,0,0,0))) +
  theme(axis.text.x = element_text(size=8, colour = "black", hjust = 1,  angle = 45, margin = margin(10,0,0,0))) +
  theme(axis.text.y = element_text(size =8, colour = "black")) +
  theme(legend.title = element_text(size = 8)) +
  theme(legend.text = element_text(size=8)) +
  theme(plot.margin= margin(0.5,0.5,0.5,0.5, "cm")) +
  theme(legend.position = "bottom")+
  scale_fill_manual(values=c("#4575b4", "#d73027")) +
  labs(fill = "Disease")
o
o + annotation_custom(tableGrob(comp_autoimmun6_groups, rows = NULL, 
                                theme = ttheme_default(base_size = 4)), xmin =12, xmax=18, ymin = 0.3, ymax=0.5)
tt = ttheme_default(colhead=list(fg_params = list(parse=TRUE)), base_size = 8)
tbl <- tableGrob(comp_autoimmun6_groups, rows = NULL, theme=tt)
grid.arrange(o,tbl,ncol =  2, heights = c(14,1,0.2))
o + annotate(geom = "table", x = 20,y = 1, label = list(comp_autoimmun6_groups))

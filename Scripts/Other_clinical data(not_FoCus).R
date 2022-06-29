############      data about patients with chronic inflammations    ############
frequencies <- fread("/mnt/nuuk/2021/HRGM/TrypCID_16S_abundancies.csv")
metainfo <- fread("/mnt/nuuk/2021/HRGM/TrypCID_16S_metaInfo.csv")
View(metainfo)
gut_dis <- merge(frequencies, metainfo, by.x= "sample", by.y="sample")
View(gut_dis)
Auxotrophy_2[,c(4:16)] <- NULL

#get nubmer of rows for every disease
disee <- unique(gut_dis$short_desc)
list <- list()
k <- 1
for (di in disee) {
  r <- nrow(metainfo[metainfo$short_desc == di])
  new <- data.frame(nrow = r,
                    disease = di)
  list[[k]] <- new
  k <- k+1
}
nrow_metainfo <- rbindlist(list)
nrow_metainfo

##loop for mergin auxotrophies of genomes with metainfo
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


###################### Kruskal Wallis test #####################################
##all diseases
Kruskall <- sumfreq_all_diseases[diseases == "Psoriasis" | diseases == "Rheumatoid arthritis with rheumatoid factor" |
                                     diseases == "Ankylosing spondylitis" | diseases == "Systemic lupus erythematosus (SLE)" |
                                     diseases == "Crohn's disease [regional enteritis]" | diseases == "Ulcerative colitis"]
describe(Kruskall)
View(Kruskall)
relAA <- unique(Kruskall$AA)
remove(w)
w <- list()
k <- 1

for (AAi in relAA) {
  print(AAi)
  kruskall_auxo <- Kruskall[AA == AAi]
  kruskal <- kruskal.test(x~diseases, data = kruskall_auxo)
  p <- kruskal$p.value
  auxo_kruskal <- data.table(pvalue = kruskal$p.value,
                             AA = AAi)
  w[[k]] <- auxo_kruskal
  k <- k+1
}
w
Kruskal_all <- rbindlist(w)
Kruskal_all$padjust <- p.adjust(Kruskal_all$pvalue, method = "fdr")
Kruskal_all$Sign <- ifelse(Kruskal_all$padjust < 0.05, "*", "NS")
Kruskal_all

##################### rheumatic diseases #######################################
Kruskall_rheu <- sumfreq_all_diseases[diseases == "Rheumatoid arthritis with rheumatoid factor" |
                                   diseases == "Ankylosing spondylitis" | diseases == "Systemic lupus erythematosus (SLE)"]
relAA <- unique(Kruskall_rheu$AA)
rheu <- list()
k <- 1

for (AAi in relAA) {
  print(AAi)
  rheu_auxo <- Kruskall_rheu[AA == AAi]
  kruskal_rheu <- kruskal.test(x~diseases, data = rheu_auxo)
  p <- kruskal_rheu$p.value
  auxo_kruskal_rheu <- data.table(pvalue = kruskal_rheu$p.value,
                             AA = AAi)
  rheu[[k]] <- auxo_kruskal_rheu
  k <- k+1
}

Kruskal_Rheu <- rbindlist(rheu)
remove(rheu, rheu_auxo, kruskal_rheu,auxo_kruskal_rheu)
Kruskal_Rheu$padjust <- p.adjust(Kruskal_Rheu$pvalue, method = "fdr")
Kruskal_Rheu$Sign <- ifelse(Kruskal_Rheu$padjust < 0.05, "*", "NS")
Kruskal_Rheu

####boxplot
x4 <- ggplot(Kruskall_rheu[AA == "Gln"], aes(diseases, x)) +
  geom_boxplot(outlier.shape = NA)+
  ylab("Frequence of gln auxotrophic bacteria") +
  theme(axis.line = element_line(size =0.4, colour = "black")) +
  theme(axis.text.y = element_text(colour = "black")) +
  theme(axis.text.x = element_text(colour = "black")) +
  theme(axis.title.x = element_blank()) +
  theme(panel.background = element_rect(fill = "white", colour = "white")) +
  annotate(geom="text", x = 3, y = 0.08, size = 2, label = "adj.pvalue < 0.05 (Kruskal-Wallis-Test)") +
  ggtitle("Rheumathic diseases") +
  scale_x_discrete("diseases", labels= c("Ankylosing spondylitis" = "Ankylosing\nspondylitis",
                                         "Rheumatoid arthritis with rheumatoid factor" = "Rheumatoid arthritis\nwith rheumatoid factor",
                                         "Systemic lupus erythematosus (SLE)" = "Systemic lupus\nerythematosus (SLE)")) 
  
x4
ggsave("output/plots/boxplot_rheumathic_Kruskal_wallis.pdf", plot = x4,
       width = 6, height = 4)


########################   pairwise wilcox test    #############################
wilcox_all <- sumfreq_all_diseases[diseases == "Psoriasis" | diseases == "Rheumatoid arthritis with rheumatoid factor" |
                                               diseases == "Ankylosing spondylitis" | diseases == "Systemic lupus erythematosus (SLE)" |
                                               diseases == "Crohn's disease [regional enteritis]" | diseases == "Ulcerative colitis"]


#testing for every auxotrophie
relAA <- unique(wilcox_all$AA)
remove(c)
c <- list()
k <- 1
for (AAi in relAA) {
  print(AAi)
  wilcox_auxo <- wilcox_all[AA == AAi]
  pair_wilcox_all <- pairwise.wilcox.test(wilcox_auxo$x, wilcox_auxo$diseases, p.adjust = "fdr", exact = FALSE)
  pairwise_wilcox <- pair_wilcox_all$p.value
  test <- data.frame(pairwise_wilcox)
  test$AA <- AAi
  test$Ankylosing.spondylitis <- ifelse(test$Ankylosing.spondylitis < 0.05, "*", "N.S")
  test$Crohn.s.disease..regional.enteritis. <- ifelse(test$Crohn.s.disease..regional.enteritis. < 0.05, "*", "N.S")
  test$Psoriasis <- ifelse(test$Psoriasis < 0.05, "*", "N.S")
  test$Rheumatoid.arthritis.with.rheumatoid.factor <- ifelse(test$Rheumatoid.arthritis.with.rheumatoid.factor < 0.05, "*", "N.S")
  test$Systemic.lupus.erythematosus..SLE. <- ifelse(test$Systemic.lupus.erythematosus..SLE. < 0.05, "*", "N.S")
  c[[k]] <- test
  k <- k+1
}
c
wilcox_auxo <- wilcox_all[AA == "Lys"]
pair_wilcox_all <- pairwise.wilcox.test(wilcox_auxo$x, wilcox_auxo$diseases, p.adjust = "fdr", exact = FALSE)
############################    visualization   ################################
cbPalette <- c("#999999", "#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7")
x1 <- ggplot(wilcox_all[AA == "Trp"], aes(diseases, x, fill =diseases)) +
          geom_boxplot(outlier.shape = NA) +
  xlab("") +
  ylab("Frequency of auxotrophic bacteria") +
  geom_signif(comparisons = list(c("Crohn's disease [regional enteritis]", "Psoriasis")),
              map_signif_level = TRUE, annotation = "*") +
  theme(axis.line = element_line(size =0.4, colour = "black")) +
  #theme(axis.text.x = element_text(colour = "black", vjust = 0.5)) +
  theme(axis.text.y = element_text(colour = "black")) +
  scale_x_discrete("diseases", labels= c("Ankylosing spondylitis" = "Ankylosing\nspondylitis",
                                         "Crohn's disease [regional enteritis]" = "	Crohn's disease\n[regional enteritis]",
                                         "Psoriasis" = "Psoriasis",
                                         "Rheumatoid arthritis with rheumatoid factor" = "Rheumatoid arthritis\nwith rheumatoid factor",
                                         "Systemic lupus erythematosus (SLE)" = "Systemic lupus\nerythematosus (SLE)",
                                         "Ulcerative colitis" = "Ulcerative colitis")) +
  #annotate(geom="text", x = 5.5, y = 0.8, size = 2, label = "adj.pvalue < 0.05 (Pairwise Wilcoxon Test)") +
  theme(panel.background = element_rect(fill = "white", colour = "white")) +
  theme(axis.title.x = element_blank())+
  theme(axis.title.y = element_text(face="bold")) +
  scale_fill_manual(values = cbPalette) +
  theme(legend.position = "none") +
  theme(axis.text.x = element_blank()) +
  facet_grid(. ~ AA, scales = "free_x", space= "free_x") +
  theme(strip.text.x = element_text(size=12, face="bold")) +
  theme(legend.title = element_text(color="black", face = "bold"))
  
x1
ggsave("output/plots/boxplot_diseases_Trp.pdf", plot = x1,
       width = 9, height = 5)


x2 <- ggplot(wilcox_all[AA == "Lys"], aes(diseases, x, fill = diseases)) +
  geom_boxplot(outlier.shape = NA) +
  xlab("Diseases") +
  ylab("") +
  theme(axis.line = element_line(size =0.4, colour = "black")) +
  geom_signif(comparisons = list(c("Crohn's disease [regional enteritis]", "Psoriasis")),
              map_signif_level = TRUE, annotation = "*") +
  #theme(axis.text.x = element_text(colour = "black", vjust = 0.5)) +
  theme(axis.text.y = element_text(colour = "black")) +
  scale_x_discrete("diseases", labels= c("Ankylosing spondylitis" = "Ankylosing\nspondylitis",
                                         "Crohn's disease [regional enteritis]" = "	Crohn's disease\n[regional enteritis]",
                                         "Psoriasis" = "Psoriasis",
                                         "Rheumatoid arthritis with rheumatoid factor" = "Rheumatoid arthritis\nwith rheumatoid factor",
                                         "Systemic lupus erythematosus (SLE)" = "Systemic lupus\nerythematosus (SLE)",
                                         "Ulcerative colitis" = "Ulcerative colitis")) +
  #annotate(geom="text", x = 5.5, y = 0.09, size = 2, label = "adj.pvalue < 0.05 (Pairwise Wilcoxon Test)") +
  theme(panel.background = element_rect(fill = "white", colour = "white")) +
  theme(axis.title.x = element_blank()) +
  scale_fill_manual(values = cbPalette) +
  theme(legend.position = "none") +
  theme(axis.text.x = element_blank()) +
  facet_grid(. ~ AA, scales = "free_x", space= "free_x")+
  theme(strip.text.x = element_text(size=12, face="bold"))+
  theme(legend.title = element_text(color="black", face = "bold"))

x2

ggsave("output/plots/boxplot_diseases_lys.pdf", plot = x2,
       width = 9, height = 5)

x3 <- ggplot(wilcox_all[AA == "Gln"], aes(diseases, x, fill=diseases)) +
  geom_boxplot(outlier.shape = NA) +
  xlab("Diseases") +
  ylab("") +
  theme(axis.line = element_line(size =0.4, colour = "black")) +
  theme(panel.background = element_rect(fill = "white", colour = "white"))+
  geom_signif(comparisons = list(c("Ankylosing spondylitis", "Systemic lupus erythematosus (SLE)"),
                                 c("Psoriasis", "Systemic lupus erythematosus (SLE)"), 
                                 c("Rheumatoid arthritis with rheumatoid factor", "Systemic lupus erythematosus (SLE)"),
                                 c("Ulcerative colitis", "Systemic lupus erythematosus (SLE)")),
              map_signif_level = TRUE, annotation = "*", step_increase = 0.12) +
  #theme(axis.text.x = element_text(colour = "black",vjust = 0.5)) +
  theme(axis.text.y = element_text(colour = "black")) +
  scale_x_discrete("diseases", labels= c("Ankylosing spondylitis" = "Ankylosing\nspondylitis",
                                      "Crohn's disease [regional enteritis]" = "	Crohn's disease\n[regional enteritis]",
                                      "Psoriasis" = "Psoriasis",
                                      "Rheumatoid arthritis with rheumatoid factor" = "Rheumatoid arthritis\nwith rheumatoid factor",
                                      "Systemic lupus erythematosus (SLE)" = "Systemic lupus\nerythematosus (SLE)",
                                      "Ulcerative colitis" = "Ulcerative colitis")) +
  #annotate(geom="text", x = 5.5, y = 0.11, size = 2, label = "adj.pvalue < 0.05 (Pairwise Wilcoxon Test)") +
  theme(panel.background = element_rect(fill = "white", colour = "white")) +
  theme(axis.title.x = element_blank()) +
  scale_fill_manual(values = cbPalette) +
  theme(legend.position = "none") +
  theme(axis.text.x = element_blank()) +
  facet_grid(. ~ AA, scales = "free_x", space= "free_x") +
  theme(strip.text.x = element_text(size=12, face="bold")) +
  theme(legend.title = element_text(color="black", face = "bold"))

x3
ggsave("output/plots/boxplot_diseases_gln.pdf", plot = x3,
       width = 9, height = 5)


x4 <- ggarrange(x1,x2,x3, labels = "A",
          ncol=3, nrow=1, common.legend = TRUE,
          legend = c("bottom"))

ggsave("output/plots/boxplot_chron_diseases_all.pdf", plot = x4,
       width =11, height = 5)






boxplot_pso_crohn <- wilcox_all[diseases == "Crohn's disease [regional enteritis]" | diseases == "Psoriasis"]
boxplot_pso_crohn_filt <- boxplot_pso_crohn[AA == "Trp" | AA == "Lys"]

trp <- ggplot(boxplot_pso_crohn_filt[AA == "Trp"], aes(diseases, x)) +
  geom_boxplot() +
  xlab("Diseases") +
  ylab("Frequence of trp auxotrophic bacteria") +
  theme_minimal() +
  geom_signif(comparisons = list(c("Crohn's disease [regional enteritis]", "Psoriasis")),
              map_signif_level = TRUE, annotation = "*") +
  theme(axis.text.x = element_text(colour = "black")) +
  theme(axis.text.y = element_text(colour = "black"))
 
 ggsave("output/plots/boxplot_Crohn_psoriasis_Trp.pdf", plot = trp,
        width = 6, height = 4)

w <- ggplot(boxplot_pso_crohn_filt[AA == "Lys"], aes(diseases, x)) +
  geom_boxplot() +
  xlab("Diseases") +
  ylab("Frequence of lys auxotrophic bacteria") +
  theme_minimal() +
  geom_signif(comparisons = list(c("Crohn's disease [regional enteritis]", "Psoriasis")),
              map_signif_level = TRUE, annotation = "*") +
  theme(axis.text.x = element_text(colour = "black")) +
  theme(axis.text.y = element_text(colour = "black"))

ggsave("output/plots/boxplot_Crohn_psoriasis.pdf_Trp", plot = w,
       width = 6, height = 4)

# ######################   Differences in the two IBD groups   ###################
# ###wilcoxon test
# 
# sumfreq_IBD_wilcox <- sumfreq_all_diseases[diseases == "Crohn's disease [regional enteritis]" | diseases == "Ulcerative colitis"]
# ggplot(sumfreq)
# AA <- unique(sumfreq_IBD_wilcox$AA)
# l <- list()
# k <- 1
# 
# for (AAi in AA) {
#   print(AAi)
#   wilcox_data <- sumfreq_IBD_wilcox[AA == AAi]
#   wilcox <- wilcox_test(x~ diseases, data = wilcox_data)
#   wilcox.pvalue <- data.table(wilcox.p = wilcox$p)
#   wilcox.pvalue$AA <- AAi
#   l[[k]] <- wilcox.pvalue
#   k <- k +1
# }
# 
# comp_IBD_groups <- rbindlist(l) 
# remove(wilcox_data, wilcox, wilcox.pvalue,l)
# 
# comp_IBD_groups$padjust = p.adjust(comp_IBD_groups$wilcox.p, method = "fdr")
# comp_IBD_groups[padjust < 0.05, sign.label := "P < 0.05"]
# comp_IBD_groups
# sumfreq_IBD
# pvalue_sumfreq_IBD <- merge(sumfreq_IBD, IBD_auxos, by.x="AA", by.y="AA")
# 
# h <- ggplot(sumfreq_IBD_wilcox, aes(AA, x, fill = diseases)) +
#   geom_boxplot(outlier.shape = NA) +
#   ylab("Frequence of auxotrophic bacteria[%]") +
#   xlab("Amino acid auxotrophies") +
#   theme(axis.line = element_line(size=0.2, colour = "black")) +
#   theme(panel.background = element_rect(fill="white", colour= "white")) +
#   theme(axis.title.y = element_text(colour = "black", size = 10, face = "bold", margin = margin(0,10,0,0)))+
#   theme(axis.title.x = element_text(colour = "black", size = 10, face = "bold", margin = margin(10,0,0,0))) +
#   theme(axis.text.x = element_text(size=8, colour = "black", hjust = 1,  angle = 45, margin = margin(10,0,0,0))) +
#   theme(axis.text.y = element_text(size =8, colour = "black")) +
#   theme(legend.title = element_text(size = 8)) +
#   theme(legend.text = element_text(size=8)) +
#   theme(plot.margin= margin(0.5,0.5,0.5,0.5, "cm")) +
#   scale_fill_manual(values=c("#4575b4", "#d73027")) +
#   labs(fill = "Disease") +
#   theme(legend.position = "bottom")
# h
# tbl <- tableGrob(comp_IBD_groups, rows = NULL, theme=tt)
# grid.arrange(h,tbl,ncol =  2, heights = c(14,1,0.2))
# #######################################################################
# autoimmun1_wilcox <- sumfreq_all_diseases[diseases == "Psoriasis" | diseases == "Ankylosing spondylitis"]
# 
# AA <- unique(autoimmun1_wilcox$AA)
# l <- list()
# k <- 1
# 
# for (AAi in AA) {
#   print(AAi)
#   wilcox_data <- autoimmun1_wilcox[AA == AAi]
#   wilcox <- wilcox_test(x~ diseases, data = wilcox_data)
#   wilcox.pvalue <- data.table(wilcox.p = wilcox$p)
#   wilcox.pvalue$AA <- AAi
#   l[[k]] <- wilcox.pvalue
#   k <- k +1
# }
# 
# comp_autoimmun1_groups <- rbindlist(l) 
# remove(wilcox_data, wilcox, wilcox.pvalue,l)
# 
# comp_autoimmun1_groups$padjust = p.adjust(comp_autoimmun1_groups$wilcox.p, method = "fdr")
# comp_autoimmun1_groups[padjust < 0.05, sign.label := "P < 0.05"]
# comp_autoimmun1_groups
# sumfreq_IBD
# 
# j <- ggplot(autoimmun1_wilcox, aes(AA, x, fill = diseases)) +
#   geom_boxplot(outlier.shape = NA) +
#   ylab("Frequence of auxotrophic bacteria[%]") +
#   xlab("Amino acid auxotrophies") +
#   theme(axis.line = element_line(size=0.2, colour = "black")) +
#   theme(panel.background = element_rect(fill="white", colour= "white")) +
#   theme(axis.title.y = element_text(colour = "black", size = 10, face = "bold", margin = margin(0,10,0,0)))+
#   theme(axis.title.x = element_text(colour = "black", size = 10, face = "bold", margin = margin(10,0,0,0))) +
#   theme(axis.text.x = element_text(size=8, colour = "black", hjust = 1,  angle = 45, margin = margin(10,0,0,0))) +
#   theme(axis.text.y = element_text(size =8, colour = "black")) +
#   theme(legend.title = element_text(size = 8)) +
#   theme(legend.text = element_text(size=8)) +
#   theme(plot.margin= margin(0.5,0.5,0.5,0.5, "cm")) +
#   scale_fill_manual(values=c("#4575b4", "#d73027")) +
#   labs(fill = "Disease") +
#   theme(legend.position = "bottom")
# j
# 
# tbl1 <- tableGrob(comp_autoimmun1_groups, rows = NULL, theme=tt)
# grid.arrange(j,tbl1,ncol =  2, heights = c(14,1,0.2))
# ########################################################################
# autoimmun2_wilcox <- sumfreq_all_diseases[diseases == "Psoriasis" | diseases == "Rheumatoid arthritis with rheumatoid factor"]
# View(autoimmun2_wilcox)
# AA <- unique(autoimmun2_wilcox$AA)
# l <- list()
# k <- 1
# 
# for (AAi in AA) {
#   print(AAi)
#   wilcox_data <- autoimmun2_wilcox[AA == AAi]
#   wilcox <- wilcox_test(x~ diseases, data = wilcox_data)
#   wilcox.pvalue <- data.table(wilcox.p = wilcox$p)
#   wilcox.pvalue$AA <- AAi
#   l[[k]] <- wilcox.pvalue
#   k <- k +1
# }
# 
# comp_autoimmun2_groups <- rbindlist(l) 
# remove(wilcox_data, wilcox, wilcox.pvalue,l)
# 
# comp_autoimmun2_groups$padjust = p.adjust(comp_autoimmun2_groups$wilcox.p, method = "fdr")
# comp_autoimmun2_groups[padjust < 0.05, sign.label := "P < 0.05"]
# comp_autoimmun2_groups
# sumfreq_IBD
# 
# k <- ggplot(autoimmun2_wilcox, aes(AA, x, fill = diseases)) +
#   geom_boxplot(outlier.shape = NA) +
#   ylab("Frequence of auxotrophic bacteria[%]") +
#   xlab("Amino acid auxotrophies") +
#   theme(axis.line = element_line(size=0.2, colour = "black")) +
#   theme(panel.background = element_rect(fill="white", colour= "white")) +
#   theme(axis.title.y = element_text(colour = "black", size = 10, face = "bold", margin = margin(0,10,0,0)))+
#   theme(axis.title.x = element_text(colour = "black", size = 10, face = "bold", margin = margin(10,0,0,0))) +
#   theme(axis.text.x = element_text(size=8, colour = "black", hjust = 1,  angle = 45, margin = margin(10,0,0,0))) +
#   theme(axis.text.y = element_text(size =8, colour = "black")) +
#   theme(legend.title = element_text(size = 8)) +
#   theme(legend.text = element_text(size=8)) +
#   theme(plot.margin= margin(0.5,0.5,0.5,0.5, "cm")) +
#   scale_fill_manual(values=c("#4575b4", "#d73027")) +
#   labs(fill = "Disease")  +
#   theme(legend.position = "bottom")
# k
# 
# tbl2 <- tableGrob(comp_autoimmun2_groups, rows = NULL, theme=tt)
# grid.arrange(k,tbl2,ncol =  2, heights = c(14,1,0.2))
# #################################################################################
# autoimmun3_wilcox <- sumfreq_all_diseases[diseases == "Psoriasis" | diseases == "Systemic lupus erythematosus (SLE)"]
# 
# AA <- unique(autoimmun3_wilcox$AA)
# l <- list()
# k <- 1
# 
# for (AAi in AA) {
#   print(AAi)
#   wilcox_data <- autoimmun3_wilcox[AA == AAi]
#   wilcox <- wilcox_test(x~ diseases, data = wilcox_data)
#   wilcox.pvalue <- data.table(wilcox.p = wilcox$p)
#   wilcox.pvalue$AA <- AAi
#   l[[k]] <- wilcox.pvalue
#   k <- k +1
# }
# 
# comp_autoimmun3_groups <- rbindlist(l) 
# remove(wilcox_data, wilcox, wilcox.pvalue,l)
# 
# comp_autoimmun3_groups$padjust = p.adjust(comp_autoimmun3_groups$wilcox.p, method = "fdr")
# comp_autoimmun3_groups[padjust < 0.05, sign.label := "P < 0.05"]
# comp_autoimmun3_groups
# 
# 
# l <- ggplot(autoimmun3_wilcox, aes(AA, x, fill = diseases)) +
#   geom_boxplot(outlier.shape = NA) +
#   ylab("Frequence of auxotrophic bacteria[%]") +
#   xlab("Amino acid auxotrophies") +
#   theme(axis.line = element_line(size=0.2, colour = "black")) +
#   theme(panel.background = element_rect(fill="white", colour= "white")) +
#   theme(axis.title.y = element_text(colour = "black", size = 10, face = "bold", margin = margin(0,10,0,0)))+
#   theme(axis.title.x = element_text(colour = "black", size = 10, face = "bold", margin = margin(10,0,0,0))) +
#   theme(axis.text.x = element_text(size=8, colour = "black", hjust = 1,  angle = 45, margin = margin(10,0,0,0))) +
#   theme(axis.text.y = element_text(size =8, colour = "black")) +
#   theme(legend.title = element_text(size = 8)) +
#   theme(legend.text = element_text(size=8)) +
#   theme(plot.margin= margin(0.5,0.5,0.5,0.5, "cm")) +
#   scale_fill_manual(values=c("#4575b4", "#d73027")) +
#   labs(fill = "Disease")  +
#   theme(legend.position = "bottom")
# l
# 
# tbl3 <- tableGrob(comp_autoimmun3_groups, rows = NULL, theme=tt)
# grid.arrange(l,tbl3,ncol =  2, heights = c(14,1,0.2))
# ################################################################################
# autoimmun4_wilcox <- sumfreq_all_diseases[diseases == "Ankylosing spondylitis" | diseases == "Rheumatoid arthritis with rheumatoid factor"]
# 
# AA <- unique(autoimmun4_wilcox$AA)
# l <- list()
# k <- 1
# 
# for (AAi in AA) {
#   print(AAi)
#   wilcox_data <- autoimmun4_wilcox[AA == AAi]
#   wilcox <- wilcox_test(x~ diseases, data = wilcox_data)
#   wilcox.pvalue <- data.table(wilcox.p = wilcox$p)
#   wilcox.pvalue$AA <- AAi
#   l[[k]] <- wilcox.pvalue
#   k <- k +1
# }
# 
# comp_autoimmun4_groups <- rbindlist(l) 
# remove(wilcox_data, wilcox, wilcox.pvalue,l)
# 
# comp_autoimmun4_groups$padjust = p.adjust(comp_autoimmun4_groups$wilcox.p, method = "fdr")
# comp_autoimmun4_groups[padjust < 0.05, sign.label := "P < 0.05"]
# comp_autoimmun4_groups
# 
# 
# m <- ggplot(autoimmun4_wilcox, aes(AA, x, fill = diseases)) +
#   geom_boxplot(outlier.shape = NA) +
#   ylab("Frequence of auxotrophic bacteria[%]") +
#   xlab("Amino acid auxotrophies") +
#   theme(axis.line = element_line(size=0.2, colour = "black")) +
#   theme(panel.background = element_rect(fill="white", colour= "white")) +
#   theme(axis.title.y = element_text(colour = "black", size = 10, face = "bold", margin = margin(0,10,0,0)))+
#   theme(axis.title.x = element_text(colour = "black", size = 10, face = "bold", margin = margin(10,0,0,0))) +
#   theme(axis.text.x = element_text(size=8, colour = "black", hjust = 1,  angle = 45, margin = margin(10,0,0,0))) +
#   theme(axis.text.y = element_text(size =8, colour = "black")) +
#   theme(legend.title = element_text(size = 8)) +
#   theme(legend.text = element_text(size=8)) +
#   theme(plot.margin= margin(0.5,0.5,0.5,0.5, "cm")) +
#   scale_fill_manual(values=c("#4575b4", "#d73027")) +
#   labs(fill = "Disease") +
#   theme(legend.position = "bottom")
# m
# tbl4 <- tableGrob(comp_autoimmun4_groups, rows = NULL, theme=tt)
# grid.arrange(m,tbl4,ncol =  2, heights = c(14,1,0.2))
# ##############################################################################
# autoimmun5_wilcox <- sumfreq_all_diseases[diseases == "Ankylosing spondylitis" | diseases == "Systemic lupus erythematosus (SLE)"]
# 
# AA <- unique(autoimmun5_wilcox$AA)
# l <- list()
# k <- 1
# 
# for (AAi in AA) {
#   print(AAi)
#   wilcox_data <- autoimmun5_wilcox[AA == AAi]
#   wilcox <- wilcox_test(x~ diseases, data = wilcox_data)
#   wilcox.pvalue <- data.table(wilcox.p = wilcox$p)
#   wilcox.pvalue$AA <- AAi
#   l[[k]] <- wilcox.pvalue
#   k <- k +1
# }
# 
# comp_autoimmun5_groups <- rbindlist(l) 
# remove(wilcox_data, wilcox, wilcox.pvalue,l)
# 
# comp_autoimmun5_groups$padjust = p.adjust(comp_autoimmun5_groups$wilcox.p, method = "fdr")
# comp_autoimmun5_groups[padjust < 0.05, sign.label := "P < 0.05"]
# comp_autoimmun5_groups
# 
# 
# 
# n <- ggplot(autoimmun5_wilcox, aes(AA, x, fill = diseases)) +
#   geom_boxplot(outlier.shape = NA) +
#   ylab("Frequence of auxotrophic bacteria[%]") +
#   xlab("Amino acid auxotrophies") +
#   theme(axis.line = element_line(size=0.2, colour = "black")) +
#   theme(panel.background = element_rect(fill="white", colour= "white")) +
#   theme(axis.title.y = element_text(colour = "black", size = 10, face = "bold", margin = margin(0,10,0,0)))+
#   theme(axis.title.x = element_text(colour = "black", size = 10, face = "bold", margin = margin(10,0,0,0))) +
#   theme(axis.text.x = element_text(size=8, colour = "black", hjust = 1,  angle = 45, margin = margin(10,0,0,0))) +
#   theme(axis.text.y = element_text(size =8, colour = "black")) +
#   theme(legend.title = element_text(size = 8)) +
#   theme(legend.text = element_text(size=8)) +
#   theme(plot.margin= margin(0.5,0.5,0.5,0.5, "cm")) +
#   scale_fill_manual(values=c("#4575b4", "#d73027")) +
#   labs(fill = "Disease") +
#   theme(legend.position = "bottom")
# n
# tbl5 <- tableGrob(comp_autoimmun5_groups, rows = NULL, theme=tt)
# grid.arrange(n,tbl5,ncol =  2, heights = c(14,1,0.2))
# ################################################################################
# autoimmun6_wilcox <- sumfreq_all_diseases[diseases == "Rheumatoid arthritis with rheumatoid factor" | diseases == "Systemic lupus erythematosus (SLE)"]
# 
# AA <- unique(autoimmun6_wilcox$AA)
# l <- list()
# k <- 1
# 
# for (AAi in AA) {
#   print(AAi)
#   wilcox_data <- autoimmun6_wilcox[AA == AAi]
#   wilcox <- wilcox_test(x~ diseases, data = wilcox_data)
#   wilcox.pvalue <- data.table(wilcox.p = wilcox$p)
#   wilcox.pvalue$AA <- AAi
#   l[[k]] <- wilcox.pvalue
#   k <- k +1
# }
# 
# comp_autoimmun6_groups <- rbindlist(l) 
# remove(wilcox_data, wilcox, wilcox.pvalue,l)
# 
# comp_autoimmun6_groups$padjust = p.adjust(comp_autoimmun6_groups$wilcox.p, method = "fdr")
# comp_autoimmun6_groups[padjust < 0.05, sign.label := "P < 0.05"]
# comp_autoimmun6_groups
# library(grid)
# library(gridExtra)
# library(scales)
# 
# o <- ggplot(autoimmun6_wilcox, aes(AA, x, fill = diseases)) +
#   geom_boxplot(outlier.shape = NA) +
#   ylab("Frequence of auxotrophic bacteria[%]") +
#   xlab("Amino acid auxotrophies") +
#   theme(axis.line = element_line(size=0.2, colour = "black")) +
#   theme(panel.background = element_rect(fill="white", colour= "white")) +
#   theme(axis.title.y = element_text(colour = "black", size = 10, face = "bold", margin = margin(0,10,0,0)))+
#   theme(axis.title.x = element_text(colour = "black", size = 10, face = "bold", margin = margin(10,0,0,0))) +
#   theme(axis.text.x = element_text(size=8, colour = "black", hjust = 1,  angle = 45, margin = margin(10,0,0,0))) +
#   theme(axis.text.y = element_text(size =8, colour = "black")) +
#   theme(legend.title = element_text(size = 8)) +
#   theme(legend.text = element_text(size=8)) +
#   theme(plot.margin= margin(0.5,0.5,0.5,0.5, "cm")) +
#   theme(legend.position = "bottom")+
#   scale_fill_manual(values=c("#4575b4", "#d73027")) +
#   labs(fill = "Disease")
# o
# o + annotation_custom(tableGrob(comp_autoimmun6_groups, rows = NULL, 
#                                 theme = ttheme_default(base_size = 4)), xmin =12, xmax=18, ymin = 0.3, ymax=0.5)
# tt = ttheme_default(colhead=list(fg_params = list(parse=TRUE)), base_size = 8)
# tbl <- tableGrob(comp_autoimmun6_groups, rows = NULL, theme=tt)
# grid.arrange(o,tbl,ncol =  2, heights = c(14,1,0.2))
# o + annotate(geom = "table", x = 20,y = 1, label = list(comp_autoimmun6_groups))
# 

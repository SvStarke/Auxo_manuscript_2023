######################### arranging the figures for publication ################


#################figure 3#################
#Load models (completeness >= 85% and contamination <=2)

source("Scripts/init_models_filtered.R")

#Predict auxotrophies

source("Scripts/predict_auxos.R")

#Add information about the genomes

source("Scripts/Rasch_sampler.R")

#Analyze and visualize the distribution of the number of auxotrophies in every phylum

source("Scripts/Occurence_Auxos_together.R")

#create figure
fi3.1 <- ggarrange(o,t,
                 labels = c("A","B"),
                 ncol=2, nrow= 1, common.legend = FALSE)
fi3.1

ggsave("output/plots/figure3_18.07.22.pdf", plot = fi3.1,
       width = 9, height = 4.5)


################# supplementary figure 3 #################
#Load models (completeness >= 85% and contamination <=2)

source("Scripts/init_models_filtered.R")

#Predict auxotrophies

source("Scripts/predict_auxos.R")

#Add information about the genomes

source("Scripts/auxotable_not_melted.R")

#Analyze and visualize the distribution of the number of auxotrophies in every phylum

source("Scripts/number_auxo_per_phylum.R")

#create figure
fi4 <- ggarrange(fi_or, ba_or, ac_or, pr_or, ba_fam,
                 labels = c("a","","","", "b"),
                 ncol=1, nrow= 5, common.legend = FALSE,
                 heights = c(1,1,1,1), widths= c(1,1,1,1))
fi4

ggsave("output/plots/figure4_28.04.22.pdf", plot = fi4,
       width = 10, height = 14)

################ supplementary figure 4 ####################################
#completeness of the pathways
#Load models (Completeness >=85%, Contamination <=2)

source("Scripts/init_models_filtered.R")

#Predict auxotrophies

source("Scripts/predict_auxos.R")

#Add information about the genomes

source("Scripts/auxotable_not_melted.R")

#Analyze and then visualize the completeness of the amino acid biosynthesis pathways

source("Scripts/Completeness_pathways.R")

#### create figure ###
fi5 <- ggarrange(tr,hi,ch,se,le,il1,il2,il3,va,il4,il5,
                 labels = c("a","b","c", "d","e"),
                 hjust = c(-0.1,0.5,0.5,0.5,-0.1),
                 ncol=4, nrow= 3, common.legend = TRUE, legend = "bottom",
                 heights = c(1,1,1,1), widths= c(1,1,1,1))
fi5

ggsave("output/plots/figure5_28.04.22_new.pdf", plot = fi5,
       width = 9, height = 9)


#################figure 4 #################
#### fermentation products
#Analyze the production of by products with statistical analysis

source("Scripts/byproduct_production.R")

#### create figure #####
fi6.1 <- ggplot(stat_BP_x_auxo[auxo.compound != "Gly"], aes(auxo.compound, by.product,
                                                        fill = -log2(fisher.or))) +
  geom_tile() +
  geom_point(aes(shape = sign.label1), size = 1, show.legend = FALSE) +
  scale_fill_gradient2(high = "#ca0020", mid = "white", low = "#0571b0") +
  scale_shape_manual(values = 8, na.translate = FALSE) +
  scale_x_discrete(expand = c(0,0)) + scale_y_discrete(expand = c(0,0)) +
  labs(x = "Auxotrophy", y = "Fermentation\nproduct", shape = "",
       fill = expression(log[2]~'(odds ratio)')) +
  theme_bw() +
  theme(legend.position = "bottom",
        legend.justification = 	1,
        axis.text.x = element_text(color = "black", angle = 45, hjust = 1, size = 9),
        axis.text.y = element_text(color = "black", size = 9)) +
  theme(axis.title.y = element_text(margin = margin(t = 0, r = 10, b= 0, l = 0))) +
  theme(axis.title.x = element_text(margin = margin(t = 10, r = 0, b= 0, l = 0))) +
  theme(panel.background = element_blank())

fi6.1 +theme(plot.margin = unit(c(1,0.5,2,0.5), "cm")) +
  theme(legend.text = element_text(size=10),
        legend.title = element_text(size=8))

fi6 <- ggarrange(fi6.1,
                 labels = c("A"),
                 ncol=1, nrow= 1, common.legend = FALSE)
fi6

ggsave("output/plots/figure6_28.04.22.pdf", plot = fi6,
       width = 6, height = 3)


################# figure 5 ###########
#Load models 

source("Scripts/DZHK_data_init.R")

#Predict auxotrophies

source("Scripts/predict_auxos.R")

#Add information about the genomes

source("Scripts/auxotable_melted_merged.R")

#Analyze associations with diseases or health factors(BMI, weight, age)

source("Scripts/DZHK_Healthmarkers.R")

#Analyze associations with diseases or health factors(BMI, weight, age)

source("Scripts/Abundancies_Gut_DZHK.R")

###diversity and frequency of auxos
source("Scripts/diversity_Auxos.R")

###number of auxos and diversity
source("Scripts/numb_auxos_div_DZHK.R")

###number of auxos and hamming distance
source("Scripts/Hamming_DZHK.R")

###metabolomics and auxos
source("Scripts/init_dzhk_metabolome_analysis.R")
source("Scripts/dzhk_metabolome.R")


###combine all figures in one figure
fi7 <- ggarrange(ü, corr_health_div_plot,div_auxos, Hamming_shannon, ncol=1,
                   nrow=4, heights = c(1,1.3,1.1,1.1), widths= c(1,1,1,1),
                 labels = c("A","B", "C", "D"), hjust = c(-0.5,-0.5, -0.5, -0.5), vjust = c(1.5,1,-0.5,1))
fi7
fi7.1 <- ggarrange(fi7,met_DZHK, ncol = 2, nrow=1, widths = c(1.2,1.2), labels = "E",hjust = c(-45))
fi7.1
ggsave("output/plots/figure7_01.06.22_new_Abund_gut_DHZK.pdf", plot = fi7.1,
      width = 12, height = 15)

######## figure 6
fi8 <- ggarrange(stability, stability_Hamming, 
                 ncol=2,
                 nrow=1, 
                 labels = c("A","B"))
fi8
ggsave("output/plots/Stability_AuxosHamming.pdf", plot = fi8,
       width = 9.5, height = 5)

#partial spearman correlation and abundance of amino acid auxotrophies in the gut 
# fi7 <- ggarrange(b, linear_health,dis_health, ncol=3,
#                    nrow=1, heights = c(1,1,1), widths= c(3,1.5,1.8),
#                  labels = "AUTO", hjust = c(-0.5,1,-3))
# fi7
# ggsave("output/plots/figure7_01.06.22_new_Abund_gut.pdf", plot = fi7,
#        width = 10, height = 5)


###### supplementary material figure 2 ######

prop_all <- ggarrange(plot1,plot2,plot3,plot4,plot5,plot6,
                           nrow=3, ncol = 2, common.legend = TRUE)
prop_all

prop_all <- ggarrange(Asp, Asn, Leu,Ile,Glu,Gln,Val,Thr,Chor,Trp,Ser,Cys,Tyr,Phe,Met,Gly,Lys,Arg, Ala,Pro,His,
                      nrow=7, ncol = 3, common.legend = TRUE)
ggsave("output/plots/suppfigure2_03.06.pdf", plot = prop_all,
       width = 10, height = 20)


###supplementary material figure 5

nutr_F1F2_freq <- ggarrange(nutrition_popgen_F1, nutrition_popgen_F2,
                            nrow =1, ncol=2, common.legend = TRUE,
                            legend = c("bottom"),
                            labels = c("a","b"))

ggsave("output/plots/suppfigure_intakeAA_freqAuxos.pdf", plot=nutr_F1F2_freq,
       width= 11, height=5)




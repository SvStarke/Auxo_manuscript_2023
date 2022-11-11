###########################     SCFA production   ##############################
#Load models (completeness >=85% and a contamination <=2)

source("Scripts/init_models_filtered.R")

#Predict auxotrophies

source("Scripts/predict_auxos.R")

#Add information about the genomes

source("Scripts/auxotable_melted_merged.R")

cutoff_prodrate <- 1 # at which mmol/gDW the rate is considered as 'real'production

exchange <- get_exchanges(models)

relCompounds <- c("Butyrate","Propionate","Acetate","DL-Lactate",
                  "Succinate")

SCFAs <- exchange[name %in% relCompounds]


#get growth rates
m_gr <- get_growth(models)

m_growth <- data.table(Genome = names(m_gr),
                       Growth = unlist(m_gr))

#merge the files
relAuxos <- unique(Auxotrophy_2$Compound)

scfa_prod1 <- merge(SCFAs, m_growth, by.x = "model",
                    by.y = "Genome")
scfa_prod1[, prod_rate := flux / Growth]


# double loop for statistics
stat_BP_x_auxo <- list()
k <- 1
for(cpdi in relCompounds) {
  print(cpdi)
  for(auxoi in relAuxos) {
    tmp_cpd <- scfa_prod1[name == cpdi & !is.nan(prod_rate)]
    tmp_axt <- Auxotrophy_2[Compound == auxoi]
    tmp_axt  <- tmp_axt[complete.cases(tmp_axt[, "Prototrophy"]),]
    tmp <- merge(tmp_cpd, tmp_axt, by.x = "model", by.y = "Genomes")
    
    # Fisher Test
    cont_tab <- table(tmp$Prototrophy, tmp$prod_rate > cutoff_prodrate)
    
    if(nrow(cont_tab) == 2 & ncol(cont_tab) == 2) {
      # Fisher test
      test_fish <- fisher.test(cont_tab)
      dttmp <- data.table(by.product = cpdi,
                          auxo.compound = auxoi,
                          fisher.p = test_fish$p.value,
                          fisher.or = test_fish$estimate)
      
      stat_BP_x_auxo[[k]] <- dttmp
      k <- k + 1
    }
    
  }
}

stat_BP_x_auxo <- rbindlist(stat_BP_x_auxo)
stat_BP_x_auxo[, fisher.padj := p.adjust(fisher.p, method = "fdr")]
stat_BP_x_auxo[, fisher.or.log2 := -log2(fisher.or)]
stat_BP_x_auxo[fisher.padj < 0.05, sign.label1 := "Padj < 0.05"]

###visualization
p <- ggplot(stat_BP_x_auxo[auxo.compound != "Gly"], aes(auxo.compound, by.product,
                                fill = -log2(fisher.or))) +
  geom_tile() +
  geom_point(aes(shape = sign.label1), size = 1) +
  scale_fill_gradient2(high = "#ca0020", mid = "white", low = "#0571b0") +
  scale_shape_manual(values = 8, na.translate = FALSE) +
  scale_x_discrete(expand = c(0,0)) + scale_y_discrete(expand = c(0,0)) +
  labs(x = "Auxotrophy", y = "Fermentation\nproduct", shape = "",
       fill = expression(log[2]~'(odds ratio)')) +
  theme_bw() +
  theme(legend.position = "bottom",
        legend.justification = 	1,
        axis.text.x = element_text(color = "black", angle = 45, hjust = 1, size = 10),
        axis.text.y = element_text(color = "black", size = 8)) +
  theme(axis.title.y = element_text(margin = margin(t = 0, r = 20, b= 0, l = 0))) +
  theme(axis.title.x = element_text(margin = margin(t = 20, r = 0, b= 0, l = 0))) +
  theme(panel.background = element_blank())

p +theme(plot.margin = unit(c(1,0.5,2,0.5), "cm")) +
  theme(legend.text = element_text(size=8),
        legend.title = element_text(size=6))


ggsave("output/plots/Heatmap_auxo-X-byproducts_FoldChange_Fisher.pdf", plot = p,
       width = 6, height = 3.0)

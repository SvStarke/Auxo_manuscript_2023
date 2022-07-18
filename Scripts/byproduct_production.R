###########################     SCFA production   ##############################
#Load models (completeness >=85% and a contamination <=2)

source("Scripts/init_models_filtered.R")

#Predict auxotrophies

source("Scripts/predict_auxos.R")

#Add information about the genomes

source("Scripts/auxotable_melted_merged.R")

cutoff_prodrate <- 1 # at which mmol/gDW the rate is considered as 'real'production

exchange <- get_exchanges(models)
#fwrite(exchange, file = "exchange.csv")

relCompounds <- c("Butyrate","Propionate","Acetate","DL-Lactate",
                  "Succinate")

SCFAs <- exchange[name %in% relCompounds]
#is.data.table(SCFAs)


#get growth rates
m_gr <- get_growth(models)
#head(m_gr)
m_growth <- data.table(Genome = names(m_gr),
                       Growth = unlist(m_gr))
#fwrite(m_growth, file = "m_grwoth.csv")

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
      
      # Wilcoxon test
      tmp_wicoxdat <- copy(tmp)
      tmp_wicoxdat <- tmp_wicoxdat[, tmp_prod_rate := prod_rate] # neue Spalte erzeugen in der negative werte mit null ersetzt werden sollen
      tmp_wicoxdat <- tmp_wicoxdat[tmp_prod_rate < 0, tmp_prod_rate := 0] # ersetzen der neg. Werte mit 0
      
      wilcox <- wilcox_test(tmp_prod_rate ~ Prototrophy, data = tmp_wicoxdat)
      mean <- tapply(tmp_wicoxdat$tmp_prod_rate, tmp_wicoxdat$Prototrophy, mean)
      median <- tapply(tmp_wicoxdat$tmp_prod_rate, tmp_wicoxdat$Prototrophy, median)
      new1 <- t(mean)
      new2 <- t(median)
      mean1 <- new1[1,1]
      mean2 <- new1[1,2]
      median1 <- new2[1,1]
      median2 <- new2[1,2]
      FC_median <- log2(median1/median2)
      FC_mean <- log2(mean1/mean2)
      dttmp <- data.table(by.product = cpdi,
                          auxo.compound = auxoi,
                          fisher.p = test_fish$p.value,
                          fisher.or = test_fish$estimate,
                          wilcox.p = wilcox$p,
                          FC.log_median = FC_median,
                          FC.log_mean = FC_mean)
      
      stat_BP_x_auxo[[k]] <- dttmp
      k <- k + 1
    }
    
  }
}

stat_BP_x_auxo <- rbindlist(stat_BP_x_auxo)
stat_BP_x_auxo[, fisher.padj := p.adjust(fisher.p, method = "fdr")]
stat_BP_x_auxo[, wilcox.padj := p.adjust(wilcox.p, method = "fdr")]
stat_BP_x_auxo[, fisher.or.log2 := -log2(fisher.or)]
stat_BP_x_auxo[fisher.padj < 0.05, sign.label1 := "Padj < 0.05"]
stat_BP_x_auxo[wilcox.padj < 0.05, sign.label2 := "Padj < 0.05"]

#export data in csv
ferm_byprod <- stat_BP_x_auxo[,!c(5,6,7,9,11,12)]
write.csv(ferm_byprod, "output/files/fermentation_byproducts.csv")

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
        axis.text.y = element_text(color = "black", size = 10)) +
  theme(axis.title.y = element_text(margin = margin(t = 0, r = 20, b= 0, l = 0))) +
  theme(axis.title.x = element_text(margin = margin(t = 20, r = 0, b= 0, l = 0))) +
  theme(panel.background = element_blank())

p +theme(plot.margin = unit(c(1,0.5,2,0.5), "cm")) +
  theme(legend.text = element_text(size=10),
        legend.title = element_text(size=8))


ggsave("output/plots/Heatmap_auxo-X-byproducts_FoldChange_Fisher.pdf", plot = p,
       width = 6, height = 3.0)


q <- ggplot(stat_BP_x_auxo[auxo.compound != "Gly"], aes(auxo.compound, by.product,
                                                        fill = FC.log_mean)) +
  geom_tile() +
  geom_point(aes(shape = sign.label2), size = 0.5) +
  scale_fill_gradient2(high = "#ca0020", mid = "white", low = "#0571b0") +
  scale_shape_manual(values = 8, na.translate = FALSE) +
  scale_x_discrete(expand = c(0,0)) + scale_y_discrete(expand = c(0,0)) +
  labs(x = "Auxotrophy", y = "Fermentation\nproduct", shape = "",
       fill = expression(log[2]~'(Fold Change)')) +
  theme_bw() +
  theme(legend.position = "bottom",
        legend.justification = 	1,
        axis.text.x = element_text(color = "black", angle = 45, hjust = 1, size = 10),
        axis.text.y = element_text(color = "black", size = 10)) +
  theme(axis.title.y = element_text(margin = margin(t = 0, r = 20, b= 0, l = 0))) +
  theme(panel.background = element_blank())

q +theme(plot.margin = unit(c(1,0.5,2,0.5), "cm")) +
  theme(legend.text = element_text(size=10),
        legend.title = element_text(size=8))
q
  ggsave("output/plots/Heatmap_auxo-X-byproducts_FoldChange_Wilcox.pdf", plot = q,
       width = 6, height = 3.0)




wicoxdata
aggregate(tmp_wicoxdat$tmp_prod_rate, list(tmp_wicoxdat$name,tmp_wicoxdat$Prototrophy), FUN = mean)

          
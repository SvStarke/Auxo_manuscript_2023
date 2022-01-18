###########################     SCFA production   ##############################
library(MicrobiomeGS2)
library(ggplot2)
library(data.table)
library(MetBrewer)

cutoff_prodrate <- 1 # at which mmol/gDW the rate is considered as 'real'production

exchange <- get_exchanges(models)
View(exchange)
fwrite(exchange, file = "exchange.csv")

relCompounds <- c("Butyrate","Propionate","Acetate","DL-Lactate",
                  "Succinate")

SCFAs <- exchange[name %in% relCompounds]
is.data.table(SCFAs)
View(SCFAs)

#get growth rates
m_gr <- lapply(models, FUN = get_growth)
head(m_gr)
m_growth <- data.table(Genome = names(m_gr),
                       Growth = unlist(m_gr))
fwrite(m_growth, file = "m_grwoth.csv")

#merge the files
relAuxos <- unique(Auxotrophy_2$Compound)

scfa_prod1 <- merge(SCFAs, m_growth, by.x = "model",
                    by.y = "Genome")
scfa_prod1[, prod_rate := flux / Growth]


# douple loop for statistics
stat_BP_x_auxo <- list()
k <- 1
for(cpdi in relCompounds) {
  print(cpdi)
  for(auxoi in relAuxos) {
    tmp_cpd <- scfa_prod1[name == cpdi & !is.nan(prod_rate)]
    tmp_axt <- Auxotrophy_2[Compound == auxoi]
    tmp <- merge(tmp_cpd, tmp_axt, by.x = "model", by.y = "Genomes")
    
    # Fisher Test
    cont_tab <- table(tmp$Prototrophy, tmp$prod_rate > cutoff_prodrate)
    
    if(nrow(cont_tab) == 2 & ncol(cont_tab) == 2) {
      # Fisher test
      test_fish <- fisher.test(cont_tab)
      
      # Wilcoxon test
      
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
stat_BP_x_auxo[, fisher.or.log10 := log10(fisher.or)]
stat_BP_x_auxo[fisher.padj < 0.05, sign.label := "Padj < 0.05"]

p <- ggplot(stat_BP_x_auxo[auxo.compound != "Gly"], aes(auxo.compound, by.product,
                                fill = log10(fisher.or))) +
  geom_tile() +
  geom_point(aes(shape = sign.label), size = 0.5) +
  scale_fill_gradient2(high = "#ca0020", mid = "white", low = "#0571b0") +
  scale_shape_manual(values = 8, na.translate = FALSE) +
  scale_x_discrete(expand = c(0,0)) + scale_y_discrete(expand = c(0,0)) +
  labs(x = "Auxotrophy", y = "Fermentation\nproduct", shape = "",
       fill = expression(log[10]~'(odds ratio)')) +
  theme_bw() +
  theme(legend.position = "bottom",
        legend.box	= "vertical",
        axis.text.x = element_text(color = "black", angle = 45, hjust = 1),
        axis.text.y = element_text(color = "black")
  )
p

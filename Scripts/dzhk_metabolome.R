#–––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––#
# Correlate auxotrophy frequencies with serum metabolite levels #
#–––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––#
library(PResiduals)
library(parallel)
library(ggplot2)

rel_mets <- colnames(dzhk_metabolome)[-1]
rel_auxo <- dzhk_auxofreq[, sd(protofreq), by = aa][V1 > 0.001, aa]

pvals <- matrix(NA_real_, nrow=length(rel_mets), ncol=length(rel_auxo))
colnames(pvals) <- rel_auxo
rownames(pvals) <- rel_mets


ncores <- detectCores()
cl <- makeCluster(ncores)
clusterExport(cl, c("dzhk_auxofreq","dzhk_spl_info","dzhk_metabolome",
                    "rel_mets"))

partSpearAuxo <- function(iaa) {
  require(PResiduals)
  require(data.table)
  
  dt_tmp <- merge(dzhk_auxofreq[aa == iaa, .(EMGE, protofreq)],
                  dzhk_metabolome)
  dt_tmp <- merge(dt_tmp, dzhk_spl_info[, .(EMGE, Age, BMI, Geschlecht)])
  dt_tmp[, age_sc := scale(Age)]
  dt_tmp[, bmi_sc := scale(BMI)]
  dt_tmp$Geschlecht <- as.factor(dt_tmp$Geschlecht)
  
  # pvals
  outp <- rep(NA_real_, length(rel_mets))
  names(outp) <- rel_mets
  
  # estimates
  outest <- rep(NA_real_, length(rel_mets))
  names(outest) <- rel_mets
  
  # estimate variance
  outvar <- rep(NA_real_, length(rel_mets))
  names(outvar) <- rel_mets
  
  # estimate SD
  
  for(imet in rel_mets) {
    cat("   ",imet)
    dt_tmp$y <- dt_tmp[[imet]]
    zut <- partial_Spearman(protofreq|y ~ bmi_sc + age_sc + Geschlecht,
                            data = dt_tmp)
    outp[imet] <- zut$TS$TB$pval
    outest[imet] <- zut$TS$TB$ts
    outvar[imet] <- zut$TS$TB$var
  }
  
  return(list(p = outp,
              est = outest,
              var = outvar))
}

lstmp <- parLapply(cl, rel_auxo, fun = partSpearAuxo)
names(lstmp) <- rel_auxo
stopCluster(cl)

lstmp2 <- lapply(lstmp, FUN = function(x) {
  dt <- data.table(met = names(x$p),
                   pval = x$p,
                   est = x$est,
                   var = x$var)
})
lstmp2 <- rbindlist(lstmp2, idcol = "auxo")

dt_auxo_serumMets <- copy(lstmp2)
dt_auxo_serumMets[, est:= -est]
rm(lstmp2)
rm(lstmp)
dt_auxo_serumMets <- merge(dt_auxo_serumMets, BCdat$mets, by.x = "met",
                           by.y = "met_name")
dt_auxo_serumMets[, padj := p.adjust(pval), by = .(auxo, met_class)]

met_classes <- unique(dt_auxo_serumMets$met_class)

ggplot(dt_auxo_serumMets[!(met_class %in% c("Triacylglycerols","Sphingolipids",
                                            "Glycerophospholipids","Ceramides",
                                            "Glycosylceramides",
                                            "Cholesterol Esters","Dihydroceramides",
                                            "Diacylglycerols","Acylcarnitines"))],
       aes(auxo, met, fill = est)) +
  geom_tile() +
  geom_point(aes(shape = padj < 0.05)) +
  scale_fill_gradient2(high = "#ca0020", mid = "#f7f7f7", low = "#0571b0") +
  scale_shape_manual(values = c(NA, 19)) +
  scale_x_discrete(expand = c(0,0)) + scale_y_discrete(expand = c(0,0)) +
  facet_grid(met_class ~ ., space = "free", scales = "free") +
  labs(x = "Auxotrophy", y = "Metabolite",
       shape = expression(p[adj] < 0.05)) +
  theme_bw() +
  theme(strip.background = element_blank(),
        strip.text.y =  element_text(angle = 0, hjust = 0),
        axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1))

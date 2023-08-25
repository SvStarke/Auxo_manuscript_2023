##Prediction of auxotrophies in Price genomes
library(MicrobiomeGS2)
library(dplyr)
library(ggplot2)

models_price <- fetch_model_collection("/mnt/nuuk/2023/Price_prototrophs/models")
model.auxo_price <- predict_auxotrophies(models_price)
Auxotrophie2 <- data.frame(model.auxo_price)
#column und rows tauschen
Auxotroph2 <- t(Auxotrophie2)
Auxotroph2[which(is.na(Auxotroph2), arr.ind = T)] <- 1
#is.matrix(Auxotroph)
#data frame erzeugen
Auxotrophy2 <- data.frame(Auxotroph2)
# is.data.frame(Auxotrophy)
# str(Auxotrophy)
# Auxotrophy
Genome <- rownames(Auxotrophy2)
Auxotrophy2$Genomes <- Genome
fwrite(Auxotrophy2, file ="output/plots/Auxo_Price_Val.csv", sep = ",")
# ----
Auxotrophy2 <- as.data.table(Auxotrophy2)
Auxotrophy2

Auxotrophy_2_price  <- melt(Auxotrophy2, id.vars = "Genomes",
                     value.name = "Prototrophy", variable.name = "Compound")

Auxotrophy_2_price$Prototrophy_exp <- 1 
cont_table_price <- table(Auxotrophy_2_price$Prototrophy_exp, Auxotrophy_2_price$Prototrophy)

Corrected_prototrophy_prediction <- cont_table_price[1,2] / sum(cont_table_price)
False_prototrophy_prediction <- cont_table_price[1,1] / sum(cont_table_price)










# AA <- unique(Auxotrophy_2$Compound)
# k <- 1
# l <- list()
# 
# for(i in AA) {
#   tmp_Auxos <- Auxotrophy_2[Compound == i]
#   Perc_Proto <- ((nrow(tmp_Auxos[tmp_Auxos$Prototrophy == 1]) / 124) * 100)
#   tmp_perc_table <- data.frame(Perc_Proto,i)
#   l[[k]] <- tmp_perc_table
#   k <- k+1
# }
# Proto_Perc <- rbindlist(l)
# Proto_Perc
# 
# fwrite(Proto_Perc, file ="output/plots/Perc_Price_Val.csv", sep = ",")
# 
# #visualization
# 
# Price_Vali <- ggplot(Proto_Perc, aes(i, Perc_Proto))+
#   geom_bar(stat = "identity") +
#   xlab("Amino Acid") +
#   ylab("Percentage of prototrophies in genomes [%]")+
#   theme_bw()
# Price_Vali
# 
# ggsave("output/plots/Price_Val.pdf", plot = Price_Vali,
#        width = 6,height =4)
# 
# ##number of organism with at least one amino acid auxotrophy
# Auxotrophy$count <- rowSums(Auxotrophy == 0)
# 
# nrow(Auxotrophy[Auxotrophy$count !=0])
# 
# #validation of loop
# (nrow(Auxotrophy_2[Auxotrophy_2$Prototrophy == 1 & Auxotrophy_2$Compound == "Trp"]))
# 
# ##combine both HRGM prediction and Price predictions
# 
# combined_HRGm_price <- merge(Freq_auxos, Proto_Perc, by.x= "AA", by.y= "i")
# combined_HRGm_price$Proto_HRGM <- (100-combined_HRGm_price$percentage)
# 
# Price_HRGM_Vali <- ggplot(combined_HRGm_price) +
#   geom_bar(data = combined_HRGm_price,aes(x = AA, y=Perc_Proto, fill =  "Prototrophic genomes"), stat = "identity", position = "identity", alpha = 0.6) +
#   geom_point(data = combined_HRGm_price, aes(x= AA, y = Proto_HRGM, colour = "HRGM genomes"))+
#   labs(fill = "", colour = "") +
#   scale_fill_manual(values = c("grey38", "#D55E00")) +
#   theme_bw()+
#   xlab("Amino acids")+
#   ylab("Predicted prototrophies [%]")
# Price_HRGM_Vali
# 
# ggsave("output/plots/Price_HRGM_Val.pdf", plot = Price_HRGM_Vali,
#        width = 7,height =4)

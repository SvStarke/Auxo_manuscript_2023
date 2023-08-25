###load scripts

source("Scripts/Popgen_data_init.R")
source("Scripts/predict_auxos.R")



##prepare data table about meta information
##load info about all samples
popgen_samples <- fread("data/mgx_abundances/atlas_samples.csv")
popgen_samples$new_name <- sub("-L001|-L002", "", popgen_samples$Full_Name)
#View(popgen_samples)
##load first timepoint
popgen_F1 <- fread("data/meta/1.metadata1.tsv")
popgen_F1$newID <- sub("$", "_F1", popgen_F1$new_id)
#head(popgen_F1)
#colnames(popgen_F1)

popgen_F1 <- popgen_F1[,c(1,2,20)]

##load second timepoint
popgen_F2 <- fread("data/meta/1.metadata2.tsv")
popgen_F2$newID <- sub("$", "_F2", popgen_F2$new_id)
colnames(popgen_F2) [2] <- "new_id"

popgen_F2 <- popgen_F2[,c(1,2,19)]

F1_F2 <- rbind(popgen_F2, popgen_F1)

#View(popgen_F1)
##merge information about time points with sample info

popgen_data <- merge(F1_F2, popgen_samples, by.x= "newID", by.y="new_name")
#View(popgen_data)
#delete sample ID with only  available information about one timepoint in original dataframe about samples
popgen_data <- popgen_data[popgen_data$atlas_name != "S26" & popgen_data$atlas_name != "S38" & popgen_data$atlas_name!= "S66"]
popgen_data
# 
###calculate bray curtis dissimilarity between two time points
popgen_hrgm_abun1 <- t(popgen_hrgm_abun1)
popgen_hrgm_abun1 <- data.frame(popgen_hrgm_abun1)
popgen_hrgm_abun1$sample <- row.names(popgen_hrgm_abun1)

#View(popgen_mags_abun1)
test <- merge(popgen_hrgm_abun1, popgen_data, by.x= "sample", by.y = "atlas_name")
test <- data.table(test)
#View(test)
# 
# ##create a loop
popgen_data_names <- unique(test$SampleID)
F1_names <- unique(popgen_F1$SampleID)
F1_names <- F1_names[F1_names %in% popgen_data_names]
F2_names <- unique(popgen_F2$SampleID)
F2_names <- F2_names[F2_names %in% popgen_data_names]


####        Bray Curtis Distance with normalised abundance data             ####
##create a loop
test2 <- merge(popgen_norm, popgen_data, by.x= "sample", by.y = "atlas_name")
test2 <- data.table(test2)
#View(test2)

k <- 1
Bray_all1 <- list()
for(F1 in F1_names) {
  print(F1)
  for(F2 in F2_names) {
    tmp3<- test2[SampleID == F1 | SampleID == F2, ]
    if(tmp3$new_id[1] == tmp3$new_id[2]) {
      tmp4 <- tmp3[,c(2:571)]
      df3 <- mutate_all(tmp4, function(x) as.numeric(as.character(x)))
      Bray_tmp <-vegdist(df3, method = "bray")
      #Bray_tmp
      Bray1 <- data.table(Sample = F1,
                         Bray_distance = Bray_tmp)
      Bray_all1[[k]] <- Bray1
      k <- k+1
    }
  }
}

Bray_samples2 <- rbindlist(Bray_all1)



###get abundance weight average of auxotrophies

##### number of auxotrophies and diversity
Auxotrophy$count <- rowSums(Auxotrophy == 0)

numb_auxos_popgen <- merge(Auxotrophy, popgen_relabun, by.x= "Genomes", by.y="model")
numb_auxos_popgen <- data.table(numb_auxos_popgen)

numb_auxos_popgen[is.na(count), count:= 0]

new <- numb_auxos_popgen[ ,sum(count*prop), by = sample]
tmp_popgen_div_numb_auxos <- merge(new, popgen_data, by.x="sample", by.y="atlas_name")

#keep only auxotrophic data for F1 
popgen_div_numb_auxos <- tmp_popgen_div_numb_auxos[tmp_popgen_div_numb_auxos$SampleID %in% F1_names]

#merge diversity files and bray curtis based on normalized abundance data

popgen_bray_numb_auxos2 <- merge(popgen_div_numb_auxos, Bray_samples2, by.x="SampleID", by.y="Sample")


##correlation analysis for stability
cor.test(popgen_bray_numb_auxos2$Bray_distance, popgen_bray_numb_auxos2$V1, method = "spearman", exact = FALSE)

##visualization
stability <- ggplot(popgen_bray_numb_auxos2, aes(V1, 1-Bray_distance)) +
  geom_smooth(method=lm, col = "black") +
  geom_point(shape = 21) +
  theme_bw() +
  xlab("Abundance-weighted average of auxotrophies") +
  ylab("1 - Bray Curtis distance") +
  theme(axis.text.x = element_text(colour="black")) +
  theme(axis.text.y = element_text(colour= "black")) +
  theme(axis.title.y = element_text(size = 12, margin = margin(r = 10))) +
  theme(axis.title.x = element_text(size = 12, margin = margin(t = 10))) +
  stat_cor(method = "spearman")
stability

ggsave("output/plots/stability_auxos.pdf", plot = stability,
       width = 4, height = 4)


# ###############           Jaccard distance              ################
# 
# ##create a loop
# test3 <- merge(popgen_norm, popgen_data, by.x= "sample", by.y = "atlas_name")
# test3 <- data.table(test3)
# 
# #View(test2)
# 
# k <- 1
# Jac_all <- list()
# for(F1 in F1_names) {
#   print(F1)
#   for(F2 in F2_names) {
#     tmp4<- test2[SampleID == F1 | SampleID == F2, ]
#     
#     if(tmp4$new_id[1] == tmp4$new_id[2]) {
#       tmp5 <- tmp4[,c(2:571)]
#       df4 <- mutate_all(tmp5, function(x) as.numeric(as.character(x)))
#       Jac_tmp <-vegdist(df4, method = "jaccard")
#       #Bray_tmp
#       Jac <- data.table(Sample = F1,
#                           Jac_distance = Jac_tmp)
#       Jac_all[[k]] <- Jac
#       k <- k+1
#     }
#   }
# }
# 
# Jac_samples <- rbindlist(Jac_all)
# 
# 
# #merge diversity files and bray curtis based on normalized abundance data
# 
# x<- merge(popgen_div_numb_auxos, Jac_samples, by.x="SampleID", by.y="Sample")
# 
# ##correlation analysis for stability
# cor.test(x$Jac_distance, x$V1, method = "spearman", exact = FALSE)
# 
# ### new column with information about timepoints
# tmp_popgen_div_numb_auxos[, Time := gsub("^.{9}", "", SampleID)]
# df <- tmp_popgen_div_numb_auxos
# 
# ###create new dataframe for visualization of F1 and F2
# F1 <- df[Time == "F1"]
# F1 <- F1$V1
# F2 <- df[Time == "F2"]
# F2 <- F2$V1
# F1F2 <- cbind(F1, F2)
# F1F2 <- data.frame(F1F2)
# 
# 
# 
# corr_F1F2 <- ggplot(F1F2, aes(F1, F2))+
#   geom_point() +
#   geom_smooth() +
#   xlab("Abundance-weighted average of auxotrophies per MAG at F1")+
#   ylab("Abundance-weighted average of auxotrophies per MAG at F2") +
#   theme(panel.background = element_rect(fill="white", colour= "white")) +
#   theme(axis.line = element_line(size=0.2))
# 
# corr_F1F2.2 <-ggplot(F1F2, aes(F1, F2))+
#   geom_point() +
#   xlab("Abundance-weighted average of auxotrophies per MAG at F1")+
#   ylab("Abundance-weighted average of auxotrophies per MAG at F2") +
#   theme(panel.background = element_rect(fill="white", colour= "white")) +
#   theme(axis.line = element_line(size=0.2))
# 
# cor.test(F1F2$F1, F1F2$F2, method = "spearman")
# 
# ggsave("output/plots/F1F2_with_line.pdf", plot = corr_F1F2,
#        width = 6, height = 5)
# 
# ggsave("output/plots/F1F2without_line.pdf", plot = corr_F1F2.2,
#        width = 6, height = 5)
# 




####       ordered F1 and F2 plot ##############
# sel_order <- 
#   df %>% 
#   filter(Time == "F1") %>% 
#   arrange(desc(V1)) %>% 
#   mutate(new_id = factor(new_id))
# 
# 
# order_stability <- df %>% 
#   mutate(new_id = factor(new_id, levels = sel_order$new_id, ordered = TRUE)) %>% 
#   ggplot(aes(x = new_id, y = V1), group = new_id) +
#   geom_point(aes(colour = Time)) +
#   xlab("Samples") +
#   ylab("Abundance-weighted average of auxotrophies per MAG") +
#   theme(axis.text.x = element_blank()) +
#   labs(colour = "Timepoint") +
#   guides(fill=guide_legend(title="Timepoints")) +
#   theme(panel.background = element_rect(fill="white", colour= "white")) +
#   theme(axis.line = element_line(size=0.2)) 
# 
# order_stability 
# ggsave("output/plots/order_stabilityx.pdf", plot = order_stability,
#        width = 10, height = 5)

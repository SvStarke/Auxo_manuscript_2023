##prepare data table about meta information

##load info about all samples
popgen_samples <- fread("/mnt/nuuk/2022/MR_popgen_MGX/atlas/atlas_samples.csv")
popgen_samples$new_name <- sub("-L001|-L002", "", popgen_samples$Full_Name)

##load first timepoint
popgen_F1 <- fread("/mnt/nuuk/2022/MR_popgen_MGX/meta/1.metadata1.tsv")
popgen_F1$newID <- sub("$", "_F1", popgen_F1$new_id)
head(popgen_F1)
colnames(popgen_F1)

popgen_F1 <- popgen_F1[,c(1,2,20)]

##load second timepoint
popgen_F2 <- fread("/mnt/nuuk/2022/MR_popgen_MGX/meta/1.metadata2.tsv")
popgen_F2$newID <- sub("$", "_F2", popgen_F2$new_id)
colnames(popgen_F2) [2] <- "new_id"

popgen_F2 <- popgen_F2[,c(1,2,19)]

F1_F2 <- rbind(popgen_F2, popgen_F1)

##merge information about time points with sample info

popgen_data <- merge(F1_F2, popgen_samples, by.x= "newID", by.y="new_name")
View(popgen_data)
#delete sample ID with only  available information about one timepoint
popgen_data <- popgen_data[popgen_data$atlas_name != "S26" & popgen_data$atlas_name != "S38" & popgen_data$atlas_name!= "S66"]
popgen_data

###calculate bray curtis dissimilarity between two time points
rownames(popgen_mags_abun1) <- popgen_mags_abun1[,1]
popgen_mags_abun1 <- t(popgen_mags_abun1)
popgen_mags_abun1 <- data.frame(popgen_mags_abun1)
popgen_mags_abun1$sample <- row.names(popgen_mags_abun1)

test <- merge(popgen_mags_abun1, popgen_data, by.x= "sample", by.y = "atlas_name")
test <- data.table(test)
View(test)

##create a loop
popgen_data_names <- unique(test$SampleID)
F1_names <- unique(popgen_F1$SampleID)
F1_names <- F1_names[F1_names %in% popgen_data_names]
F2_names <- unique(popgen_F2$SampleID)
F2_names <- F2_names[F2_names %in% popgen_data_names]

View(test)
k <- 1
Bray_all <- list()
for(F1 in F1_names) {
  print(F1)
  for(F2 in F2_names) {
    tmp<- test[SampleID == F1 | SampleID == F2, ]

    if(tmp$new_id[1] == tmp$new_id[2]) {
      tmp2 <- tmp[,c(2:571)]
      df2 <- mutate_all(tmp2, function(x) as.numeric(as.character(x)))
      Bray_tmp <-vegdist(df2, method = "bray")
      #Bray_tmp
      Bray <- data.table(Sample = F1,
                         Bray_dist = Bray_tmp)
      Bray_all[[k]] <- Bray
      k <- k+1
     }
  }
}

Bray_samples <- rbindlist(Bray_all)
Bray_all
View(Bray_samples)

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

popgen_bray_numb_auxos <- merge(popgen_div_numb_auxos, Bray_samples, by.x="SampleID", by.y="Sample")


##correlation analysis for stability
cor.test(popgen_bray_numb_auxos$Bray_dist, popgen_bray_numb_auxos$V1, method = "spearman", exact = FALSE)


##visualization
stability <- ggplot(popgen_bray_numb_auxos, aes(V1, Bray_dist,)) +
  geom_point() +
  geom_smooth(method=lm) +
  theme_bw() +
  xlab("Abundance-weighted average of auxotrophies per MAG") +
  ylab("Bray distance") +
  theme(axis.text.x = element_text(colour="black")) +
  theme(axis.text.y = element_text(colour= "black")) +
  theme(axis.title.y = element_text(size = 10, margin = margin(r = 10))) +
  theme(axis.title.x = element_text(size = 10, margin = margin(t = 10))) 
stability

ggsave("output/plots/stability_auxos.pdf", plot = stability,
       width = 5, height = 5)



####

###visualization of number of average

### new column with information about timepoints
tmp_popgen_div_numb_auxos[, Time := gsub("^.{9}", "", SampleID)]
df <- tmp_popgen_div_numb_auxos
View(tmp_popgen_div_numb_auxos)




sel_order <- 
  df %>% 
  filter(Time == "F1") %>% 
  arrange(desc(V1)) %>% 
  mutate(new_id = factor(new_id))


order_stability <- df %>% 
  mutate(new_id = factor(new_id, levels = sel_order$new_id, ordered = TRUE)) %>% 
  ggplot(aes(x = new_id, y = V1), group = new_id) +
  geom_point(aes(colour = Time)) +
  xlab("Samples") +
  ylab("Abundance-weighted average of auxotrophies per MAG") +
  theme(axis.text.x = element_blank()) +
  labs(colour = "Timepoint") +
  guides(fill=guide_legend(title="Timepoints")) +
  theme(panel.background = element_rect(fill="white", colour= "white")) +
  theme(axis.line = element_line(size=0.2))

order_stability 
ggsave("output/plots/order_stabilityx.pdf", plot = order_stability,
       width = 10, height = 5)

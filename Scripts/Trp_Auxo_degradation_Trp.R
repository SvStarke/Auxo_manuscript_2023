###############trp auxo vs proto degradation of trp 

Metadata <- fread("/mnt/nuuk/2021/HRGM/REPR_Genomes_metadata.tsv")
View(Metadata)
relGenomes <- Metadata[`Completeness (%)`>= 85 & `Contamination (%)` <=2 & !grepl("^d__Archaea", `GTDB Taxonomy`), `HRGM name`]

models <- fetch_model_collection("/mnt/nuuk/2021/HRGM/models/", IDs = relGenomes)
rxns <- fetch_model_collection("/mnt/nuuk/2021/HRGM/models", file.type = "reactions", IDs = relGenomes)
pwys <- fetch_model_collection("/mnt/nuuk/2021/HRGM/models", file.type = "pathways", IDs = relGenomes)
pwys_cov_Trpcat <- get_pathway_coverage(models,rxns,pwys, pathways.of.interest = c("TRPCAT-PWY", 
                                                                                   "TRYPDEG-PWY",
                                                                                   "TRPKYNCAT-PWY",
                                                                                   "PWY-3181",
                                                                                   "PWY-5081",
                                                                                   "PWY-5655",
                                                                                   "PWY-5651"))

Trp_auxo <- Auxotrophy[,c(7,22)]

trp_auxo <- merge(pwys_cov_Trpcat,Trp_auxo, by.x="model", by.y="Genomes")                                                                          

###statistical analysis
relEC <- unique(trp_auxo$ec)
View(trp_auxo)
relPathways <- unique(trp_auxo$pathway)
remove(trp_auxo_x)
trp_auxo_x <- list()
k <- 1

for(pi in relPathways) {
  print(pi)
  for (eci in relEC){
    tmp_trp_x <- trp_auxo[pathway == pi]
    tmp_trp <- tmp_trp_x[ec == eci]
    
    # Fisher Test
    cont_tab <- table(tmp_trp$Trp, tmp_trp$prediction)
    if(nrow(cont_tab) == 2 & ncol(cont_tab) == 2) {
    col_order <- c("TRUE", "FALSE")
    cont_tab <- cont_tab[, col_order]
    test_fish <- fisher.test(cont_tab)
    dttmp <- data.table(EC = eci,
                        fisher.p = test_fish$p.value,
                        fisher.or = test_fish$estimate,
                        Pathway = pi,
                        AA = "Trp")
      
      trp_auxo_x[[k]] <- dttmp
      k <- k + 1
    }
  }
}

trp_auxos <- rbindlist(trp_auxo_x)
trp_auxos[, fisher.padj := p.adjust(fisher.p, method = "fdr")]
trp_auxos[, fisher.or.log2 := -log2(fisher.or)]
trp_auxos[fisher.padj < 0.05, sign.label1 := "Padj < 0.05"]

trp_auxos_mod <- trp_auxos %>% 
  mutate(Pathway = recode(Pathway, "PWY-3162" = "PWY-3162","TRPKYNCAT-PWY" = "TRPKYN\nCAT-PWY", "TRYPDEG-PWY" = "TRYP\nDEG\n-PWY"))

p <- ggplot(trp_auxos_mod, aes(AA, EC,fill = -log2(fisher.or))) +
  geom_tile() +
  geom_point(aes(shape = sign.label1), size = 1) +
  scale_fill_gradient2(high = "#ca0020", mid = "white", low = "#0571b0") +
  scale_shape_manual(values = 8, na.translate = FALSE) +
  facet_grid(Pathway ~., scales = "free_y", space= "free_y") +
  labs(x = "", y = "Enzymatic reaction", shape = "",
       fill = expression(log[2]~'(odds ratio)')) +
  theme_bw() +
  theme(legend.position = "right",
        legend.justification = 	1,
        axis.text.x = element_blank(),
        axis.text.y = element_text(color = "black", size = 8)) +
  theme(axis.title.y=  element_text(size=9, margin = margin(t = 0, r = 20, b= 0, l = 0), face = "bold")) +
  theme(panel.background = element_blank()) +
  theme(strip.text.y = element_text(size=6))

p

ggsave("output/plots/TRP_auxo_Trp_degradation.pdf", plot = p,
       width = , height = 5)

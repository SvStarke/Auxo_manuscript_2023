#bile acid metabolism
sybil::SYBIL_SETTINGS("SOLVER","cplexAPI")
library(MicrobiomeGS2)
models <- readRDS("/home/svenja/workspace/2022/auxotrophies_hrgm/data/models.RDS")
flux <- lapply(models, FUN = get_flux_distribution, exclude.unused = F)
saveRDS(flux, "flux.RDS")
models <- readRDS("models.RDS")
fwrite(flux, file = "flux.csv")
View(flux)

flux_i <- rbindlist(flux, idcol = "model")
View(flux_i)
fwrite(flux_i, file = "flux_i.csv")
flux_i <- read.csv("/Users/svenjabusche/Desktop/flux_i.csv")
remove(flux_i)

#bei produced metabolites schauen ob Produkte entstanden sind
#bile acid deconjugation
#colate, chenodeoxycholate
flu_deconjug_BA <- flux_i$ec == "3.5.1.24" | flux_i$ec == "3.5.1.74" 
flu_deconjug_BA <- flux_i[flu_deconjug_BA, ]
View(flu_deconjug_BA)

fwrite(flu_deconjug_BA, file="flu_deconjug_BA.csv")

#delete archaes
Auxotrophy_2 <- Auxotrophy_2[!(Auxotrophy_2$phylum == "Euryarchaeota" | Auxotrophy_2$phylum == "Thermoplasmatota"), ]

#put all Firmicutes phyla to one Firmicutes phylum
Auxotrophy_2$phylum[Auxotrophy_2$phylum == "Firmicutes_A"] <- "Firmicutes"
Auxotrophy_2$phylum[Auxotrophy_2$phylum == "Firmicutes_B"] <- "Firmicutes"
Auxotrophy_2$phylum[Auxotrophy_2$phylum== "Firmicutes_G"] <- "Firmicutes"
Auxotrophy_2$phylum[Auxotrophy_2$phylum == "Firmicutes_I"] <- "Firmicutes"

#mergen mit tryptophan auxotrophie daten
flux_deconjug_BA_all <- merge(Auxotrophy_2, flu_deconjug_BA, by.x = "Genomes", by.y = "model", allow.cartesian = TRUE)
View(flux_deconjug_BA_all)

#trp auxotrophic microbiota
fluxes_deconjug_BA_trp_auxo <- flux_deconjug_BA_all$Prototrophy == 0 &  flux_deconjug_BA_all$Compound == "Trp" &  flux_deconjug_BA_all$status == "good_blast"
fluxes_deconjug_BA__trp_auxo <- flux_deconjug_BA_all[fluxes_deconjug_BA_trp_auxo, ]
View(fluxes_deconjug_BA__trp_auxo)

fwrite(fluxes_deconjug_BA__trp_auxo, file="fluxes_deconjug_BA__trp_auxo.csv")

ggplot(data=fluxes_deconjug_BA__trp_auxo, aes(seedID, fill = phylum)) +
  geom_bar() +
  ggtitle("Reactions with a good blast for bile acid deconjugation \nby tryptophan auxotrophic microbiota") +
  ylab("number of genomes") +
  xlab("seedID")

ggplot(data=fluxes_deconjug_BA__trp_auxo, aes(seedID, fill = ec)) +
  geom_bar() +
  ggtitle("Enzymes with a good blast for bile acid deconjugation \n by tryptophan auxotrophic microbiota") +
  ylab("number of genomes") +
  xlab("seedID")

ggplot(data=fluxes_deconjug_BA__trp_auxo, aes(name, fill = ec)) +
  geom_bar() +
  ggtitle("Enzymes with a good blast for bile acid deconjugation \n by tryptophan auxotrophic microbiota") +
  ylab("number of genomes") +
  xlab("seedID")


#compare to tryptophan prototrophic microbiota
fluxes_deconjug_BA_trp_proto <- flux_deconjug_BA_all$Compound == "Trp" & flux_deconjug_BA_all$status == "good_blast"
fluxes_deconjug_BA_trp_proto <- flux_deconjug_BA_all[fluxes_deconjug_BA_trp_proto, ]
View(fluxes_deconjug_BA_trp_proto)

ggplot(data=fluxes_deconjug_BA_trp_proto, aes(Prototrophy, fill = seedID)) +
  geom_bar() 


#fisher test for rxn02015
#create a dataframe 
#number of genomes who are tryptophan auxotrophic but not H2S producer
#filter for tryptophan auxotrophic bacteria
numb_auxo <- Auxotrophy_2$Compound == "Trp" & Auxotrophy_2$Prototrophy == 0
numb_auxo <- Auxotrophy_2[numb_auxo, ]
View(numb_B_auxo)
#delete archaen
numb_auxo<- numb_auxo[!(numb_auxo$phylum == "Euryarchaeota" | numb_auxo$phylum == "Thermoplasmatota"), ]
View(numb_auxo)
#calculate number for contingency table
#number of auxotrophic bacteria rxn02015 (ec 3.5.1.74)
fluxes_deconjug_BA_trp_rxn02015_auxo <- flux_deconjug_BA_all$Compound == "Trp" & flux_deconjug_BA_all$status == "good_blast" &
  flux_deconjug_BA_all$seedID == "rxn02015" & flux_deconjug_BA_all$Prototrophy == 0 
fluxes_deconjug_BA_trp_rxn02015_auxo <- flux_deconjug_BA_all[fluxes_deconjug_BA_trp_rxn02015_auxo, ]
nrow(fluxes_deconjug_BA_trp_rxn02015_auxo)
#number prototrophic bacteria rxn02015
fluxes_deconjug_BA_trp_rxn02015_proto <- flux_deconjug_BA_all$Compound == "Trp" & flux_deconjug_BA_all$status == "good_blast" &
  flux_deconjug_BA_all$seedID == "rxn02015" & flux_deconjug_BA_all$Prototrophy == 1 
fluxes_deconjug_BA_trp_rxn02015_proto <- flux_deconjug_BA_all[fluxes_deconjug_BA_trp_rxn02015_proto, ]
nrow(fluxes_deconjug_BA_trp_rxn02015_proto)
#number trp auxotrophic and non producer
sum_no_rxn02015_auxo <- nrow(numb_auxo) - nrow(fluxes_deconjug_BA_trp_rxn02015_auxo)
sum_no_rxn02015_auxo 
#number trp prototrophic and non producer
sum_no_rxn02015_proto <- 5416 - (nrow(fluxes_deconjug_BA_trp_rxn02015_auxo) + nrow(fluxes_deconjug_BA_trp_rxn02015_proto) + sum_no_rxn02015_auxo)
sum_no_rxn02015_proto 

#create contigency table
fisher_rxn02015 <- data.frame("Trp-Auxotrophy" = c(nrow(fluxes_deconjug_BA_trp_rxn02015_auxo), sum_no_rxn02015_auxo), 
                       "Trp-Prototrophy"= c(nrow(fluxes_deconjug_BA_trp_rxn02015_proto), sum_no_rxn02015_proto), 
                       row.names = c("rxn02015: yes", "rxn02015: no"))

fisher_rxn02015
#control numbers
nrow(H2S_prod_auxo) + nrow(H2S_prod_proto) + sum_no_H2S_auxo + sum_no_H2S_proto == 5416

#create mosaic plot
mosaicplot(fisher_rxn02015, color = TRUE, main="rxn02015 - trp auxotrophic and prototrophic microbiota")
#fisher test
fisher.test(fisher_rxn02015)

#for rxn02795
#number of auxotrophic bacteria  rxn02795 (ec 3.5.1.24)
fluxes_deconjug_BA_trp_rxn02795_auxo <- flux_deconjug_BA_all$Compound == "Trp" & flux_deconjug_BA_all$status == "good_blast" &
  flux_deconjug_BA_all$seedID == "rxn02795" & flux_deconjug_BA_all$Prototrophy == 0 
fluxes_deconjug_BA_trp_rxn02795_auxo <- flux_deconjug_BA_all[fluxes_deconjug_BA_trp_rxn02795_auxo, ]
nrow(fluxes_deconjug_BA_trp_rxn02795_auxo)
#number prototrophic bacteria  rxn02795
fluxes_deconjug_BA_trp_rxn02795_proto <- flux_deconjug_BA_all$Compound == "Trp" & flux_deconjug_BA_all$status == "good_blast" &
  flux_deconjug_BA_all$seedID == "rxn02795" & flux_deconjug_BA_all$Prototrophy == 1 
fluxes_deconjug_BA_trp_rxn02795_proto <- flux_deconjug_BA_all[fluxes_deconjug_BA_trp_rxn02795_proto, ]
nrow(fluxes_deconjug_BA_trp_rxn02795_proto)
#number trp auxotrophic and non producer
sum_no_rxn02795_auxo <- nrow(numb_auxo) - nrow(fluxes_deconjug_BA_trp_rxn02795_auxo)
sum_no_rxn02795_auxo 
#number trp prototrophic and non producer
sum_no_rxn02795_proto <- 5416 - (nrow(fluxes_deconjug_BA_trp_rxn02795_auxo) + nrow(fluxes_deconjug_BA_trp_rxn02795_proto) + sum_no_rxn02795_auxo)
sum_no_rxn02795_proto 

#create contigency table
fisher_rxn02795 <- data.frame("Trp-Auxotrophy" = c(nrow(fluxes_deconjug_BA_trp_rxn02795_auxo), sum_no_rxn02795_auxo), 
                              "Trp-Prototrophy"= c(nrow(fluxes_deconjug_BA_trp_rxn02795_proto), sum_no_rxn02795_proto), 
                              row.names = c(" rxn02795: yes", " rxn02795: no"))

fisher_rxn02795
#control numbers
nrow(fluxes_deconjug_BA_trp_rxn02795_auxo) + sum_no_rxn02795_auxo + nrow(fluxes_deconjug_BA_trp_rxn02795_proto) + sum_no_rxn02795_proto == 5416

#create mosaic plot
mosaicplot(fisher_rxn02795, color = TRUE, main=" rxn02795 - trp auxotrophic and prototrophic microbiota")
#fisher test
fisher.test(fisher_rxn02795)


#for rxn02796
#number of auxotrophic bacteria  rxn02796 (3.5.1.74)
fluxes_deconjug_BA_trp_rxn02796_auxo <- flux_deconjug_BA_all$Compound == "Trp" & flux_deconjug_BA_all$status == "good_blast" &
  flux_deconjug_BA_all$seedID == "rxn02796" & flux_deconjug_BA_all$Prototrophy == 0 
fluxes_deconjug_BA_trp_rxn02796_auxo <- flux_deconjug_BA_all[fluxes_deconjug_BA_trp_rxn02796_auxo, ]
nrow(fluxes_deconjug_BA_trp_rxn02796_auxo)
#number prototrophic bacteria  rxn02796
fluxes_deconjug_BA_trp_rxn02796_proto <- flux_deconjug_BA_all$Compound == "Trp" & flux_deconjug_BA_all$status == "good_blast" &
  flux_deconjug_BA_all$seedID == "rxn02796" & flux_deconjug_BA_all$Prototrophy == 1 
fluxes_deconjug_BA_trp_rxn02796_proto <- flux_deconjug_BA_all[fluxes_deconjug_BA_trp_rxn02796_proto, ]
nrow(fluxes_deconjug_BA_trp_rxn02796_proto)
#number trp auxotrophic and non producer
sum_no_rxn02796_auxo <- nrow(numb_auxo) - nrow(fluxes_deconjug_BA_trp_rxn02796_auxo)
sum_no_rxn02796_auxo 
#number trp prototrophic and non producer
sum_no_rxn02796_proto <- 5416 - (nrow(fluxes_deconjug_BA_trp_rxn02796_auxo) + nrow(fluxes_deconjug_BA_trp_rxn02796_proto) + sum_no_rxn02796_auxo)
sum_no_rxn02796_proto 

#create contigency table
fisher_rxn02796 <- data.frame("Trp-Auxotrophy" = c(nrow(fluxes_deconjug_BA_trp_rxn02796_auxo), sum_no_rxn02796_auxo), 
                              "Trp-Prototrophy"= c(nrow(fluxes_deconjug_BA_trp_rxn02796_proto), sum_no_rxn02796_proto), 
                              row.names = c(" rxn02796: yes", " rxn02796: no"))

fisher_rxn02796
#control numbers
nrow(fluxes_deconjug_BA_trp_rxn02796_auxo) + sum_no_rxn02796_auxo + nrow(fluxes_deconjug_BA_trp_rxn02796_proto) + sum_no_rxn02796_proto == 5416

#create mosaic plot
mosaicplot(fisher_rxn02796, color = TRUE, main=" rxn02796 - trp auxotrophic and prototrophic microbiota")
#fisher test
fisher.test(fisher_rxn02796)

#for rxn03094
#number of auxotrophic bacteria  rxn03094 (3.5.1.24)
fluxes_deconjug_BA_trp_rxn03094_auxo <- flux_deconjug_BA_all$Compound == "Trp" & flux_deconjug_BA_all$status == "good_blast" &
  flux_deconjug_BA_all$seedID == "rxn03094" & flux_deconjug_BA_all$Prototrophy == 0 
fluxes_deconjug_BA_trp_rxn03094_auxo <- flux_deconjug_BA_all[fluxes_deconjug_BA_trp_rxn03094_auxo, ]
nrow(fluxes_deconjug_BA_trp_rxn03094_auxo)
#number prototrophic bacteria  rxn03094
fluxes_deconjug_BA_trp_rxn03094_proto <- flux_deconjug_BA_all$Compound == "Trp" & flux_deconjug_BA_all$status == "good_blast" &
  flux_deconjug_BA_all$seedID == "rxn03094" & flux_deconjug_BA_all$Prototrophy == 1 
fluxes_deconjug_BA_trp_rxn03094_proto <- flux_deconjug_BA_all[fluxes_deconjug_BA_trp_rxn03094_proto, ]
nrow(fluxes_deconjug_BA_trp_rxn03094_proto)
#number trp auxotrophic and non producer
sum_no_rxn03094_auxo <- nrow(numb_auxo) - nrow(fluxes_deconjug_BA_trp_rxn03094_auxo)
sum_no_rxn03094_auxo 
#number trp prototrophic and non producer
sum_no_rxn03094_proto <- 5416 - (nrow(fluxes_deconjug_BA_trp_rxn03094_auxo) + nrow(fluxes_deconjug_BA_trp_rxn03094_proto) + sum_no_rxn03094_auxo)
sum_no_rxn03094_proto 

#create contigency table
fisher_rxn03094 <- data.frame("Trp-Auxotrophy" = c(nrow(fluxes_deconjug_BA_trp_rxn03094_auxo), sum_no_rxn03094_auxo), 
                              "Trp-Prototrophy"= c(nrow(fluxes_deconjug_BA_trp_rxn03094_proto), sum_no_rxn03094_proto), 
                              row.names = c(" rxn03094: yes", " rxn03094: no"))

fisher_rxn03094
#control numbers
nrow(fluxes_deconjug_BA_trp_rxn03094_auxo) + sum_no_rxn03094_auxo + nrow(fluxes_deconjug_BA_trp_rxn03094_proto) + sum_no_rxn03094_proto == 5416

#create mosaic plot
mosaicplot(fisher_rxn03094, color = TRUE, main=" rxn03094 - trp auxotrophic and prototrophic microbiota")
#fisher test
fisher.test(fisher_rxn03094)

#number of auxotrophic bacteria  rxn03095, ec 3.5.1.24
fluxes_deconjug_BA_trp_rxn03095_auxo_ec3.5.1.24 <- flux_deconjug_BA_all$Compound == "Trp" & flux_deconjug_BA_all$status == "good_blast" &
  flux_deconjug_BA_all$seedID == "rxn03095" & flux_deconjug_BA_all$Prototrophy == 0 & flux_deconjug_BA_all$ec == "3.5.1.74"
fluxes_deconjug_BA_trp_rxn03095_auxo_ec3.5.1.24 <- flux_deconjug_BA_all[fluxes_deconjug_BA_trp_rxn03095_auxo_ec3.5.1.24 , ]
nrow(fluxes_deconjug_BA_trp_rxn03095_auxo_ec3.5.1.24 )
#number prototrophic bacteria  rxn03095, , ec 3.5.1.24
fluxes_deconjug_BA_trp_rxn03095_proto_ec3.5.1.24 <- flux_deconjug_BA_all$Compound == "Trp" & flux_deconjug_BA_all$status == "good_blast" &
  flux_deconjug_BA_all$seedID == "rxn03095" & flux_deconjug_BA_all$Prototrophy == 1 & flux_deconjug_BA_all$ec == "3.5.1.74"
fluxes_deconjug_BA_trp_rxn03095_proto_ec3.5.1.24 <- flux_deconjug_BA_all[fluxes_deconjug_BA_trp_rxn03095_proto_ec3.5.1.24, ]
nrow(fluxes_deconjug_BA_trp_rxn03095_proto_ec3.5.1.24)
#number trp auxotrophic and non producer
sum_no_rxn03095_auxo_ec3.5.1.24 <- nrow(numb_auxo) - nrow(fluxes_deconjug_BA_trp_rxn03095_auxo_ec3.5.1.24)
sum_no_rxn03095_auxo_ec3.5.1.24 
#number trp prototrophic and non producer
sum_no_rxn03095_proto_ec3.5.1.24 <- 5416 - (nrow(fluxes_deconjug_BA_trp_rxn03095_auxo_ec3.5.1.24) + nrow(fluxes_deconjug_BA_trp_rxn03095_proto_ec3.5.1.24) + sum_no_rxn03095_auxo_ec3.5.1.24)
sum_no_rxn03095_proto_ec3.5.1.24

#create contigency table
fisher_rxn03095_ec_3.5.1.24 <- data.frame("Trp-Auxotrophy" = c(nrow(fluxes_deconjug_BA_trp_rxn03095_auxo_ec3.5.1.24), sum_no_rxn03095_auxo_ec3.5.1.24), 
                              "Trp-Prototrophy"= c(nrow(fluxes_deconjug_BA_trp_rxn03095_proto_ec3.5.1.24), sum_no_rxn03095_proto_ec3.5.1.24), 
                              row.names = c(" rxn03095: yes", " rxn03095: no"))

fisher_rxn03095_ec_3.5.1.24
#control numbers
nrow(fluxes_deconjug_BA_trp_rxn03095_auxo_ec3.5.1.24) + sum_no_rxn03095_auxo_ec3.5.1.24 + nrow(fluxes_deconjug_BA_trp_rxn03095_proto_ec3.5.1.24) + sum_no_rxn03095_proto_ec3.5.1.24 == 5416

#create mosaic plot
mosaicplot(fisher_rxn03095_ec_3.5.1.24, color = TRUE, main=" rxn03095 (ec 3.5.1.24) - trp auxotrophic and prototrophic microbiota")
#fisher test
fisher.test(fisher_rxn03095_ec_3.5.1.24)

#number of auxotrophic bacteria  rxn03095, ec 3.5.1.74
fluxes_deconjug_BA_trp_rxn03095_auxo_ec3.5.1.74 <- flux_deconjug_BA_all$Compound == "Trp" & flux_deconjug_BA_all$status == "good_blast" &
  flux_deconjug_BA_all$seedID == "rxn03095" & flux_deconjug_BA_all$Prototrophy == 0 & flux_deconjug_BA_all$ec == "3.5.1.74"
fluxes_deconjug_BA_trp_rxn03095_auxo_ec3.5.1.74 <- flux_deconjug_BA_all[fluxes_deconjug_BA_trp_rxn03095_auxo_ec3.5.1.74 , ]
nrow(fluxes_deconjug_BA_trp_rxn03095_auxo_ec3.5.1.74 )
#number prototrophic bacteria  rxn03095, , ec 3.5.1.74
fluxes_deconjug_BA_trp_rxn03095_proto_ec3.5.1.74 <- flux_deconjug_BA_all$Compound == "Trp" & flux_deconjug_BA_all$status == "good_blast" &
  flux_deconjug_BA_all$seedID == "rxn03095" & flux_deconjug_BA_all$Prototrophy == 1 & flux_deconjug_BA_all$ec == "3.5.1.74"
fluxes_deconjug_BA_trp_rxn03095_proto_ec3.5.1.74 <- flux_deconjug_BA_all[fluxes_deconjug_BA_trp_rxn03095_proto_ec3.5.1.74, ]
nrow(fluxes_deconjug_BA_trp_rxn03095_proto_ec3.5.1.74)
#number trp auxotrophic and non producer
sum_no_rxn03095_auxo_ec3.5.1.74 <- nrow(numb_auxo) - nrow(fluxes_deconjug_BA_trp_rxn03095_auxo_ec3.5.1.74)
sum_no_rxn03095_auxo_ec3.5.1.74 
#number trp prototrophic and non producer
sum_no_rxn03095_proto_ec3.5.1.74 <- 5416 - (nrow(fluxes_deconjug_BA_trp_rxn03095_auxo_ec3.5.1.74) + nrow(fluxes_deconjug_BA_trp_rxn03095_proto_ec3.5.1.74) + sum_no_rxn03095_auxo_ec3.5.1.74)
sum_no_rxn03095_proto_ec3.5.1.74

#create contigency table
fisher_rxn03095_ec3.5.1.74 <- data.frame("Trp-Auxotrophy" = c(nrow(fluxes_deconjug_BA_trp_rxn03095_auxo_ec3.5.1.74), sum_no_rxn03095_auxo_ec3.5.1.74), 
                              "Trp-Prototrophy"= c(nrow(fluxes_deconjug_BA_trp_rxn03095_proto_ec3.5.1.74), sum_no_rxn03095_proto_ec3.5.1.74), 
                              row.names = c(" rxn03095: yes", " rxn03095: no"))

fisher_rxn03095_ec3.5.1.74
#control numbers
nrow(fluxes_deconjug_BA_trp_rxn03095_auxo_ec3.5.1.74) + sum_no_rxn03095_auxo_ec3.5.1.74 + nrow(fluxes_deconjug_BA_trp_rxn03095_proto_ec3.5.1.74) + sum_no_rxn03095_proto_ec3.5.1.74 == 5416

#create mosaic plot
mosaicplot(fisher_rxn03095_ec3.5.1.74, color = TRUE, main=" rxn03095 (ec 3.5.1.74) - trp auxotrophic and prototrophic microbiota")
#fisher test
fisher.test(fisher_rxn03095_ec3.5.1.74)

#fisher test
#number of auxotrophic bacteria  rxn04068 (3.5.1.24)
fluxes_deconjug_BA_trp_rxn04068_auxo <- flux_deconjug_BA_all$Compound == "Trp" & flux_deconjug_BA_all$status == "good_blast" &
  flux_deconjug_BA_all$seedID == "rxn04068" & flux_deconjug_BA_all$Prototrophy == 0 
fluxes_deconjug_BA_trp_rxn04068_auxo <- flux_deconjug_BA_all[fluxes_deconjug_BA_trp_rxn04068_auxo, ]
nrow(fluxes_deconjug_BA_trp_rxn04068_auxo)
#number prototrophic bacteria  rxn04068
fluxes_deconjug_BA_trp_rxn04068_proto <- flux_deconjug_BA_all$Compound == "Trp" & flux_deconjug_BA_all$status == "good_blast" &
  flux_deconjug_BA_all$seedID == "rxn04068" & flux_deconjug_BA_all$Prototrophy == 1 
fluxes_deconjug_BA_trp_rxn04068_proto <- flux_deconjug_BA_all[fluxes_deconjug_BA_trp_rxn04068_proto, ]
nrow(fluxes_deconjug_BA_trp_rxn04068_proto)
#number trp auxotrophic and non producer
sum_no_rxn04068_auxo <- nrow(numb_auxo) - nrow(fluxes_deconjug_BA_trp_rxn04068_auxo)
sum_no_rxn04068_auxo 
#number trp prototrophic and non producer
sum_no_rxn04068_proto <- 5416 - (nrow(fluxes_deconjug_BA_trp_rxn04068_auxo) + nrow(fluxes_deconjug_BA_trp_rxn04068_proto) + sum_no_rxn04068_auxo)
sum_no_rxn04068_proto 

#create contigency table
fisher_rxn04068 <- data.frame("Trp-Auxotrophy" = c(nrow(fluxes_deconjug_BA_trp_rxn04068_auxo), sum_no_rxn04068_auxo), 
                              "Trp-Prototrophy"= c(nrow(fluxes_deconjug_BA_trp_rxn04068_proto), sum_no_rxn04068_proto), 
                              row.names = c(" rxn04068: yes", " rxn04068: no"))

fisher_rxn04068
#control numbers
nrow(fluxes_deconjug_BA_trp_rxn04068_auxo) + sum_no_rxn04068_auxo + nrow(fluxes_deconjug_BA_trp_rxn04068_proto) + sum_no_rxn04068_proto == 5416

#create mosaic plot
mosaicplot(fisher_rxn04068, color = TRUE, main=" rxn04068 (3.5.1.24) - trp auxotrophic and prototrophic microbiota")
#fisher test
fisher.test(fisher_rxn04068)


#percentage completeness pathways
#ec 3.5.1.24: rxn02795, rxn03094, rxn04068 --> same number of genomes with good blast --> same number of genomes no blast
#sum_no_rxn02795_auxo/ nrow(numb_auxo) * 100
#sum_no_rxn03094_auxo/ nrow(numb_auxo) * 100
#found enzymes
rxn04068_ec3.5.1.24_2 <- nrow(fluxes_deconjug_BA_trp_rxn04068_auxo)/ nrow(numb_auxo) * 100
rxn03095_ec3.5.1.24_2 <- nrow(fluxes_deconjug_BA_trp_rxn03095_auxo_ec3.5.1.24)/ nrow(numb_auxo) * 100
ec3.5.1.24_2 <- (rxn04068_ec3.5.1.24_2 + rxn03095_ec3.5.1.24_2) / 2

#ec 3.5.1.74: rxn02015, rxn02796, rxn03095
#w5 <- sum_no_rxn02015_auxo/ nrow(numb_auxo) * 100
rxn02796_ec3.5.1.74_2 <-  nrow(fluxes_deconjug_BA_trp_rxn02796_auxo)/ nrow(numb_auxo) * 100
rxn03095_ec3.5.1.74_2 <-  nrow(fluxes_deconjug_BA_trp_rxn03095_auxo_ec3.5.1.74)/ nrow(numb_auxo) * 100
ec3.5.1.74_2 <- (rxn02796_ec3.5.1.74_2 +rxn03095_ec3.5.1.74_2) / 2

new2 <- data.frame(rxn04068_ec3.5.1.24_2, rxn03095_ec3.5.1.24_2, rxn02796_ec3.5.1.74_2,rxn03095_ec3.5.1.74_2)
new2 <- t(new2)
colnames(new2) <- "percentage"
new2 <- data.frame(new2)
new2$new <- rownames(new2)

ggplot (new2,aes(new, percentage)) +
  geom_bar(stat = "identity") +
  ggtitle("Number of percentage of enzymes found in the pathway of bile acid deconjugation(BSH)\nof trp auxotrophic microbiota")

bsh2 <- data.frame(ec3.5.1.24_2,ec3.5.1.74_2)
bsh2 <- t(bsh2)
colnames(bsh2) <- "percentage"
bsh2 <- data.frame(bsh2)
bsh2$ec <- rownames(bsh2)

ggplot (bsh2,aes(ec, percentage)) +
  geom_bar(stat = "identity") +
  ggtitle("Completeness of the bile acid deconjugation(BSH) pathway in tryptophan auxotrophic microbiota") +
  xlab("enzymes") +
  ylab("pathway completeness [%]") +
  ylim(0,100)


#not found enzymes
rxn04068_ec3.5.1.24 <- sum_no_rxn04068_auxo/ nrow(numb_auxo) * 100
rxn03095_ec3.5.1.24 <- sum_no_rxn03095_auxo_ec3.5.1.24/ nrow(numb_auxo) * 100
ec3.5.1.24 <- (rxn04068_ec3.5.1.24 + rxn03095_ec3.5.1.24) / 2

#ec 3.5.1.74: rxn02015, rxn02796, rxn03095
#w5 <- sum_no_rxn02015_auxo/ nrow(numb_auxo) * 100
rxn02796_ec3.5.1.74 <- sum_no_rxn02796_auxo/ nrow(numb_auxo) * 100
rxn03095_ec3.5.1.74 <- sum_no_rxn03095_auxo_ec3.5.1.74/ nrow(numb_auxo) * 100
ec3.5.1.74 <- (rxn02796_ec3.5.1.74 +rxn03095_ec3.5.1.74) / 2

new <- data.frame(rxn04068_ec3.5.1.24, rxn03095_ec3.5.1.24, rxn02796_ec3.5.1.74, rxn03095_ec3.5.1.74)
new1 <- t(new)
colnames(new1) <- "percentage"
new1 <- data.frame(new1)
new1$new <- rownames(new1)


bsh <- data.frame(ec3.5.1.24,ec3.5.1.74)
bsh <- t(bsh)
colnames(bsh) <- "percentage"
bsh <- data.frame(bsh)
bsh$ec <- rownames(bsh)

ggplot (bsh,aes(ec, percentage)) +
  geom_bar(stat = "identity") +
  ggtitle("Incompleteness of the bile acid deconjugation (BSH) pathway in tryptophan auxotrophic microbiota")+
  xlab("enzymes") +
  ylab("pathway incompleteness [%]") +
  ylim(0,100)


#bile acid beta dehydroxylation
flu_beta_dehydroxy_BA <- flux_i$ec == "6.2.1.7" | flux_i$ec == "1.1.1.395" | flux_i$ec == "1.3.1.116" | 
  flux_i$ec == "4.2.1.M3" | flux_i$ec == "2.8.3.25" | flux_i$ec == "1.3.1.114" | flux_i$ec == "1.1.1.52"   

flu_beta_dehydroxy_BA<- flux_i[  flu_beta_dehydroxy_BA, ]
View(flu_beta_dehydroxy_BA)

fwrite(flu_beta_dehydroxy_BA, file="flu_beta_dehydroxy_BA.csv")

#mergen mit tryptophan auxotrophie daten
flux_beta_dehydroxy_BA_all <- merge(Auxotrophy_2, flu_beta_dehydroxy_BA, by.x = "Genomes", by.y = "model", allow.cartesian = TRUE)
View(flux_beta_dehydroxy_BA_all)

fluxes_beta_dehydroxy_BA_trp_auxo <- flux_beta_dehydroxy_BA_all$Prototrophy == 0 & 
  flux_beta_dehydroxy_BA_all$Compound == "Trp" &  flux_beta_dehydroxy_BA_all$status == "good_blast"
fluxes_beta_dehydroxy_BA_trp_auxo <-  flux_beta_dehydroxy_BA_all[fluxes_beta_dehydroxy_BA_trp_auxo, ]
View(fluxes_beta_dehydroxy_BA_trp_auxo)

fwrite(fluxes_beta_dehydroxy_BA_trp_auxo, file="fluxes_beta_dehydroxy_BA_trp_auxo.csv")
View(  fluxes_beta_dehydroxy_BA_trp_auxo)

#visualization
ggplot(data=  fluxes_beta_dehydroxy_BA_trp_auxo, aes(rxn, fill = phylum)) +
  geom_bar() +
  ggtitle("Reactions with a good blast for bile acid beta dehydroxylation\nby tryptophan auxotrophic microbiota") +
  ylab("number of genomes") +
  xlab("rxn") +
  theme(axis.text.x = element_text (angle=90))

ggplot(data=  fluxes_beta_dehydroxy_BA_trp_auxo, aes(seed, fill = ec)) +
  geom_bar() +
  ggtitle("Enzymes with a good blast for bile acid beta dehydroxylation\nby tryptophan auxotrophic microbiota") +
  ylab("number of genomes") +
  xlab("seed") +
  theme(axis.text.x = element_text (angle=90))


#fisher test/comparison trp auxotrophic and prototrophic
#CAUTION: due to the fact that the number of genomes are the same for every reactions of the enzymes(see to the last ggplot)
#for one reaction of every enzyme the fisher t-test will be done because the results would be the same
#number of auxotrophic bacteria  rxn02890 (ec 1.1.1.52)
fluxes_beta_dehydroxy_BA_trp_rxn02890_auxo <- flux_beta_dehydroxy_BA_all$Compound == "Trp" & flux_beta_dehydroxy_BA_all $status == "good_blast" &
  flux_beta_dehydroxy_BA_all$seed == "rxn02890" & flux_beta_dehydroxy_BA_all$Prototrophy == 0 
fluxes_beta_dehydroxy_BA_trp_rxn02890_auxo <- flux_beta_dehydroxy_BA_all[fluxes_beta_dehydroxy_BA_trp_rxn02890_auxo, ]
nrow(fluxes_beta_dehydroxy_BA_trp_rxn02890_auxo)
#number prototrophic bacteria  rxn02890
fluxes_beta_dehydroxy_BA_trp_rxn02890_proto <- flux_beta_dehydroxy_BA_all$Compound == "Trp" & flux_beta_dehydroxy_BA_all$status == "good_blast" &
  flux_beta_dehydroxy_BA_all$seed == "rxn02890" & flux_beta_dehydroxy_BA_all$Prototrophy == 1 
fluxes_beta_dehydroxy_BA_trp_rxn02890_proto <- flux_beta_dehydroxy_BA_all[fluxes_beta_dehydroxy_BA_trp_rxn02890_proto, ]
nrow(fluxes_beta_dehydroxy_BA_trp_rxn02890_proto)
#number trp auxotrophic and non producer
sum_no_rxn02890_auxo <- nrow(numb_auxo) - nrow(fluxes_beta_dehydroxy_BA_trp_rxn02890_auxo)
sum_no_rxn02890_auxo 
#number trp prototrophic and non producer
sum_no_rxn02890_proto <- 5416 - (nrow(fluxes_beta_dehydroxy_BA_trp_rxn02890_auxo) + nrow(fluxes_beta_dehydroxy_BA_trp_rxn02890_proto) + sum_no_rxn02890_auxo)
sum_no_rxn02890_proto 

#create contigency table
fisher_ec1.1.1.52 <- data.frame("Trp-Auxotrophy" = c(nrow(fluxes_beta_dehydroxy_BA_trp_rxn02890_auxo), sum_no_rxn02890_auxo), 
                              "Trp-Prototrophy"= c(nrow(fluxes_beta_dehydroxy_BA_trp_rxn02890_proto), sum_no_rxn02890_proto), 
                              row.names = c(" rxn02890: yes", " rxn02890: no"))

fisher_ec1.1.1.52
#control numbers
nrow(fluxes_beta_dehydroxy_BA_trp_rxn02890_auxo) + sum_no_rxn02890_auxo + nrow(fluxes_beta_dehydroxy_BA_trp_rxn02890_proto) + sum_no_rxn02890_proto == 5416

#create mosaic plot
mosaicplot(fisher_ec1.1.1.52, color = TRUE, main=" rxn02890 (ec 1.1.1.52) - trp auxotrophic and prototrophic microbiota")
#fisher test
f1 <- fisher.test(fisher_ec1.1.1.52)
pvalue_ec1.1.1.52 <- f1$p.value


#fisher test/comparison trp auxotrophic and prototrophic
#CAUTION: due to the fact that the number of genomes are the same for every reactions of the enzymes(see to the last ggplot)
#for one reaction of every enzyme the fisher t-test will be done because the results would be the same
#number of auxotrophic bacteria  rxn36658 (ec 1.1.1.395)
fluxes_beta_dehydroxy_BA_trp_rxn36658_auxo <- flux_beta_dehydroxy_BA_all$Compound == "Trp" & flux_beta_dehydroxy_BA_all $status == "good_blast" &
  flux_beta_dehydroxy_BA_all$seed == "rxn36658" & flux_beta_dehydroxy_BA_all$Prototrophy == 0 
fluxes_beta_dehydroxy_BA_trp_rxn36658_auxo <- flux_beta_dehydroxy_BA_all[fluxes_beta_dehydroxy_BA_trp_rxn36658_auxo, ]
nrow(fluxes_beta_dehydroxy_BA_trp_rxn36658_auxo)
#number prototrophic bacteria  rxn36658 
fluxes_beta_dehydroxy_BA_trp_rxn36658_proto <- flux_beta_dehydroxy_BA_all$Compound == "Trp" & flux_beta_dehydroxy_BA_all$status == "good_blast" &
  flux_beta_dehydroxy_BA_all$seed == "rxn36658" & flux_beta_dehydroxy_BA_all$Prototrophy == 1 
fluxes_beta_dehydroxy_BA_trp_rxn36658_proto <- flux_beta_dehydroxy_BA_all[fluxes_beta_dehydroxy_BA_trp_rxn36658_proto, ]
nrow(fluxes_beta_dehydroxy_BA_trp_rxn36658_proto)
#number trp auxotrophic and non producer
sum_no_rxn36658_auxo <- nrow(numb_auxo) - nrow(fluxes_beta_dehydroxy_BA_trp_rxn36658_auxo)
sum_no_rxn36658_auxo 
#number trp prototrophic and non producer
sum_no_rxn36658_proto <- 5416 - (nrow(fluxes_beta_dehydroxy_BA_trp_rxn36658_auxo) + nrow(fluxes_beta_dehydroxy_BA_trp_rxn36658_proto) + sum_no_rxn36658_auxo)
sum_no_rxn36658_proto 

#create contigency table
fisher_ec1.1.1.395  <- data.frame("Trp-Auxotrophy" = c(nrow(fluxes_beta_dehydroxy_BA_trp_rxn36658_auxo), sum_no_rxn36658_auxo), 
                              "Trp-Prototrophy"= c(nrow(fluxes_beta_dehydroxy_BA_trp_rxn36658_proto), sum_no_rxn36658_proto), 
                              row.names = c(" rxn36658 : yes", " rxn36658 : no"))

fisher_ec1.1.1.395
#control numbers
nrow(fluxes_beta_dehydroxy_BA_trp_rxn36658_auxo) + sum_no_rxn36658_auxo + nrow(fluxes_beta_dehydroxy_BA_trp_rxn36658_proto) + sum_no_rxn36658_proto == 5416

#create mosaic plot
mosaicplot(fisher_ec1.1.1.395 , color = TRUE, main=" rxn36658  (ec 1.1.1.395) - trp auxotrophic and prototrophic microbiota")
#fisher test
fisher.test(fisher_ec1.1.1.395)
#fisher test
f2 <- fisher.test(fisher_ec1.1.1.395)
pvalue_ec1.1.1.395 <- f2$p.value 

if(pvalue_ec1.1.1.395 < 0.05){ 
print("p <0.05")
} else {
  "p > 0.05"
}


#fisher test/comparison trp auxotrophic and prototrophic
#CAUTION: due to the fact that the number of genomes are the same for every reactions of the enzymes(see to the last ggplot)
#for one reaction of every enzyme the fisher t-test will be done because the results would be the same
#number of auxotrophic bacteria  rxnrxn41554 (ec 1.3.1.116)
fluxes_beta_dehydroxy_BA_trp_rxn41554_auxo <- flux_beta_dehydroxy_BA_all$Compound == "Trp" & flux_beta_dehydroxy_BA_all $status == "good_blast" &
  flux_beta_dehydroxy_BA_all$seed == "rxn41554" & flux_beta_dehydroxy_BA_all$Prototrophy == 0 
fluxes_beta_dehydroxy_BA_trp_rxn41554_auxo <- flux_beta_dehydroxy_BA_all[fluxes_beta_dehydroxy_BA_trp_rxn41554_auxo, ]
nrow(fluxes_beta_dehydroxy_BA_trp_rxn41554_auxo)
#number prototrophic bacteria  rxn41554  
fluxes_beta_dehydroxy_BA_trp_rxn41554_proto <- flux_beta_dehydroxy_BA_all$Compound == "Trp" & flux_beta_dehydroxy_BA_all$status == "good_blast" &
  flux_beta_dehydroxy_BA_all$seed == "rxn41554" & flux_beta_dehydroxy_BA_all$Prototrophy == 1 
fluxes_beta_dehydroxy_BA_trp_rxn41554_proto <- flux_beta_dehydroxy_BA_all[fluxes_beta_dehydroxy_BA_trp_rxn41554_proto, ]
nrow(fluxes_beta_dehydroxy_BA_trp_rxn41554_proto)
#number trp auxotrophic and non producer
sum_no_rxn41554_auxo <- nrow(numb_auxo) - nrow(fluxes_beta_dehydroxy_BA_trp_rxn41554_auxo)
sum_no_rxn41554_auxo 
#number trp prototrophic and non producer
sum_no_rxn41554_proto <- 5416 - (nrow(fluxes_beta_dehydroxy_BA_trp_rxn41554_auxo) + nrow(fluxes_beta_dehydroxy_BA_trp_rxn41554_proto) + sum_no_rxn41554_auxo)
sum_no_rxn41554_proto 

#create contigency table
fisher_ec1.3.1.116 <- data.frame("Trp-Auxotrophy" = c(nrow(fluxes_beta_dehydroxy_BA_trp_rxn41554_auxo), sum_no_rxn41554_auxo), 
                               "Trp-Prototrophy"= c(nrow(fluxes_beta_dehydroxy_BA_trp_rxn41554_proto), sum_no_rxn41554_proto), 
                               row.names = c(" rxn41554  : yes", " rxn41554  : no"))

fisher_ec1.3.1.116
#control numbers
nrow(fluxes_beta_dehydroxy_BA_trp_rxn41554_auxo) + sum_no_rxn41554_auxo + nrow(fluxes_beta_dehydroxy_BA_trp_rxn41554_proto) + sum_no_rxn41554_proto == 5416

#create mosaic plot
mosaicplot(fisher_ec1.3.1.116  , color = TRUE, main=" rxn41554   (ec 1.3.1.116) - trp auxotrophic and prototrophic microbiota")
#fisher test
fisher.test(fisher_ec1.3.1.116)
#pvalue
f3 <- fisher.test(fisher_ec1.3.1.116)
pvalue_ec1.3.1.116 <- f3$p.value

#fisher test/comparison trp auxotrophic and prototrophic
#CAUTION: due to the fact that the number of genomes are the same for every reactions of the enzymes(see to the last ggplot)
#for one reaction of every enzyme the fisher t-test will be done because the results would be the same
#number of auxotrophic bacteria  rxn43871 (ec 1.3.1.114)
fluxes_beta_dehydroxy_BA_trp_rxn43871_auxo <- flux_beta_dehydroxy_BA_all$Compound == "Trp" & flux_beta_dehydroxy_BA_all $status == "good_blast" &
  flux_beta_dehydroxy_BA_all$seed == "rxn43871" & flux_beta_dehydroxy_BA_all$Prototrophy == 0 
fluxes_beta_dehydroxy_BA_trp_rxn43871_auxo <- flux_beta_dehydroxy_BA_all[fluxes_beta_dehydroxy_BA_trp_rxn43871_auxo, ]
nrow(fluxes_beta_dehydroxy_BA_trp_rxn43871_auxo)
#number prototrophic bacteria  rxn43871   
fluxes_beta_dehydroxy_BA_trp_rxn43871_proto <- flux_beta_dehydroxy_BA_all$Compound == "Trp" & flux_beta_dehydroxy_BA_all$status == "good_blast" &
  flux_beta_dehydroxy_BA_all$seed == "rxn43871" & flux_beta_dehydroxy_BA_all$Prototrophy == 1 
fluxes_beta_dehydroxy_BA_trp_rxn43871_proto <- flux_beta_dehydroxy_BA_all[fluxes_beta_dehydroxy_BA_trp_rxn43871_proto, ]
nrow(fluxes_beta_dehydroxy_BA_trp_rxn43871_proto)
#number trp auxotrophic and non producer
sum_no_rxn43871_auxo <- nrow(numb_auxo) - nrow(fluxes_beta_dehydroxy_BA_trp_rxn43871_auxo)
sum_no_rxn43871_auxo 
#number trp prototrophic and non producer
sum_no_rxn43871_proto <- 5416 - (nrow(fluxes_beta_dehydroxy_BA_trp_rxn43871_auxo) + nrow(fluxes_beta_dehydroxy_BA_trp_rxn43871_proto) + sum_no_rxn43871_auxo)
sum_no_rxn43871_proto 

#create contigency table
fisher_ec1.3.1.114 <- data.frame("Trp-Auxotrophy" = c(nrow(fluxes_beta_dehydroxy_BA_trp_rxn43871_auxo), sum_no_rxn43871_auxo), 
                                "Trp-Prototrophy"= c(nrow(fluxes_beta_dehydroxy_BA_trp_rxn43871_proto), sum_no_rxn43871_proto), 
                                row.names = c(" rxn43871   : yes", " rxn43871   : no"))

fisher_ec1.3.1.114  
#control numbers
nrow(fluxes_beta_dehydroxy_BA_trp_rxn43871_auxo) + sum_no_rxn43871_auxo + nrow(fluxes_beta_dehydroxy_BA_trp_rxn43871_proto) + sum_no_rxn43871_proto == 5416

#create mosaic plot
mosaic5 <- mosaicplot(fisher_rxn43871, color = TRUE, main=" rxn43871 (ec 1.3.1.114) - trp auxotrophic and prototrophic microbiota")
mosaic6 <- mosaicplot(fisher_rxn43871, color = TRUE, main=" rxn43871 (ec 1.3.1.114) - trp auxotrophic and prototrophic microbiota")


#fisher test
fisher.test(fisher_rxn43871)
#pvalue
f4 <- fisher.test(fisher_ec1.3.1.114)
pvalue_ec1.3.1.114 <- f4$p.value


#completeness pathways
#percentage completeness pathways
#found enzymes
ec1.1.1.52 <- nrow(fluxes_beta_dehydroxy_BA_trp_rxn02890_auxo)/ nrow(numb_auxo) * 100
ec1.1.1.395 <- nrow(fluxes_beta_dehydroxy_BA_trp_rxn36658_auxo)/ nrow(numb_auxo) * 100
ec1.3.1.114 <-  nrow(fluxes_beta_dehydroxy_BA_trp_rxn41554_auxo)/ nrow(numb_auxo) * 100
ec1.3.1.116 <-  nrow(fluxes_beta_dehydroxy_BA_trp_rxn43871_auxo)/ nrow(numb_auxo) * 100

#originally the pathway exists of 7 enzymes, only four are shown in the plot so they are also missing but test
flux_i$ec == "6.2.1.7" | flux_i$ec == "1.1.1.395" | flux_i$ec == "1.3.1.116" | 
  flux_i$ec == "4.2.1.M3" | flux_i$ec == "2.8.3.25" | flux_i$ec == "1.3.1.114" | flux_i$ec == "1.1.1.52"

ec6.2.1.7 <-flux_beta_dehydroxy_BA_all$Compound == "Trp" & flux_beta_dehydroxy_BA_all $status == "good_blast" &
  flux_beta_dehydroxy_BA_all$ec == "6.2.1.7" & flux_beta_dehydroxy_BA_all$Prototrophy == 0 
ec6.2.1.7 <- flux_beta_dehydroxy_BA_all[ec6.2.1.7, ]
ec6.2.1.7 <- nrow(ec6.2.1.7)
ec4.2.1.M3 <-flux_beta_dehydroxy_BA_all$Compound == "Trp" & flux_beta_dehydroxy_BA_all $status == "good_blast" &
  flux_beta_dehydroxy_BA_all$ec == "4.2.1.M3" & flux_beta_dehydroxy_BA_all$Prototrophy == 0 
ec4.2.1.M3 <- flux_beta_dehydroxy_BA_all[ec4.2.1.M3, ]
ec4.2.1.M3  <- nrow(ec4.2.1.M3)
ec2.8.3.25 <-flux_beta_dehydroxy_BA_all$Compound == "Trp" & flux_beta_dehydroxy_BA_all $status == "good_blast" &
  flux_beta_dehydroxy_BA_all$ec == "2.8.3.25" & flux_beta_dehydroxy_BA_all$Prototrophy == 0 
ec2.8.3.25 <- flux_beta_dehydroxy_BA_all[ec2.8.3.25, ]
ec2.8.3.25 <- nrow(ec2.8.3.25)

dehydroxy <- data.frame(ec1.1.1.52,ec1.1.1.395,ec1.3.1.114,ec1.3.1.116,ec6.2.1.7, ec4.2.1.M3, ec2.8.3.25)
dehydroxy<- t(dehydroxy)
colnames(dehydroxy) <- "percentage"
dehydroxy<- data.frame(dehydroxy)
dehydroxy$new <- rownames(dehydroxy)
pvalue <- c(pvalue_ec1.1.1.52, pvalue_ec1.1.1.395, pvalue_ec1.3.1.114, pvalue_ec1.3.1.116, 0.05, 0.05, 0.05)
dehydroxy$pvalues <- pvalue

ggplot (dehydroxy,aes(new, percentage, label = ifelse(pvalues < 0.05, "*","NS"))) +
  geom_bar(stat = "identity") +
  ggtitle(expression(atop("Completeness of the beta dehydroxylation pathway in tryptophan auxotrophic microbiota", atop("added results of Fisher exact t-test for comparison to tryptophan prototrophic microbiota")))) +
  xlab("enzymes") +
  ylab("pathway completeness [%]") +
  ylim(0,100) +
  geom_text(vjust = -1) +
  labs(caption = "* - statistically significant to trp prototrophy microbiota\n\n NS - not statistically significant to trp protrotrophy microbiota")


#percentage completeness pathways
#found enzymes
ec1.1.1.52_no <- sum_no_rxn02890_auxo/ nrow(numb_auxo) * 100
ec1.1.1.395_no <- sum_no_rxn43871_auxo / nrow(numb_auxo) * 100
ec1.3.1.114_no <-  sum_no_rxn41554_auxo / nrow(numb_auxo) * 100
ec1.3.1.116_no <-   sum_no_rxn43871_auxo / nrow(numb_auxo) * 100

#originally the pathway exists of 7 enzymes, only four are shown in the plot so they are also missing but test
flux_i$ec == "6.2.1.7" | flux_i$ec == "1.1.1.395" | flux_i$ec == "1.3.1.116" | 
  flux_i$ec == "4.2.1.M3" | flux_i$ec == "2.8.3.25" | flux_i$ec == "1.3.1.114" | flux_i$ec == "1.1.1.52"

ec6.2.1.7 <-flux_beta_dehydroxy_BA_all$Compound == "Trp" & flux_beta_dehydroxy_BA_all $status == "good_blast" &
  flux_beta_dehydroxy_BA_all$ec == "6.2.1.7" & flux_beta_dehydroxy_BA_all$Prototrophy == 0 
ec6.2.1.7 <- flux_beta_dehydroxy_BA_all[ec6.2.1.7, ]
ec6.2.1.7 <- nrow(ec6.2.1.7)
ec4.2.1.M3 <-flux_beta_dehydroxy_BA_all$Compound == "Trp" & flux_beta_dehydroxy_BA_all $status == "good_blast" &
  flux_beta_dehydroxy_BA_all$ec == "4.2.1.M3" & flux_beta_dehydroxy_BA_all$Prototrophy == 0 
ec4.2.1.M3 <- flux_beta_dehydroxy_BA_all[ec4.2.1.M3, ]
ec4.2.1.M3  <- nrow(ec4.2.1.M3)
ec2.8.3.25 <-flux_beta_dehydroxy_BA_all$Compound == "Trp" & flux_beta_dehydroxy_BA_all $status == "good_blast" &
  flux_beta_dehydroxy_BA_all$ec == "2.8.3.25" & flux_beta_dehydroxy_BA_all$Prototrophy == 0 
ec2.8.3.25 <- flux_beta_dehydroxy_BA_all[ec2.8.3.25, ]
ec2.8.3.25 <- nrow(ec2.8.3.25)

beta_dehydroxy_no <- data.frame(ec1.1.1.52_no,ec1.1.1.395_no,ec1.3.1.114_no,ec1.3.1.116_no,100,100,100)
beta_dehydroxy_no <- t(beta_dehydroxy_no)
colnames(beta_dehydroxy_no) <- "percentage"
beta_dehydroxy_no<- data.frame(beta_dehydroxy_no)
beta_dehydroxy_no$new <- rownames(beta_dehydroxy_no)

ggplot (beta_dehydroxy_no,aes(new, percentage)) +
  geom_bar(stat = "identity") +
  ggtitle("Incompleteness of the beta dehydroxylation pathway in tryptophan auxotrophic microbiota") +
  scale_x_discrete(breaks=c("ec1.1.1.395_no","ec1.1.1.52_no","ec1.3.1.114_no","ec1.3.1.116_no","X100", "X100.1", "X100.2"), 
                   labels=c("ec1.1.1.395","ec1.1.1.52","ec1.3.1.114","ec1.3.1.116","ec2.8.3.25", "ec4.2.1.M3", "ec6.2.1.7")) +
  xlab("enzymes") +
  ylab("pathway incompleteness [%]") +
  ylim(0,100)

#bile acid epimerization
flu_epimerization_BA <- flux_i$ec == "1.1.1.159" | flux_i$ec == "1.1.1.176" | flux_i$ec == "1.1.1.238" | 
  flux_i$ec == "1.1.1.201" | flux_i$ec == "1.1.1.52" | flux_i$ec == "1.1.1.391" 

flu_epimerization_BA<- flux_i[flu_epimerization_BA, ]
View(flu_epimerization_BA)

fwrite(flu_epimerization_BA, file="flu_epimerization_BA.csv")

#mergen mit tryptophan auxotrophie daten
flux_epimerization_BA_all <- merge(Auxotrophy_2, flu_epimerization_BA, by.x = "Genomes", by.y = "model", allow.cartesian = TRUE)
View(flux_epimerization_BA_all)

fluxes_epimerization_BA_trp_auxo <- flux_epimerization_BA_all$Prototrophy == 0 & 
  flux_epimerization_BA_all$Compound == "Trp" &  flux_epimerization_BA_all$status == "good_blast"
fluxes_epimerization_BA_trp_auxo <-  flux_epimerization_BA_all[fluxes_epimerization_BA_trp_auxo, ]
View(fluxes_epimerization_BA_trp_auxo)

fwrite(fluxes_epimerization_BA_trp_auxo, file="fluxes_epimerization_BA_trp_auxo.csv")
View(  fluxes_epimerization_BA_trp_auxo)

#visualization
ggplot(data=  fluxes_epimerization_BA_trp_auxo, aes(rxn, fill = phylum)) +
  geom_bar() +
  ggtitle("Reactions with a good blast for bile acid beta epimerization\nby tryptophan auxotrophic microbiota") +
  ylab("number of genomes") +
  xlab("rxn") +
  theme(axis.text.x = element_text (angle=90))

ggplot(data=  fluxes_epimerization_BA_trp_auxo, aes(rxn, fill = ec)) +
  geom_bar() +
  ggtitle("Enzymes with a good blast for bile acid epimerization\n by tryptophan auxotrophic microbiota") +
  ylab("number of genomes") +
  xlab("rxn") +
  theme(axis.text.x = element_text (angle=90))

#fisher test/comparison trp auxotrophic and prototrophic
#CAUTION: due to the fact that the number of genomes are the same for every reactions of the enzymes(see to the last ggplot)
#for one reaction of every enzyme the fisher t-test will be done because the results would be the same
#number of auxotrophic bacteria  rxn02012 (ec 1.1.1.159)
fluxes_epimerization_BA_trp_rxn02012_auxo <- flux_epimerization_BA_all$Compound == "Trp" & flux_epimerization_BA_all$status == "good_blast" &
  flux_epimerization_BA_all$seed == "rxn02012" & flux_epimerization_BA_all$Prototrophy == 0 
fluxes_epimerization_BA_trp_rxn02012_auxo <- flux_epimerization_BA_all[fluxes_epimerization_BA_trp_rxn02012_auxo, ]
nrow(fluxes_epimerization_BA_trp_rxn02012_auxo)
#number prototrophic bacteria  rxn02012   
fluxes_epimerization_BA_trp_rxn02012_proto <- flux_epimerization_BA_all$Compound == "Trp" & flux_epimerization_BA_all$status == "good_blast" &
  flux_epimerization_BA_all$seed == "rxn02012" & flux_epimerization_BA_all$Prototrophy == 1 
fluxes_epimerization_BA_trp_rxn02012_proto <- flux_epimerization_BA_all[fluxes_epimerization_BA_trp_rxn02012_proto, ]
nrow(fluxes_epimerization_BA_trp_rxn02012_proto)
#number trp auxotrophic and non producer
sum_no_rxn02012_auxo <- nrow(numb_auxo) - nrow(fluxes_epimerization_BA_trp_rxn02012_auxo)
sum_no_rxn02012_auxo 
#number trp prototrophic and non producer
sum_no_rxn02012_proto <- 5416 - (nrow(fluxes_epimerization_BA_trp_rxn02012_auxo) + nrow(fluxes_epimerization_BA_trp_rxn02012_proto) + sum_no_rxn02012_auxo)
sum_no_rxn02012_proto 

#create contigency table
fisher_ec1.1.1.159 <- data.frame("Trp-Auxotrophy" = c(nrow(fluxes_epimerization_BA_trp_rxn02012_auxo), sum_no_rxn02012_auxo), 
                                 "Trp-Prototrophy"= c(nrow(fluxes_epimerization_BA_trp_rxn02012_proto), sum_no_rxn02012_proto), 
                                 row.names = c(" rxn02012   : yes", " rxn02012   : no"))

fisher_ec1.1.1.159
#pvalue
f5 <- fisher.test(fisher_ec1.1.1.159)
pvalue_ec1.1.1.159 <- f5$p.value

#fisher test/comparison trp auxotrophic and prototrophic
#CAUTION: due to the fact that the number of genomes are the same for every reactions of the enzymes(see to the last ggplot)
#for one reaction of every enzyme the fisher t-test will be done because the results would be the same
#number of auxotrophic bacteria  rxn02013 (ec 1.1.1.176)
fluxes_epimerization_BA_trp_rxn02013_auxo <- flux_epimerization_BA_all$Compound == "Trp" & flux_epimerization_BA_all$status == "good_blast" &
  flux_epimerization_BA_all$seed == "rxn02013" & flux_epimerization_BA_all$Prototrophy == 0 
fluxes_epimerization_BA_trp_rxn02013_auxo <- flux_epimerization_BA_all[fluxes_epimerization_BA_trp_rxn02013_auxo, ]
nrow(fluxes_epimerization_BA_trp_rxn02013_auxo)
#number prototrophic bacteria  rxn02013   
fluxes_epimerization_BA_trp_rxn02013_proto <- flux_epimerization_BA_all$Compound == "Trp" & flux_epimerization_BA_all$status == "good_blast" &
  flux_epimerization_BA_all$seed == "rxn02013" & flux_epimerization_BA_all$Prototrophy == 1 
fluxes_epimerization_BA_trp_rxn02013_proto <- flux_epimerization_BA_all[fluxes_epimerization_BA_trp_rxn02013_proto, ]
nrow(fluxes_epimerization_BA_trp_rxn02013_proto)
#number trp auxotrophic and non producer
sum_no_rxn02013_auxo <- nrow(numb_auxo) - nrow(fluxes_epimerization_BA_trp_rxn02013_auxo)
sum_no_rxn02013_auxo 
#number trp prototrophic and non producer
sum_no_rxn02013_proto <- 5416 - (nrow(fluxes_epimerization_BA_trp_rxn02013_auxo) + nrow(fluxes_epimerization_BA_trp_rxn02013_proto) + sum_no_rxn02013_auxo)
sum_no_rxn02013_proto 

#create contigency table
fisher_ec1.1.1.176   <- data.frame("Trp-Auxotrophy" = c(nrow(fluxes_epimerization_BA_trp_rxn02013_auxo), sum_no_rxn02013_auxo), 
                                 "Trp-Prototrophy"= c(nrow(fluxes_epimerization_BA_trp_rxn02013_proto), sum_no_rxn02013_proto), 
                                 row.names = c(" rxn02013   : yes", " rxn02013   : no"))

fisher_ec1.1.1.176 
#control numbers
nrow(fluxes_epimerization_BA_trp_rxn02013_auxo) + sum_no_rxn02013_auxo + nrow(fluxes_epimerization_BA_trp_rxn02013_proto) + sum_no_rxn02013_proto == 5416

#create mosaic plot
mosaic5 <- mosaicplot(fisher_rxn02013, color = TRUE, main=" rxn02013 (ec 1.1.1.176) - trp auxotrophic and prototrophic microbiota")

#fisher test
fisher.test(fisher_ec1.1.1.176)
#pvalue
f6 <- fisher.test(fisher_ec1.1.1.176)
pvalue_ec1.1.1.176<- f6$p.value

#fisher test/comparison trp auxotrophic and prototrophic
#CAUTION: due to the fact that the number of genomes are the same for every reactions of the enzymes(see to the last ggplot)
#for one reaction of every enzyme the fisher t-test will be done because the results would be the same
#number of auxotrophic bacteria  rxn02665 (ec 1.1.1.201)
fluxes_epimerization_BA_trp_rxn02665_auxo <- flux_epimerization_BA_all$Compound == "Trp" & flux_epimerization_BA_all$status == "good_blast" &
  flux_epimerization_BA_all$seed == "rxn02665" & flux_epimerization_BA_all$Prototrophy == 0 
fluxes_epimerization_BA_trp_rxn02665_auxo <- flux_epimerization_BA_all[fluxes_epimerization_BA_trp_rxn02665_auxo, ]
nrow(fluxes_epimerization_BA_trp_rxn02665_auxo)
#number prototrophic bacteria  rxn02665   
fluxes_epimerization_BA_trp_rxn02665_proto <- flux_epimerization_BA_all$Compound == "Trp" & flux_epimerization_BA_all$status == "good_blast" &
  flux_epimerization_BA_all$seed == "rxn02665" & flux_epimerization_BA_all$Prototrophy == 1 
fluxes_epimerization_BA_trp_rxn02665_proto <- flux_epimerization_BA_all[fluxes_epimerization_BA_trp_rxn02665_proto, ]
nrow(fluxes_epimerization_BA_trp_rxn02665_proto)
#number trp auxotrophic and non producer
sum_no_rxn02665_auxo <- nrow(numb_auxo) - nrow(fluxes_epimerization_BA_trp_rxn02665_auxo)
sum_no_rxn02665_auxo 
#number trp prototrophic and non producer
sum_no_rxn02665_proto <- 5416 - (nrow(fluxes_epimerization_BA_trp_rxn02665_auxo) + nrow(fluxes_epimerization_BA_trp_rxn02665_proto) + sum_no_rxn02665_auxo)
sum_no_rxn02665_proto 

#create contigency table
fisher_ec1.1.1.201 <- data.frame("Trp-Auxotrophy" = c(nrow(fluxes_epimerization_BA_trp_rxn02665_auxo), sum_no_rxn02665_auxo), 
                                 "Trp-Prototrophy"= c(nrow(fluxes_epimerization_BA_trp_rxn02665_proto), sum_no_rxn02665_proto), 
                                 row.names = c(" rxn02665   : yes", " rxn02665   : no"))

fisher_ec1.1.1.201
#control numbers
nrow(fluxes_epimerization_BA_trp_rxn02665_auxo) + sum_no_rxn02665_auxo + nrow(fluxes_epimerization_BA_trp_rxn02665_proto) + sum_no_rxn02665_proto == 5416

#create mosaic plot
mosaic5 <- mosaicplot(fisher_ec1.1.1.201, color = TRUE, main=" rxn02665 (ec 1.1.1.201) - trp auxotrophic and prototrophic microbiota")

#fisher test
fisher.test(fisher_ec1.1.1.201)
#pvalue
f7 <- fisher.test(fisher_ec1.1.1.201)
pvalue_ec1.1.1.201 <- f7$p.value

#fisher test/comparison trp auxotrophic and prototrophic
#CAUTION: due to the fact that the number of genomes are the same for every reactions of the enzymes(see to the last ggplot)
#for one reaction of every enzyme the fisher t-test will be done because the results would be the same
#number of auxotrophic bacteria  rxn02890 (ec 1.1.1.52)
fluxes_epimerization_BA_trp_rxn02890_auxo <- flux_epimerization_BA_all$Compound == "Trp" & flux_epimerization_BA_all $status == "good_blast" &
  flux_epimerization_BA_all$seed == "rxn02890" & flux_epimerization_BA_all$Prototrophy == 0 
fluxes_epimerization_BA_trp_rxn02890_auxo <- flux_epimerization_BA_all[fluxes_epimerization_BA_trp_rxn02890_auxo, ]
nrow(fluxes_epimerization_BA_trp_rxn02890_auxo)
#number prototrophic bacteria  rxn02890   
fluxes_epimerization_BA_trp_rxn02890_proto <- flux_epimerization_BA_all$Compound == "Trp" & flux_epimerization_BA_all$status == "good_blast" &
  flux_epimerization_BA_all$seed == "rxn02890" & flux_epimerization_BA_all$Prototrophy == 1 
fluxes_epimerization_BA_trp_rxn02890_proto <- flux_epimerization_BA_all[fluxes_epimerization_BA_trp_rxn02890_proto, ]
nrow(fluxes_epimerization_BA_trp_rxn02890_proto)
#number trp auxotrophic and non producer
sum_no_rxn02890_auxo <- nrow(numb_auxo) - nrow(fluxes_epimerization_BA_trp_rxn02890_auxo)
sum_no_rxn02890_auxo 
#number trp prototrophic and non producer
sum_no_rxn02890_proto <- 5416 - (nrow(fluxes_epimerization_BA_trp_rxn02890_auxo) + nrow(fluxes_epimerization_BA_trp_rxn02890_proto) + sum_no_rxn02890_auxo)
sum_no_rxn02890_proto 

#create contigency table
fisher_ec1.1.1.52   <- data.frame("Trp-Auxotrophy" = c(nrow(fluxes_epimerization_BA_trp_rxn02890_auxo), sum_no_rxn02890_auxo), 
                                 "Trp-Prototrophy"= c(nrow(fluxes_epimerization_BA_trp_rxn02890_proto), sum_no_rxn02890_proto), 
                                 row.names = c(" rxn02890   : yes", " rxn02890   : no"))

fisher_ec1.1.1.52   
#control numbers
nrow(fluxes_epimerization_BA_trp_rxn02890_auxo) + sum_no_rxn02890_auxo + nrow(fluxes_epimerization_BA_trp_rxn02890_proto) + sum_no_rxn02890_proto == 5416

#create mosaic plot
mosaic5 <- mosaicplot(fisher_ec1.1.1.52, color = TRUE, main=" rxn02890 (ec 1.1.1.52) - trp auxotrophic and prototrophic microbiota")

#fisher test
fisher.test(fisher_ec1.1.1.52)
#pvalue
f8 <- fisher.test(fisher_ec1.1.1.52)
pvalue_ec1.1.1.52<- f8$p.value

#fisher test/comparison trp auxotrophic and prototrophic
#CAUTION: due to the fact that the number of genomes are the same for every reactions of the enzymes(see to the last ggplot)
#for one reaction of every enzyme the fisher t-test will be done because the results would be the same
#number of auxotrophic bacteria  rxn21855 (ec 1.1.1.391)
fluxes_epimerization_BA_trp_rxn21855_auxo <- flux_epimerization_BA_all$Compound == "Trp" & flux_epimerization_BA_all $status == "good_blast" &
  flux_epimerization_BA_all$seed == "rxn21855" & flux_epimerization_BA_all$Prototrophy == 0 
fluxes_epimerization_BA_trp_rxn21855_auxo <- flux_epimerization_BA_all[fluxes_epimerization_BA_trp_rxn21855_auxo, ]
nrow(fluxes_epimerization_BA_trp_rxn21855_auxo)
#number prototrophic bacteria  rxn21855   
fluxes_epimerization_BA_trp_rxn21855_proto <- flux_epimerization_BA_all$Compound == "Trp" & flux_epimerization_BA_all$status == "good_blast" &
  flux_epimerization_BA_all$seed == "rxn21855" & flux_epimerization_BA_all$Prototrophy == 1 
fluxes_epimerization_BA_trp_rxn21855_proto <- flux_epimerization_BA_all[fluxes_epimerization_BA_trp_rxn21855_proto, ]
nrow(fluxes_epimerization_BA_trp_rxn21855_proto)
#number trp auxotrophic and non producer
sum_no_rxn21855_auxo <- nrow(numb_auxo) - nrow(fluxes_epimerization_BA_trp_rxn21855_auxo)
sum_no_rxn21855_auxo 
#number trp prototrophic and non producer
sum_no_rxn21855_proto <- 5416 - (nrow(fluxes_epimerization_BA_trp_rxn21855_auxo) + nrow(fluxes_epimerization_BA_trp_rxn21855_proto) + sum_no_rxn21855_auxo)
sum_no_rxn21855_proto 

#create contigency table
fisher_ec1.1.1.391 <- data.frame("Trp-Auxotrophy" = c(nrow(fluxes_epimerization_BA_trp_rxn21855_auxo), sum_no_rxn21855_auxo), 
                                 "Trp-Prototrophy"= c(nrow(fluxes_epimerization_BA_trp_rxn21855_proto), sum_no_rxn21855_proto), 
                                 row.names = c(" rxn21855   : yes", " rxn21855   : no"))

fisher_ec1.1.1.391

#control numbers
nrow(fluxes_epimerization_BA_trp_rxn21855_auxo) + sum_no_rxn21855_auxo + nrow(fluxes_epimerization_BA_trp_rxn21855_proto) + sum_no_rxn21855_proto == 5416

#create mosaic plot
mosaic5 <- mosaicplot(fisher_ec1.1.1.391, color = TRUE, main=" rxn21855 (ec 1.1.1.391) - trp auxotrophic and prototrophic microbiota")

#fisher test
fisher.test(fisher_ec1.1.1.391)
#pvalue
f9 <- fisher.test(fisher_ec1.1.1.391)
pvalue_ec1.1.1.391 <- f9$p.value


#completeness pathways
#percentage completeness pathways
#found enzymes
ec1.1.1.159 <- nrow(fluxes_epimerization_BA_trp_rxn02012_auxo)/ nrow(numb_auxo) * 100
ec1.1.1.176 <- nrow(fluxes_epimerization_BA_trp_rxn02013_auxo)/ nrow(numb_auxo) * 100
ec1.1.1.201 <- nrow(fluxes_epimerization_BA_trp_rxn02665_auxo)/ nrow(numb_auxo) * 100
ec1.1.1.391 <-  nrow(fluxes_epimerization_BA_trp_rxn21855_auxo)/ nrow(numb_auxo) * 100
ec1.1.1.52 <- nrow(fluxes_epimerization_BA_trp_rxn02890_auxo)/ nrow(numb_auxo) * 100

#originally the pathway exists of 6 enzymes, only four are shown in the plot so they are also missing but test
flu_epimerization_BA <- flux_i$ec == "1.1.1.159" | flux_i$ec == "1.1.1.176" | flux_i$ec == "1.1.1.238" | 
  flux_i$ec == "1.1.1.201" | flux_i$ec == "1.1.1.52" | flux_i$ec == "1.1.1.391" 

ec1.1.1.238 <- flux_epimerization_BA_all$Compound == "Trp" & flux_epimerization_BA_all$status == "good_blast" &
  flux_epimerization_BA_all$ec == "1.1.1.238" & flux_epimerization_BA_all$Prototrophy == 0 
ec1.1.1.238<- flux_epimerization_BA_all[ec1.1.1.238, ]
ec1.1.1.238 <- nrow(ec1.1.1.238)


epimerization <- data.frame(ec1.1.1.159,ec1.1.1.176,ec1.3.1.201,ec1.1.1.391,ec1.1.1.52, ec1.1.1.238)
epimerization <- t(epimerization)
colnames(epimerization) <- "percentage"
epimerization  <- data.frame(epimerization)
epimerization $new <- rownames(epimerization)
pvalue <- c(pvalue_ec1.1.1.159,pvalue_ec1.1.1.176,pvalue_ec1.1.1.201,pvalue_ec1.1.1.391,pvalue_ec1.1.1.52, 0.05)
epimerization$pvalues <- pvalue

ggplot (epimerization,aes(new, percentage, label = ifelse(pvalues < 0.05, "*","NS"))) +
  geom_bar(stat = "identity") +
  ggtitle(expression(atop("Completeness of the bile acid epimerization pathway in tryptophan auxotrophic microbiota", atop("added results of Fisher exact t-test for comparison to tryptophan prototrophic microbiota")))) +
  xlab("enzymes") +
  ylab("pathway completeness [%]") +
  ylim(0,100) +
  geom_text(vjust = -1) +
  labs(caption = "* - statistically significant to trp prototrophy microbiota\n\n NS - not statistically significant to trp protrotrophy microbiota")


#found enzymes
ec1.1.1.159_no <- sum_no_rxn21855_auxo/ nrow(numb_auxo) * 100
ec1.1.1.176_no <- sum_no_rxn02013_auxo / nrow(numb_auxo) * 100
ec1.3.1.201_no <- sum_no_rxn02665_auxo / nrow(numb_auxo) * 100
ec1.3.1.391_no <- sum_no_rxn21855_auxo / nrow(numb_auxo) * 100
ec1.3.1.52_no <- sum_no_rxn02890_auxo / nrow(numb_auxo) * 100

#originally the pathway exists of 6 enzymes, only four are shown in the plot so they are also missing but test
flu_epimerization_BA <- flux_i$ec == "1.1.1.159" | flux_i$ec == "1.1.1.176" | flux_i$ec == "1.1.1.238" | 
  flux_i$ec == "1.1.1.201" | flux_i$ec == "1.1.1.52" | flux_i$ec == "1.1.1.391" 

ec1.1.1.238 <- flux_epimerization_BA_all$Compound == "Trp" & flux_epimerization_BA_all$status == "good_blast" &
  flux_epimerization_BA_all$ec == "1.1.1.238" & flux_epimerization_BA_all$Prototrophy == 0 
ec1.1.1.238<- flux_epimerization_BA_all[ec1.1.1.238, ]
ec1.1.1.238 <- nrow(ec1.1.1.238)


epimerization_no <- data.frame(ec1.1.1.159_no,ec1.1.1.176_no,ec1.3.1.201_no,ec1.3.1.391_no,ec1.3.1.52_no, 100)
epimerization_no <- t(epimerization_no)
colnames(epimerization_no) <- "percentage"
epimerization_no  <- data.frame(epimerization_no)
epimerization_no$new <- rownames(epimerization_no)

ggplot (epimerization_no,aes(new, percentage)) +
  geom_bar(stat = "identity") +
ggtitle("Incompleteness of the bile acid epimerization pathway by tryptophan auxotrophic microbiota") +
  scale_x_discrete(breaks=c("ec1.1.1.159_no","ec1.1.1.176_no","ec1.3.1.201_no","ec1.3.1.391_no","ec1.3.1.52_no", "X100"), 
                   labels=c("ec1.1.1.159","ec1.1.1.176","ec1.3.1.201","ec1.3.1.391","ec1.3.1.52", "ec1.1.1.238")) +
  xlab("enzymes") +
  ylab("pathway incompleteness [%]") +
  ylim(0,100)

#ES FEHLT ALPHA dehyroxylation!!!!
flux_i <- read.csv("/Users/svenjabusche/Desktop/flux_i.csv")

flu_alpha_dehydroxy_BA <- flux_i$ec == "6.2.1.7" | flux_i$ec == "1.1.1.395" | flux_i$ec == "1.3.1.115" | 
  flux_i$ec == "2.8.3.25" | flux_i$ec == "4.2.1.106" | flux_i$ec == "1.3.1.114"  | flux_i$ec == "1.1.1.52" 

flu_alpha_dehydroxy_BA<- flux_i[flu_alpha_dehydroxy_BA, ]
View(flu_alpha_dehydroxy_BA)

ggplot(flu_alpha_dehydroxy_BA, aes(rxn, fill = ec))+
  geom_bar()

#mergen mit tryptophan auxotrophie daten
flux_alpha_dehydroxy_BA_all <- merge(Auxotrophy_2, flu_alpha_dehydroxy_BA, by.x = "Genomes", by.y = "model", allow.cartesian = TRUE)


fluxes_alpha_dehydroxy_BA_trp_auxo <- flux_alpha_dehydroxy_BA_all$Prototrophy == 0 & 
  flux_alpha_dehydroxy_BA_all$Compound == "Trp" &  flux_alpha_dehydroxy_BA_all$status == "good_blast"
fluxes_alpha_dehydroxy_BA_trp_auxo <-  flux_alpha_dehydroxy_BA_all[fluxes_alpha_dehydroxy_BA_trp_auxo, ]
View(fluxes_alpha_dehydroxy_BA_trp_auxo)

View(  fluxes_alpha_dehydroxy_BA_trp_auxo)

#visualization
ggplot(data=  fluxes_alpha_dehydroxy_BA_trp_auxo, aes(rxn, fill = phylum)) +
  geom_bar() +
  ggtitle("Reactions with a good blast for bile acid alpha dehydroxylation \nby tryptophan auxotrophic microbiota") +
  ylab("number of genomes") +
  xlab("rxn") +
  theme(axis.text.x = element_text (angle=90)) +
  theme(legend.text = element_text (size=8)) +
  theme(legend.position = "right") 


ggplot(data=  fluxes_alpha_dehydroxy_BA_trp_auxo, aes(seed, fill = ec)) +
  geom_bar() +
  ggtitle("Enzymes with a good blast for bile acid alpha dehydroxylation \nby tryptophan auxotrophic microbiota") +
  ylab("number of genomes") +
  xlab("rxn") +
  theme(axis.text.x = element_text (angle=90))


#fisher test/comparison trp auxotrophic and prototrophic
#CAUTION: due to the fact that the number of genomes are the same for every reactions of the enzymes(see to the last ggplot)
#for one reaction of every enzyme the fisher t-test will be done because the results would be the same
#number of auxotrophic bacteria  rxn02890 (ec 1.1.1.52)
flux_alpha_dehydroxy_BA_all_trp_rxn02890_auxo <- flux_alpha_dehydroxy_BA_all$Compound == "Trp" & flux_alpha_dehydroxy_BA_all$status == "good_blast" &
  flux_alpha_dehydroxy_BA_all$seed == "rxn02890" & flux_alpha_dehydroxy_BA_all$Prototrophy == 0 
flux_alpha_dehydroxy_BA_all_trp_rxn02890_auxo <- flux_alpha_dehydroxy_BA_all[flux_alpha_dehydroxy_BA_all_trp_rxn02890_auxo, ]
nrow(flux_alpha_dehydroxy_BA_all_trp_rxn02890_auxo)
#number prototrophic bacteria  rxn02890   
flux_alpha_dehydroxy_BA_all_trp_rxn02890_proto <- flux_alpha_dehydroxy_BA_all$Compound == "Trp" & flux_alpha_dehydroxy_BA_all$status == "good_blast" &
  flux_alpha_dehydroxy_BA_all$seed == "rxn02890" & flux_alpha_dehydroxy_BA_all$Prototrophy == 1 
flux_alpha_dehydroxy_BA_all_trp_rxn02890_proto <- flux_alpha_dehydroxy_BA_all[flux_alpha_dehydroxy_BA_all_trp_rxn02890_proto, ]
nrow(flux_alpha_dehydroxy_BA_all_trp_rxn02890_proto)
#number trp auxotrophic and non producer
sum_no_rxn02890_auxo <- nrow(numb_auxo) - nrow(flux_alpha_dehydroxy_BA_all_trp_rxn02890_auxo)
sum_no_rxn02890_auxo 
#number trp prototrophic and non producer
sum_no_rxn02890_proto <- 5416 - (nrow(flux_alpha_dehydroxy_BA_all_trp_rxn02890_auxo) + nrow(flux_alpha_dehydroxy_BA_all_trp_rxn02890_proto) + sum_no_rxn02890_auxo)
sum_no_rxn02890_proto 

#create contigency table
fisher_ec1.1.1.52  <- data.frame("Trp-Auxotrophy" = c(nrow(flux_alpha_dehydroxy_BA_all_trp_rxn02890_auxo), sum_no_rxn02890_auxo), 
                                 "Trp-Prototrophy"= c(nrow(flux_alpha_dehydroxy_BA_all_trp_rxn02890_proto), sum_no_rxn02890_proto), 
                                 row.names = c(" rxn02890   : yes", " rxn02890   : no"))

fisher_ec1.1.1.52  
#control numbers
nrow(flux_alpha_dehydroxy_BA_all_trp_rxn02890_auxo) + sum_no_rxn02890_auxo + nrow(flux_alpha_dehydroxy_BA_all_trp_rxn02890_proto) + sum_no_rxn02890_proto == 5416

#create mosaic plot
mosaic5 <- mosaicplot(fisher_ec1.1.1.52, color = TRUE, main=" rxn02890 (ec 1.1.1.52) - trp auxotrophic and prototrophic microbiota")

#fisher test
fisher.test(fisher_ec1.1.1.52)
#pvalue
f10 <- fisher.test(fisher_ec1.1.1.52)
pvalue_ec1.1.1.52 <- f10$p.value

#fisher test/comparison trp auxotrophic and prototrophic
#CAUTION: due to the fact that the number of genomes are the same for every reactions of the enzymes(see to the last ggplot)
#for one reaction of every enzyme the fisher t-test will be done because the results would be the same
#number of auxotrophic bacteria  rxn05066 (ec 4.2.1.106)
flux_alpha_dehydroxy_BA_all_trp_rxn05066_auxo <- flux_alpha_dehydroxy_BA_all$Compound == "Trp" & flux_alpha_dehydroxy_BA_all$status == "good_blast" &
  flux_alpha_dehydroxy_BA_all$seed == "rxn05066" & flux_alpha_dehydroxy_BA_all$Prototrophy == 0 
flux_alpha_dehydroxy_BA_all_trp_rxn05066_auxo <- flux_alpha_dehydroxy_BA_all[flux_alpha_dehydroxy_BA_all_trp_rxn05066_auxo, ]
nrow(flux_alpha_dehydroxy_BA_all_trp_rxn05066_auxo)
#number prototrophic bacteria  rxn05066   
flux_alpha_dehydroxy_BA_all_trp_rxn05066_proto <- flux_alpha_dehydroxy_BA_all$Compound == "Trp" & flux_alpha_dehydroxy_BA_all$status == "good_blast" &
  flux_alpha_dehydroxy_BA_all$seed == "rxn05066" & flux_alpha_dehydroxy_BA_all$Prototrophy == 1 
flux_alpha_dehydroxy_BA_all_trp_rxn05066_proto <- flux_alpha_dehydroxy_BA_all[flux_alpha_dehydroxy_BA_all_trp_rxn05066_proto, ]
nrow(flux_alpha_dehydroxy_BA_all_trp_rxn05066_proto)
#number trp auxotrophic and non producer
sum_no_rxn05066_auxo <- nrow(numb_auxo) - nrow(flux_alpha_dehydroxy_BA_all_trp_rxn05066_auxo)
sum_no_rxn05066_auxo 
#number trp prototrophic and non producer
sum_no_rxn05066_proto <- 5416 - (nrow(flux_alpha_dehydroxy_BA_all_trp_rxn05066_auxo) + nrow(flux_alpha_dehydroxy_BA_all_trp_rxn05066_proto) + sum_no_rxn05066_auxo)
sum_no_rxn05066_proto 

#create contigency table
fisher_ec4.2.1.106   <- data.frame("Trp-Auxotrophy" = c(nrow(flux_alpha_dehydroxy_BA_all_trp_rxn05066_auxo), sum_no_rxn05066_auxo), 
                                 "Trp-Prototrophy"= c(nrow(flux_alpha_dehydroxy_BA_all_trp_rxn05066_proto), sum_no_rxn05066_proto), 
                                 row.names = c(" rxn05066   : yes", " rxn05066   : no"))

fisher_ec4.2.1.106 
#control numbers
nrow(flux_alpha_dehydroxy_BA_all_trp_rxn05066_auxo) + sum_no_rxn05066_auxo + nrow(flux_alpha_dehydroxy_BA_all_trp_rxn05066_proto) + sum_no_rxn05066_proto == 5416

#create mosaic plot
mosaic5 <- mosaicplot(fisher_ec4.2.1.106, color = TRUE, main=" rxn05066 (ec 4.2.1.106) - trp auxotrophic and prototrophic microbiota")

#fisher test
fisher.test(fisher_ec4.2.1.106)

#pvalue
f11 <- fisher.test(fisher_ec4.2.1.106)
pvalue_ec4.2.1.106 <- f11$p.value

#fisher test/comparison trp auxotrophic and prototrophic
#CAUTION: due to the fact that the number of genomes are the same for every reactions of the enzymes(see to the last ggplot)
#for one reaction of every enzyme the fisher t-test will be done because the results would be the same
#number of auxotrophic bacteria  rxn19220 (ec 1.3.1.115)
flux_alpha_dehydroxy_BA_all_trp_rxn19220_auxo <- flux_alpha_dehydroxy_BA_all$Compound == "Trp" & flux_alpha_dehydroxy_BA_all$status == "good_blast" &
  flux_alpha_dehydroxy_BA_all$seed == "rxn19220" & flux_alpha_dehydroxy_BA_all$Prototrophy == 0 
flux_alpha_dehydroxy_BA_all_trp_rxn19220_auxo <- flux_alpha_dehydroxy_BA_all[flux_alpha_dehydroxy_BA_all_trp_rxn19220_auxo, ]
nrow(flux_alpha_dehydroxy_BA_all_trp_rxn19220_auxo)
#number prototrophic bacteria  rxn19220   
flux_alpha_dehydroxy_BA_all_trp_rxn19220_proto <- flux_alpha_dehydroxy_BA_all$Compound == "Trp" & flux_alpha_dehydroxy_BA_all$status == "good_blast" &
  flux_alpha_dehydroxy_BA_all$seed == "rxn19220" & flux_alpha_dehydroxy_BA_all$Prototrophy == 1 
flux_alpha_dehydroxy_BA_all_trp_rxn19220_proto <- flux_alpha_dehydroxy_BA_all[flux_alpha_dehydroxy_BA_all_trp_rxn19220_proto, ]
nrow(flux_alpha_dehydroxy_BA_all_trp_rxn19220_proto)
#number trp auxotrophic and non producer
sum_no_rxn19220_auxo <- nrow(numb_auxo) - nrow(flux_alpha_dehydroxy_BA_all_trp_rxn19220_auxo)
sum_no_rxn19220_auxo 
#number trp prototrophic and non producer
sum_no_rxn19220_proto <- 5416 - (nrow(flux_alpha_dehydroxy_BA_all_trp_rxn19220_auxo) + nrow(flux_alpha_dehydroxy_BA_all_trp_rxn19220_proto) + sum_no_rxn19220_auxo)
sum_no_rxn19220_proto 

#create contigency table
fisher_ec1.3.1.115   <- data.frame("Trp-Auxotrophy" = c(nrow(flux_alpha_dehydroxy_BA_all_trp_rxn19220_auxo), sum_no_rxn19220_auxo), 
                                 "Trp-Prototrophy"= c(nrow(flux_alpha_dehydroxy_BA_all_trp_rxn19220_proto), sum_no_rxn19220_proto), 
                                 row.names = c(" rxn19220   : yes", " rxn19220   : no"))

fisher_ec1.3.1.115  
#control numbers
nrow(flux_alpha_dehydroxy_BA_all_trp_rxn19220_auxo) + sum_no_rxn19220_auxo + nrow(flux_alpha_dehydroxy_BA_all_trp_rxn19220_proto) + sum_no_rxn19220_proto == 5416

#create mosaic plot
mosaic5 <- mosaicplot(fisher_ec1.3.1.115, color = TRUE, main=" rxn19220 (ec 1.3.1.115) - trp auxotrophic and prototrophic microbiota")

#fisher test
fisher.test(fisher_ec1.3.1.115)
#pvalue
f12 <- fisher.test(fisher_ec1.3.1.115)
pvalue_ec1.3.1.115 <- f12$p.value

#fisher test/comparison trp auxotrophic and prototrophic
#CAUTION: due to the fact that the number of genomes are the same for every reactions of the enzymes(see to the last ggplot)
#for one reaction of every enzyme the fisher t-test will be done because the results would be the same
#number of auxotrophic bacteria  rxn36658 (ec 1.1.1.395)
flux_alpha_dehydroxy_BA_all_trp_rxn36658_auxo <- flux_alpha_dehydroxy_BA_all$Compound == "Trp" &flux_alpha_dehydroxy_BA_all$status == "good_blast" &
  flux_alpha_dehydroxy_BA_all$seed == "rxn36658" & flux_alpha_dehydroxy_BA_all$Prototrophy == 0 
flux_alpha_dehydroxy_BA_all_trp_rxn36658_auxo <- flux_alpha_dehydroxy_BA_all[flux_alpha_dehydroxy_BA_all_trp_rxn36658_auxo, ]
nrow(flux_alpha_dehydroxy_BA_all_trp_rxn36658_auxo)

#number prototrophic bacteria  rxn36658   
flux_alpha_dehydroxy_BA_all_trp_rxn36658_proto <- flux_alpha_dehydroxy_BA_all$Compound == "Trp" & flux_alpha_dehydroxy_BA_all$status == "good_blast" &
  flux_alpha_dehydroxy_BA_all$seed == "rxn36658" & flux_alpha_dehydroxy_BA_all$Prototrophy == 1 
flux_alpha_dehydroxy_BA_all_trp_rxn36658_proto <- flux_alpha_dehydroxy_BA_all[flux_alpha_dehydroxy_BA_all_trp_rxn36658_proto, ]
nrow(flux_alpha_dehydroxy_BA_all_trp_rxn36658_proto)
#number trp auxotrophic and non producer
sum_no_rxn36658_auxo <- nrow(numb_auxo) - nrow(flux_alpha_dehydroxy_BA_all_trp_rxn36658_auxo)
sum_no_rxn36658_auxo 
#number trp prototrophic and non producer
sum_no_rxn36658_proto <- 5416 - (nrow(flux_alpha_dehydroxy_BA_all_trp_rxn36658_auxo) + nrow(flux_alpha_dehydroxy_BA_all_trp_rxn36658_proto) + sum_no_rxn36658_auxo)
sum_no_rxn36658_proto 

#create contigency table
fisher_ec1.1.1.395   <- data.frame("Trp-Auxotrophy" = c(nrow(flux_alpha_dehydroxy_BA_all_trp_rxn36658_auxo), sum_no_rxn36658_auxo), 
                                 "Trp-Prototrophy"= c(nrow(flux_alpha_dehydroxy_BA_all_trp_rxn36658_proto), sum_no_rxn36658_proto), 
                                 row.names = c(" rxn36658   : yes", " rxn36658   : no"))
fisher_ec1.1.1.395 
#control numbers
nrow(flux_alpha_dehydroxy_BA_all_trp_rxn36658_auxo) + sum_no_rxn36658_auxo + nrow(flux_alpha_dehydroxy_BA_all_trp_rxn36658_proto) + sum_no_rxn36658_proto == 5416

#create mosaic plot
mosaic5 <- mosaicplot(fisher_ec1.1.1.395, color = TRUE, main=" rxn36658 (ec 1.1.1.395) - trp auxotrophic and prototrophic microbiota")

#fisher test
fisher.test(fisher_ec1.1.1.395)
#pvalue
f13 <- fisher.test(fisher_ec1.1.1.395)
fisher_ec1.1.1.395 <- f13$p.value

#fisher test/comparison trp auxotrophic and prototrophic
#CAUTION: due to the fact that the number of genomes are the same for every reactions of the enzymes(see to the last ggplot)
#for one reaction of every enzyme the fisher t-test will be done because the results would be the same
#number of auxotrophic bacteria  rxn43871 (ec 1.3.1.114)
flux_alpha_dehydroxy_BA_all_trp_rxn43871_auxo <- flux_alpha_dehydroxy_BA_all$Compound == "Trp" &flux_alpha_dehydroxy_BA_all$status == "good_blast" &
  flux_alpha_dehydroxy_BA_all$seed == "rxn43871" & flux_alpha_dehydroxy_BA_all$Prototrophy == 0 
flux_alpha_dehydroxy_BA_all_trp_rxn43871_auxo <- flux_alpha_dehydroxy_BA_all[flux_alpha_dehydroxy_BA_all_trp_rxn43871_auxo, ]
nrow(flux_alpha_dehydroxy_BA_all_trp_rxn43871_auxo)

#number prototrophic bacteria  rxn43871   
flux_alpha_dehydroxy_BA_all_trp_rxn43871_proto <- flux_alpha_dehydroxy_BA_all$Compound == "Trp" & flux_alpha_dehydroxy_BA_all$status == "good_blast" &
  flux_alpha_dehydroxy_BA_all$seed == "rxn43871" & flux_alpha_dehydroxy_BA_all$Prototrophy == 1 
flux_alpha_dehydroxy_BA_all_trp_rxn43871_proto <- flux_alpha_dehydroxy_BA_all[flux_alpha_dehydroxy_BA_all_trp_rxn43871_proto, ]
nrow(flux_alpha_dehydroxy_BA_all_trp_rxn43871_proto)
#number trp auxotrophic and non producer
sum_no_rxn43871_auxo <- nrow(numb_auxo) - nrow(flux_alpha_dehydroxy_BA_all_trp_rxn43871_auxo)
sum_no_rxn43871_auxo 
#number trp prototrophic and non producer
sum_no_rxn43871_proto <- 5416 - (nrow(flux_alpha_dehydroxy_BA_all_trp_rxn43871_auxo) + nrow(flux_alpha_dehydroxy_BA_all_trp_rxn43871_proto) + sum_no_rxn43871_auxo)
sum_no_rxn43871_proto 

#create contigency table
fisher_ec1.3.1.114<- data.frame("Trp-Auxotrophy" = c(nrow(flux_alpha_dehydroxy_BA_all_trp_rxn43871_auxo), sum_no_rxn43871_auxo), 
                                 "Trp-Prototrophy"= c(nrow(flux_alpha_dehydroxy_BA_all_trp_rxn43871_proto), sum_no_rxn43871_proto), 
                                 row.names = c(" rxn43871   : yes", " rxn43871   : no"))

fisher_ec1.3.1.114
#control numbers
nrow(flux_alpha_dehydroxy_BA_all_trp_rxn43871_auxo) + sum_no_rxn43871_auxo + nrow(flux_alpha_dehydroxy_BA_all_trp_rxn43871_proto) + sum_no_rxn43871_proto == 5416

#create mosaic plot
mosaic5 <- mosaicplot(fisher_ec1.3.1.114, color = TRUE, main=" rxn43871 (ec 1.3.1.114) - trp auxotrophic and prototrophic microbiota")

#fisher test
fisher.test(fisher_ec1.3.1.114)
#pvalue
f14 <- fisher.test(fisher_ec1.3.1.114)
fisher_ec1.3.1.114 <- f14$p.value

##completeness pathways
#percentage completeness pathways
#found enzymes
ec1.3.1.114 <- nrow(flux_alpha_dehydroxy_BA_all_trp_rxn43871_auxo)/ nrow(numb_auxo) * 100
ec1.1.1.395 <- nrow(flux_alpha_dehydroxy_BA_all_trp_rxn36658_auxo)/ nrow(numb_auxo) * 100
ec1.3.1.115 <- nrow(flux_alpha_dehydroxy_BA_all_trp_rxn19220_auxo)/ nrow(numb_auxo) * 100
ec1.1.1.52 <-  nrow(flux_alpha_dehydroxy_BA_all_trp_rxn02890_auxo)/ nrow(numb_auxo) * 100
ec4.2.1.106 <- nrow(flux_alpha_dehydroxy_BA_all_trp_rxn05066_auxo)/ nrow(numb_auxo) * 100

#originally the pathway exists of 6 enzymes, only four are shown in the plot so they are also missing but test
flu_alpha_dehydroxy_BA <- flux_i$ec == "6.2.1.7" | flux_i$ec == "1.1.1.395" | flux_i$ec == "1.3.1.115" | 
  flux_i$ec == "2.8.3.25" | flux_i$ec == "4.2.1.106" | flux_i$ec == "1.3.1.114"  | flux_i$ec == "1.1.1.52" 

ec6.2.1.7 <- flux_alpha_dehydroxy_BA_all$Compound == "Trp" & flux_alpha_dehydroxy_BA_all$status == "good_blast" &
  flux_alpha_dehydroxy_BA_all$ec == "6.2.1.7" & flux_alpha_dehydroxy_BA_all$Prototrophy == 0 
ec6.2.1.7<- flux_alpha_dehydroxy_BA_all[ec6.2.1.7, ]
ec6.2.1.7 <- nrow(ec6.2.1.7)
ec2.8.3.25 <- flux_alpha_dehydroxy_BA_all$Compound == "Trp" & flux_alpha_dehydroxy_BA_all$status == "good_blast" &
  flux_alpha_dehydroxy_BA_all$ec == "2.8.3.25" & flux_alpha_dehydroxy_BA_all$Prototrophy == 0 
ec2.8.3.25 <- flux_alpha_dehydroxy_BA_all[ec2.8.3.25, ]
ec2.8.3.25 <- nrow(ec2.8.3.25)

alpha_dehydroxy <- data.frame(ec1.3.1.114,ec1.1.1.395,ec1.3.1.115,ec1.1.1.52,ec4.2.1.106, ec6.2.1.7, ec2.8.3.25)
alpha_dehydroxy <- t(alpha_dehydroxy)
colnames(alpha_dehydroxy) <- "percentage"
alpha_dehydroxy <- data.frame(alpha_dehydroxy)
alpha_dehydroxy$new <- rownames(alpha_dehydroxy)
pvalue <- c(pvalue_ec1.3.1.114,pvalue_ec1.1.1.395,pvalue_ec1.3.1.115,pvalue_ec1.1.1.52,pvalue_ec4.2.1.106, 0.05, 0.05)
alpha_dehydroxy$pvalues <- pvalue

ggplot (alpha_dehydroxy,aes(new, percentage, label = ifelse(pvalues < 0.05, "*","NS"))) +
  geom_bar(stat = "identity") +
  ggtitle(expression(atop("Completeness of the bile acid alpha dehydroxylation pathway in tryptophan auxotrophic microbiota", atop("added results of Fisher exact t-test for comparison to tryptophan prototrophic microbiota")))) +
  xlab("enzymes") +
  ylab("pathway completeness [%]") +
  ylim(0,100) +
  geom_text(vjust = -1) +
  labs(caption = "* - statistically significant to trp prototrophy microbiota\n\n NS - not statistically significant to trp protrotrophy microbiota")

#not found enzymes
ec1.3.1.114_no <- sum_no_rxn43871_auxo/ nrow(numb_auxo) * 100
ec1.1.1.395_no <- sum_no_rxn36658_auxo/ nrow(numb_auxo) * 100
ec1.3.1.115_no <-  sum_no_rxn19220_auxo/ nrow(numb_auxo) * 100
ec1.1.1.52_no <-  sum_no_rxn02890_auxo/ nrow(numb_auxo) * 100
ec4.2.1.106_no <- sum_no_rxn05066_auxo/ nrow(numb_auxo) * 100

#originally the pathway exists of 6 enzymes, only four are shown in the plot so they are also missing but test
flu_alpha_dehydroxy_BA <- flux_i$ec == "6.2.1.7" | flux_i$ec == "1.1.1.395" | flux_i$ec == "1.3.1.115" | 
  flux_i$ec == "2.8.3.25" | flux_i$ec == "4.2.1.106" | flux_i$ec == "1.3.1.114"  | flux_i$ec == "1.1.1.52" 

ec6.2.1.7 <- flux_alpha_dehydroxy_BA_all$Compound == "Trp" & flux_alpha_dehydroxy_BA_all$status == "good_blast" &
  flux_alpha_dehydroxy_BA_all$ec == "6.2.1.7" & flux_alpha_dehydroxy_BA_all$Prototrophy == 0 
ec6.2.1.7<- flux_alpha_dehydroxy_BA_all[ec6.2.1.7, ]
ec6.2.1.7 <- numb(auxo) - nrow(ec6.2.1.7)
ec2.8.3.25 <- flux_alpha_dehydroxy_BA_all$Compound == "Trp" & flux_alpha_dehydroxy_BA_all$status == "good_blast" &
  flux_alpha_dehydroxy_BA_all$ec == "2.8.3.25" & flux_alpha_dehydroxy_BA_all$Prototrophy == 0 
ec2.8.3.25 <- flux_alpha_dehydroxy_BA_all[ec2.8.3.25, ]
ec2.8.3.25 <- nrow(ec2.8.3.25)

alpha_dehydroxy_no <- data.frame(ec1.3.1.114_no,ec1.1.1.395_no,ec1.3.1.115_no,ec1.1.1.52_no,ec4.2.1.106_no, 100, 100)
alpha_dehydroxy_no <- t(alpha_dehydroxy_no)
colnames(alpha_dehydroxy_no) <- "percentage"
alpha_dehydroxy_no <- data.frame(alpha_dehydroxy_no)
alpha_dehydroxy_no$new <- rownames(alpha_dehydroxy_no)

ggplot (alpha_dehydroxy_no,aes(new, percentage)) +
  geom_bar(stat = "identity") +
  ggtitle("Incompleteness of the alpha dehydroxylation pathway in trp auxotrophic microbiota") +
  scale_x_discrete(breaks=c("ec1.1.1.395_no", "ec1.1.1.52_no", "ec1.3.1.114_no","ec1.3.1.115_no", "ec4.2.1.106_no", "X100", "X100.1"), 
                   labels=c("ec1.1.1.395", "ec1.1.1.1.52", "ec1.3.1.114", "ec1.3.1.115", "ec4.2.1.106", "ec6.2.1.7", "ec2.8.3.25")) +
  xlab("enzymes") +
  ylab("pathway incompleteness [%]") +
  ylim(0,100)
  

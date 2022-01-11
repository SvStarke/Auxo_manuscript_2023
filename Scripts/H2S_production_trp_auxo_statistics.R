############    H2S production by trp auxotrophic microbiota   #################
 
#get production of metabolites
model.production <- lapply(models, FUN = get_produced_metabolites)

head(model.production)

modprod_DT <- rbindlist(model.production, idcol = "model")
View(modprod_DT)
fwrite(modprod_DT, file="modprod_DT.csv")

is.data.frame(modprod_DT)

#get growth rates of bacteria
m_gr <- lapply(models, FUN = get_growth)
head(m_gr)
m_growth <- data.table(Genome = names(m_gr),
                       Growth = unlist(m_gr))
fwrite(m_growth, file = "m_growth.csv")
is.list(m_gr)
View(m_growth)

#merge the files
auxo_growth <- merge(Auxotrophy_2, m_growth, by.x = "Genomes",
                     by.y = "Genome", allow.cartesian=TRUE)
fwrite(auxo_growth, file = "auxo_growth.csv")
View(auxo_growth)

auxo_Prod1 <- merge(auxo_growth, modprod_DT, by.x = "Genomes",
                    by.y = "model", allow.cartesian=TRUE)
fwrite(auxo_Prod1, file= "all_Genomes_Production.csv")
auxo_Prod1 <- read.csv("/Users/svenjabusche/Desktop/all_Genomes_Production.csv")
is.data.frame(auxo_Prod1)
is.data.table(auxo_Prod1)
View(auxo_Prod1)

#delete archaes
auxo_Prod1<- auxo_Prod1[!(auxo_Prod1$phylum == "Euryarchaeota" | auxo_Prod1$phylum == "Thermoplasmatota"), ]
View(auxo_Prod1)


#put all Firmicutes phyla to one Firmicutes phylum
auxo_Prod1$phylum[auxo_Prod1$phylum == "Firmicutes_A"] <- "Firmicutes"
auxo_Prod1$phylum[auxo_Prod1$phylum == "Firmicutes_B"] <- "Firmicutes"
auxo_Prod1$phylum[auxo_Prod1$phylum== "Firmicutes_G"] <- "Firmicutes"
auxo_Prod1$phylum[auxo_Prod1$phylum == "Firmicutes_I"] <- "Firmicutes"

#H2S production for all auxotrophic microbiota
H2S_prod <- auxo_Prod1[grepl("H2S-e0 Exchange", rxn.name) & mtf.flux > 0 & Prototrophy == 0][order(-mtf.flux)]
View(prod_auxo_sulfate)

prod_auxo_sulfate<- auxo_Prod1$rxn.name == "H2S-e0 Exchange" & mtf.flux > 0 & Prototrophy == 0

prod_auxo_sulfate$Production <- prod_auxo_sulfate$mtf.flux/prod_auxo_sulfate$Growth

ggplot(data=prod_auxo_sulfate, aes(x = Compound, y = Production, fill = phylum)) +
  geom_bar(stat="identity") +
  ggtitle("Hydrogen-Sulfide Production by auxotrophic microbiota") +
  xlab("Amino Acid Auxotrophy") +
  ylab("Production rate [mmol/gDW]") +
  scale_fill_manual(values = caPalette) +
  scale_color_manual(values = caPalette)

#H2S production all tryptophan auxotrophic
H2S_prod_auxo <- prod_auxo_sulfate[grepl("H2S-e0 Exchange", rxn.name) & mtf.flux > 0 & Prototrophy == 0 & Compound == "Trp"][order(-mtf.flux)]

ggplot(data=prod_auxotrp_sulfate, aes(x = phylum, y = Production, fill = phylum)) +
  geom_bar(stat="identity") +
  ggtitle("Hydrogen-Sulfide Production by tryptophan auxotrophic microbiota") +
  xlab("Phylum") +
  ylab("Production rate [mmol/gDW]") +
  scale_fill_manual(values = caPalette) +
  scale_color_manual(values = caPalette) +
  theme(axis.text.x = element_text (angle=90, size=10)) +
  theme(legend.position = "none")


ggplot(data=prod_auxotrp_sulfate, aes(x = Production, y = phylum, fill = phylum, colour = phylum)) +
  geom_point() +
  ggtitle("Hydrogen-Sulfide Production by tryptophan auxotrophic microbiota") +
  ylab("Phylum") +
  xlab("Production rate [mmol/gDW]") +
  scale_fill_manual(values = caPalette) +
  scale_color_manual(values = caPalette) +
  theme(axis.text.x = element_text (angle=90, size=10))+
  theme(legend.position = "none")

#number of genomes
ggplot(data=prod_auxotrp_sulfate, aes(rxn.name, fill = phylum, colour = phylum)) +
  geom_bar() +
  ggtitle("Number of genomes with Hydrogen-Sulfide Production by\ntryptophan auxotrophic microbiota") +
  xlab("H2S") +
  ylab("number of genomes") +
  scale_fill_manual(values = caPalette) +
  scale_color_manual(values = caPalette) +
  theme(axis.text.x = element_blank()) +
  theme(legend.position ="right")

#Wilcoxon Test
#get the H2S production of tryptophan prototrophic microbiota
H2S_prod_proto <- auxo_Prod1[grepl("H2S-e0 Exchange", rxn.name) & mtf.flux > 0 & Prototrophy == 1 & Compound == "Trp"][order(-mtf.flux)]
View(H2S_prod_proto)

#H2S production prototrophic and auxotrophic
H2S_prod <- auxo_Prod1$rxn.name == "H2S-e0 Exchange" & auxo_Prod1$mtf.flux > 0 &  auxo_Prod1$Compound == "Trp"
H2S_prod <- auxo_Prod1[H2S_prod,]
View(prod_proto_sulfate)

H2S_prod$Production <- H2S_prod$mtf.flux/H2S_prod$Growth

#Test and Visualization
install.packages("ggstatsplot")
library(ggstatsplot)

# plot with statistical results --> unabnhängige Stichprobe (n ist unterschiedlich und auxotroph/prototrophe Stichprobe)
p1 <- ggbetweenstats( # independent samples
  data = H2S_prod,
  x = Prototrophy,
  y = Production,
  plot.type = "violin", # for boxplot
  type = "nonparametric", # for wilcoxon
  centrality.plotting = TRUE,
  pairwise.display ="significant",
  title = "All phyla",
  ylab = "Productionrates [mmol/gDW]",
  xlab = "",
  results.subtitle = TRUE,
  conf.level = 0.95,
  caption = "p > 0.05",
  point.args = list(
    position = ggplot2::position_jitterdodge(dodge.width = 0.20),
    alpha = 0.2,
    size = 1,
    stroke = 0)
)

p1.1 <- p1 + theme(text = element_text(size= 10)) +
  theme(plot.subtitle = element_text(size=6)) +
  theme(plot.title = element_text(size=10))

#single plot
ggbetweenstats( # independent samples
  data = H2S_prod,
  x = Prototrophy,
  y = Production,
  plot.type = "violin", # for boxplot
  type = "nonparametric", # for wilcoxon
  centrality.plotting = TRUE,
  pairwise.display ="significant",
  title = "Comparison of H2S production between tryptophan auxotrophic\nand prototrophic microbiota ",
  ylab = "Productionrates [mmol/gDW]",
  xlab = "",
  results.subtitle = TRUE,
  conf.level = 0.95,
  caption = "Legend: 0 - Auxotrophy, 1 - Prototrophy, Test: Mann-Whitney U Test (CI = 0.95%) p > 0.05",
  point.args = list(
    position = ggplot2::position_jitterdodge(dodge.width = 0.20),
    alpha = 0.2,
    size = 1,
    stroke = 0)
)
#filter Bacteroidetes
H2S_prod_Bact <- auxo_Prod1$rxn.name == "H2S-e0 Exchange" & auxo_Prod1$mtf.flux > 0 &  auxo_Prod1$Compound == "Trp" & auxo_Prod1$phylum == "Bacteroidota"
H2S_prod_Bact <- auxo_Prod1[H2S_prod_Bact, ]

H2S_prod_Bact$Production <- H2S_prod_Bact$mtf.flux/H2S_prod_Bact$Growth

#plot
p2 <- ggbetweenstats( # independent samples
  data = H2S_prod_Bact,
  x = Prototrophy,
  y = Production,
  plot.type = "violin", # for boxplot
  type = "nonparametric", # for wilcoxon
  centrality.plotting = TRUE,
  pairwise.display ="significant",
  title = "Bacteroidota",
  ylab = "Productionrates [mmol/gDW]",
  xlab = "",
  results.subtitle = TRUE,
  conf.level = 0.95,
  caption = "",
  point.args = list(
    position = ggplot2::position_jitterdodge(dodge.width = 0.20),
    alpha = 0.2,
    size = 1,
    stroke = 0)
)

p2
p2.1 <- p2 + theme(text = element_text(size= 10)) +
  theme(plot.subtitle = element_text(size=6)) +
  theme(plot.title = element_text(size=10))


ggbetweenstats( # independent samples
  data = H2S_prod_Bact,
  x = Prototrophy,
  y = Production,
  plot.type = "violin", # for boxplot
  type = "nonparametric", # for wilcoxon
  centrality.plotting = TRUE,
  pairwise.display ="significant",
  title = "Comparison of H2S production between tryptophan auxotrophic\nand prototrophic Bacteroidetes ",
  ylab = "Productionrates [mmol/gDW]",
  xlab = "",
  results.subtitle = TRUE,
  conf.level = 0.95,
  caption = "Legend: 0 - Auxotrophy, 1 - Prototrophy, Test: Mann-Whitney U Test (CI = 0.95%)",
  point.args = list(
    position = ggplot2::position_jitterdodge(dodge.width = 0.20),
    alpha = 0.2,
    size = 1,
    stroke = 0)
)

#filter Firmicutes
H2S_prod_Firm <- auxo_Prod1$rxn.name == "H2S-e0 Exchange" & auxo_Prod1$mtf.flux > 0 &  auxo_Prod1$Compound == "Trp" & auxo_Prod1$phylum == "Firmicutes"
H2S_prod_Firm  <- auxo_Prod1[H2S_prod_Firm , ]

H2S_prod_Firm $Production <- H2S_prod_Firm $mtf.flux/H2S_prod_Firm $Growth

p3 <- ggbetweenstats( # independent samples
  data = H2S_prod_Firm ,
  x = Prototrophy,
  y = Production,
  plot.type = "violin", # for boxplot
  type = "nonparametric", # for wilcoxon
  centrality.plotting = TRUE,
  pairwise.display ="significant",
  title = "Firmicutes",
  ylab = "",
  xlab = "",
  results.subtitle = TRUE,
  conf.level = 0.95,
  caption = "",
  point.args = list(
    position = ggplot2::position_jitterdodge(dodge.width = 0.20),
    alpha = 0.2,
    size = 1,
    stroke = 0)
)
p3
p3.1 <- p3 + theme(text = element_text(size= 10)) +
  theme(plot.subtitle = element_text(size=6)) +
  theme(plot.title = element_text(size=10))

ggbetweenstats( 
  data = H2S_prod_Firm ,
  x = Prototrophy,
  y = Production,
  plot.type = "violin", 
  type = "nonparametric", 
  centrality.plotting = TRUE,
  pairwise.display ="significant",
  title = "Comparison of H2S production between tryptophan auxotrophic\nand prototrophic Firmicutes",
  ylab = "Productionrates [mmol/gDW]",
  xlab = "",
  results.subtitle = TRUE,
  conf.level = 0.95,
  caption = "Legend: 0 - Auxotrophy, 1 - Prototrophy, Test: Mann-Whitney U Test (CI = 0.95%)",
  point.args = list(
    position = ggplot2::position_jitterdodge(dodge.width = 0.20),
    alpha = 0.2,
    size = 1,
    stroke = 0)
)

#filter Proteobacteria
H2S_prod_Proteo  <- auxo_Prod1$rxn.name == "H2S-e0 Exchange" & auxo_Prod1$mtf.flux > 0 &  auxo_Prod1$Compound == "Trp" & auxo_Prod1$phylum == "Proteobacteria"
H2S_prod_Proteo <- auxo_Prod1[H2S_prod_Proteo , ]

H2S_prod_Proteo $Production <- H2S_prod_Proteo $mtf.flux/H2S_prod_Proteo $Growth

p4 <- ggbetweenstats( # independent samples
  data = H2S_prod_Proteo,
  x = Prototrophy,
  y = Production,
  plot.type = "violin", 
  type = "nonparametric", # for wilcoxon
  centrality.plotting = TRUE,
  pairwise.display ="significant",
  title = "Proteobacteria",
  ylab = "Productionrates [mmol/gDW]",
  xlab = "",
  results.subtitle = TRUE,
  conf.level = 0.95,
  caption = "",
  point.args = list(
    position = ggplot2::position_jitterdodge(dodge.width = 0.20),
    alpha = 0.2,
    size = 1,
    stroke = 0),
)
p4
p4.1 <- p4 + theme(text = element_text(size= 10)) +
  theme(plot.subtitle = element_text(size=6)) +
  theme(plot.title = element_text(size=10))

ggbetweenstats( 
  data = H2S_prod_Proteo ,
  x = Prototrophy,
  y = Production,
  plot.type = "violin", 
  type = "nonparametric", 
  centrality.plotting = TRUE,
  pairwise.display ="significant",
  title = "Comparison of H2S production between tryptophan auxotrophic\nand prototrophic Proteobacteria",
  ylab = "Productionrates [mmol/gDW]",
  xlab = "",
  results.subtitle = TRUE,
  conf.level = 0.95,
  caption = "Legend: 0 - Auxotrophy, 1 - Prototrophy, Test: Mann-Whitney U Test (CI = 0.95%)",
  point.args = list(
    position = ggplot2::position_jitterdodge(dodge.width = 0.20),
    alpha = 0.2,
    size = 1,
    stroke = 0)
)

#filter Actinobacteria
H2S_prod_Actino <- auxo_Prod1$rxn.name == "H2S-e0 Exchange" & auxo_Prod1$mtf.flux > 0 &  auxo_Prod1$Compound == "Trp" & auxo_Prod1$phylum == "Actinobacteriota"
H2S_prod_Actino<- auxo_Prod1[H2S_prod_Actino, ]

H2S_prod_Actino$Production <- H2S_prod_Actino$mtf.flux/H2S_prod_Actino$Growth

p5 <- ggbetweenstats( # independent samples
  data = H2S_prod_Actino,
  x = Prototrophy,
  y = Production,
  plot.type = "violin", # for boxplot
  type = "nonparametric", # for wilcoxon
  centrality.plotting = TRUE,
  pairwise.display ="significant",
  title = "Actinobacteria",
  ylab = "",
  xlab = "",
  results.subtitle = TRUE,
  conf.level = 0.95,
  caption = "",
  point.args = list(
    position = ggplot2::position_jitterdodge(dodge.width = 0.20),
    alpha = 0.2,
    size = 1,
    stroke = 0),
)
p5
p5.1 <- p5 + theme(text = element_text(size= 10)) +
  theme(plot.subtitle = element_text(size=6)) +
  theme(plot.title = element_text(size=10))

ggbetweenstats( 
  data = H2S_prod_Actino,
  x = Prototrophy,
  y = Production,
  plot.type = "violin", 
  type = "nonparametric", 
  centrality.plotting = TRUE,
  pairwise.display ="significant",
  title = "Comparison of H2S production between tryptophan auxotrophic\nand prototrophic Actinobacteria",
  ylab = "Productionrates [mmol/gDW]",
  xlab = "",
  results.subtitle = TRUE,
  conf.level = 0.95,
  caption = "Legend: 0 - Auxotrophy, 1 - Prototrophy, Test: Mann-Whitney U Test (CI = 0.95%)",
  point.args = list(
    position = ggplot2::position_jitterdodge(dodge.width = 0.20),
    alpha = 0.2,
    size = 1,
    stroke = 0)
) 

#combine all plots
combine_plots(list(p2.1, p3.1, p4.1, p5.1),
              plotgrid.args = list(nrow=2),
              annotation.args = list(title = "Comparison of H2S productionrates by tryptophan auxotrophic and prototrophic microbiota",
                                     caption = "Legend: 0 - Auxotrophy, 1 - Prototrophy, Test: Mann-Whitney U Test (CI = 0.95%)",
                                     theme = ggplot2::theme(
                                       plot.subtitle = element_text(size = 8),
                                       axis.title.y = element_text(face = "plain"),
                                       text = element_text(size=8, face = "plain"),
                                       plot.title = element_text(size = 11, face = "bold"))))


#fisher test
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
#Propionate production
nrow(H2S_prod_auxo)
nrow(H2S_prod_proto)
#number trp auxotrophic and non Butyrate producer
sum_no_H2S_auxo <- nrow(numb_auxo) - nrow(H2S_prod_auxo)
sum_no_H2S_auxo
#number trp prototrophic and non butyrate producer
sum_no_H2S_proto <- 5416 -(nrow(H2S_prod_auxo) + nrow(H2S_prod_proto) + sum_no_H2S_auxo)
sum_no_H2S_proto 

#create contigency table
fisher_P <- data.frame("Trp-Auxotrophy" = c(nrow(H2S_prod_auxo), sum_no_H2S_auxo), 
                       "Trp-Prototrophy"= c(nrow(H2S_prod_proto),sum_no_H2S_proto), 
                       row.names = c("H2S production: yes", "H2S production: no"))

fisher_P
#control numbers
nrow(H2S_prod_auxo) + nrow(H2S_prod_proto) + sum_no_H2S_auxo + sum_no_H2S_proto == 5416

#create mosaic plot
mosaicplot(fisher_P, color = TRUE, main="H2S production of trp auxotrophic and prototrophic microbiota")
#fisher test
fisher.test(fisher_P)

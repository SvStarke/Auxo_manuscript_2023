#SCFA production by tryptophan auxotrophic microbiota
#statistical analysis(Wilcox, Fisher)

model.production <- lapply(models, FUN = get_produced_metabolites)

head(model.production)


modprod_DT <- rbindlist(model.production, idcol = "model")
View(modprod_DT)
fwrite(modprod_DT, file="modprod_DT.csv")

is.data.frame(modprod_DT)

m_gr <- lapply(models, FUN = get_growth)
head(m_gr)
m_growth <- data.table(Genome = names(m_gr),
                       Growth = unlist(m_gr))
fwrite(m_growth, file = "m_growth.csv")
is.list(m_gr)
View(m_growth)

# Anwendung (Hier Auxotrophies vorhersagen)
model.auxo <- lapply(models, FUN = predict_auxotrohies)

head(model.auxo)

Metadata <- fread("/mnt/nuuk/2021/HRGM/REPR_Genomes_metadata.tsv")

#create the first data.frame
Auxotrophie <- data.frame(model.auxo)
head(Auxotrophie) 
summary(Auxotrophie)
str(Auxotrophie)
is.data.frame(Auxotrophie)

#column und rows tauschen
Auxotroph <- t(Auxotrophie)

is.matrix(Auxotroph)
#data frame erzeugen
Auxotrophy <- data.frame(Auxotroph)
is.data.frame(Auxotrophy)
str(Auxotrophy)
Auxotrophy
Genome <- rownames(Auxotrophy)
Auxotrophy$Genomes <- Genome
# ----
Auxotrophy_2 <- as.data.table(Auxotrophy)
Auxotrophy_2 <- melt(Auxotrophy_2, id.vars = "Genomes",
                     value.name = "Prototrophy", variable.name = "Compound")
fwrite(Auxotrophy_2, file = "Info_allGenomes_Auxotrophy_Protrophy_Metadata.csv")
Auxotrophy_2 <- read.csv("/Users/svenjabusche/Desktop/Info_allGenomes_Auxotrophy_Protrophy_Metadata.csv")
head(Auxotrophy_2)
Metadata <- fread("/mnt/nuuk/2021/HRGM/REPR_Genomes_metadata.tsv")

Auxotrophy_2 <- merge(Auxotrophy_2, Metadata, by.x = "Genomes",
                      by.y = "HRGM name")
Auxotrophy_2[, phylum := str_match(`GTDB Taxonomy`, "p__.*;c__")[,1]]
Auxotrophy_2[, phylum := gsub("p__|;c__","", phylum)]
Auxotrophy_2[, phylum := gsub("_C$","", phylum)]
#merge the files
auxo_growth <- merge(Auxotrophy_2, m_growth, by.x = "Genomes",
                     by.y = "Genome", allow.cartesian=TRUE)
fwrite(auxo_growth, file = "auxo_growth.csv")
View(auxo_growth)

auxo_Prod1 <- merge(auxo_growth, modprod_DT, by.x = "Genomes",
                    by.y = "model", allow.cartesian=TRUE)
fwrite(auxo_Prod1, file= "all_Genomes_Production.csv")
auxo_Prod1 <- read.csv("/Users/svenjabusche/Desktop/all_Genomes_Production.csv")
View(auxo_Prod1)

#delete archaes
auxo_Prod1<- auxo_Prod1[!(auxo_Prod1$phylum == "Euryarchaeota" | auxo_Prod1$phylum == "Thermoplasmatota"), ]
View(auxo_Prod1)


#put all Firmicutes phyla to one Firmicutes phylum
auxo_Prod1$phylum[auxo_Prod1$phylum == "Firmicutes_A"] <- "Firmicutes"
auxo_Prod1$phylum[auxo_Prod1$phylum == "Firmicutes_B"] <- "Firmicutes"
auxo_Prod1$phylum[auxo_Prod1$phylum== "Firmicutes_G"] <- "Firmicutes"
auxo_Prod1$phylum[auxo_Prod1$phylum == "Firmicutes_I"] <- "Firmicutes"

#Ab der Stelle Filtern nach tryptophan auxotropher Microbiota
# bsp: Filter nach Butyrate und sortiert
B_trp_auxo <- auxo_Prod1[grepl("Butyrate", rxn.name) & mtf.flux > 0 & Prototrophy == 0 & Compound == "Trp"][order(-mtf.flux)]
P_trp_auxo <- auxo_Prod1[grepl("Propionate", rxn.name) & mtf.flux > 0 & Prototrophy == 0 & Compound == "Trp"][order(-mtf.flux)]
A_trp_auxo <- auxo_Prod1[grepl("Acetate", rxn.name) & mtf.flux > 0 & Prototrophy == 0 & Compound == "Trp"][order(-mtf.flux)]

View(B_trp2)
B_trp
P_trp
A_trp

#create one dataframe with the production of all SCFA
SCFA_trp_auxo <- bind_rows(B_trp, P_trp, A_trp)
is.data.frame(SCFA_trp_auxo)

#Production of organic acids
#Involvement of auxotrophic bacteria in the production of organic acids

Production <- SCFA_trp_auxo$mtf.flux / SCFA_trp_auxo$Growth
SCFA_trp_auxo$prod <- SCFA_trp_auxo$mtf.flux/SCFA_trp_auxo$Growth
View(SCFA_trp_auxo)
fwrite(SCFA_trp_auxo, file = "SCFA_trp_auxo.csv")

#Visualization barplot
ggplot(data= SCFA_trp_auxo, aes(x= rxn.name, y= mtf.flux, fill= phylum)) +
  geom_bar(stat= "identity")

#colour for disabled people
caPalette <- c("#560133", "#C7007C", "#DA00FD", "#7CFFFA", "#005745", "#00306F", 
               "#00C2F9", "#004002", "#00B408", "#F60239", "#CD022D","#E69F00", 
               "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7",
               "#450270", "#FFCCFE", "#5A000F", "#FF5AAF", "#003D30", "#FFAC3B",
               "#00735C")

#visualization
ggplot(data= SCFA_trp_auxo, aes(x=rxn.name, y= prod, fill = phylum, col = phylum)) +
  geom_bar(stat= "identity") +
  theme(panel.background = element_rect(fill="white", colour= "black")) +
  scale_fill_manual(values = caPalette) +
  scale_color_manual(values = caPalette) +
  xlab ("SCFA") +
  ylab("Production rate [mmol/gDW]")  + 
  scale_x_discrete(breaks=c("Acetate-e0 Exchange", "Butyrate-e0 Exchange", "Propionate-e0 Exchange"), 
                   labels=c("Acetate", "Butyrate", "Propionate")) +
  ggtitle("SCFA Production of tryptophan auxotrophic microbiota from HRGM") 

#visualization boxplot
ggplot(data=SCFA_trp_auxo, aes(x=rxn.name, y=prod)) +
  geom_boxplot()

#visualization number of genomes
View(SCFA_trp_auxo)
ggplot(data= SCFA_trp_auxo, aes(rxn.name, fill = phylum)) +
  geom_bar () +
  xlab("SCFA") +
  ylab("number of genomes")+
  ggtitle("Number of genomes from tryptophan auxotrophic bacteria producing SCFA") +
  theme(panel.background = element_rect(fill="white", colour= "black")) +
  scale_fill_manual(values = caPalette) +
  scale_color_manual(values = caPalette) +
  scale_x_discrete(breaks=c("Acetate-e0 Exchange", "Butyrate-e0 Exchange", "Propionate-e0 Exchange"), 
                   labels=c("Acetate", "Butyrate", "Propionate"))

#Wilcoxtest (statistical analysis of trp auxotrophic and prototrophic microbiota)
#first get productionrates of trp auxotrophic microbiota in the already existing dataframe
B_trp$prod <- B_trp$mtf.flux/B_trp$Growth

#tryptophan prototrophic
B_trp_Pr <-auxo_Prod1$rxn.name == "Butyrate-e0 Exchange" & auxo_Prod1$mtf.flux > 0 & auxo_Prod1$Prototrophy == 1 & auxo_Prod1$Compound == "Trp"
B_trp_Pr <-auxo_Prod1[B_trp_Pr, ]
View(B_trp_Proto)

#new column with productionrates
B_trp_Pr$prod <- B_trp_Pr$mtf.flux/B_trp_Pr$Growth

ggplot(data=B_trp_Proto, aes(rxn.name, prod, fill = phylum))+
  geom_bar(stat = "identity") +
  ggtitle("Butyrate production of tryptophan prototrophic microbiota")

#wilcox test
test <- wilcox.test(B_trp$prod, B_trp_Pr$prod,
                    paired = FALSE, conf.level = 0.95)
test

#Protrophy and Auxotorphy data in one dataset for a visualization of the statistical results
B_trp_auxo_proto <-auxo_Prod1$rxn.name == "Butyrate-e0 Exchange" & auxo_Prod1$mtf.flux > 0 & auxo_Prod1$Compound == "Trp"
B_trp_auxo_proto <-auxo_Prod1[B_trp_auxo_proto, ]

#new column with productionrates
B_trp_auxo_proto$prod <- B_trp_auxo_proto$mtf.flux/B_trp_auxo_proto$Growth
B_trp_auxo_proto <-B_trp_auxo_proto[order(B_trp_auxo_proto$prod),]


# load package
install.packages("ggstatsplot")
library(ggstatsplot)

# plot with statistical results --> independent 
ggbetweenstats( # independent samples
  data = B_trp_auxo_proto,
  x = Prototrophy,
  y = prod,
  plot.type = "violin",
  type = "nonparametric", # for wilcoxon
  centrality.plotting = TRUE,
  pairwise.display ="significant",
  title = "Comparison of the butyrate productionrates by tryptophan auxotrophic \nand prototrophic microbiota",
  ylab = "Productionrates [mmol/gDW]",
  xlab = "",
  results.subtitle = TRUE,
  conf.level = 0.95,
  caption = "0 - Auxotrophy, 1 - Prototrophy\nMann-Whitney U Test, p < 0.05, CI = 0.95%",
  point.args = list(
    position = ggplot2::position_jitterdodge(dodge.width = 0.20),
    alpha = 0.2,
    size = 1,
    stroke = 0)
)

p_B1 <- ggbetweenstats( # independent samples
  data = B_trp_auxo_proto,
  x = Prototrophy,
  y = prod,
  plot.type = "violin",
  type = "nonparametric", # for wilcoxon
  centrality.plotting = TRUE,
  pairwise.display ="significant",
  title = "All phyla",
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
)  + theme(text = element_text(size= 10)) +
  theme(plot.subtitle = element_text(size=6)) +
  theme(plot.title = element_text(size=10))

#filter phyla
#Bacteroidetes
B_trp_Bact <-auxo_Prod1$rxn.name == "Butyrate-e0 Exchange" & auxo_Prod1$mtf.flux > 0 & auxo_Prod1$Compound == "Trp" & auxo_Prod1$phylum == "Bacteroidota"
B_trp_Bact <-auxo_Prod1[B_trp_Bact, ]
View(B_trp_Bact)

#new column with productionrates
B_trp_Bact$prod <- B_trp_Bact$mtf.flux/B_trp_Bact$Growth
B_trp_Bact <-B_trp_Bact[order(B_trp_Bact$prod),]

#visualization
ggbetweenstats(
  data = B_trp_Bact,
  x = Prototrophy,
  y = prod,
  plot.type = "violin", # for boxplot
  type = "nonparametric", # for wilcoxon
  centrality.plotting = TRUE,
  pairwise.display ="significant",
  title = "Comparison of the butyrate productionrates by tryptophan auxotrophic \nand prototrophic Bacteroidetes",
  ylab = "Productionrates [mmol/gDW]",
  xlab = "",
  results.subtitle = TRUE,
  conf.level = 0.95,
  caption = "0 - Auxotrophy, 1 - Prototrophy\nMann-Whitney U Test (CI = 0.95%)",
  point.args = list(
    position = ggplot2::position_jitterdodge(dodge.width = 0.20),
    alpha = 0.2,
    size = 1,
    stroke = 0)
)

p_B2 <- ggbetweenstats(
  data = B_trp_Bact,
  x = Prototrophy,
  y = prod,
  plot.type = "violin", # for boxplot
  type = "nonparametric", # for wilcoxon
  centrality.plotting = TRUE,
  pairwise.display ="significant",
  title = "Bacteroidetes",
  ylab = "Productionrates [mmol/g DW]",
  xlab = "",
  results.subtitle = TRUE,
  conf.level = 0.95,
  caption = "",
  point.args = list(
    position = ggplot2::position_jitterdodge(dodge.width = 0.20),
    alpha = 0.2,
    size = 1,
    stroke = 0)
)  + theme(text = element_text(size= 10)) +
  theme(plot.subtitle = element_text(size=6)) +
  theme(plot.title = element_text(size=10))


#Firmicutes#
B_trp_Firm <-auxo_Prod1$rxn.name == "Butyrate-e0 Exchange" & auxo_Prod1$mtf.flux > 0 & auxo_Prod1$Compound == "Trp" & auxo_Prod1$phylum == "Firmicutes"
B_trp_Firm <-auxo_Prod1[B_trp_Firm, ]
View(B_trp_Firm)

#new column with productionrates
B_trp_Firm$prod <- B_trp_Firm$mtf.flux/B_trp_Firm$Growth
B_trp_Firm <-B_trp_Firm[order(B_trp_Firm$prod),]

#visualization
ggbetweenstats(
  data = B_trp_Firm,
  x = Prototrophy,
  y = prod,
  plot.type = "violin", 
  type = "nonparametric", # for wilcoxon
  centrality.plotting = TRUE,
  pairwise.display ="significant",
  title = "Comparison of the butyrate productionrates by tryptophan auxotrophic \nand prototrophic Firmicutes",
  ylab = "Productionrates [mmol/gDW]",
  xlab = "",
  results.subtitle = TRUE,
  conf.level = 0.95,
  caption = "0 - Auxotrophy, 1 - Prototrophy\nMann-Whitney U Test (CI = 0.95%)",
  point.args = list(
    position = ggplot2::position_jitterdodge(dodge.width = 0.20),
    alpha = 0.2,
    size = 1,
    stroke = 0)
)
p_B3 <- ggbetweenstats(
  data = B_trp_Firm,
  x = Prototrophy,
  y = prod,
  plot.type = "violin",
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
)  + theme(text = element_text(size= 10)) +
  theme(plot.subtitle = element_text(size=5.5)) +
  theme(plot.title = element_text(size=10))

#Proteobacteria
B_trp_Proteo <-auxo_Prod1$rxn.name == "Butyrate-e0 Exchange" & auxo_Prod1$mtf.flux > 0 & auxo_Prod1$Compound == "Trp" & auxo_Prod1$phylum == "Proteobacteria"
B_trp_Proteo <-auxo_Prod1[B_trp_Proteo, ]
View(B_trp_Proteo)

#new column with productionrates
B_trp_Proteo$prod <- B_trp_Proteo$mtf.flux/B_trp_Proteo$Growth
B_trp_Proteo <-B_trp_Proteo[order(B_trp_Proteo$prod),]

#visualization
ggbetweenstats(
  data = B_trp_Proteo,
  x = Prototrophy,
  y = prod,
  plot.type = "violin", 
  type = "nonparametric", # for wilcoxon
  centrality.plotting = TRUE,
  pairwise.display ="significant",
  title = "Comparison of the butyrate productionrates by tryptophan auxotrophic \nand prototrophic Proteobacteria",
  ylab = "Productionrates [mmol/gDW]",
  xlab = "",
  results.subtitle = TRUE,
  conf.level = 0.95,
  caption = "0 - Auxotrophy, 1 - Prototrophy\nMann-Whitney U Test (CI = 0.95%)",
  point.args = list(
    position = ggplot2::position_jitterdodge(dodge.width = 0.20),
    alpha = 0.2,
    size = 1,
    stroke = 0)
)

p_B4 <- ggbetweenstats(
  data = B_trp_Proteo,
  x = Prototrophy,
  y = prod,
  plot.type = "violin", 
  type = "nonparametric", # for wilcoxon
  centrality.plotting = TRUE,
  pairwise.display ="significant",
  title = "Proteobacteria",
  ylab = "Productionrates [mmol/g DW]",
  xlab = "",
  results.subtitle = TRUE,
  conf.level = 0.95,
  caption = "",
  point.args = list(
    position = ggplot2::position_jitterdodge(dodge.width = 0.20),
    alpha = 0.2,
    size = 1,
    stroke = 0)
)  + theme(text = element_text(size= 10)) +
  theme(plot.subtitle = element_text(size=6)) +
  theme(plot.title = element_text(size=10))

#Actinobacteria
B_trp_Actino <-auxo_Prod1$rxn.name == "Butyrate-e0 Exchange" & auxo_Prod1$mtf.flux > 0 & auxo_Prod1$Compound == "Trp" & auxo_Prod1$phylum == "Actinobacteriota"
B_trp_Actino <-auxo_Prod1[B_trp_Actino, ]
View(B_trp_Actino)

#new column with productionrates
B_trp_Actino$prod <- B_trp_Actino$mtf.flux/B_trp_Actino$Growth
B_trp_Actino <-B_trp_Actino[order(B_trp_Actino$prod),]

#visualization
ggbetweenstats(
  data = B_trp_Actino,
  x = Prototrophy,
  y = prod,
  plot.type = "violin",
  type = "nonparametric", # for wilcoxon
  centrality.plotting = TRUE,
  pairwise.display ="significant",
  title = "Comparison of the butyrate productionrates by tryptophan auxotrophic \nand prototrophic Actinobacteria",
  ylab = "Productionrates [mmol/gDW]",
  xlab = "",
  results.subtitle = TRUE,
  conf.level = 0.95,
  caption = "0 - Auxotrophy, 1 - Prototrophy\nMann-Whitney U Test (CI = 0.95%)",
  point.args = list(
    position = ggplot2::position_jitterdodge(dodge.width = 0.20),
    alpha = 0.2,
    size = 1,
    stroke = 0)
)
p_B5 <- ggbetweenstats(
  data = B_trp_Actino,
  x = Prototrophy,
  y = prod,
  plot.type = "violin", 
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
    stroke = 0)
)  + theme(text = element_text(size= 10)) +
  theme(plot.subtitle = element_text(size=6)) +
  theme(plot.title = element_text(size=10))

#combine all plots
combine_plots(list(p_B2, p_B3, p_B4, p_B5),
              plotgrid.args = list(nrow=2),
              annotation.args = list(title = "Comparison of butyrate productionrates by tryptophan auxotrophic and\nprototrophic microbiota",
                                     caption = "Legend: 0 - Auxotrophy, 1 - Prototrophy, Test: Mann-Whitney U Test (CI = 0.95%)",
                                     theme = ggplot2::theme(
                                       plot.subtitle = element_text(size = 6),
                                       axis.title.y = element_text(face = "plain"),
                                       text = element_text(size=8, face = "plain"),
                                       plot.title = element_text(size = 12, face = "bold"))))

#fisher test
#create a dataframe for Butyrate production
#number of genomes who are tryptophan auxotrophic but not butyrate producer
#filter for tryptophan auxotrophic bacteria
Auxotrophy_2 <- read.csv ("/Users/svenjabusche/Desktop/Info_allGenomes_Auxotrophy_Protrophy_Metadata.csv")
View(Auxotrophy_2)
numb_B_auxo <- Auxotrophy_2$Compound == "Trp" & Auxotrophy_2$Prototrophy == 0
numb_B_auxo <- Auxotrophy_2[numb_B_auxo, ]
View(numb_B_auxo)
#delete archaen
numb_B_auxo<- numb_B_auxo[!(numb_B_auxo$phylum == "Euryarchaeota" | numb_B_auxo$phylum == "Thermoplasmatota"), ]
View(numb_B_auxo)
#calculate number for contingency table
#Butyrate production
nrow(B_trp_A)
nrow(B_trp_Pr)
#number trp auxotrophic and non Butyrate producer
sum_no_B <- nrow(numb_B_auxo) - nrow(B_trp_A)
sum_no_B
#number trp prototrophic and non butyrate producer
sum <- nrow(B_trp_A) + nrow(B_trp_Pr) + sum_no_B
5416 - sum 

#filter for non Butyrate producer
numb_trp_auxo_not_B <- numb_B_auxo[!(numb_B_auxo$rxn.name == "Butyrate-e0 Exchange") , ]
View(numb_trp_auxo_not_B)


5416 - 4079
fisher_B <- data.frame("Trp-Auxotrophy" = c(704, 3073), "Trp-Prototrophy"= c(302,1367), 
                       row.names = c("Butyrate production: yes", "Butyrate production: no"))

fisher_B
#control numbers
704 + 302 + 3073 + 1337 == 5416

#mosaicplot
mosaicplot(fisher_B, color = TRUE)

#fisher test
fisher.test(fisher_B)

#Acetate
#tryptophan auxotrophic
A_trp_auxo <-auxo_Prod1$rxn.name == "Acetate-e0 Exchange" & auxo_Prod1$mtf.flux > 0 & auxo_Prod1$Prototrophy == 0 & auxo_Prod1$Compound == "Trp"
A_trp_auxo  <-auxo_Prod1[A_trp_auxo , ]

#new column with the production rates
A_trp_auxo$prod <- A_trp_auxo$mtf.flux/A_trp_auxo$Growth

#visualization
ggplot(data=A_trp_auxo, aes(rxn.name, prod, fill = phylum))+
  geom_bar(stat = "identity") +
  ggtitle("Acetate production of Tryptophan auxotrophy microbiota")

#tryptophan prototrophic
A_trp_Pr <-auxo_Prod1$rxn.name == "Acetate-e0 Exchange" & auxo_Prod1$mtf.flux > 0 & auxo_Prod1$Prototrophy == 1 & auxo_Prod1$Compound == "Trp"
A_trp_Pr <-auxo_Prod1[A_trp_Pr, ]

#new column with productionrates
A_trp_Pr$prod <- A_trp_Pr$mtf.flux/A_trp_Pr$Growth

#visualization
ggplot(data=A_trp_Proto, aes(rxn.name, prod, fill = phylum))+
  geom_bar(stat = "identity") +
  ggtitle("Acetate production of Tryptophan prototrophy microbiota")

#
test <- wilcox.test(A_trp_auxo$prod, A_trp_Pr$prod,
                    paired = FALSE, conf.level = 0.95)
test

#Protrophy and Auxotorphy data in one dataset
A_trp_auxo_proto <-auxo_Prod1$rxn.name == "Acetate-e0 Exchange" & auxo_Prod1$mtf.flux > 0 & auxo_Prod1$Compound == "Trp"
A_trp_auxo_proto <-auxo_Prod1[A_trp_auxo_proto, ]

#new column with productionrates
A_trp_auxo_proto$prod <- A_trp_auxo_proto$mtf.flux/A_trp_auxo_proto$Growth
A_trp_auxo_proto <-A_trp_auxo_proto[order(A_trp_auxo_proto$prod),]

#
test <- wilcox.test(A_trp_auxo$prod, A_trp_Pr$prod,
                    paired = FALSE)
test

# plot with statistical results --> independent (n ist unterschiedlich und auxotroph/prototrophe Stichprobe)
ggbetweenstats( # independent samples
  data = A_trp_auxo_proto,
  x = Prototrophy,
  y = prod,
  plot.type = "boxviolin", # for boxplot
  type = "nonparametric", # for wilcoxon
  centrality.plotting = TRUE,
  pairwise.display ="significant",
  title = "Comparison of the acetate productionrates by tryptophan auxotrophic \nand prototrophic microbiota",
  ylab = "Productionrates[mmol/gDW]",
  xlab = "",
  results.subtitle = TRUE,
  conf.level = 0.95,
  caption = "0 - Auxotrophy, 1 - Prototrophy\nMann-Whitney U Test, p < 0.05, CI = 0.95%"
)

p_A1 <- ggbetweenstats(
  data = A_trp_auxo_proto,
  x = Prototrophy,
  y = prod,
  plot.type = "violin", # for boxplot
  type = "nonparametric", # for wilcoxon
  centrality.plotting = TRUE,
  pairwise.display ="significant",
  title = "All phyla",
  ylab = "Productionrates [mmol/g DW]",
  xlab = "",
  results.subtitle = TRUE,
  conf.level = 0.95,
  caption = "",
  point.args = list(
    position = ggplot2::position_jitterdodge(dodge.width = 0.20),
    alpha = 0.2,
    size = 1,
    stroke = 0)
)  + theme(text = element_text(size= 10)) +
  theme(plot.subtitle = element_text(size=6)) +
  theme(plot.title = element_text(size=10))

#filter phyla
#Bacteroidetes
A_trp_Bact <-auxo_Prod1$rxn.name == "Acetate-e0 Exchange" & auxo_Prod1$mtf.flux > 0 & auxo_Prod1$Compound == "Trp" & auxo_Prod1$phylum == "Bacteroidota"
A_trp_Bact <-auxo_Prod1[A_trp_Bact, ]
View(A_trp_Bact)



#new column with productionrates
A_trp_Bact$prod <- A_trp_Bact$mtf.flux/A_trp_Bact$Growth
A_trp_Bact <-A_trp_Bact[order(A_trp_Bact$prod),]

#remove all NA values
B_Trp <- B_trp[complete.cases(B_trp), ]

#visualization
ggbetweenstats(
  data = A_trp_Bact,
  x = Prototrophy,
  y = prod,
  plot.type = "violin", # for boxplot
  type = "nonparametric", # for wilcoxon
  centrality.plotting = TRUE,
  pairwise.display ="significant",
  title = "Comparison of the acetate productionrates by tryptophan auxotrophic \nand prototrophic Bacteroidetes",
  ylab = "Productionrates [mmol/gDW]",
  xlab = "",
  results.subtitle = TRUE,
  conf.level = 0.95,
  caption = "0 - Auxotrophy, 1 - Prototrophy\nMann-Whitney U Test (CI = 0.95%)",
  point.args = list(
    position = ggplot2::position_jitterdodge(dodge.width = 0.20),
    alpha = 0.2,
    size = 1,
    stroke = 0)
)

p_A2 <- ggbetweenstats(
  data = A_trp_Bact,
  x = Prototrophy,
  y = prod,
  plot.type = "violin", # for boxplot
  type = "nonparametric", # for wilcoxon
  centrality.plotting = TRUE,
  pairwise.display ="significant",
  title = "Bacteroidetes",
  ylab = "Productionrates [mmol/g DW]",
  xlab = "",
  results.subtitle = TRUE,
  conf.level = 0.95,
  caption = "",
  point.args = list(
    position = ggplot2::position_jitterdodge(dodge.width = 0.20),
    alpha = 0.2,
    size = 1,
    stroke = 0)
)  + theme(text = element_text(size= 10)) +
  theme(plot.subtitle = element_text(size=6)) +
  theme(plot.title = element_text(size=10))


#Firmicutes#
A_trp_Firm <-auxo_Prod1$rxn.name == "Acetate-e0 Exchange" & auxo_Prod1$mtf.flux > 0 & auxo_Prod1$Compound == "Trp" & auxo_Prod1$phylum == "Firmicutes"
A_trp_Firm <-auxo_Prod1[A_trp_Firm, ]
View(A_trp_Firm)

#new column with productionrates
A_trp_Firm$prod <- A_trp_Firm$mtf.flux/A_trp_Firm$Growth
A_trp_Firm <-A_trp_Firm[order(A_trp_Firm$prod),]

#remove all NA values
B_Trp <- B_trp[complete.cases(B_trp), ]

#visualization
ggbetweenstats(
  data = A_trp_Firm,
  x = Prototrophy,
  y = prod,
  plot.type = "violin", # for boxplot
  type = "nonparametric", # for wilcoxon
  centrality.plotting = TRUE,
  pairwise.display ="significant",
  title = "Comparison of the acetate productionrates by tryptophan auxotrophic \nand prototrophic Firmicutes",
  ylab = "Productionrates [mmol/gDW]",
  xlab = "",
  results.subtitle = TRUE,
  conf.level = 0.95,
  caption = "0 - Auxotrophy, 1 - Prototrophy\nMann-Whitney U Test (CI = 0.95%)",
  point.args = list(
    position = ggplot2::position_jitterdodge(dodge.width = 0.20),
    alpha = 0.2,
    size = 1,
    stroke = 0)
)
p_A3 <- ggbetweenstats(
  data = A_trp_Firm,
  x = Prototrophy,
  y = prod,
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
)  + theme(text = element_text(size= 10)) +
  theme(plot.subtitle = element_text(size=5.5)) +
  theme(plot.title = element_text(size=10))

#Proteobacteria
A_trp_Proteo <-auxo_Prod1$rxn.name == "Acetate-e0 Exchange" & auxo_Prod1$mtf.flux > 0 & auxo_Prod1$Compound == "Trp" & auxo_Prod1$phylum == "Proteobacteria"
A_trp_Proteo <-auxo_Prod1[A_trp_Proteo, ]
View(A_trp_Proteo)

#new column with productionrates
A_trp_Proteo$prod <- A_trp_Proteo$mtf.flux/A_trp_Proteo$Growth
A_trp_Proteo <-A_trp_Proteo[order(A_trp_Proteo$prod),]

#remove all NA values
B_Trp <- B_trp[complete.cases(B_trp), ]

#visualization
ggbetweenstats(
  data = A_trp_Proteo,
  x = Prototrophy,
  y = prod,
  plot.type = "violin", # for boxplot
  type = "nonparametric", # for wilcoxon
  centrality.plotting = TRUE,
  pairwise.display ="significant",
  title = "Comparison of the acetate productionrates by tryptophan auxotrophic \nand prototrophic Proteobacteria",
  ylab = "Productionrates [mmol/gDW]",
  xlab = "",
  results.subtitle = TRUE,
  conf.level = 0.95,
  caption = "0 - Auxotrophy, 1 - Prototrophy\nMann-Whitney U Test (CI = 0.95%)",
  point.args = list(
    position = ggplot2::position_jitterdodge(dodge.width = 0.20),
    alpha = 0.2,
    size = 1,
    stroke = 0)
)
p_A4 <- ggbetweenstats(
  data = A_trp_Proteo,
  x = Prototrophy,
  y = prod,
  plot.type = "violin", # for boxplot
  type = "nonparametric", # for wilcoxon
  centrality.plotting = TRUE,
  pairwise.display ="significant",
  title = "Proteobacteria",
  ylab = "Productionrates [mmol/g DW]",
  xlab = "",
  results.subtitle = TRUE,
  conf.level = 0.95,
  caption = "",
  point.args = list(
    position = ggplot2::position_jitterdodge(dodge.width = 0.20),
    alpha = 0.2,
    size = 1,
    stroke = 0)
)  + theme(text = element_text(size= 10)) +
  theme(plot.subtitle = element_text(size=6)) +
  theme(plot.title = element_text(size=10))

#Actinobacteria
A_trp_Actino <-auxo_Prod1$rxn.name == "Acetate-e0 Exchange" & auxo_Prod1$mtf.flux > 0 & auxo_Prod1$Compound == "Trp" & auxo_Prod1$phylum == "Actinobacteriota"
A_trp_Actino <-auxo_Prod1[A_trp_Actino, ]
View(A_trp_Actino)

remove(A_trp_Actino)
#new column with productionrates
A_trp_Actino$prod <- A_trp_Actino$mtf.flux/A_trp_Actino$Growth
A_trp_Actino <-A_trp_Actino[order(A_trp_Actino$prod),]

#remove all NA values
B_Trp <- B_trp[complete.cases(B_trp), ]

#visualization
ggbetweenstats(
  data = A_trp_Actino,
  x = Prototrophy,
  y = prod,
  plot.type = "violin", # for boxplot
  type = "nonparametric", # for wilcoxon
  centrality.plotting = TRUE,
  pairwise.display ="significant",
  title = "Comparison of the acetate productionrates by tryptophan auxotrophic \nand prototrophic Actinobacteria",
  ylab = "Productionrates [mmol/gDW]",
  xlab = "",
  results.subtitle = TRUE,
  conf.level = 0.95,
  caption = "0 - Auxotrophy, 1 - Prototrophy\nMann-Whitney U Test (CI = 0.95%)",
  point.args = list(
    position = ggplot2::position_jitterdodge(dodge.width = 0.20),
    alpha = 0.2,
    size = 1,
    stroke = 0)
)
p_A5 <- ggbetweenstats(
  data = A_trp_Actino,
  x = Prototrophy,
  y = prod,
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
    stroke = 0)
)  + theme(text = element_text(size= 10)) +
  theme(plot.subtitle = element_text(size=6)) +
  theme(plot.title = element_text(size=10))

#combine all plots
combine_plots(list(p_A2, p_A3, p_A4, p_A5),
              plotgrid.args = list(nrow=2),
              annotation.args = list(title = "Comparison of acetate productionrates by tryptophan auxotrophic and\nprototrophic microbiota",
                                     caption = "Legend: 0 - Auxotrophy, 1 - Prototrophy, Test: Mann-Whitney U Test (CI = 0.95%)",
                                     theme = ggplot2::theme(
                                       plot.subtitle = element_text(size = 6),
                                       axis.title.y = element_text(face = "plain"),
                                       text = element_text(size=8, face = "plain"),
                                       plot.title = element_text(size = 12, face = "bold"))))

#fisher test
#create a dataframe for Acetate production
#number of genomes who are tryptophan auxotrophic but not butyrate producer
#filter for tryptophan auxotrophic bacteria
numb_auxo <- Auxotrophy_2$Compound == "Trp" & Auxotrophy_2$Prototrophy == 0
numb_auxo <- Auxotrophy_2[numb_auxo, ]
View(numb_B_auxo)
#delete archaen
numb_auxo<- numb_auxo[!(numb_auxo$phylum == "Euryarchaeota" | numb_auxo$phylum == "Thermoplasmatota"), ]
View(numb_auxo)
#calculate number for contingency table
#Acetate production
nrow(A_trp_A)
nrow(A_trp_Pr)
#number trp auxotrophic and non Butyrate producer
sum_no_A_auxo <- nrow(numb_auxo) - nrow(A_trp_A)
sum_no_A_auxo
#number trp prototrophic and non butyrate producer
sum_no_A_proto <- 5416 -(nrow(A_trp_A) + nrow(A_trp_Pr) + sum_no_A_auxo)
sum_no_A_proto 

#create contigency table
fisher_A <- data.frame("Trp-Auxotrophy" = c(nrow(A_trp_A), sum_no_A_auxo), "Trp-Prototrophy"= c(nrow(A_trp_Pr),sum_no_A_proto), 
                       row.names = c("Acetate production: yes", "Acetate production: no"))

fisher_A
#control numbers
nrow(A_trp_A) + nrow(A_trp_Pr) + sum_no_A_auxo + sum_no_A_proto == 5416

#create mosaic plot
mosaicplot(fisher_A, color = TRUE, main="Acetate production of trp-auxotrophic and prototrophic microbiota")
#fisher test
fisher.test(fisher_A)

#Propionate
#tryptophan auxotrophic
P_trp_A <-auxo_Prod1$rxn.name == "Propionate-e0 Exchange" & auxo_Prod1$mtf.flux > 0 & auxo_Prod1$Prototrophy == 0 & auxo_Prod1$Compound == "Trp"
P_trp_A <-auxo_Prod1[P_trp_A, ]

#new column with the production rates
P_trp_A$prod <- P_trp_A$mtf.flux/ P_trp_A$Growth

ggplot(data=P_trp_A, aes(rxn.name, prod, fill = phylum))+
  geom_bar(stat = "identity") +
  ggtitle("Propionate production of Tryptophan auxotrophy microbiota")

#tryptophan prototrophic
P_trp_Pr <-auxo_Prod1$rxn.name == "Propionate-e0 Exchange" & auxo_Prod1$mtf.flux > 0 & auxo_Prod1$Prototrophy == 1 & auxo_Prod1$Compound == "Trp"
P_trp_Pr <-auxo_Prod1[P_trp_Pr, ]

View(P_trp_Proto)
#change class from character to numeric for columns mtf.flux and Growth
class(P_trp_A$mtf.flux)

P_trp_Pr$mtf.flux <- as.numeric(as.character(P_trp_Pr$mtf.flux))
P_trp_Pr$Growth <- as.numeric(as.character(P_trp_Pr$Growth))

#new column with productionrates
P_trp_Pr$prod <- P_trp_Pr$mtf.flux/P_trp_Pr$Growth

#remove all NA values
P_trp_Proto <- P_trp_Pr[complete.cases(P_trp_Pr), ]

ggplot(data=P_trp_Proto, aes(rxn.name, prod, fill = phylum))+
  geom_bar(stat = "identity") +
  ggtitle("Propionate production of Tryptophan prototrophy microbiota")

#
test <- wilcox.test(P_trp_Auxo$prod, P_trp_Proto$prod,
                    paired = FALSE)
test

#Protrophy and Auxotorphy data in one dataset
P_trp <-auxo_Prod1$rxn.name == "Propionate-e0 Exchange" & auxo_Prod1$mtf.flux > 0 & auxo_Prod1$Compound == "Trp"
P_trp <-auxo_Prod1[P_trp, ]

#new column with productionrates
P_trp$prod <- P_trp$mtf.flux/P_trp$Growth
P_trp <-P_trp[order(P_trp$prod),]

# plot with statistical results --> independent (n ist unterschiedlich und auxotroph/prototrophe Stichprobe)
ggbetweenstats( # independent samples
  data = P_trp,
  x = Prototrophy,
  y = prod,
  plot.type = "boxviolin", 
  type = "nonparametric", # for wilcoxon
  pairwise.display ="significant",
  title = "Comparison of the propionate productionrates by tryptophan auxotrophic \nand prototrophic microbiota", 
  ylab = "Productionrates [mmol/gDW]",
  xlab = "",
  results.subtitle = TRUE,
  conf.level = 0.95,
  caption = "0 - Auxotrophy, 1 - Prototrophy\nMann-Whitney U Test, p < 0.05, CI = 0.95%"
)

p_P1 <- ggbetweenstats(
  data = P_trp,
  x = Prototrophy,
  y = prod,
  plot.type = "violin",
  type = "nonparametric", # for wilcoxon
  centrality.plotting = TRUE,
  pairwise.display ="significant",
  title = "All phyla",
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
)  + theme(text = element_text(size= 10)) +
  theme(plot.subtitle = element_text(size=6)) +
  theme(plot.title = element_text(size=10))
ggbetweenstats(data=P_trp, x= Prototrophy, y= prod,
               type=)

#filter phyla
#Bacteroidetes
P_trp_Bact <-auxo_Prod1$rxn.name == "Propionate-e0 Exchange" & auxo_Prod1$mtf.flux > 0 & auxo_Prod1$Compound == "Trp" & auxo_Prod1$phylum == "Bacteroidota"
P_trp_Bact <-auxo_Prod1[P_trp_Bact, ]
View(P_trp_Bact)

#new column with productionrates
P_trp_Bact$prod <- P_trp_Bact$mtf.flux/P_trp_Bact$Growth
P_trp_Bact <-P_trp_Bact[order(P_trp_Bact$prod),]

#visualization
ggbetweenstats(
  data = P_trp_Bact,
  x = Prototrophy,
  y = prod,
  plot.type = "violin",
  type = "nonparametric", # for wilcoxon
  centrality.plotting = TRUE,
  pairwise.display ="significant",
  title = "Comparison of the propionate productionrates by tryptophan auxotrophic \nand prototrophic Bacteroidetes",
  ylab = "Productionrates [mmol/gDW]",
  xlab = "",
  results.subtitle = TRUE,
  conf.level = 0.95,
  caption = "0 - Auxotrophy, 1 - Prototrophy\nMann-Whitney U Test (CI = 0.95%)",
  point.args = list(
    position = ggplot2::position_jitterdodge(dodge.width = 0.20),
    alpha = 0.2,
    size = 1,
    stroke = 0)
)

p_P2 <- ggbetweenstats(
  data = P_trp_Bact,
  x = Prototrophy,
  y = prod,
  plot.type = "violin", 
  type = "nonparametric", # for wilcoxon
  centrality.plotting = TRUE,
  pairwise.display ="significant",
  title = "Bacteroidetes",
  ylab = "Productionrates [mmol/g DW]",
  xlab = "",
  results.subtitle = TRUE,
  conf.level = 0.95,
  caption = "",
  point.args = list(
    position = ggplot2::position_jitterdodge(dodge.width = 0.20),
    alpha = 0.2,
    size = 1,
    stroke = 0)
)  + theme(text = element_text(size= 10)) +
  theme(plot.subtitle = element_text(size=6)) +
  theme(plot.title = element_text(size=10))


#Firmicutes
P_trp_Firm <-auxo_Prod1$rxn.name == "Propionate-e0 Exchange" & auxo_Prod1$mtf.flux > 0 & auxo_Prod1$Compound == "Trp" & auxo_Prod1$phylum == "Firmicutes"
P_trp_Firm <-auxo_Prod1[P_trp_Firm, ]
View(P_trp_Firm)

#new column with productionrates
P_trp_Firm$prod <- P_trp_Firm$mtf.flux/P_trp_Firm$Growth
P_trp_Firm <-P_trp_Firm[order(P_trp_Firm$prod),]

#visualization
ggbetweenstats(
  data = P_trp_Firm,
  x = Prototrophy,
  y = prod,
  plot.type = "violin", 
  type = "nonparametric", # for wilcoxon
  centrality.plotting = TRUE,
  pairwise.display ="significant",
  title = "Comparison of the propionate productionrates by tryptophan auxotrophic \nand prototrophic Firmicutes",
  ylab = "Productionrates [mmol/gDW]",
  xlab = "",
  results.subtitle = TRUE,
  conf.level = 0.95,
  caption = "0 - Auxotrophy, 1 - Prototrophy\nMann-Whitney U Test (CI = 0.95%)",
  point.args = list(
    position = ggplot2::position_jitterdodge(dodge.width = 0.20),
    alpha = 0.2,
    size = 1,
    stroke = 0)
)
p_P3 <- ggbetweenstats(
  data = P_trp_Firm,
  x = Prototrophy,
  y = prod,
  plot.type = "violin",
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
)  + theme(text = element_text(size= 10)) +
  theme(plot.subtitle = element_text(size=5.5)) +
  theme(plot.title = element_text(size=10))

#Proteobacteria
P_trp_Proteo <-auxo_Prod1$rxn.name == "Propionate-e0 Exchange" & auxo_Prod1$mtf.flux > 0 & auxo_Prod1$Compound == "Trp" & auxo_Prod1$phylum == "Proteobacteria"
P_trp_Proteo <-auxo_Prod1[P_trp_Proteo, ]
View(P_trp_Proteo)

#new column with productionrates
P_trp_Proteo$prod <- P_trp_Proteo$mtf.flux/P_trp_Proteo$Growth
P_trp_Proteo <-P_trp_Proteo[order(P_trp_Proteo$prod),]

#remove all NA values
P_Trp <- P_trp[complete.cases(P_trp), ]

#visualization
ggbetweenstats(
  data = P_trp_Proteo,
  x = Prototrophy,
  y = prod,
  plot.type = "violin", # for boxplot
  type = "nonparametric", # for wilcoxon
  centrality.plotting = TRUE,
  pairwise.display ="significant",
  title = "Comparison of the propionate productionrates by tryptophan auxotrophic \nand prototrophic Proteobacteria",
  ylab = "Productionrates [mmol/gDW]",
  xlab = "",
  results.subtitle = TRUE,
  conf.level = 0.95,
  caption = "0 - Auxotrophy, 1 - Prototrophy\nMann-Whitney U Test (CI = 0.95%)",
  point.args = list(
    position = ggplot2::position_jitterdodge(dodge.width = 0.20),
    alpha = 0.2,
    size = 1,
    stroke = 0)
)
p_P4 <- ggbetweenstats(
  data = P_trp_Proteo,
  x = Prototrophy,
  y = prod,
  plot.type = "violin", # for boxplot
  type = "nonparametric", # for wilcoxon
  centrality.plotting = TRUE,
  pairwise.display ="significant",
  title = "Proteobacteria",
  ylab = "Productionrates [mmol/g DW]",
  xlab = "",
  results.subtitle = TRUE,
  conf.level = 0.95,
  caption = "",
  point.args = list(
    position = ggplot2::position_jitterdodge(dodge.width = 0.20),
    alpha = 0.2,
    size = 1,
    stroke = 0)
)  + theme(text = element_text(size= 10)) +
  theme(plot.subtitle = element_text(size=6)) +
  theme(plot.title = element_text(size=10))

#Actinobacteria
P_trp_Actino <-auxo_Prod1$rxn.name == "Propionate-e0 Exchange" & auxo_Prod1$mtf.flux > 0 & auxo_Prod1$Compound == "Trp" & auxo_Prod1$phylum == "Actinobacteriota"
P_trp_Actino <-auxo_Prod1[P_trp_Actino, ]
View(P_trp_Proteo)

#new column with productionrates
P_trp_Actino$prod <- P_trp_Actino$mtf.flux/P_trp_Actino$Growth
P_trp_Actino <-P_trp_Actino[order(P_trp_Actino$prod),]

#remove all NA values
P_trp_Actino <- P_trp_Actino[complete.cases(P_trp_Actino), ]

#visualization
ggbetweenstats(
  data = P_trp_Actino,
  x = Prototrophy,
  y = prod,
  plot.type = "violin",
  type = "nonparametric", # for wilcoxon
  centrality.plotting = TRUE,
  pairwise.display ="significant",
  title = "Comparison of the propionate productionrates by tryptophan auxotrophic \nand prototrophic Actinobacteria",
  ylab = "Productionrates [mmol/gDW]",
  xlab = "",
  results.subtitle = TRUE,
  conf.level = 0.95,
  caption = "0 - Auxotrophy, 1 - Prototrophy\nMann-Whitney U Test (CI = 0.95%)",
  point.args = list(
    position = ggplot2::position_jitterdodge(dodge.width = 0.20),
    alpha = 0.2,
    size = 1,
    stroke = 0)
)
p_P5 <- ggbetweenstats(
  data = P_trp_Actino,
  x = Prototrophy,
  y = prod,
  plot.type = "violin",
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
    stroke = 0)
)  + theme(text = element_text(size= 10)) +
  theme(plot.subtitle = element_text(size=6)) +
  theme(plot.title = element_text(size=10))

#combine all plots
combine_plots(list(p_P2, p_P3, p_P4, p_P5),
              plotgrid.args = list(nrow=2),
              annotation.args = list(title = "Comparison of propionate productionrates by tryptophan auxotrophic and\nprototrophic microbiota",
                                     caption = "Legend: 0 - Auxotrophy, 1 - Prototrophy, Test: Mann-Whitney U Test (CI = 0.95%)",
                                     theme = ggplot2::theme(
                                       plot.subtitle = element_text(size = 6),
                                       axis.title.y = element_text(face = "plain"),
                                       text = element_text(size=8, face = "plain"),
                                       plot.title = element_text(size = 12, face = "bold"))))

#fisher test
#create a dataframe for propionate production
#number of genomes who are tryptophan auxotrophic but not propioante producer
#filter for tryptophan auxotrophic bacteria
numb_auxo <- Auxotrophy_2$Compound == "Trp" & Auxotrophy_2$Prototrophy == 0
numb_auxo <- Auxotrophy_2[numb_auxo, ]
View(numb_B_auxo)
#delete archaen
numb_auxo<- numb_auxo[!(numb_auxo$phylum == "Euryarchaeota" | numb_auxo$phylum == "Thermoplasmatota"), ]
View(numb_auxo)
#calculate number for contingency table
#Propionate production
nrow(P_trp_A)
nrow(P_trp_Pr)
#number trp auxotrophic and non Butyrate producer
sum_no_P_auxo <- nrow(numb_auxo) - nrow(P_trp_A)
sum_no_P_auxo
#number trp prototrophic and non butyrate producer
sum_no_P_proto <- 5416 -(nrow(P_trp_A) + nrow(P_trp_Pr) + sum_no_P_auxo)
sum_no_P_proto 

#create contigency table
fisher_P <- data.frame("Trp-Auxotrophy" = c(nrow(P_trp_A), sum_no_P_auxo), 
                       "Trp-Prototrophy"= c(nrow(P_trp_Pr),sum_no_P_proto), 
                       row.names = c("Propionate production: yes", "Propionate production: no"))

fisher_P
#control numbers
nrow(P_trp_A) + nrow(P_trp_Pr) + sum_no_P_auxo + sum_no_P_proto == 5416

#create mosaic plot
mosaicplot(fisher_P, color = TRUE, main="Propionate production of trp auxotrophic and prototrophic microbiota")
#fisher test
fisher.test(fisher_P)

#Lactate
#L_trp_A <- auxo_Prod1[grepl("L-Lactate", rxn.name) & mtf.flux > 0 & Prototrophy == 0 & Compound == "Trp"][order(-mtf.flux)]

#tryptophan auxotrophic
L_trp_A <-auxo_Prod1$rxn.name == "L-Lactate-e0 Exchange" & auxo_Prod1$mtf.flux > 0 & auxo_Prod1$Prototrophy == 0 & auxo_Prod1$Compound == "Trp"
L_trp_A <-auxo_Prod1[L_trp_A, ]
View(L_trp_A)

#new column with the production rates
L_trp_A$prod <- L_trp_A$mtf.flux/ L_trp_A$Growth

ggplot(data=L_trp_A, aes(rxn.name, prod, fill = phylum))+
  geom_bar(stat = "identity") +
  ggtitle("L-Lactateproduction of Tryptophan auxotrophy microbiota")

#tryptophan prototrophic
L_trp_Pr <-auxo_Prod1$rxn.name == "L-Lactate-e0 Exchange" & auxo_Prod1$mtf.flux > 0 & auxo_Prod1$Prototrophy == 1 & auxo_Prod1$Compound == "Trp"
L_trp_Pr <-auxo_Prod1[L_trp_Pr, ]

View(P_trp_Proto)
#new column with productionrates
L_trp_Pr$prod <- L_trp_Pr$mtf.flux/L_trp_Pr$Growth

ggplot(data=L_trp_Proto, aes(rxn.name, prod, fill = phylum))+
  geom_bar(stat = "identity") +
  ggtitle("Propionate production of Tryptophan prototrophy microbiota")

#
test <- wilcox.test(L_trp_Auxo$prod, L_trp_Proto$prod,
                    paired = FALSE)
test

#Protrophy and Auxotorphy data in one dataset
L_trp <-auxo_Prod1$rxn.name == "L-Lactate-e0 Exchange" & auxo_Prod1$mtf.flux > 0 & auxo_Prod1$Compound == "Trp"
L_trp <-auxo_Prod1[L_trp, ]

#new column with productionrates
L_trp$prod <- L_trp$mtf.flux/L_trp$Growth
L_trp <-L_trp[order(L_trp$prod),]

# plot with statistical results --> independent (n ist unterschiedlich und auxotroph/prototrophe Stichprobe)
ggbetweenstats( # independent samples
  data = L_trp,
  x = Prototrophy,
  y = prod,
  plot.type = "boxviolin", 
  type = "nonparametric", # for wilcoxon
  pairwise.display ="significant",
  title = "Comparison of the L-Lactate productionrates by tryptophan auxotrophic \nand prototrophic microbiota", 
  ylab = "Productionrates [mmol/gDW]",
  xlab = "",
  results.subtitle = TRUE,
  conf.level = 0.95,
  caption = "0 - Auxotrophy, 1 - Prototrophy\nMann-Whitney U Test, p < 0.05, CI = 0.95%"
)

p_L1 <- ggbetweenstats(
  data = L_trp,
  x = Prototrophy,
  y = prod,
  plot.type = "violin", 
  type = "nonparametric", # for wilcoxon
  centrality.plotting = TRUE,
  pairwise.display ="significant",
  title = "All phyla",
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
)  + theme(text = element_text(size= 10)) +
  theme(plot.subtitle = element_text(size=6)) +
  theme(plot.title = element_text(size=10))
ggbetweenstats(data=P_trp, x= Prototrophy, y= prod,
               type=)


#filter phyla
#Bacteroidetes
L_trp_Bact <-auxo_Prod1$rxn.name == "L-Lactate-e0 Exchange" & auxo_Prod1$mtf.flux > 0 & auxo_Prod1$Compound == "Trp" & auxo_Prod1$phylum == "Bacteroidota"
L_trp_Bact <-auxo_Prod1[L_trp_Bact, ]
View(L_trp_Bact)

#new column with productionrates
L_trp_Bact$prod <- L_trp_Bact$mtf.flux/L_trp_Bact$Growth
L_trp_Bact <-L_trp_Bact[order(L_trp_Bact$prod),]

#visualization
ggbetweenstats(
  data = L_trp_Bact,
  x = Prototrophy,
  y = prod,
  plot.type = "violin", 
  type = "nonparametric", # for wilcoxon
  centrality.plotting = TRUE,
  pairwise.display ="significant",
  title = "Comparison of the L-Lactate productionrates by tryptophan auxotrophic \nand prototrophic Bacteroidetes",
  ylab = "Productionrates [mmol/gDW]",
  xlab = "",
  results.subtitle = TRUE,
  conf.level = 0.95,
  caption = "0 - Auxotrophy, 1 - Prototrophy\nMann-Whitney U Test (CI = 0.95%)",
  point.args = list(
    position = ggplot2::position_jitterdodge(dodge.width = 0.20),
    alpha = 0.2,
    size = 1,
    stroke = 0)
)

p_L2 <- ggbetweenstats(
  data = L_trp_Bact,
  x = Prototrophy,
  y = prod,
  plot.type = "violin", 
  type = "nonparametric", # for wilcoxon
  centrality.plotting = TRUE,
  pairwise.display ="significant",
  title = "Bacteroidetes",
  ylab = "Productionrates [mmol/g DW]",
  xlab = "",
  results.subtitle = TRUE,
  conf.level = 0.95,
  caption = "",
  point.args = list(
    position = ggplot2::position_jitterdodge(dodge.width = 0.20),
    alpha = 0.2,
    size = 1,
    stroke = 0)
)  + theme(text = element_text(size= 10)) +
  theme(plot.subtitle = element_text(size=6)) +
  theme(plot.title = element_text(size=10))


#Firmicutes#
L_trp_Firm <-auxo_Prod1$rxn.name == "L-Lactate-e0 Exchange" & auxo_Prod1$mtf.flux > 0 & auxo_Prod1$Compound == "Trp" & auxo_Prod1$phylum == "Firmicutes"
L_trp_Firm <-auxo_Prod1[L_trp_Firm, ]
View(P_trp_Firm)

#new column with productionrates
L_trp_Firm$prod <- L_trp_Firm$mtf.flux/L_trp_Firm$Growth
L_trp_Firm <-L_trp_Firm[order(L_trp_Firm$prod),]

#visualization
ggbetweenstats(
  data = L_trp_Firm,
  x = Prototrophy,
  y = prod,
  plot.type = "violin", # for boxplot
  type = "nonparametric", # for wilcoxon
  centrality.plotting = TRUE,
  pairwise.display ="significant",
  title = "Comparison of the L-Lactate productionrates by tryptophan auxotrophic \nand prototrophic Firmicutes",
  ylab = "Productionrates [mmol/gDW]",
  xlab = "",
  results.subtitle = TRUE,
  conf.level = 0.95,
  caption = "0 - Auxotrophy, 1 - Prototrophy\nMann-Whitney U Test (CI = 0.95%)",
  point.args = list(
    position = ggplot2::position_jitterdodge(dodge.width = 0.20),
    alpha = 0.2,
    size = 1,
    stroke = 0)
)
p_L3 <- ggbetweenstats(
  data = L_trp_Firm,
  x = Prototrophy,
  y = prod,
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
)  + theme(text = element_text(size= 10)) +
  theme(plot.subtitle = element_text(size=5.5)) +
  theme(plot.title = element_text(size=10))

#Proteobacteria
L_trp_Proteo <-auxo_Prod1$rxn.name == "L-Lactate-e0 Exchange" & auxo_Prod1$mtf.flux > 0 & auxo_Prod1$Compound == "Trp" & auxo_Prod1$phylum == "Proteobacteria"
L_trp_Proteo <-auxo_Prod1[L_trp_Proteo, ]
View(L_trp_Proteo)

#new column with productionrates
L_trp_Proteo$prod <- L_trp_Proteo$mtf.flux/L_trp_Proteo$Growth
L_trp_Proteo <-L_trp_Proteo[order(L_trp_Proteo$prod),]

#visualization
ggbetweenstats(
  data = L_trp_Proteo,
  x = Prototrophy,
  y = prod,
  plot.type = "violin", 
  type = "nonparametric", # for wilcoxon
  centrality.plotting = TRUE,
  pairwise.display ="significant",
  title = "Comparison of the L-Lactate productionrates by tryptophan auxotrophic \nand prototrophic Proteobacteria",
  ylab = "Productionrates [mmol/gDW]",
  xlab = "",
  results.subtitle = TRUE,
  conf.level = 0.95,
  caption = "0 - Auxotrophy, 1 - Prototrophy\nMann-Whitney U Test (CI = 0.95%)",
  point.args = list(
    position = ggplot2::position_jitterdodge(dodge.width = 0.20),
    alpha = 0.2,
    size = 1,
    stroke = 0)
)
p_L4 <- ggbetweenstats(
  data = L_trp_Proteo,
  x = Prototrophy,
  y = prod,
  plot.type = "violin", 
  type = "nonparametric", # for wilcoxon
  centrality.plotting = TRUE,
  pairwise.display ="significant",
  title = "Proteobacteria",
  ylab = "Productionrates [mmol/g DW]",
  xlab = "",
  results.subtitle = TRUE,
  conf.level = 0.95,
  caption = "",
  point.args = list(
    position = ggplot2::position_jitterdodge(dodge.width = 0.20),
    alpha = 0.2,
    size = 1,
    stroke = 0)
)  + theme(text = element_text(size= 10)) +
  theme(plot.subtitle = element_text(size=6)) +
  theme(plot.title = element_text(size=10))

#Actinobacteria
L_trp_Actino <-auxo_Prod1$rxn.name == "L-Lactate-e0 Exchange" & auxo_Prod1$mtf.flux > 0 & auxo_Prod1$Compound == "Trp" & auxo_Prod1$phylum == "Actinobacteriota"
L_trp_Actino <-auxo_Prod1[L_trp_Actino, ]
View(L_trp_Actino)

#new column with productionrates
L_trp_Actino$prod <- L_trp_Actino$mtf.flux/L_trp_Actino$Growth
L_trp_Actino <-L_trp_Actino[order(L_trp_Actino$prod),]

#visualization
ggbetweenstats(
  data = L_trp_Actino,
  x = Prototrophy,
  y = prod,
  plot.type = "violin",
  type = "nonparametric", # for wilcoxon
  centrality.plotting = TRUE,
  pairwise.display ="significant",
  title = "Comparison of the L-Lactate productionrates by tryptophan auxotrophic \nand prototrophic Actinobacteria",
  ylab = "Productionrates [mmol/gDW]",
  xlab = "",
  results.subtitle = TRUE,
  conf.level = 0.95,
  caption = "0 - Auxotrophy, 1 - Prototrophy\nMann-Whitney U Test (CI = 0.95%)",
  point.args = list(
    position = ggplot2::position_jitterdodge(dodge.width = 0.20),
    alpha = 0.2,
    size = 1,
    stroke = 0)
)
p_L5 <- ggbetweenstats(
  data = L_trp_Actino,
  x = Prototrophy,
  y = prod,
  plot.type = "violin",
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
    stroke = 0)
)  + theme(text = element_text(size= 10)) +
  theme(plot.subtitle = element_text(size=6)) +
  theme(plot.title = element_text(size=10))

#combine all plots
combine_plots(list(p_L3, p_L4, p_L5),
              plotgrid.args = list(nrow=2),
              annotation.args = list(title = "Comparison of L-Lactate productionrates by tryptophan auxotrophic and\nprototrophic microbiota",
                                     caption = "Legend: 0 - Auxotrophy, 1 - Prototrophy, Test: Mann-Whitney U Test (CI = 0.95%)",
                                     theme = ggplot2::theme(
                                       plot.subtitle = element_text(size = 6),
                                       axis.title.y = element_text(face = "plain"),
                                       text = element_text(size=8, face = "plain"),
                                       plot.title = element_text(size = 12, face = "bold"))))

#fisher test
#create a dataframe for lactate production
#number of genomes who are tryptophan auxotrophic but not lactate producer
#filter for tryptophan auxotrophic bacteria
numb_auxo <- Auxotrophy_2$Compound == "Trp" & Auxotrophy_2$Prototrophy == 0
numb_auxo <- Auxotrophy_2[numb_auxo, ]
View(numb_B_auxo)
#delete archaen
numb_auxo<- numb_auxo[!(numb_auxo$phylum == "Euryarchaeota" | numb_auxo$phylum == "Thermoplasmatota"), ]
View(numb_auxo)
#calculate number for contingency table
#Propionate production
nrow(L_trp_A)
nrow(L_trp_Pr)
#number trp auxotrophic and non Butyrate producer
sum_no_L_auxo <- nrow(numb_auxo) - nrow(L_trp_A)
sum_no_L_auxo
#number trp prototrophic and non butyrate producer
sum_no_L_proto <- 5416 -(nrow(L_trp_A) + nrow(L_trp_Pr) + sum_no_L_auxo)
sum_no_L_proto 

#create contigency table
fisher_P <- data.frame("Trp-Auxotrophy" = c(nrow(L_trp_A), sum_no_L_auxo), 
                       "Trp-Prototrophy"= c(nrow(L_trp_Pr),sum_no_L_proto), 
                       row.names = c("Propionate production: yes", "Propionate production: no"))

fisher_P
#control numbers
nrow(P_trp_A) + nrow(P_trp_Pr) + sum_no_P_auxo + sum_no_P_proto == 5416

#create mosaic plot
mosaicplot(fisher_P, color = TRUE, main="Lactate production of trp auxotrophic and prototrophic microbiota")
#fisher test
fisher.test(fisher_P)

#Formate
#tryptophan auxotrophic
F_trp_A <-auxo_Prod1$rxn.name == "Formate-e0 Exchange" & auxo_Prod1$mtf.flux > 0 & auxo_Prod1$Prototrophy == 0 & auxo_Prod1$Compound == "Trp"
F_trp_A <-auxo_Prod1[F_trp_A, ]
View(F_trp_A)
#new column with the production rates
F_trp_A$prod <- F_trp_A$mtf.flux/ F_trp_A$Growth
#visualization
ggplot(data=F_trp_A, aes(rxn.name, prod, fill = phylum))+
  geom_bar(stat = "identity") +
  ggtitle("Formateproduction of Tryptophan auxotrophy microbiota")

#tryptophan prototrophic
F_trp_Pr <-auxo_Prod1$rxn.name == "Formate-e0 Exchange" & auxo_Prod1$mtf.flux > 0 & auxo_Prod1$Prototrophy == 1 & auxo_Prod1$Compound == "Trp"
F_trp_Pr <-auxo_Prod1[F_trp_Pr, ]
View(F_trp_Proto)

#new column with productionrates
F_trp_Pr$prod <- F_trp_Pr$mtf.flux/F_trp_Pr$Growth

ggplot(data=F_trp_Proto, aes(rxn.name, prod, fill = phylum))+
  geom_bar(stat = "identity") +
  ggtitle("Propionate production of Tryptophan prototrophy microbiota")

#
test <- wilcox.test(F_trp_Auxo$prod, F_trp_Proto$prod,
                    paired = FALSE)
test

#Protrophy and Auxotorphy data in one dataset
F_trp <-auxo_Prod1$rxn.name == "Formate-e0 Exchange" & auxo_Prod1$mtf.flux > 0 & auxo_Prod1$Compound == "Trp"
F_trp <-auxo_Prod1[F_trp, ]

#new column with productionrates
F_trp$prod <- F_trp$mtf.flux/F_trp$Growth
F_trp <-F_trp[order(F_trp$prod),]

# plot with statistical results --> unabnhängige Stichprobe (n ist unterschiedlich und auxotroph/prototrophe Stichprobe)
ggbetweenstats( # independent samples
  data = F_trp,
  x = Prototrophy,
  y = prod,
  plot.type = "boxviolin", 
  type = "nonparametric", # for wilcoxon
  pairwise.display ="significant",
  title = "Comparison of the formate productionrates by tryptophan auxotrophic \nand prototrophic microbiota", 
  ylab = "Productionrates [mmol/gDW]",
  xlab = "",
  results.subtitle = TRUE,
  conf.level = 0.95,
  caption = "0 - Auxotrophy, 1 - Prototrophy\nMann-Whitney U Test, p < 0.05, CI = 0.95%"
)

p_F1 <- ggbetweenstats(
  data = F_trp,
  x = Prototrophy,
  y = prod,
  plot.type = "violin",
  type = "nonparametric", # for wilcoxon
  centrality.plotting = TRUE,
  pairwise.display ="significant",
  title = "All phyla",
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
)  + theme(text = element_text(size= 10)) +
  theme(plot.subtitle = element_text(size=6)) +
  theme(plot.title = element_text(size=10))
ggbetweenstats(data=P_trp, x= Prototrophy, y= prod,
               type=)


#filter phyla
#Bacteroidetes
F_trp_Bact <-auxo_Prod1$rxn.name == "Formate-e0 Exchange" & auxo_Prod1$mtf.flux > 0 & auxo_Prod1$Compound == "Trp" & auxo_Prod1$phylum == "Bacteroidota"
F_trp_Bact <-auxo_Prod1[F_trp_Bact, ]
View(F_trp_Bact)

#new column with productionrates
F_trp_Bact$prod <- F_trp_Bact$mtf.flux/F_trp_Bact$Growth
F_trp_Bact <-F_trp_Bact[order(F_trp_Bact$prod),]

#visualization
ggbetweenstats(
  data = F_trp_Bact,
  x = Prototrophy,
  y = prod,
  plot.type = "violin", 
  type = "nonparametric", # for wilcoxon
  centrality.plotting = TRUE,
  pairwise.display ="significant",
  title = "Comparison of formate productionrates by tryptophan auxotrophic \nand prototrophic Bacteroidetes",
  ylab = "Productionrates [mmol/gDW]",
  xlab = "",
  results.subtitle = TRUE,
  conf.level = 0.95,
  caption = "0 - Auxotrophy, 1 - Prototrophy\nMann-Whitney U Test (CI = 0.95%)",
  point.args = list(
    position = ggplot2::position_jitterdodge(dodge.width = 0.20),
    alpha = 0.2,
    size = 1,
    stroke = 0)
)

p_F2 <- ggbetweenstats(
  data = F_trp_Bact,
  x = Prototrophy,
  y = prod,
  plot.type = "violin",
  type = "nonparametric", # for wilcoxon
  centrality.plotting = TRUE,
  pairwise.display ="significant",
  title = "Bacteroidetes",
  ylab = "Productionrates [mmol/g DW]",
  xlab = "",
  results.subtitle = TRUE,
  conf.level = 0.95,
  caption = "",
  point.args = list(
    position = ggplot2::position_jitterdodge(dodge.width = 0.20),
    alpha = 0.2,
    size = 1,
    stroke = 0)
)  + theme(text = element_text(size= 10)) +
  theme(plot.subtitle = element_text(size=6)) +
  theme(plot.title = element_text(size=10))


#Firmicutes
F_trp_Firm <-auxo_Prod1$rxn.name == "Formate-e0 Exchange" & auxo_Prod1$mtf.flux > 0 & auxo_Prod1$Compound == "Trp" & auxo_Prod1$phylum == "Firmicutes"
F_trp_Firm <-auxo_Prod1[F_trp_Firm, ]
View(F_trp_Firm)

#new column with productionrates
F_trp_Firm$prod <- F_trp_Firm$mtf.flux/F_trp_Firm$Growth
F_trp_Firm <-F_trp_Firm[order(F_trp_Firm$prod),]

#visualization
ggbetweenstats(
  data = F_trp_Firm,
  x = Prototrophy,
  y = prod,
  plot.type = "violin", 
  type = "nonparametric", # for wilcoxon
  centrality.plotting = TRUE,
  pairwise.display ="significant",
  title = "Comparison of the Formate productionrates by tryptophan auxotrophic \nand prototrophic Firmicutes",
  ylab = "Productionrates [mmol/gDW]",
  xlab = "",
  results.subtitle = TRUE,
  conf.level = 0.95,
  caption = "0 - Auxotrophy, 1 - Prototrophy\nMann-Whitney U Test (CI = 0.95%)",
  point.args = list(
    position = ggplot2::position_jitterdodge(dodge.width = 0.20),
    alpha = 0.2,
    size = 1,
    stroke = 0)
)
p_F3 <- ggbetweenstats(
  data = F_trp_Firm,
  x = Prototrophy,
  y = prod,
  plot.type = "violin", 
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
)  + theme(text = element_text(size= 10)) +
  theme(plot.subtitle = element_text(size=5.5)) +
  theme(plot.title = element_text(size=10))

#Proteobacteria
F_trp_Proteo <-auxo_Prod1$rxn.name == "Formate-e0 Exchange" & auxo_Prod1$mtf.flux > 0 & auxo_Prod1$Compound == "Trp" & auxo_Prod1$phylum == "Proteobacteria"
F_trp_Proteo <-auxo_Prod1[F_trp_Proteo, ]

#new column with productionrates
F_trp_Proteo$prod <- F_trp_Proteo$mtf.flux/F_trp_Proteo$Growth
F_trp_Proteo <-F_trp_Proteo[order(F_trp_Proteo$prod),]

#visualization
ggbetweenstats(
  data = F_trp_Proteo,
  x = Prototrophy,
  y = prod,
  plot.type = "violin",
  type = "nonparametric", # for wilcoxon
  centrality.plotting = TRUE,
  pairwise.display ="significant",
  title = "Comparison of formate productionrates by tryptophan auxotrophic \nand prototrophic Proteobacteria",
  ylab = "Productionrates [mmol/gDW]",
  xlab = "",
  results.subtitle = TRUE,
  conf.level = 0.95,
  caption = "0 - Auxotrophy, 1 - Prototrophy\nMann-Whitney U Test (CI = 0.95%)",
  point.args = list(
    position = ggplot2::position_jitterdodge(dodge.width = 0.20),
    alpha = 0.2,
    size = 1,
    stroke = 0)
)
p_F4 <- ggbetweenstats(
  data = F_trp_Proteo,
  x = Prototrophy,
  y = prod,
  plot.type = "violin", 
  type = "nonparametric", # for wilcoxon
  centrality.plotting = TRUE,
  pairwise.display ="significant",
  title = "Proteobacteria",
  ylab = "Productionrates [mmol/g DW]",
  xlab = "",
  results.subtitle = TRUE,
  conf.level = 0.95,
  caption = "",
  point.args = list(
    position = ggplot2::position_jitterdodge(dodge.width = 0.20),
    alpha = 0.2,
    size = 1,
    stroke = 0)
)  + theme(text = element_text(size= 10)) +
  theme(plot.subtitle = element_text(size=6)) +
  theme(plot.title = element_text(size=10))

#Actinobacteria
F_trp_Actino <-auxo_Prod1$rxn.name == "Formate-e0 Exchange" & auxo_Prod1$mtf.flux > 0 & auxo_Prod1$Compound == "Trp" & auxo_Prod1$phylum == "Actinobacteriota"
F_trp_Actino <-auxo_Prod1[F_trp_Actino, ]
View(F_trp_Actino)

#new column with productionrates
F_trp_Actino$prod <- F_trp_Actino$mtf.flux/F_trp_Actino$Growth
F_trp_Actino <-F_trp_Actino[order(F_trp_Actino$prod),]

#visualization
ggbetweenstats(
  data = F_trp_Actino,
  x = Prototrophy,
  y = prod,
  plot.type = "violin", 
  type = "nonparametric", # for wilcoxon
  centrality.plotting = TRUE,
  pairwise.display ="significant",
  title = "Comparison of the formate productionrates by tryptophan auxotrophic \nand prototrophic Actinobacteria",
  ylab = "Productionrates [mmol/gDW]",
  xlab = "",
  results.subtitle = TRUE,
  conf.level = 0.95,
  caption = "0 - Auxotrophy, 1 - Prototrophy\nMann-Whitney U Test (CI = 0.95%)",
  point.args = list(
    position = ggplot2::position_jitterdodge(dodge.width = 0.20),
    alpha = 0.2,
    size = 1,
    stroke = 0)
)
p_F5 <- ggbetweenstats(
  data = F_trp_Actino,
  x = Prototrophy,
  y = prod,
  plot.type = "violin",
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
    stroke = 0)
)  + theme(text = element_text(size= 10)) +
  theme(plot.subtitle = element_text(size=6)) +
  theme(plot.title = element_text(size=10))

#combine all plots
combine_plots(list(p_F2,p_F3, p_F4, p_F5),
              plotgrid.args = list(nrow=2),
              annotation.args = list(title = "Comparison of Formate productionrates by tryptophan auxotrophic and\nprototrophic microbiota",
                                     caption = "Legend: 0 - Auxotrophy, 1 - Prototrophy, Test: Mann-Whitney U Test (CI = 0.95%)",
                                     theme = ggplot2::theme(
                                       plot.subtitle = element_text(size = 6),
                                       axis.title.y = element_text(face = "plain"),
                                       text = element_text(size=8, face = "plain"),
                                       plot.title = element_text(size = 12, face = "bold"))))

#fisher test
#create a dataframe for lactate production
#number of genomes who are tryptophan auxotrophic but not lactate producer
#filter for tryptophan auxotrophic bacteria
numb_auxo <- Auxotrophy_2$Compound == "Trp" & Auxotrophy_2$Prototrophy == 0
numb_auxo <- Auxotrophy_2[numb_auxo, ]
View(numb_B_auxo)
#delete archaen
numb_auxo<- numb_auxo[!(numb_auxo$phylum == "Euryarchaeota" | numb_auxo$phylum == "Thermoplasmatota"), ]
View(numb_auxo)
#calculate number for contingency table
#Propionate production
nrow(F_trp_A)
nrow(F_trp_Pr)
#number trp auxotrophic and non Butyrate producer
sum_no_F_auxo <- nrow(numb_auxo) - nrow(F_trp_A)
sum_no_F_auxo
#number trp prototrophic and non butyrate producer
sum_no_F_proto <- 5416 -(nrow(F_trp_A) + nrow(F_trp_Pr) + sum_no_F_auxo)
sum_no_F_proto 

#create contigency table
fisher_P <- data.frame("Trp-Auxotrophy" = c(nrow(F_trp_A), sum_no_F_auxo), 
                       "Trp-Prototrophy"= c(nrow(F_trp_Pr),sum_no_F_proto), 
                       row.names = c("Propionate production: yes", "Propionate production: no"))

fisher_P
#control numbers
nrow(P_trp_A) + nrow(P_trp_Pr) + sum_no_P_auxo + sum_no_P_proto == 5416

#create mosaic plot
mosaicplot(fisher_P, color = TRUE, main="Formate production of trp auxotrophic and prototrophic microbiota")
#fisher test
fisher.test(fisher_P)

##############      Abundancies of amino acid auxotrophies     ################

#get the number of genomes
numbauxogenomes <- nrow(Auxotrophy)


#first visualization to give an overview
ggplot(Auxotrophy_2[Prototrophy == 0],aes(Compound, fill = phylum)) +
  geom_bar()


# get abundancies
#number of genomes with a zero (find number of genomes with an auxotrophy)
nauxoAla <- nrow(Auxotrophy_2[c(Auxotrophy_2$Compound == "Ala" & Auxotrophy_2$Prototrophy == 0)])
nauxoVal <- nrow(Auxotrophy_2[c(Auxotrophy_2$Compound == "Val" & Auxotrophy_2$Prototrophy == 0)])
nauxoMet <- nrow(Auxotrophy_2[c(Auxotrophy_2$Compound == "Met" & Auxotrophy_2$Prototrophy == 0)])
nauxoLeu <- nrow(Auxotrophy_2[c(Auxotrophy_2$Compound == "Leu" & Auxotrophy_2$Prototrophy == 0)])
nauxoIle <- nrow(Auxotrophy_2[c(Auxotrophy_2$Compound == "Ile" & Auxotrophy_2$Prototrophy == 0)])
nauxoPro <- nrow(Auxotrophy_2[c(Auxotrophy_2$Compound == "Pro" & Auxotrophy_2$Prototrophy == 0)])
nauxoTrp <- nrow(Auxotrophy_2[c(Auxotrophy_2$Compound == "Trp" & Auxotrophy_2$Prototrophy == 0)])
nauxoPhe <- nrow(Auxotrophy_2[c(Auxotrophy_2$Compound == "Phe" & Auxotrophy_2$Prototrophy == 0)])
nauxoLys <- nrow(Auxotrophy_2[c(Auxotrophy_2$Compound == "Lys" & Auxotrophy_2$Prototrophy == 0)])
nauxoArg <- nrow(Auxotrophy_2[c(Auxotrophy_2$Compound == "Arg" & Auxotrophy_2$Prototrophy == 0)])
nauxoHis <- nrow(Auxotrophy_2[c(Auxotrophy_2$Compound == "His" & Auxotrophy_2$Prototrophy == 0)])
nauxoTyr <- nrow(Auxotrophy_2[c(Auxotrophy_2$Compound == "Tyr" & Auxotrophy_2$Prototrophy == 0)])
nauxoThr <- nrow(Auxotrophy_2[c(Auxotrophy_2$Compound == "Thr" & Auxotrophy_2$Prototrophy == 0)])
nauxoGln <- nrow(Auxotrophy_2[c(Auxotrophy_2$Compound == "Gln" & Auxotrophy_2$Prototrophy == 0)])
nauxoGlu <- nrow(Auxotrophy_2[c(Auxotrophy_2$Compound == "Glu" & Auxotrophy_2$Prototrophy == 0)])
nauxoGly <- nrow(Auxotrophy_2[c(Auxotrophy_2$Compound == "Gly" & Auxotrophy_2$Prototrophy == 0)])
nauxoSer <- nrow(Auxotrophy_2[c(Auxotrophy_2$Compound == "Ser" & Auxotrophy_2$Prototrophy == 0)])
nauxoCys <- nrow(Auxotrophy_2[c(Auxotrophy_2$Compound == "Cys" & Auxotrophy_2$Prototrophy == 0)])
nauxoAsp <- nrow(Auxotrophy_2[c(Auxotrophy_2$Compound == "Asp" & Auxotrophy_2$Prototrophy == 0)])
nauxoAsn <- nrow(Auxotrophy_2[c(Auxotrophy_2$Compound == "Asn" & Auxotrophy_2$Prototrophy == 0)])
nauxoChor <- nrow(Auxotrophy_2[c(Auxotrophy_2$Compound == "Chor" & Auxotrophy_2$Prototrophy == 0)])


#create a dataframe
NumberAuxamino <- data.frame(nauxoAla, nauxoVal, nauxoMet, nauxoLeu, nauxoIle, 
                             nauxoPro, nauxoTrp, nauxoPhe, nauxoLys, nauxoArg, nauxoHis, nauxoTyr, 
                             nauxoThr, nauxoGlu, nauxoGln, 
                             nauxoGly, nauxoSer, nauxoCys, nauxoAsp, nauxoAsn, nauxoChor)
rownames(NumberAuxamino) <- c("Genomes")
NumberAuxamino

#Calculate the % of the amino acid auxotrophies from all analzyed genomes
abunamino <- NumberAuxamino[1, ] /numbauxogenomes *100
abunaminot <- t(abunamino)
abunaminos <- data.frame(abunaminot)
is.data.frame(abunaminos)
abunaminos

#add a second column with the names of amino acids
amino <- c("Ala", "Val", "Met", "Leu", "Ile", "Pro", "Trp", "Phe", "Lys", "Arg", "His", "Tyr", "Thr", "Glu", "Gln",
           "Gly", "Ser", "Cys", "Asp", "Asn", "Chor")
amino
abunaminos$aminoacids <- amino
abunaminos



#if necessary round the values f.e. when creating also a table
round(abunamino)

#adding a third column with the species/genera whatever is available
#Visualization
library(ggplot2)

#palette for colour blindness
cbPalette <- c("#560133", "#C7007C", "#9F0162", "#C7007C", "#005745", "#00306F", 
               "#00C2F9", "#004002", "#00B408", "#F60239", "#CD022D","#E69F00", 
               "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7",
               "#450270", "#FFCCFE", "#5A000F")

r <- ggplot(data=abunaminos, aes(x=aminoacids, y=Genomes)) +
  geom_bar(stat="identity") + 
  coord_cartesian(ylim=c(0,100)) +
  ylab("Auxotrophies in HRGM Genomes[%]") +
  xlab("Amino acids") +
  ggtitle("Abundance of amino acid auxotrophies in HRGM Genomes") +
  theme(panel.background = element_rect(fill="white", colour= "black")) +
  scale_fill_manual(values = cbPalette)

r

q <- ggplot(data = Auxotrophy_2, aes(x=aminoacids, y= Complet))
q

#Alanine
nauxoAlaAct <- nrow(Auxotrophy_2[c(Auxotrophy_2$Compound == "Ala" & 
                                     Auxotrophy_2$Prototrophy == 0 & Auxotrophy_2$phylum == "Actinobacteriota")]) / numbauxogenomes  *100
nauxoAlaBac <- nrow(Auxotrophy_2[c(Auxotrophy_2$Compound == "Ala" & 
                                     Auxotrophy_2$Prototrophy == 0 & Auxotrophy_2$phylum == "Bacteroidota")]) / numbauxogenomes *100
nauxoAlaBdel <- nrow(Auxotrophy_2[c(Auxotrophy_2$Compound == "Ala" & 
                                      Auxotrophy_2$Prototrophy == 0 & Auxotrophy_2$phylum == "Bdellovibrionota")]) /numbauxogenomes  *100
nauxoAlaCamp <- nrow(Auxotrophy_2[c(Auxotrophy_2$Compound == "Ala" & 
                                      Auxotrophy_2$Prototrophy == 0 & Auxotrophy_2$phylum == "Campylobacterota")]) / numbauxogenomes  *100
nauxoAlaCyan <- nrow(Auxotrophy_2[c(Auxotrophy_2$Compound == "Ala" & 
                                      Auxotrophy_2$Prototrophy == 0 & Auxotrophy_2$phylum == "Cyanobacteria")]) / numbauxogenomes  *100
nauxoAlaDesulf <- nrow(Auxotrophy_2[c(Auxotrophy_2$Compound == "Ala" & 
                                        Auxotrophy_2$Prototrophy == 0 & Auxotrophy_2$phylum == "Desulfobacterota_A")]) / numbauxogenomes  *100
nauxoAlaElusi <- nrow(Auxotrophy_2[c(Auxotrophy_2$Compound == "Ala" & 
                                       Auxotrophy_2$Prototrophy == 0 & Auxotrophy_2$phylum == "Elusimicrobiota")]) / numbauxogenomes  *100
nauxoAlaEremi <- nrow(Auxotrophy_2[c(Auxotrophy_2$Compound == "Ala" & 
                                       Auxotrophy_2$Prototrophy == 0 & Auxotrophy_2$phylum == "Eremiobacterota")]) / numbauxogenomes  *100
nauxoAlaEury <- nrow(Auxotrophy_2[c(Auxotrophy_2$Compound == "Ala" & 
                                      Auxotrophy_2$Prototrophy == 0 & Auxotrophy_2$phylum == "Euryarchaeota")]) / numbauxogenomes  *100
nauxoAlaFibro <- nrow(Auxotrophy_2[c(Auxotrophy_2$Compound == "Ala" & 
                                       Auxotrophy_2$Prototrophy == 0 & Auxotrophy_2$phylum == "Fibrobacterota")]) / numbauxogenomes  *100
nauxoAlaFirm <- nrow(Auxotrophy_2[c(Auxotrophy_2$Compound == "Ala" & 
                                      Auxotrophy_2$Prototrophy == 0 & Auxotrophy_2$phylum == "Firmicutes")]) / numbauxogenomes  *100
nauxoAlaFirmA <- nrow(Auxotrophy_2[c(Auxotrophy_2$Compound == "Ala" & 
                                       Auxotrophy_2$Prototrophy == 0 & Auxotrophy_2$phylum == "Firmicutes_A")]) / numbauxogenomes  *100
nauxoAlaFirmB <- nrow(Auxotrophy_2[c(Auxotrophy_2$Compound == "Ala" & 
                                       Auxotrophy_2$Prototrophy == 0 & Auxotrophy_2$phylum == "Firmicutes_B")]) / numbauxogenomes  *100
nauxoAlaFirmG <- nrow(Auxotrophy_2[c(Auxotrophy_2$Compound == "Ala" & 
                                       Auxotrophy_2$Prototrophy == 0 & Auxotrophy_2$phylum == "Firmicutes_G")]) / numbauxogenomes  *100
nauxoAlaFirmI <- nrow(Auxotrophy_2[c(Auxotrophy_2$Compound == "Ala" & 
                                       Auxotrophy_2$Prototrophy == 0 & Auxotrophy_2$phylum == "Firmicutes_I")]) / numbauxogenomes  *100
nauxoAlaFuso <- nrow(Auxotrophy_2[c(Auxotrophy_2$Compound == "Ala" & 
                                      Auxotrophy_2$Prototrophy == 0 & Auxotrophy_2$phylum == "Fusobacteriota")]) / numbauxogenomes  *100
nauxoAlaHalo <- nrow(Auxotrophy_2[c(Auxotrophy_2$Compound == "Ala" & 
                                      Auxotrophy_2$Prototrophy == 0 & Auxotrophy_2$phylum == "Halobacterota")]) / numbauxogenomes  *100
nauxoAlaMyxo <- nrow(Auxotrophy_2[c(Auxotrophy_2$Compound == "Ala" & 
                                      Auxotrophy_2$Prototrophy == 0 & Auxotrophy_2$phylum == "Myxococcota")]) / numbauxogenomes  *100
nauxoAlaPates <- nrow(Auxotrophy_2[c(Auxotrophy_2$Compound == "Ala" & 
                                       Auxotrophy_2$Prototrophy == 0 & Auxotrophy_2$phylum == "Patescibacteria")]) / numbauxogenomes  *100
nauxoAlaProteo <- nrow(Auxotrophy_2[c(Auxotrophy_2$Compound == "Ala" & 
                                        Auxotrophy_2$Prototrophy == 0 & Auxotrophy_2$phylum == "Proteobacteria")]) / numbauxogenomes  *100
nauxoAlaSpiro <- nrow(Auxotrophy_2[c(Auxotrophy_2$Compound == "Ala" & 
                                       Auxotrophy_2$Prototrophy == 0 & Auxotrophy_2$phylum == "Spirochaetota")]) / numbauxogenomes  *100
nauxoAlaSyner <- nrow(Auxotrophy_2[c(Auxotrophy_2$Compound == "Ala" & 
                                       Auxotrophy_2$Prototrophy == 0 & Auxotrophy_2$phylum == "Synergistota")]) / numbauxogenomes  *100
nauxoAlaTherm <- nrow(Auxotrophy_2[c(Auxotrophy_2$Compound == "Ala" & 
                                       Auxotrophy_2$Prototrophy == 0 & Auxotrophy_2$phylum == "Thermoplasmatota")]) / numbauxogenomes  *100
nauxoAlaVerru <- nrow(Auxotrophy_2[c(Auxotrophy_2$Compound == "Ala" & 
                                       Auxotrophy_2$Prototrophy == 0 & Auxotrophy_2$phylum == "Verrucomicrobiota")]) / numbauxogenomes  *100
nauxoAlaNa <- nrow(Auxotrophy_2[c(Auxotrophy_2$Compound == "Ala" & 
                                    Auxotrophy_2$Prototrophy == 0 & Auxotrophy_2$phylum == "NA")]) /numbauxogenomes  *100
#dataframe for auxotrophies (Alanine):

auxoAla <- data.frame(nauxoAlaAct, nauxoAlaBac, nauxoAlaBdel, nauxoAlaCamp, nauxoAlaCyan, nauxoAlaDesulf, nauxoAlaElusi, nauxoAlaEremi,
                      nauxoAlaEury, nauxoAlaFibro, nauxoAlaFirm, nauxoAlaFirmA, nauxoAlaFirmB, nauxoAlaFirmG, nauxoAlaFirmI, nauxoAlaFuso,
                      nauxoAlaHalo, nauxoAlaMyxo, nauxoAlaNa, nauxoAlaPates, nauxoAlaProteo, nauxoAlaSpiro, nauxoAlaSyner, nauxoAlaTherm,
                      nauxoAlaVerru) 
auxoAlaphylum <-  t(auxoAla)
auxoAlaphyla <- data.frame(auxoAlaphylum)
colnames(auxoAlaphyla) <- c("Auxo")

View(Auxotrophy_2)
#Valine
nauxoValAct <- nrow(Auxotrophy_2[c(Auxotrophy_2$Compound == "Val" & 
                                     Auxotrophy_2$Prototrophy == 0 & Auxotrophy_2$phylum == "Actinobacteriota")]) / numbauxogenomes  *100
nauxoValBac <- nrow(Auxotrophy_2[c(Auxotrophy_2$Compound == "Val" & 
                                     Auxotrophy_2$Prototrophy == 0 & Auxotrophy_2$phylum == "Bacteroidota")]) / numbauxogenomes  *100
nauxoValBdel <- nrow(Auxotrophy_2[c(Auxotrophy_2$Compound == "Val" & 
                                      Auxotrophy_2$Prototrophy == 0 & Auxotrophy_2$phylum == "Bdellovibrionota")]) / numbauxogenomes  *100
nauxoValCamp <- nrow(Auxotrophy_2[c(Auxotrophy_2$Compound == "Val" & 
                                      Auxotrophy_2$Prototrophy == 0 & Auxotrophy_2$phylum == "Campylobacterota")]) / numbauxogenomes  *100
nauxoValCyan <- nrow(Auxotrophy_2[c(Auxotrophy_2$Compound == "Val" & 
                                      Auxotrophy_2$Prototrophy == 0 & Auxotrophy_2$phylum == "Cyanobacteria")]) / numbauxogenomes  *100
nauxoValDesulf <- nrow(Auxotrophy_2[c(Auxotrophy_2$Compound == "Val" & 
                                        Auxotrophy_2$Prototrophy == 0 & Auxotrophy_2$phylum == "Desulfobacterota_A")]) / numbauxogenomes  *100
nauxoValElusi <- nrow(Auxotrophy_2[c(Auxotrophy_2$Compound == "Val" & 
                                       Auxotrophy_2$Prototrophy == 0 & Auxotrophy_2$phylum == "Elusimicrobiota")]) / numbauxogenomes  *100
nauxoValEremi <- nrow(Auxotrophy_2[c(Auxotrophy_2$Compound == "Val" & 
                                       Auxotrophy_2$Prototrophy == 0 & Auxotrophy_2$phylum == "Eremiobacterota")]) / numbauxogenomes  *100
nauxoValEury <- nrow(Auxotrophy_2[c(Auxotrophy_2$Compound == "Val" & 
                                      Auxotrophy_2$Prototrophy == 0 & Auxotrophy_2$phylum == "Euryarchaeota")]) / numbauxogenomes  *100
nauxoValFibro <- nrow(Auxotrophy_2[c(Auxotrophy_2$Compound == "Val" & 
                                       Auxotrophy_2$Prototrophy == 0 & Auxotrophy_2$phylum == "Fibrobacterota")]) / numbauxogenomes  *100
nauxoValFirm <- nrow(Auxotrophy_2[c(Auxotrophy_2$Compound == "Val" & 
                                      Auxotrophy_2$Prototrophy == 0 & Auxotrophy_2$phylum == "Firmicutes")]) / numbauxogenomes  *100
nauxoValFirmA <- nrow(Auxotrophy_2[c(Auxotrophy_2$Compound == "Val" & 
                                       Auxotrophy_2$Prototrophy == 0 & Auxotrophy_2$phylum == "Firmicutes_A")]) / numbauxogenomes  *100
nauxoValFirmB <- nrow(Auxotrophy_2[c(Auxotrophy_2$Compound == "Val" & 
                                       Auxotrophy_2$Prototrophy == 0 & Auxotrophy_2$phylum == "Firmicutes_B")]) / numbauxogenomes  *100
nauxoValFirmG <- nrow(Auxotrophy_2[c(Auxotrophy_2$Compound == "Val" & 
                                       Auxotrophy_2$Prototrophy == 0 & Auxotrophy_2$phylum == "Firmicutes_G")]) / numbauxogenomes  *100
nauxoValFirmI <- nrow(Auxotrophy_2[c(Auxotrophy_2$Compound == "Val" & 
                                       Auxotrophy_2$Prototrophy == 0 & Auxotrophy_2$phylum == "Firmicutes_I")]) / numbauxogenomes  *100
nauxoValFuso <- nrow(Auxotrophy_2[c(Auxotrophy_2$Compound == "Val" & 
                                      Auxotrophy_2$Prototrophy == 0 & Auxotrophy_2$phylum == "Fusobacteriota")]) / numbauxogenomes  *100
nauxoValHalo <- nrow(Auxotrophy_2[c(Auxotrophy_2$Compound == "Val" & 
                                      Auxotrophy_2$Prototrophy == 0 & Auxotrophy_2$phylum == "Halobacterota")]) / numbauxogenomes  *100
nauxoValMyxo <- nrow(Auxotrophy_2[c(Auxotrophy_2$Compound == "Val" & 
                                      Auxotrophy_2$Prototrophy == 0 & Auxotrophy_2$phylum == "Myxococcota")]) / numbauxogenomes  *100
nauxoValPates <- nrow(Auxotrophy_2[c(Auxotrophy_2$Compound == "Val" & 
                                       Auxotrophy_2$Prototrophy == 0 & Auxotrophy_2$phylum == "Patescibacteria")]) / numbauxogenomes  *100
nauxoValProteo <- nrow(Auxotrophy_2[c(Auxotrophy_2$Compound == "Val" & 
                                        Auxotrophy_2$Prototrophy == 0 & Auxotrophy_2$phylum == "Proteobacteria")]) / numbauxogenomes  *100
nauxoValSpiro <- nrow(Auxotrophy_2[c(Auxotrophy_2$Compound == "Val" & 
                                       Auxotrophy_2$Prototrophy == 0 & Auxotrophy_2$phylum == "Spirochaetota")]) / numbauxogenomes  *100
nauxoValSyner <- nrow(Auxotrophy_2[c(Auxotrophy_2$Compound == "Val" & 
                                       Auxotrophy_2$Prototrophy == 0 & Auxotrophy_2$phylum == "Synergistota")]) / numbauxogenomes  *100
nauxoValTherm <- nrow(Auxotrophy_2[c(Auxotrophy_2$Compound == "Val" & 
                                       Auxotrophy_2$Prototrophy == 0 & Auxotrophy_2$phylum == "Thermoplasmatota")]) / numbauxogenomes  *100
nauxoValVerru <- nrow(Auxotrophy_2[c(Auxotrophy_2$Compound == "Val" & 
                                       Auxotrophy_2$Prototrophy == 0 & Auxotrophy_2$phylum == "Verrucomicrobiota")]) / numbauxogenomes  *100
nauxoValNa <- nrow(Auxotrophy_2[c(Auxotrophy_2$Compound == "Val" & 
                                    Auxotrophy_2$Prototrophy == 0 & Auxotrophy_2$phylum == "NA")]) / numbauxogenomes  *100

#dataframe for auxotrophies Valine:
auxoVal <- data.frame(nauxoValAct, nauxoValBac, nauxoValBdel, nauxoValCamp, nauxoValCyan, nauxoValDesulf, nauxoValElusi, nauxoValEremi,
                      nauxoValEury, nauxoValFibro, nauxoValFirm, nauxoValFirmA, nauxoValFirmB, nauxoValFirmG, nauxoValFirmI, nauxoValFuso,
                      nauxoValHalo, nauxoValMyxo, nauxoValNa, nauxoValPates, nauxoValProteo, nauxoValSpiro, nauxoValSyner, nauxoValTherm,
                      nauxoValVerru) 
auxoValphylum <-  t(auxoVal)
auxoValphyla <- data.frame(auxoValphylum)
colnames(auxoValphyla) <- c("Auxo")

#Methionine
nauxoMetAct <- nrow(Auxotrophy_2[c(Auxotrophy_2$Compound == "Met" & 
                                     Auxotrophy_2$Prototrophy == 0 & Auxotrophy_2$phylum == "Actinobacteriota")]) / numbauxogenomes  *100
nauxoMetBac <- nrow(Auxotrophy_2[c(Auxotrophy_2$Compound == "Met" & 
                                     Auxotrophy_2$Prototrophy == 0 & Auxotrophy_2$phylum == "Bacteroidota")]) / numbauxogenomes  *100
nauxoMetBdel <- nrow(Auxotrophy_2[c(Auxotrophy_2$Compound == "Met" & 
                                      Auxotrophy_2$Prototrophy == 0 & Auxotrophy_2$phylum == "Bdellovibrionota")]) / numbauxogenomes  *100
nauxoMetCamp <- nrow(Auxotrophy_2[c(Auxotrophy_2$Compound == "Met" & 
                                      Auxotrophy_2$Prototrophy == 0 & Auxotrophy_2$phylum == "Campylobacterota")]) / numbauxogenomes  *100
nauxoMetCyan <- nrow(Auxotrophy_2[c(Auxotrophy_2$Compound == "Met" & 
                                      Auxotrophy_2$Prototrophy == 0 & Auxotrophy_2$phylum == "Cyanobacteria")]) / numbauxogenomes  *100
nauxoMetDesulf <- nrow(Auxotrophy_2[c(Auxotrophy_2$Compound == "Met" & 
                                        Auxotrophy_2$Prototrophy == 0 & Auxotrophy_2$phylum == "Desulfobacterota_A")]) / numbauxogenomes  *100
nauxoMetElusi <- nrow(Auxotrophy_2[c(Auxotrophy_2$Compound == "Met" & 
                                       Auxotrophy_2$Prototrophy == 0 & Auxotrophy_2$phylum == "Elusimicrobiota")]) / numbauxogenomes  *100
nauxoMetEremi <- nrow(Auxotrophy_2[c(Auxotrophy_2$Compound == "Met" & 
                                       Auxotrophy_2$Prototrophy == 0 & Auxotrophy_2$phylum == "Eremiobacterota")]) / numbauxogenomes  *100
nauxoMetEury <- nrow(Auxotrophy_2[c(Auxotrophy_2$Compound == "Met" & 
                                      Auxotrophy_2$Prototrophy == 0 & Auxotrophy_2$phylum == "Euryarchaeota")]) / numbauxogenomes  *100
nauxoMetFibro <- nrow(Auxotrophy_2[c(Auxotrophy_2$Compound == "Met" & 
                                       Auxotrophy_2$Prototrophy == 0 & Auxotrophy_2$phylum == "Fibrobacterota")]) / numbauxogenomes  *100
nauxoMetFirm <- nrow(Auxotrophy_2[c(Auxotrophy_2$Compound == "Met" & 
                                      Auxotrophy_2$Prototrophy == 0 & Auxotrophy_2$phylum == "Firmicutes")]) / numbauxogenomes  *100
nauxoMetFirmA <- nrow(Auxotrophy_2[c(Auxotrophy_2$Compound == "Met" & 
                                       Auxotrophy_2$Prototrophy == 0 & Auxotrophy_2$phylum == "Firmicutes_A")]) / numbauxogenomes  *100
nauxoMetFirmB <- nrow(Auxotrophy_2[c(Auxotrophy_2$Compound == "Met" & 
                                       Auxotrophy_2$Prototrophy == 0 & Auxotrophy_2$phylum == "Firmicutes_B")]) / numbauxogenomes  *100
nauxoMetFirmG <- nrow(Auxotrophy_2[c(Auxotrophy_2$Compound == "Met" & 
                                       Auxotrophy_2$Prototrophy == 0 & Auxotrophy_2$phylum == "Firmicutes_G")]) / numbauxogenomes  *100
nauxoMetFirmI <- nrow(Auxotrophy_2[c(Auxotrophy_2$Compound == "Met" & 
                                       Auxotrophy_2$Prototrophy == 0 & Auxotrophy_2$phylum == "Firmicutes_I")]) / numbauxogenomes  *100
nauxoMetFuso <- nrow(Auxotrophy_2[c(Auxotrophy_2$Compound == "Met" & 
                                      Auxotrophy_2$Prototrophy == 0 & Auxotrophy_2$phylum == "Fusobacteriota")]) / numbauxogenomes  *100
nauxoMetHalo <- nrow(Auxotrophy_2[c(Auxotrophy_2$Compound == "Met" & 
                                      Auxotrophy_2$Prototrophy == 0 & Auxotrophy_2$phylum == "Halobacterota")]) / numbauxogenomes  *100
nauxoMetMyxo <- nrow(Auxotrophy_2[c(Auxotrophy_2$Compound == "Met" & 
                                      Auxotrophy_2$Prototrophy == 0 & Auxotrophy_2$phylum == "Myxococcota")]) / numbauxogenomes  *100
nauxoMetPates <- nrow(Auxotrophy_2[c(Auxotrophy_2$Compound == "Met" & 
                                       Auxotrophy_2$Prototrophy == 0 & Auxotrophy_2$phylum == "Patescibacteria")]) / numbauxogenomes  *100
nauxoMetProteo <- nrow(Auxotrophy_2[c(Auxotrophy_2$Compound == "Met" & 
                                        Auxotrophy_2$Prototrophy == 0 & Auxotrophy_2$phylum == "Proteobacteria")]) / numbauxogenomes  *100
nauxoMetSpiro <- nrow(Auxotrophy_2[c(Auxotrophy_2$Compound == "Met" & 
                                       Auxotrophy_2$Prototrophy == 0 & Auxotrophy_2$phylum == "Spirochaetota")]) / numbauxogenomes  *100
nauxoMetSyner <- nrow(Auxotrophy_2[c(Auxotrophy_2$Compound == "Met" & 
                                       Auxotrophy_2$Prototrophy == 0 & Auxotrophy_2$phylum == "Synergistota")]) / numbauxogenomes  *100
nauxoMetTherm <- nrow(Auxotrophy_2[c(Auxotrophy_2$Compound == "Met" & 
                                       Auxotrophy_2$Prototrophy == 0 & Auxotrophy_2$phylum == "Thermoplasmatota")]) / numbauxogenomes  *100
nauxoMetVerru <- nrow(Auxotrophy_2[c(Auxotrophy_2$Compound == "Met" & 
                                       Auxotrophy_2$Prototrophy == 0 & Auxotrophy_2$phylum == "Verrucomicrobiota")]) / numbauxogenomes  *100
nauxoMetNa <- nrow(Auxotrophy_2[c(Auxotrophy_2$Compound == "Met" & 
                                    Auxotrophy_2$Prototrophy == 0 & Auxotrophy_2$phylum == "NA")]) / numbauxogenomes  *100

#dataframe for auxotrophies Methionine:
auxoMet <- data.frame(nauxoMetAct, nauxoMetBac, nauxoMetBdel, nauxoMetCamp, nauxoMetCyan, nauxoMetDesulf, nauxoMetElusi, nauxoMetEremi,
                      nauxoMetEury, nauxoMetFibro, nauxoMetFirm, nauxoMetFirmA, nauxoMetFirmB, nauxoMetFirmG, nauxoMetFirmI, nauxoMetFuso,
                      nauxoMetHalo, nauxoMetMyxo, nauxoMetNa, nauxoMetPates, nauxoMetProteo, nauxoMetSpiro, nauxoMetSyner, nauxoMetTherm,
                      nauxoMetVerru) 
auxoMetphylum <-  t(auxoMet)
auxoMetphyla <- data.frame(auxoMetphylum)
colnames(auxoMetphyla) <- c("Auxo")

#Leucine
nauxoLeuAct <- nrow(Auxotrophy_2[c(Auxotrophy_2$Compound == "Leu" & 
                                     Auxotrophy_2$Prototrophy == 0 & Auxotrophy_2$phylum == "Actinobacteriota")]) / numbauxogenomes  *100
nauxoLeuBac <- nrow(Auxotrophy_2[c(Auxotrophy_2$Compound == "Leu" & 
                                     Auxotrophy_2$Prototrophy == 0 & Auxotrophy_2$phylum == "Bacteroidota")]) / numbauxogenomes  *100
nauxoLeuBdel <- nrow(Auxotrophy_2[c(Auxotrophy_2$Compound == "Leu" & 
                                      Auxotrophy_2$Prototrophy == 0 & Auxotrophy_2$phylum == "Bdellovibrionota")]) / numbauxogenomes  *100
nauxoLeuCamp <- nrow(Auxotrophy_2[c(Auxotrophy_2$Compound == "Leu" & 
                                      Auxotrophy_2$Prototrophy == 0 & Auxotrophy_2$phylum == "Campylobacterota")]) / numbauxogenomes  *100
nauxoLeuCyan <- nrow(Auxotrophy_2[c(Auxotrophy_2$Compound == "Leu" & 
                                      Auxotrophy_2$Prototrophy == 0 & Auxotrophy_2$phylum == "Cyanobacteria")]) / numbauxogenomes  *100
nauxoLeuDesulf <- nrow(Auxotrophy_2[c(Auxotrophy_2$Compound == "Leu" & 
                                        Auxotrophy_2$Prototrophy == 0 & Auxotrophy_2$phylum == "Desulfobacterota_A")]) / numbauxogenomes  *100
nauxoLeuElusi <- nrow(Auxotrophy_2[c(Auxotrophy_2$Compound == "Leu" & 
                                       Auxotrophy_2$Prototrophy == 0 & Auxotrophy_2$phylum == "Elusimicrobiota")]) / numbauxogenomes  *100
nauxoLeuEremi <- nrow(Auxotrophy_2[c(Auxotrophy_2$Compound == "Leu" & 
                                       Auxotrophy_2$Prototrophy == 0 & Auxotrophy_2$phylum == "Eremiobacterota")]) / numbauxogenomes  *100
nauxoLeuEury <- nrow(Auxotrophy_2[c(Auxotrophy_2$Compound == "Leu" & 
                                      Auxotrophy_2$Prototrophy == 0 & Auxotrophy_2$phylum == "Euryarchaeota")]) / numbauxogenomes  *100
nauxoLeuFibro <- nrow(Auxotrophy_2[c(Auxotrophy_2$Compound == "Leu" & 
                                       Auxotrophy_2$Prototrophy == 0 & Auxotrophy_2$phylum == "Fibrobacterota")]) / numbauxogenomes  *100
nauxoLeuFirm <- nrow(Auxotrophy_2[c(Auxotrophy_2$Compound == "Leu" & 
                                      Auxotrophy_2$Prototrophy == 0 & Auxotrophy_2$phylum == "Firmicutes")]) / numbauxogenomes  *100
nauxoLeuFirmA <- nrow(Auxotrophy_2[c(Auxotrophy_2$Compound == "Leu" & 
                                       Auxotrophy_2$Prototrophy == 0 & Auxotrophy_2$phylum == "Firmicutes_A")]) / numbauxogenomes  *100
nauxoLeuFirmB <- nrow(Auxotrophy_2[c(Auxotrophy_2$Compound == "Leu" & 
                                       Auxotrophy_2$Prototrophy == 0 & Auxotrophy_2$phylum == "Firmicutes_B")]) / numbauxogenomes  *100
nauxoLeuFirmG <- nrow(Auxotrophy_2[c(Auxotrophy_2$Compound == "Leu" & 
                                       Auxotrophy_2$Prototrophy == 0 & Auxotrophy_2$phylum == "Firmicutes_G")]) / numbauxogenomes  *100
nauxoLeuFirmI <- nrow(Auxotrophy_2[c(Auxotrophy_2$Compound == "Leu" & 
                                       Auxotrophy_2$Prototrophy == 0 & Auxotrophy_2$phylum == "Firmicutes_I")]) / numbauxogenomes  *100
nauxoLeuFuso <- nrow(Auxotrophy_2[c(Auxotrophy_2$Compound == "Leu" & 
                                      Auxotrophy_2$Prototrophy == 0 & Auxotrophy_2$phylum == "Fusobacteriota")]) / numbauxogenomes  *100
nauxoLeuHalo <- nrow(Auxotrophy_2[c(Auxotrophy_2$Compound == "Leu" & 
                                      Auxotrophy_2$Prototrophy == 0 & Auxotrophy_2$phylum == "Halobacterota")]) / numbauxogenomes  *100
nauxoLeuMyxo <- nrow(Auxotrophy_2[c(Auxotrophy_2$Compound == "Leu" & 
                                      Auxotrophy_2$Prototrophy == 0 & Auxotrophy_2$phylum == "Myxococcota")]) / numbauxogenomes  *100
nauxoLeuPates <- nrow(Auxotrophy_2[c(Auxotrophy_2$Compound == "Leu" & 
                                       Auxotrophy_2$Prototrophy == 0 & Auxotrophy_2$phylum == "Patescibacteria")]) / numbauxogenomes  *100
nauxoLeuProteo <- nrow(Auxotrophy_2[c(Auxotrophy_2$Compound == "Leu" & 
                                        Auxotrophy_2$Prototrophy == 0 & Auxotrophy_2$phylum == "Proteobacteria")]) / numbauxogenomes  *100
nauxoLeuSpiro <- nrow(Auxotrophy_2[c(Auxotrophy_2$Compound == "Leu" & 
                                       Auxotrophy_2$Prototrophy == 0 & Auxotrophy_2$phylum == "Spirochaetota")]) / numbauxogenomes  *100
nauxoLeuSyner <- nrow(Auxotrophy_2[c(Auxotrophy_2$Compound == "Leu" & 
                                       Auxotrophy_2$Prototrophy == 0 & Auxotrophy_2$phylum == "Synergistota")]) / numbauxogenomes  *100
nauxoLeuTherm <- nrow(Auxotrophy_2[c(Auxotrophy_2$Compound == "Leu" & 
                                       Auxotrophy_2$Prototrophy == 0 & Auxotrophy_2$phylum == "Thermoplasmatota")]) / numbauxogenomes  *100
nauxoLeuVerru <- nrow(Auxotrophy_2[c(Auxotrophy_2$Compound == "Leu" & 
                                       Auxotrophy_2$Prototrophy == 0 & Auxotrophy_2$phylum == "Verrucomicrobiota")]) / numbauxogenomes  *100
nauxoLeuNa <- nrow(Auxotrophy_2[c(Auxotrophy_2$Compound == "Leu" & 
                                    Auxotrophy_2$Prototrophy == 0 & Auxotrophy_2$phylum == "NA")]) / numbauxogenomes  *100

#dataframe for auxotrophies Leucine:
auxoLeu <- data.frame(nauxoLeuAct, nauxoLeuBac, nauxoLeuBdel, nauxoLeuCamp, nauxoLeuCyan, nauxoLeuDesulf, nauxoLeuElusi, nauxoLeuEremi,
                      nauxoLeuEury, nauxoLeuFibro, nauxoLeuFirm, nauxoLeuFirmA, nauxoLeuFirmB, nauxoLeuFirmG, nauxoLeuFirmI, nauxoLeuFuso,
                      nauxoLeuHalo, nauxoLeuMyxo, nauxoLeuNa, nauxoLeuPates, nauxoLeuProteo, nauxoLeuSpiro, nauxoLeuSyner, nauxoLeuTherm,
                      nauxoLeuVerru) 
auxoLeuphylum <-  t(auxoLeu)
auxoLeuphyla <- data.frame(auxoLeuphylum)
colnames(auxoLeuphyla) <- c("Auxo")


#Isoleucine
nauxoIleAct <- nrow(Auxotrophy_2[c(Auxotrophy_2$Compound == "Ile" & 
                                     Auxotrophy_2$Prototrophy == 0 & Auxotrophy_2$phylum == "Actinobacteriota")]) / numbauxogenomes  *100
nauxoIleBac <- nrow(Auxotrophy_2[c(Auxotrophy_2$Compound == "Ile" & 
                                     Auxotrophy_2$Prototrophy == 0 & Auxotrophy_2$phylum == "Bacteroidota")]) / numbauxogenomes  *100
nauxoIleBdel <- nrow(Auxotrophy_2[c(Auxotrophy_2$Compound == "Ile" & 
                                      Auxotrophy_2$Prototrophy == 0 & Auxotrophy_2$phylum == "Bdellovibrionota")]) / numbauxogenomes  *100
nauxoIleCamp <- nrow(Auxotrophy_2[c(Auxotrophy_2$Compound == "Ile" & 
                                      Auxotrophy_2$Prototrophy == 0 & Auxotrophy_2$phylum == "Campylobacterota")]) / numbauxogenomes  *100
nauxoIleCyan <- nrow(Auxotrophy_2[c(Auxotrophy_2$Compound == "Ile" & 
                                      Auxotrophy_2$Prototrophy == 0 & Auxotrophy_2$phylum == "Cyanobacteria")]) / numbauxogenomes  *100
nauxoIleDesulf <- nrow(Auxotrophy_2[c(Auxotrophy_2$Compound == "Ile" & 
                                        Auxotrophy_2$Prototrophy == 0 & Auxotrophy_2$phylum == "Desulfobacterota_A")]) / numbauxogenomes  *100
nauxoIleElusi <- nrow(Auxotrophy_2[c(Auxotrophy_2$Compound == "Ile" & 
                                       Auxotrophy_2$Prototrophy == 0 & Auxotrophy_2$phylum == "Elusimicrobiota")]) / numbauxogenomes  *100
nauxoIleEremi <- nrow(Auxotrophy_2[c(Auxotrophy_2$Compound == "Ile" & 
                                       Auxotrophy_2$Prototrophy == 0 & Auxotrophy_2$phylum == "Eremiobacterota")]) / numbauxogenomes  *100
nauxoIleEury <- nrow(Auxotrophy_2[c(Auxotrophy_2$Compound == "Ile" & 
                                      Auxotrophy_2$Prototrophy == 0 & Auxotrophy_2$phylum == "Euryarchaeota")]) / numbauxogenomes  *100
nauxoIleFibro <- nrow(Auxotrophy_2[c(Auxotrophy_2$Compound == "Ile" & 
                                       Auxotrophy_2$Prototrophy == 0 & Auxotrophy_2$phylum == "Fibrobacterota")]) / numbauxogenomes  *100
nauxoIleFirm <- nrow(Auxotrophy_2[c(Auxotrophy_2$Compound == "Ile" & 
                                      Auxotrophy_2$Prototrophy == 0 & Auxotrophy_2$phylum == "Firmicutes")]) / numbauxogenomes  *100
nauxoIleFirmA <- nrow(Auxotrophy_2[c(Auxotrophy_2$Compound == "Ile" & 
                                       Auxotrophy_2$Prototrophy == 0 & Auxotrophy_2$phylum == "Firmicutes_A")]) / numbauxogenomes  *100
nauxoIleFirmB <- nrow(Auxotrophy_2[c(Auxotrophy_2$Compound == "Ile" & 
                                       Auxotrophy_2$Prototrophy == 0 & Auxotrophy_2$phylum == "Firmicutes_B")]) / numbauxogenomes  *100
nauxoIleFirmG <- nrow(Auxotrophy_2[c(Auxotrophy_2$Compound == "Ile" & 
                                       Auxotrophy_2$Prototrophy == 0 & Auxotrophy_2$phylum == "Firmicutes_G")]) / numbauxogenomes  *100
nauxoIleFirmI <- nrow(Auxotrophy_2[c(Auxotrophy_2$Compound == "Ile" & 
                                       Auxotrophy_2$Prototrophy == 0 & Auxotrophy_2$phylum == "Firmicutes_I")]) / numbauxogenomes  *100
nauxoIleFuso <- nrow(Auxotrophy_2[c(Auxotrophy_2$Compound == "Ile" & 
                                      Auxotrophy_2$Prototrophy == 0 & Auxotrophy_2$phylum == "Fusobacteriota")]) / numbauxogenomes  *100
nauxoIleHalo <- nrow(Auxotrophy_2[c(Auxotrophy_2$Compound == "Ile" & 
                                      Auxotrophy_2$Prototrophy == 0 & Auxotrophy_2$phylum == "Halobacterota")]) / numbauxogenomes  *100
nauxoIleMyxo <- nrow(Auxotrophy_2[c(Auxotrophy_2$Compound == "Ile" & 
                                      Auxotrophy_2$Prototrophy == 0 & Auxotrophy_2$phylum == "Myxococcota")]) / numbauxogenomes  *100
nauxoIlePates <- nrow(Auxotrophy_2[c(Auxotrophy_2$Compound == "Ile" & 
                                       Auxotrophy_2$Prototrophy == 0 & Auxotrophy_2$phylum == "Patescibacteria")]) / numbauxogenomes  *100
nauxoIleProteo <- nrow(Auxotrophy_2[c(Auxotrophy_2$Compound == "Ile" & 
                                        Auxotrophy_2$Prototrophy == 0 & Auxotrophy_2$phylum == "Proteobacteria")]) / numbauxogenomes  *100
nauxoIleSpiro <- nrow(Auxotrophy_2[c(Auxotrophy_2$Compound == "Ile" & 
                                       Auxotrophy_2$Prototrophy == 0 & Auxotrophy_2$phylum == "Spirochaetota")]) / numbauxogenomes  *100
nauxoIleSyner <- nrow(Auxotrophy_2[c(Auxotrophy_2$Compound == "Ile" & 
                                       Auxotrophy_2$Prototrophy == 0 & Auxotrophy_2$phylum == "Synergistota")]) / numbauxogenomes  *100
nauxoIleTherm <- nrow(Auxotrophy_2[c(Auxotrophy_2$Compound == "Ile" & 
                                       Auxotrophy_2$Prototrophy == 0 & Auxotrophy_2$phylum == "Thermoplasmatota")]) / numbauxogenomes  *100
nauxoIleVerru <- nrow(Auxotrophy_2[c(Auxotrophy_2$Compound == "Ile" & 
                                       Auxotrophy_2$Prototrophy == 0 & Auxotrophy_2$phylum == "Verrucomicrobiota")]) / numbauxogenomes  *100
nauxoIleNa <- nrow(Auxotrophy_2[c(Auxotrophy_2$Compound == "Ile" & 
                                    Auxotrophy_2$Prototrophy == 0 & Auxotrophy_2$phylum == "NA")]) / numbauxogenomes  *100

#dataframe for auxotrophies IleIlecine:
auxoIle <- data.frame(nauxoIleAct, nauxoIleBac, nauxoIleBdel, nauxoIleCamp, nauxoIleCyan, nauxoIleDesulf, nauxoIleElusi, nauxoIleEremi,
                      nauxoIleEury, nauxoIleFibro, nauxoIleFirm, nauxoIleFirmA, nauxoIleFirmB, nauxoIleFirmG, nauxoIleFirmI, nauxoIleFuso,
                      nauxoIleHalo, nauxoIleMyxo, nauxoIleNa, nauxoIlePates, nauxoIleProteo, nauxoIleSpiro, nauxoIleSyner, nauxoIleTherm,
                      nauxoIleVerru) 
auxoIlephylum <-  t(auxoIle)
auxoIlephyla <- data.frame(auxoIlephylum)
colnames(auxoIlephyla) <- "Auxo"

#Proline
nauxoProAct <- nrow(Auxotrophy_2[c(Auxotrophy_2$Compound == "Pro" & 
                                     Auxotrophy_2$Prototrophy == 0 & Auxotrophy_2$phylum == "Actinobacteriota")]) / numbauxogenomes  *100
nauxoProBac <- nrow(Auxotrophy_2[c(Auxotrophy_2$Compound == "Pro" & 
                                     Auxotrophy_2$Prototrophy == 0 & Auxotrophy_2$phylum == "Bacteroidota")]) / numbauxogenomes  *100
nauxoProBdel <- nrow(Auxotrophy_2[c(Auxotrophy_2$Compound == "Pro" & 
                                      Auxotrophy_2$Prototrophy == 0 & Auxotrophy_2$phylum == "Bdellovibrionota")]) / numbauxogenomes  *100
nauxoProCamp <- nrow(Auxotrophy_2[c(Auxotrophy_2$Compound == "Pro" & 
                                      Auxotrophy_2$Prototrophy == 0 & Auxotrophy_2$phylum == "Campylobacterota")]) / numbauxogenomes  *100
nauxoProCyan <- nrow(Auxotrophy_2[c(Auxotrophy_2$Compound == "Pro" & 
                                      Auxotrophy_2$Prototrophy == 0 & Auxotrophy_2$phylum == "Cyanobacteria")]) / numbauxogenomes  *100
nauxoProDesulf <- nrow(Auxotrophy_2[c(Auxotrophy_2$Compound == "Pro" & 
                                        Auxotrophy_2$Prototrophy == 0 & Auxotrophy_2$phylum == "Desulfobacterota_A")]) / numbauxogenomes  *100
nauxoProElusi <- nrow(Auxotrophy_2[c(Auxotrophy_2$Compound == "Pro" & 
                                       Auxotrophy_2$Prototrophy == 0 & Auxotrophy_2$phylum == "Elusimicrobiota")]) / numbauxogenomes  *100
nauxoProEremi <- nrow(Auxotrophy_2[c(Auxotrophy_2$Compound == "Pro" & 
                                       Auxotrophy_2$Prototrophy == 0 & Auxotrophy_2$phylum == "Eremiobacterota")]) / numbauxogenomes  *100
nauxoProEury <- nrow(Auxotrophy_2[c(Auxotrophy_2$Compound == "Pro" & 
                                      Auxotrophy_2$Prototrophy == 0 & Auxotrophy_2$phylum == "Euryarchaeota")]) / numbauxogenomes  *100
nauxoProFibro <- nrow(Auxotrophy_2[c(Auxotrophy_2$Compound == "Pro" & 
                                       Auxotrophy_2$Prototrophy == 0 & Auxotrophy_2$phylum == "Fibrobacterota")]) / numbauxogenomes  *100
nauxoProFirm <- nrow(Auxotrophy_2[c(Auxotrophy_2$Compound == "Pro" & 
                                      Auxotrophy_2$Prototrophy == 0 & Auxotrophy_2$phylum == "Firmicutes")]) / numbauxogenomes  *100
nauxoProFirmA <- nrow(Auxotrophy_2[c(Auxotrophy_2$Compound == "Pro" & 
                                       Auxotrophy_2$Prototrophy == 0 & Auxotrophy_2$phylum == "Firmicutes_A")]) / numbauxogenomes  *100
nauxoProFirmB <- nrow(Auxotrophy_2[c(Auxotrophy_2$Compound == "Pro" & 
                                       Auxotrophy_2$Prototrophy == 0 & Auxotrophy_2$phylum == "Firmicutes_B")]) / numbauxogenomes  *100
nauxoProFirmG <- nrow(Auxotrophy_2[c(Auxotrophy_2$Compound == "Pro" & 
                                       Auxotrophy_2$Prototrophy == 0 & Auxotrophy_2$phylum == "Firmicutes_G")]) / numbauxogenomes  *100
nauxoProFirmI <- nrow(Auxotrophy_2[c(Auxotrophy_2$Compound == "Pro" & 
                                       Auxotrophy_2$Prototrophy == 0 & Auxotrophy_2$phylum == "Firmicutes_I")]) / numbauxogenomes  *100
nauxoProFuso <- nrow(Auxotrophy_2[c(Auxotrophy_2$Compound == "Pro" & 
                                      Auxotrophy_2$Prototrophy == 0 & Auxotrophy_2$phylum == "Fusobacteriota")]) / numbauxogenomes  *100
nauxoProHalo <- nrow(Auxotrophy_2[c(Auxotrophy_2$Compound == "Pro" & 
                                      Auxotrophy_2$Prototrophy == 0 & Auxotrophy_2$phylum == "Halobacterota")]) / numbauxogenomes  *100
nauxoProMyxo <- nrow(Auxotrophy_2[c(Auxotrophy_2$Compound == "Pro" & 
                                      Auxotrophy_2$Prototrophy == 0 & Auxotrophy_2$phylum == "Myxococcota")]) / numbauxogenomes  *100
nauxoProPates <- nrow(Auxotrophy_2[c(Auxotrophy_2$Compound == "Pro" & 
                                       Auxotrophy_2$Prototrophy == 0 & Auxotrophy_2$phylum == "Patescibacteria")]) / numbauxogenomes  *100
nauxoProProteo <- nrow(Auxotrophy_2[c(Auxotrophy_2$Compound == "Pro" & 
                                        Auxotrophy_2$Prototrophy == 0 & Auxotrophy_2$phylum == "Proteobacteria")]) / numbauxogenomes  *100
nauxoProSpiro <- nrow(Auxotrophy_2[c(Auxotrophy_2$Compound == "Pro" & 
                                       Auxotrophy_2$Prototrophy == 0 & Auxotrophy_2$phylum == "Spirochaetota")]) / numbauxogenomes  *100
nauxoProSyner <- nrow(Auxotrophy_2[c(Auxotrophy_2$Compound == "Pro" & 
                                       Auxotrophy_2$Prototrophy == 0 & Auxotrophy_2$phylum == "Synergistota")]) / numbauxogenomes  *100
nauxoProTherm <- nrow(Auxotrophy_2[c(Auxotrophy_2$Compound == "Pro" & 
                                       Auxotrophy_2$Prototrophy == 0 & Auxotrophy_2$phylum == "Thermoplasmatota")]) / numbauxogenomes  *100
nauxoProVerru <- nrow(Auxotrophy_2[c(Auxotrophy_2$Compound == "Pro" & 
                                       Auxotrophy_2$Prototrophy == 0 & Auxotrophy_2$phylum == "Verrucomicrobiota")]) / numbauxogenomes  *100
nauxoProNa <- nrow(Auxotrophy_2[c(Auxotrophy_2$Compound == "Pro" & 
                                    Auxotrophy_2$Prototrophy == 0 & Auxotrophy_2$phylum == "NA")]) / numbauxogenomes  *100

#dataframe for auxotrophies ProProcine:
auxoPro <- data.frame(nauxoProAct, nauxoProBac, nauxoProBdel, nauxoProCamp, nauxoProCyan, nauxoProDesulf, nauxoProElusi, nauxoProEremi,
                      nauxoProEury, nauxoProFibro, nauxoProFirm, nauxoProFirmA, nauxoProFirmB, nauxoProFirmG, nauxoProFirmI, nauxoProFuso,
                      nauxoProHalo, nauxoProMyxo, nauxoProNa, nauxoProPates, nauxoProProteo, nauxoProSpiro, nauxoProSyner, nauxoProTherm,
                      nauxoProVerru) 
auxoProphylum <-  t(auxoPro)
auxoProphyla <- data.frame(auxoProphylum)
colnames(auxoProphyla) <- "Auxo"

#Tryptophan
NumberAuxamino
nauxoTrpAct <- nrow(Auxotrophy_2[c(Auxotrophy_2$Compound == "Trp" & 
                                     Auxotrophy_2$Prototrophy == 0 & Auxotrophy_2$phylum == "Actinobacteriota")]) / numbauxogenomes  *100
nauxoTrpBac <- nrow(Auxotrophy_2[c(Auxotrophy_2$Compound == "Trp" & 
                                     Auxotrophy_2$Prototrophy == 0 & Auxotrophy_2$phylum == "Bacteroidota")]) / numbauxogenomes  *100
nauxoTrpBdel <- nrow(Auxotrophy_2[c(Auxotrophy_2$Compound == "Trp" & 
                                      Auxotrophy_2$Prototrophy == 0 & Auxotrophy_2$phylum == "Bdellovibrionota")]) / numbauxogenomes  *100
nauxoTrpCamp <- nrow(Auxotrophy_2[c(Auxotrophy_2$Compound == "Trp" & 
                                      Auxotrophy_2$Prototrophy == 0 & Auxotrophy_2$phylum == "Campylobacterota")]) / numbauxogenomes  *100
nauxoTrpCyan <- nrow(Auxotrophy_2[c(Auxotrophy_2$Compound == "Trp" & 
                                      Auxotrophy_2$Prototrophy == 0 & Auxotrophy_2$phylum == "Cyanobacteria")]) / numbauxogenomes  *100
nauxoTrpDesulf <- nrow(Auxotrophy_2[c(Auxotrophy_2$Compound == "Trp" & 
                                        Auxotrophy_2$Prototrophy == 0 & Auxotrophy_2$phylum == "Desulfobacterota_A")]) / numbauxogenomes  *100
nauxoTrpElusi <- nrow(Auxotrophy_2[c(Auxotrophy_2$Compound == "Trp" & 
                                       Auxotrophy_2$Prototrophy == 0 & Auxotrophy_2$phylum == "Elusimicrobiota")]) / numbauxogenomes  *100
nauxoTrpEremi <- nrow(Auxotrophy_2[c(Auxotrophy_2$Compound == "Trp" & 
                                       Auxotrophy_2$Prototrophy == 0 & Auxotrophy_2$phylum == "Eremiobacterota")]) / numbauxogenomes  *100
nauxoTrpEury <- nrow(Auxotrophy_2[c(Auxotrophy_2$Compound == "Trp" & 
                                      Auxotrophy_2$Prototrophy == 0 & Auxotrophy_2$phylum == "Euryarchaeota")]) / numbauxogenomes  *100
nauxoTrpFibro <- nrow(Auxotrophy_2[c(Auxotrophy_2$Compound == "Trp" & 
                                       Auxotrophy_2$Prototrophy == 0 & Auxotrophy_2$phylum == "Fibrobacterota")]) / numbauxogenomes  *100
nauxoTrpFirm <- nrow(Auxotrophy_2[c(Auxotrophy_2$Compound == "Trp" & 
                                      Auxotrophy_2$Prototrophy == 0 & Auxotrophy_2$phylum == "Firmicutes")]) / numbauxogenomes  *100
nauxoTrpFirmA <- nrow(Auxotrophy_2[c(Auxotrophy_2$Compound == "Trp" & 
                                       Auxotrophy_2$Prototrophy == 0 & Auxotrophy_2$phylum == "Firmicutes_A")]) / numbauxogenomes  *100
nauxoTrpFirmB <- nrow(Auxotrophy_2[c(Auxotrophy_2$Compound == "Trp" & 
                                       Auxotrophy_2$Prototrophy == 0 & Auxotrophy_2$phylum == "Firmicutes_B")]) / numbauxogenomes  *100
nauxoTrpFirmG <- nrow(Auxotrophy_2[c(Auxotrophy_2$Compound == "Trp" & 
                                       Auxotrophy_2$Prototrophy == 0 & Auxotrophy_2$phylum == "Firmicutes_G")]) / numbauxogenomes  *100
nauxoTrpFirmI <- nrow(Auxotrophy_2[c(Auxotrophy_2$Compound == "Trp" & 
                                       Auxotrophy_2$Prototrophy == 0 & Auxotrophy_2$phylum == "Firmicutes_I")]) / numbauxogenomes  *100
nauxoTrpFuso <- nrow(Auxotrophy_2[c(Auxotrophy_2$Compound == "Trp" & 
                                      Auxotrophy_2$Prototrophy == 0 & Auxotrophy_2$phylum == "Fusobacteriota")]) / numbauxogenomes  *100
nauxoTrpHalo <- nrow(Auxotrophy_2[c(Auxotrophy_2$Compound == "Trp" & 
                                      Auxotrophy_2$Prototrophy == 0 & Auxotrophy_2$phylum == "Halobacterota")]) / numbauxogenomes  *100
nauxoTrpMyxo <- nrow(Auxotrophy_2[c(Auxotrophy_2$Compound == "Trp" & 
                                      Auxotrophy_2$Prototrophy == 0 & Auxotrophy_2$phylum == "Myxococcota")]) / numbauxogenomes  *100
nauxoTrpPates <- nrow(Auxotrophy_2[c(Auxotrophy_2$Compound == "Trp" & 
                                       Auxotrophy_2$Prototrophy == 0 & Auxotrophy_2$phylum == "Patescibacteria")]) / numbauxogenomes  *100
nauxoTrpProteo <- nrow(Auxotrophy_2[c(Auxotrophy_2$Compound == "Trp" & 
                                        Auxotrophy_2$Prototrophy == 0 & Auxotrophy_2$phylum == "Proteobacteria")]) / numbauxogenomes  *100
nauxoTrpSpiro <- nrow(Auxotrophy_2[c(Auxotrophy_2$Compound == "Trp" & 
                                       Auxotrophy_2$Prototrophy == 0 & Auxotrophy_2$phylum == "Spirochaetota")]) / numbauxogenomes  *100
nauxoTrpSyner <- nrow(Auxotrophy_2[c(Auxotrophy_2$Compound == "Trp" & 
                                       Auxotrophy_2$Prototrophy == 0 & Auxotrophy_2$phylum == "Synergistota")]) / numbauxogenomes  *100
nauxoTrpTherm <- nrow(Auxotrophy_2[c(Auxotrophy_2$Compound == "Trp" & 
                                       Auxotrophy_2$Prototrophy == 0 & Auxotrophy_2$phylum == "Thermoplasmatota")]) / numbauxogenomes  *100
nauxoTrpVerru <- nrow(Auxotrophy_2[c(Auxotrophy_2$Compound == "Trp" & 
                                       Auxotrophy_2$Prototrophy == 0 & Auxotrophy_2$phylum == "Verrucomicrobiota")]) / numbauxogenomes  *100
nauxoTrpNa <- nrow(Auxotrophy_2[c(Auxotrophy_2$Compound == "Trp" & 
                                    Auxotrophy_2$Prototrophy == 0 & Auxotrophy_2$phylum == "NA")]) / numbauxogenomes  *100

#dataframe for auxotrophies TrpTrpcine:
auxoTrp <- data.frame(nauxoTrpAct, nauxoTrpBac, nauxoTrpBdel, nauxoTrpCamp, nauxoTrpCyan, nauxoTrpDesulf, nauxoTrpElusi, nauxoTrpEremi,
                      nauxoTrpEury, nauxoTrpFibro, nauxoTrpFirm, nauxoTrpFirmA, nauxoTrpFirmB, nauxoTrpFirmG, nauxoTrpFirmI, nauxoTrpFuso,
                      nauxoTrpHalo, nauxoTrpMyxo, nauxoTrpNa, nauxoTrpPates, nauxoTrpProteo, nauxoTrpSpiro, nauxoTrpSyner, nauxoTrpTherm,
                      nauxoTrpVerru) 
auxoTrpphylum <-  t(auxoTrp)
auxoTrpphyla <- data.frame(auxoTrpphylum)
colnames(auxoTrpphyla) <- "Auxo"

#Phenylalanine
nauxoPheAct <- nrow(Auxotrophy_2[c(Auxotrophy_2$Compound == "Phe" & 
                                     Auxotrophy_2$Prototrophy == 0 & Auxotrophy_2$phylum == "Actinobacteriota")]) / numbauxogenomes  *100
nauxoPheBac <- nrow(Auxotrophy_2[c(Auxotrophy_2$Compound == "Phe" & 
                                     Auxotrophy_2$Prototrophy == 0 & Auxotrophy_2$phylum == "Bacteroidota")]) / numbauxogenomes  *100
nauxoPheBdel <- nrow(Auxotrophy_2[c(Auxotrophy_2$Compound == "Phe" & 
                                      Auxotrophy_2$Prototrophy == 0 & Auxotrophy_2$phylum == "Bdellovibrionota")]) / numbauxogenomes  *100
nauxoPheCamp <- nrow(Auxotrophy_2[c(Auxotrophy_2$Compound == "Phe" & 
                                      Auxotrophy_2$Prototrophy == 0 & Auxotrophy_2$phylum == "Campylobacterota")]) / numbauxogenomes  *100
nauxoPheCyan <- nrow(Auxotrophy_2[c(Auxotrophy_2$Compound == "Phe" & 
                                      Auxotrophy_2$Prototrophy == 0 & Auxotrophy_2$phylum == "Cyanobacteria")]) / numbauxogenomes  *100
nauxoPheDesulf <- nrow(Auxotrophy_2[c(Auxotrophy_2$Compound == "Phe" & 
                                        Auxotrophy_2$Prototrophy == 0 & Auxotrophy_2$phylum == "Desulfobacterota_A")]) / numbauxogenomes  *100
nauxoPheElusi <- nrow(Auxotrophy_2[c(Auxotrophy_2$Compound == "Phe" & 
                                       Auxotrophy_2$Prototrophy == 0 & Auxotrophy_2$phylum == "Elusimicrobiota")]) / numbauxogenomes  *100
nauxoPheEremi <- nrow(Auxotrophy_2[c(Auxotrophy_2$Compound == "Phe" & 
                                       Auxotrophy_2$Prototrophy == 0 & Auxotrophy_2$phylum == "Eremiobacterota")]) / numbauxogenomes  *100
nauxoPheEury <- nrow(Auxotrophy_2[c(Auxotrophy_2$Compound == "Phe" & 
                                      Auxotrophy_2$Prototrophy == 0 & Auxotrophy_2$phylum == "Euryarchaeota")]) / numbauxogenomes  *100
nauxoPheFibro <- nrow(Auxotrophy_2[c(Auxotrophy_2$Compound == "Phe" & 
                                       Auxotrophy_2$Prototrophy == 0 & Auxotrophy_2$phylum == "Fibrobacterota")]) / numbauxogenomes  *100
nauxoPheFirm <- nrow(Auxotrophy_2[c(Auxotrophy_2$Compound == "Phe" & 
                                      Auxotrophy_2$Prototrophy == 0 & Auxotrophy_2$phylum == "Firmicutes")]) / numbauxogenomes  *100
nauxoPheFirmA <- nrow(Auxotrophy_2[c(Auxotrophy_2$Compound == "Phe" & 
                                       Auxotrophy_2$Prototrophy == 0 & Auxotrophy_2$phylum == "Firmicutes_A")]) / numbauxogenomes  *100
nauxoPheFirmB <- nrow(Auxotrophy_2[c(Auxotrophy_2$Compound == "Phe" & 
                                       Auxotrophy_2$Prototrophy == 0 & Auxotrophy_2$phylum == "Firmicutes_B")]) / numbauxogenomes  *100
nauxoPheFirmG <- nrow(Auxotrophy_2[c(Auxotrophy_2$Compound == "Phe" & 
                                       Auxotrophy_2$Prototrophy == 0 & Auxotrophy_2$phylum == "Firmicutes_G")]) / numbauxogenomes  *100
nauxoPheFirmI <- nrow(Auxotrophy_2[c(Auxotrophy_2$Compound == "Phe" & 
                                       Auxotrophy_2$Prototrophy == 0 & Auxotrophy_2$phylum == "Firmicutes_I")]) / numbauxogenomes  *100
nauxoPheFuso <- nrow(Auxotrophy_2[c(Auxotrophy_2$Compound == "Phe" & 
                                      Auxotrophy_2$Prototrophy == 0 & Auxotrophy_2$phylum == "Fusobacteriota")]) / numbauxogenomes  *100
nauxoPheHalo <- nrow(Auxotrophy_2[c(Auxotrophy_2$Compound == "Phe" & 
                                      Auxotrophy_2$Prototrophy == 0 & Auxotrophy_2$phylum == "Halobacterota")]) / numbauxogenomes  *100
nauxoPheMyxo <- nrow(Auxotrophy_2[c(Auxotrophy_2$Compound == "Phe" & 
                                      Auxotrophy_2$Prototrophy == 0 & Auxotrophy_2$phylum == "Myxococcota")]) / numbauxogenomes  *100
nauxoPhePates <- nrow(Auxotrophy_2[c(Auxotrophy_2$Compound == "Phe" & 
                                       Auxotrophy_2$Prototrophy == 0 & Auxotrophy_2$phylum == "Patescibacteria")]) / numbauxogenomes  *100
nauxoPheProteo <- nrow(Auxotrophy_2[c(Auxotrophy_2$Compound == "Phe" & 
                                        Auxotrophy_2$Prototrophy == 0 & Auxotrophy_2$phylum == "Proteobacteria")]) / numbauxogenomes  *100
nauxoPheSpiro <- nrow(Auxotrophy_2[c(Auxotrophy_2$Compound == "Phe" & 
                                       Auxotrophy_2$Prototrophy == 0 & Auxotrophy_2$phylum == "Spirochaetota")]) / numbauxogenomes  *100
nauxoPheSyner <- nrow(Auxotrophy_2[c(Auxotrophy_2$Compound == "Phe" & 
                                       Auxotrophy_2$Prototrophy == 0 & Auxotrophy_2$phylum == "Synergistota")]) / numbauxogenomes  *100
nauxoPheTherm <- nrow(Auxotrophy_2[c(Auxotrophy_2$Compound == "Phe" & 
                                       Auxotrophy_2$Prototrophy == 0 & Auxotrophy_2$phylum == "Thermoplasmatota")]) / numbauxogenomes  *100
nauxoPheVerru <- nrow(Auxotrophy_2[c(Auxotrophy_2$Compound == "Phe" & 
                                       Auxotrophy_2$Prototrophy == 0 & Auxotrophy_2$phylum == "Verrucomicrobiota")]) / numbauxogenomes  *100
nauxoPheNa <- nrow(Auxotrophy_2[c(Auxotrophy_2$Compound == "Phe" & 
                                    Auxotrophy_2$Prototrophy == 0 & Auxotrophy_2$phylum == "NA")]) / numbauxogenomes  *100

#dataframe for auxotrophies Phenylalanine:
auxoPhe <- data.frame(nauxoPheAct, nauxoPheBac, nauxoPheBdel, nauxoPheCamp, nauxoPheCyan, nauxoPheDesulf, nauxoPheElusi, nauxoPheEremi,
                      nauxoPheEury, nauxoPheFibro, nauxoPheFirm, nauxoPheFirmA, nauxoPheFirmB, nauxoPheFirmG, nauxoPheFirmI, nauxoPheFuso,
                      nauxoPheHalo, nauxoPheMyxo, nauxoPheNa, nauxoPhePates, nauxoPheProteo, nauxoPheSpiro, nauxoPheSyner, nauxoPheTherm,
                      nauxoPheVerru) 
auxoPhephylum <-  t(auxoPhe)
auxoPhephyla <- data.frame(auxoPhephylum)
colnames(auxoPhephyla) <- "Auxo"


#Lysine
nauxoLysAct <- nrow(Auxotrophy_2[c(Auxotrophy_2$Compound == "Lys" & 
                                     Auxotrophy_2$Prototrophy == 0 & Auxotrophy_2$phylum == "Actinobacteriota")]) / numbauxogenomes  *100
nauxoLysBac <- nrow(Auxotrophy_2[c(Auxotrophy_2$Compound == "Lys" & 
                                     Auxotrophy_2$Prototrophy == 0 & Auxotrophy_2$phylum == "Bacteroidota")]) / numbauxogenomes  *100
nauxoLysBdel <- nrow(Auxotrophy_2[c(Auxotrophy_2$Compound == "Lys" & 
                                      Auxotrophy_2$Prototrophy == 0 & Auxotrophy_2$phylum == "Bdellovibrionota")]) / numbauxogenomes  *100
nauxoLysCamp <- nrow(Auxotrophy_2[c(Auxotrophy_2$Compound == "Lys" & 
                                      Auxotrophy_2$Prototrophy == 0 & Auxotrophy_2$phylum == "Campylobacterota")]) / numbauxogenomes  *100
nauxoLysCyan <- nrow(Auxotrophy_2[c(Auxotrophy_2$Compound == "Lys" & 
                                      Auxotrophy_2$Prototrophy == 0 & Auxotrophy_2$phylum == "Cyanobacteria")]) / numbauxogenomes  *100
nauxoLysDesulf <- nrow(Auxotrophy_2[c(Auxotrophy_2$Compound == "Lys" & 
                                        Auxotrophy_2$Prototrophy == 0 & Auxotrophy_2$phylum == "Desulfobacterota_A")]) / numbauxogenomes  *100
nauxoLysElusi <- nrow(Auxotrophy_2[c(Auxotrophy_2$Compound == "Lys" & 
                                       Auxotrophy_2$Prototrophy == 0 & Auxotrophy_2$phylum == "Elusimicrobiota")]) / numbauxogenomes  *100
nauxoLysEremi <- nrow(Auxotrophy_2[c(Auxotrophy_2$Compound == "Lys" & 
                                       Auxotrophy_2$Prototrophy == 0 & Auxotrophy_2$phylum == "Eremiobacterota")]) / numbauxogenomes  *100
nauxoLysEury <- nrow(Auxotrophy_2[c(Auxotrophy_2$Compound == "Lys" & 
                                      Auxotrophy_2$Prototrophy == 0 & Auxotrophy_2$phylum == "Euryarchaeota")]) / numbauxogenomes  *100
nauxoLysFibro <- nrow(Auxotrophy_2[c(Auxotrophy_2$Compound == "Lys" & 
                                       Auxotrophy_2$Prototrophy == 0 & Auxotrophy_2$phylum == "Fibrobacterota")]) / numbauxogenomes  *100
nauxoLysFirm <- nrow(Auxotrophy_2[c(Auxotrophy_2$Compound == "Lys" & 
                                      Auxotrophy_2$Prototrophy == 0 & Auxotrophy_2$phylum == "Firmicutes")]) / numbauxogenomes  *100
nauxoLysFirmA <- nrow(Auxotrophy_2[c(Auxotrophy_2$Compound == "Lys" & 
                                       Auxotrophy_2$Prototrophy == 0 & Auxotrophy_2$phylum == "Firmicutes_A")]) / numbauxogenomes  *100
nauxoLysFirmB <- nrow(Auxotrophy_2[c(Auxotrophy_2$Compound == "Lys" & 
                                       Auxotrophy_2$Prototrophy == 0 & Auxotrophy_2$phylum == "Firmicutes_B")]) / numbauxogenomes  *100
nauxoLysFirmG <- nrow(Auxotrophy_2[c(Auxotrophy_2$Compound == "Lys" & 
                                       Auxotrophy_2$Prototrophy == 0 & Auxotrophy_2$phylum == "Firmicutes_G")]) / numbauxogenomes  *100
nauxoLysFirmI <- nrow(Auxotrophy_2[c(Auxotrophy_2$Compound == "Lys" & 
                                       Auxotrophy_2$Prototrophy == 0 & Auxotrophy_2$phylum == "Firmicutes_I")]) / numbauxogenomes  *100
nauxoLysFuso <- nrow(Auxotrophy_2[c(Auxotrophy_2$Compound == "Lys" & 
                                      Auxotrophy_2$Prototrophy == 0 & Auxotrophy_2$phylum == "Fusobacteriota")]) / numbauxogenomes  *100
nauxoLysHalo <- nrow(Auxotrophy_2[c(Auxotrophy_2$Compound == "Lys" & 
                                      Auxotrophy_2$Prototrophy == 0 & Auxotrophy_2$phylum == "Halobacterota")]) / numbauxogenomes  *100
nauxoLysMyxo <- nrow(Auxotrophy_2[c(Auxotrophy_2$Compound == "Lys" & 
                                      Auxotrophy_2$Prototrophy == 0 & Auxotrophy_2$phylum == "Myxococcota")]) / numbauxogenomes  *100
nauxoLysPates <- nrow(Auxotrophy_2[c(Auxotrophy_2$Compound == "Lys" & 
                                       Auxotrophy_2$Prototrophy == 0 & Auxotrophy_2$phylum == "Patescibacteria")]) / numbauxogenomes  *100
nauxoLysProteo <- nrow(Auxotrophy_2[c(Auxotrophy_2$Compound == "Lys" & 
                                        Auxotrophy_2$Prototrophy == 0 & Auxotrophy_2$phylum == "Proteobacteria")]) / numbauxogenomes  *100
nauxoLysSpiro <- nrow(Auxotrophy_2[c(Auxotrophy_2$Compound == "Lys" & 
                                       Auxotrophy_2$Prototrophy == 0 & Auxotrophy_2$phylum == "Spirochaetota")]) / numbauxogenomes  *100
nauxoLysSyner <- nrow(Auxotrophy_2[c(Auxotrophy_2$Compound == "Lys" & 
                                       Auxotrophy_2$Prototrophy == 0 & Auxotrophy_2$phylum == "Synergistota")]) / numbauxogenomes  *100
nauxoLysTherm <- nrow(Auxotrophy_2[c(Auxotrophy_2$Compound == "Lys" & 
                                       Auxotrophy_2$Prototrophy == 0 & Auxotrophy_2$phylum == "Thermoplasmatota")]) / numbauxogenomes  *100
nauxoLysVerru <- nrow(Auxotrophy_2[c(Auxotrophy_2$Compound == "Lys" & 
                                       Auxotrophy_2$Prototrophy == 0 & Auxotrophy_2$phylum == "Verrucomicrobiota")]) / numbauxogenomes  *100
nauxoLysNa <- nrow(Auxotrophy_2[c(Auxotrophy_2$Compound == "Lys" & 
                                    Auxotrophy_2$Prototrophy == 0 & Auxotrophy_2$phylum == "NA")]) / numbauxogenomes  *100

#dataframe for auxotrophies Lysnylalanine:
auxoLys <- data.frame(nauxoLysAct, nauxoLysBac, nauxoLysBdel, nauxoLysCamp, nauxoLysCyan, nauxoLysDesulf, nauxoLysElusi, nauxoLysEremi,
                      nauxoLysEury, nauxoLysFibro, nauxoLysFirm, nauxoLysFirmA, nauxoLysFirmB, nauxoLysFirmG, nauxoLysFirmI, nauxoLysFuso,
                      nauxoLysHalo, nauxoLysMyxo, nauxoLysNa, nauxoLysPates, nauxoLysProteo, nauxoLysSpiro, nauxoLysSyner, nauxoLysTherm,
                      nauxoLysVerru) 
auxoLysphylum <-  t(auxoLys)
auxoLysphyla <- data.frame(auxoLysphylum)
colnames(auxoLysphyla) <- "Auxo"

#Arginine
nauxoArgAct <- nrow(Auxotrophy_2[c(Auxotrophy_2$Compound == "Arg" & 
                                     Auxotrophy_2$Prototrophy == 0 & Auxotrophy_2$phylum == "Actinobacteriota")]) / numbauxogenomes  *100
nauxoArgBac <- nrow(Auxotrophy_2[c(Auxotrophy_2$Compound == "Arg" & 
                                     Auxotrophy_2$Prototrophy == 0 & Auxotrophy_2$phylum == "Bacteroidota")]) / numbauxogenomes  *100
nauxoArgBdel <- nrow(Auxotrophy_2[c(Auxotrophy_2$Compound == "Arg" & 
                                      Auxotrophy_2$Prototrophy == 0 & Auxotrophy_2$phylum == "Bdellovibrionota")]) / numbauxogenomes  *100
nauxoArgCamp <- nrow(Auxotrophy_2[c(Auxotrophy_2$Compound == "Arg" & 
                                      Auxotrophy_2$Prototrophy == 0 & Auxotrophy_2$phylum == "Campylobacterota")]) / numbauxogenomes  *100
nauxoArgCyan <- nrow(Auxotrophy_2[c(Auxotrophy_2$Compound == "Arg" & 
                                      Auxotrophy_2$Prototrophy == 0 & Auxotrophy_2$phylum == "Cyanobacteria")]) / numbauxogenomes  *100
nauxoArgDesulf <- nrow(Auxotrophy_2[c(Auxotrophy_2$Compound == "Arg" & 
                                        Auxotrophy_2$Prototrophy == 0 & Auxotrophy_2$phylum == "Desulfobacterota_A")]) / numbauxogenomes  *100
nauxoArgElusi <- nrow(Auxotrophy_2[c(Auxotrophy_2$Compound == "Arg" & 
                                       Auxotrophy_2$Prototrophy == 0 & Auxotrophy_2$phylum == "Elusimicrobiota")]) / numbauxogenomes  *100
nauxoArgEremi <- nrow(Auxotrophy_2[c(Auxotrophy_2$Compound == "Arg" & 
                                       Auxotrophy_2$Prototrophy == 0 & Auxotrophy_2$phylum == "Eremiobacterota")]) / numbauxogenomes  *100
nauxoArgEury <- nrow(Auxotrophy_2[c(Auxotrophy_2$Compound == "Arg" & 
                                      Auxotrophy_2$Prototrophy == 0 & Auxotrophy_2$phylum == "Euryarchaeota")]) / numbauxogenomes  *100
nauxoArgFibro <- nrow(Auxotrophy_2[c(Auxotrophy_2$Compound == "Arg" & 
                                       Auxotrophy_2$Prototrophy == 0 & Auxotrophy_2$phylum == "Fibrobacterota")]) / numbauxogenomes  *100
nauxoArgFirm <- nrow(Auxotrophy_2[c(Auxotrophy_2$Compound == "Arg" & 
                                      Auxotrophy_2$Prototrophy == 0 & Auxotrophy_2$phylum == "Firmicutes")]) / numbauxogenomes  *100
nauxoArgFirmA <- nrow(Auxotrophy_2[c(Auxotrophy_2$Compound == "Arg" & 
                                       Auxotrophy_2$Prototrophy == 0 & Auxotrophy_2$phylum == "Firmicutes_A")]) / numbauxogenomes  *100
nauxoArgFirmB <- nrow(Auxotrophy_2[c(Auxotrophy_2$Compound == "Arg" & 
                                       Auxotrophy_2$Prototrophy == 0 & Auxotrophy_2$phylum == "Firmicutes_B")]) / numbauxogenomes  *100
nauxoArgFirmG <- nrow(Auxotrophy_2[c(Auxotrophy_2$Compound == "Arg" & 
                                       Auxotrophy_2$Prototrophy == 0 & Auxotrophy_2$phylum == "Firmicutes_G")]) / numbauxogenomes  *100
nauxoArgFirmI <- nrow(Auxotrophy_2[c(Auxotrophy_2$Compound == "Arg" & 
                                       Auxotrophy_2$Prototrophy == 0 & Auxotrophy_2$phylum == "Firmicutes_I")]) / numbauxogenomes  *100
nauxoArgFuso <- nrow(Auxotrophy_2[c(Auxotrophy_2$Compound == "Arg" & 
                                      Auxotrophy_2$Prototrophy == 0 & Auxotrophy_2$phylum == "Fusobacteriota")]) / numbauxogenomes  *100
nauxoArgHalo <- nrow(Auxotrophy_2[c(Auxotrophy_2$Compound == "Arg" & 
                                      Auxotrophy_2$Prototrophy == 0 & Auxotrophy_2$phylum == "Halobacterota")]) / numbauxogenomes  *100
nauxoArgMyxo <- nrow(Auxotrophy_2[c(Auxotrophy_2$Compound == "Arg" & 
                                      Auxotrophy_2$Prototrophy == 0 & Auxotrophy_2$phylum == "Myxococcota")]) / numbauxogenomes  *100
nauxoArgPates <- nrow(Auxotrophy_2[c(Auxotrophy_2$Compound == "Arg" & 
                                       Auxotrophy_2$Prototrophy == 0 & Auxotrophy_2$phylum == "Patescibacteria")]) / numbauxogenomes  *100
nauxoArgProteo <- nrow(Auxotrophy_2[c(Auxotrophy_2$Compound == "Arg" & 
                                        Auxotrophy_2$Prototrophy == 0 & Auxotrophy_2$phylum == "Proteobacteria")]) / numbauxogenomes  *100
nauxoArgSpiro <- nrow(Auxotrophy_2[c(Auxotrophy_2$Compound == "Arg" & 
                                       Auxotrophy_2$Prototrophy == 0 & Auxotrophy_2$phylum == "Spirochaetota")]) / numbauxogenomes  *100
nauxoArgSyner <- nrow(Auxotrophy_2[c(Auxotrophy_2$Compound == "Arg" & 
                                       Auxotrophy_2$Prototrophy == 0 & Auxotrophy_2$phylum == "Synergistota")]) / numbauxogenomes  *100
nauxoArgTherm <- nrow(Auxotrophy_2[c(Auxotrophy_2$Compound == "Arg" & 
                                       Auxotrophy_2$Prototrophy == 0 & Auxotrophy_2$phylum == "Thermoplasmatota")]) / numbauxogenomes  *100
nauxoArgVerru <- nrow(Auxotrophy_2[c(Auxotrophy_2$Compound == "Arg" & 
                                       Auxotrophy_2$Prototrophy == 0 & Auxotrophy_2$phylum == "Verrucomicrobiota")]) / numbauxogenomes  *100
nauxoArgNa <- nrow(Auxotrophy_2[c(Auxotrophy_2$Compound == "Arg" & 
                                    Auxotrophy_2$Prototrophy == 0 & Auxotrophy_2$phylum == "NA")]) / numbauxogenomes  *100

#dataframe for auxotrophies Arginine:
auxoArg <- data.frame(nauxoArgAct, nauxoArgBac, nauxoArgBdel, nauxoArgCamp, nauxoArgCyan, nauxoArgDesulf, nauxoArgElusi, nauxoArgEremi,
                      nauxoArgEury, nauxoArgFibro, nauxoArgFirm, nauxoArgFirmA, nauxoArgFirmB, nauxoArgFirmG, nauxoArgFirmI, nauxoArgFuso,
                      nauxoArgHalo, nauxoArgMyxo, nauxoArgNa, nauxoArgPates, nauxoArgProteo, nauxoArgSpiro, nauxoArgSyner, nauxoArgTherm,
                      nauxoArgVerru) 
auxoArgphylum <-  t(auxoArg)
auxoArgphyla <- data.frame(auxoArgphylum)
colnames(auxoArgphyla) <- "Auxo"

#Histidine
nauxoHisAct <- nrow(Auxotrophy_2[c(Auxotrophy_2$Compound == "His" & 
                                     Auxotrophy_2$Prototrophy == 0 & Auxotrophy_2$phylum == "Actinobacteriota")]) / numbauxogenomes  *100
nauxoHisBac <- nrow(Auxotrophy_2[c(Auxotrophy_2$Compound == "His" & 
                                     Auxotrophy_2$Prototrophy == 0 & Auxotrophy_2$phylum == "Bacteroidota")]) / numbauxogenomes  *100
nauxoHisBdel <- nrow(Auxotrophy_2[c(Auxotrophy_2$Compound == "His" & 
                                      Auxotrophy_2$Prototrophy == 0 & Auxotrophy_2$phylum == "Bdellovibrionota")]) / numbauxogenomes  *100
nauxoHisCamp <- nrow(Auxotrophy_2[c(Auxotrophy_2$Compound == "His" & 
                                      Auxotrophy_2$Prototrophy == 0 & Auxotrophy_2$phylum == "Campylobacterota")]) / numbauxogenomes  *100
nauxoHisCyan <- nrow(Auxotrophy_2[c(Auxotrophy_2$Compound == "His" & 
                                      Auxotrophy_2$Prototrophy == 0 & Auxotrophy_2$phylum == "Cyanobacteria")]) / numbauxogenomes  *100
nauxoHisDesulf <- nrow(Auxotrophy_2[c(Auxotrophy_2$Compound == "His" & 
                                        Auxotrophy_2$Prototrophy == 0 & Auxotrophy_2$phylum == "Desulfobacterota_A")]) / numbauxogenomes  *100
nauxoHisElusi <- nrow(Auxotrophy_2[c(Auxotrophy_2$Compound == "His" & 
                                       Auxotrophy_2$Prototrophy == 0 & Auxotrophy_2$phylum == "Elusimicrobiota")]) / numbauxogenomes  *100
nauxoHisEremi <- nrow(Auxotrophy_2[c(Auxotrophy_2$Compound == "His" & 
                                       Auxotrophy_2$Prototrophy == 0 & Auxotrophy_2$phylum == "Eremiobacterota")]) / numbauxogenomes  *100
nauxoHisEury <- nrow(Auxotrophy_2[c(Auxotrophy_2$Compound == "His" & 
                                      Auxotrophy_2$Prototrophy == 0 & Auxotrophy_2$phylum == "Euryarchaeota")]) / numbauxogenomes  *100
nauxoHisFibro <- nrow(Auxotrophy_2[c(Auxotrophy_2$Compound == "His" & 
                                       Auxotrophy_2$Prototrophy == 0 & Auxotrophy_2$phylum == "Fibrobacterota")]) / numbauxogenomes  *100
nauxoHisFirm <- nrow(Auxotrophy_2[c(Auxotrophy_2$Compound == "His" & 
                                      Auxotrophy_2$Prototrophy == 0 & Auxotrophy_2$phylum == "Firmicutes")]) / numbauxogenomes  *100
nauxoHisFirmA <- nrow(Auxotrophy_2[c(Auxotrophy_2$Compound == "His" & 
                                       Auxotrophy_2$Prototrophy == 0 & Auxotrophy_2$phylum == "Firmicutes_A")]) / numbauxogenomes  *100
nauxoHisFirmB <- nrow(Auxotrophy_2[c(Auxotrophy_2$Compound == "His" & 
                                       Auxotrophy_2$Prototrophy == 0 & Auxotrophy_2$phylum == "Firmicutes_B")]) / numbauxogenomes  *100
nauxoHisFirmG <- nrow(Auxotrophy_2[c(Auxotrophy_2$Compound == "His" & 
                                       Auxotrophy_2$Prototrophy == 0 & Auxotrophy_2$phylum == "Firmicutes_G")]) / numbauxogenomes  *100
nauxoHisFirmI <- nrow(Auxotrophy_2[c(Auxotrophy_2$Compound == "His" & 
                                       Auxotrophy_2$Prototrophy == 0 & Auxotrophy_2$phylum == "Firmicutes_I")]) / numbauxogenomes  *100
nauxoHisFuso <- nrow(Auxotrophy_2[c(Auxotrophy_2$Compound == "His" & 
                                      Auxotrophy_2$Prototrophy == 0 & Auxotrophy_2$phylum == "Fusobacteriota")]) / numbauxogenomes  *100
nauxoHisHalo <- nrow(Auxotrophy_2[c(Auxotrophy_2$Compound == "His" & 
                                      Auxotrophy_2$Prototrophy == 0 & Auxotrophy_2$phylum == "Halobacterota")]) / numbauxogenomes  *100
nauxoHisMyxo <- nrow(Auxotrophy_2[c(Auxotrophy_2$Compound == "His" & 
                                      Auxotrophy_2$Prototrophy == 0 & Auxotrophy_2$phylum == "Myxococcota")]) / numbauxogenomes  *100
nauxoHisPates <- nrow(Auxotrophy_2[c(Auxotrophy_2$Compound == "His" & 
                                       Auxotrophy_2$Prototrophy == 0 & Auxotrophy_2$phylum == "Patescibacteria")]) / numbauxogenomes  *100
nauxoHisProteo <- nrow(Auxotrophy_2[c(Auxotrophy_2$Compound == "His" & 
                                        Auxotrophy_2$Prototrophy == 0 & Auxotrophy_2$phylum == "Proteobacteria")]) / numbauxogenomes  *100
nauxoHisSpiro <- nrow(Auxotrophy_2[c(Auxotrophy_2$Compound == "His" & 
                                       Auxotrophy_2$Prototrophy == 0 & Auxotrophy_2$phylum == "Spirochaetota")]) / numbauxogenomes  *100
nauxoHisSyner <- nrow(Auxotrophy_2[c(Auxotrophy_2$Compound == "His" & 
                                       Auxotrophy_2$Prototrophy == 0 & Auxotrophy_2$phylum == "Synergistota")]) / numbauxogenomes  *100
nauxoHisTherm <- nrow(Auxotrophy_2[c(Auxotrophy_2$Compound == "His" & 
                                       Auxotrophy_2$Prototrophy == 0 & Auxotrophy_2$phylum == "Thermoplasmatota")]) / numbauxogenomes  *100
nauxoHisVerru <- nrow(Auxotrophy_2[c(Auxotrophy_2$Compound == "His" & 
                                       Auxotrophy_2$Prototrophy == 0 & Auxotrophy_2$phylum == "Verrucomicrobiota")]) / numbauxogenomes  *100
nauxoHisNa <- nrow(Auxotrophy_2[c(Auxotrophy_2$Compound == "His" & 
                                    Auxotrophy_2$Prototrophy == 0 & Auxotrophy_2$phylum == "NA")]) / numbauxogenomes  *100

#dataframe for auxotrophies Hisinine:
auxoHis <- data.frame(nauxoHisAct, nauxoHisBac, nauxoHisBdel, nauxoHisCamp, nauxoHisCyan, nauxoHisDesulf, nauxoHisElusi, nauxoHisEremi,
                      nauxoHisEury, nauxoHisFibro, nauxoHisFirm, nauxoHisFirmA, nauxoHisFirmB, nauxoHisFirmG, nauxoHisFirmI, nauxoHisFuso,
                      nauxoHisHalo, nauxoHisMyxo, nauxoHisNa, nauxoHisPates, nauxoHisProteo, nauxoHisSpiro, nauxoHisSyner, nauxoHisTherm,
                      nauxoHisVerru) 
auxoHisphylum <-  t(auxoHis)
auxoHisphyla <- data.frame(auxoHisphylum)
colnames(auxoHisphyla) <- "Auxo"

#Tyrosine
nauxoTyrAct <- nrow(Auxotrophy_2[c(Auxotrophy_2$Compound == "Tyr" & 
                                     Auxotrophy_2$Prototrophy == 0 & Auxotrophy_2$phylum == "Actinobacteriota")]) / numbauxogenomes  *100
nauxoTyrBac <- nrow(Auxotrophy_2[c(Auxotrophy_2$Compound == "Tyr" & 
                                     Auxotrophy_2$Prototrophy == 0 & Auxotrophy_2$phylum == "Bacteroidota")]) / numbauxogenomes  *100
nauxoTyrBdel <- nrow(Auxotrophy_2[c(Auxotrophy_2$Compound == "Tyr" & 
                                      Auxotrophy_2$Prototrophy == 0 & Auxotrophy_2$phylum == "Bdellovibrionota")]) / numbauxogenomes  *100
nauxoTyrCamp <- nrow(Auxotrophy_2[c(Auxotrophy_2$Compound == "Tyr" & 
                                      Auxotrophy_2$Prototrophy == 0 & Auxotrophy_2$phylum == "Campylobacterota")]) / numbauxogenomes  *100
nauxoTyrCyan <- nrow(Auxotrophy_2[c(Auxotrophy_2$Compound == "Tyr" & 
                                      Auxotrophy_2$Prototrophy == 0 & Auxotrophy_2$phylum == "Cyanobacteria")]) / numbauxogenomes  *100
nauxoTyrDesulf <- nrow(Auxotrophy_2[c(Auxotrophy_2$Compound == "Tyr" & 
                                        Auxotrophy_2$Prototrophy == 0 & Auxotrophy_2$phylum == "Desulfobacterota_A")]) / numbauxogenomes  *100
nauxoTyrElusi <- nrow(Auxotrophy_2[c(Auxotrophy_2$Compound == "Tyr" & 
                                       Auxotrophy_2$Prototrophy == 0 & Auxotrophy_2$phylum == "Elusimicrobiota")]) / numbauxogenomes  *100
nauxoTyrEremi <- nrow(Auxotrophy_2[c(Auxotrophy_2$Compound == "Tyr" & 
                                       Auxotrophy_2$Prototrophy == 0 & Auxotrophy_2$phylum == "Eremiobacterota")]) / numbauxogenomes  *100
nauxoTyrEury <- nrow(Auxotrophy_2[c(Auxotrophy_2$Compound == "Tyr" & 
                                      Auxotrophy_2$Prototrophy == 0 & Auxotrophy_2$phylum == "Euryarchaeota")]) / numbauxogenomes  *100
nauxoTyrFibro <- nrow(Auxotrophy_2[c(Auxotrophy_2$Compound == "Tyr" & 
                                       Auxotrophy_2$Prototrophy == 0 & Auxotrophy_2$phylum == "Fibrobacterota")]) / numbauxogenomes  *100
nauxoTyrFirm <- nrow(Auxotrophy_2[c(Auxotrophy_2$Compound == "Tyr" & 
                                      Auxotrophy_2$Prototrophy == 0 & Auxotrophy_2$phylum == "Firmicutes")]) / numbauxogenomes  *100
nauxoTyrFirmA <- nrow(Auxotrophy_2[c(Auxotrophy_2$Compound == "Tyr" & 
                                       Auxotrophy_2$Prototrophy == 0 & Auxotrophy_2$phylum == "Firmicutes_A")]) / numbauxogenomes  *100
nauxoTyrFirmB <- nrow(Auxotrophy_2[c(Auxotrophy_2$Compound == "Tyr" & 
                                       Auxotrophy_2$Prototrophy == 0 & Auxotrophy_2$phylum == "Firmicutes_B")]) / numbauxogenomes  *100
nauxoTyrFirmG <- nrow(Auxotrophy_2[c(Auxotrophy_2$Compound == "Tyr" & 
                                       Auxotrophy_2$Prototrophy == 0 & Auxotrophy_2$phylum == "Firmicutes_G")]) / numbauxogenomes  *100
nauxoTyrFirmI <- nrow(Auxotrophy_2[c(Auxotrophy_2$Compound == "Tyr" & 
                                       Auxotrophy_2$Prototrophy == 0 & Auxotrophy_2$phylum == "Firmicutes_I")]) / numbauxogenomes  *100
nauxoTyrFuso <- nrow(Auxotrophy_2[c(Auxotrophy_2$Compound == "Tyr" & 
                                      Auxotrophy_2$Prototrophy == 0 & Auxotrophy_2$phylum == "Fusobacteriota")]) / numbauxogenomes  *100
nauxoTyrHalo <- nrow(Auxotrophy_2[c(Auxotrophy_2$Compound == "Tyr" & 
                                      Auxotrophy_2$Prototrophy == 0 & Auxotrophy_2$phylum == "Halobacterota")]) / numbauxogenomes  *100
nauxoTyrMyxo <- nrow(Auxotrophy_2[c(Auxotrophy_2$Compound == "Tyr" & 
                                      Auxotrophy_2$Prototrophy == 0 & Auxotrophy_2$phylum == "Myxococcota")]) / numbauxogenomes  *100
nauxoTyrPates <- nrow(Auxotrophy_2[c(Auxotrophy_2$Compound == "Tyr" & 
                                       Auxotrophy_2$Prototrophy == 0 & Auxotrophy_2$phylum == "Patescibacteria")]) / numbauxogenomes  *100
nauxoTyrProteo <- nrow(Auxotrophy_2[c(Auxotrophy_2$Compound == "Tyr" & 
                                        Auxotrophy_2$Prototrophy == 0 & Auxotrophy_2$phylum == "Proteobacteria")]) / numbauxogenomes  *100
nauxoTyrSpiro <- nrow(Auxotrophy_2[c(Auxotrophy_2$Compound == "Tyr" & 
                                       Auxotrophy_2$Prototrophy == 0 & Auxotrophy_2$phylum == "Spirochaetota")]) / numbauxogenomes  *100
nauxoTyrSyner <- nrow(Auxotrophy_2[c(Auxotrophy_2$Compound == "Tyr" & 
                                       Auxotrophy_2$Prototrophy == 0 & Auxotrophy_2$phylum == "Synergistota")]) / numbauxogenomes  *100
nauxoTyrTherm <- nrow(Auxotrophy_2[c(Auxotrophy_2$Compound == "Tyr" & 
                                       Auxotrophy_2$Prototrophy == 0 & Auxotrophy_2$phylum == "Thermoplasmatota")]) / numbauxogenomes  *100
nauxoTyrVerru <- nrow(Auxotrophy_2[c(Auxotrophy_2$Compound == "Tyr" & 
                                       Auxotrophy_2$Prototrophy == 0 & Auxotrophy_2$phylum == "Verrucomicrobiota")]) / numbauxogenomes  *100
nauxoTyrNa <- nrow(Auxotrophy_2[c(Auxotrophy_2$Compound == "Tyr" & 
                                    Auxotrophy_2$Prototrophy == 0 & Auxotrophy_2$phylum == "NA")]) / numbauxogenomes  *100

#dataframe for auxotrophies Tyrosine:
auxoTyr <- data.frame(nauxoTyrAct, nauxoTyrBac, nauxoTyrBdel, nauxoTyrCamp, nauxoTyrCyan, nauxoTyrDesulf, nauxoTyrElusi, nauxoTyrEremi,
                      nauxoTyrEury, nauxoTyrFibro, nauxoTyrFirm, nauxoTyrFirmA, nauxoTyrFirmB, nauxoTyrFirmG, nauxoTyrFirmI, nauxoTyrFuso,
                      nauxoTyrHalo, nauxoTyrMyxo, nauxoTyrNa, nauxoTyrPates, nauxoTyrProteo, nauxoTyrSpiro, nauxoTyrSyner, nauxoTyrTherm,
                      nauxoTyrVerru) 
auxoTyrphylum <-  t(auxoTyr)
auxoTyrphyla <- data.frame(auxoTyrphylum)
colnames(auxoTyrphyla) <- "Auxo"

#Threonine
nauxoThrAct <- nrow(Auxotrophy_2[c(Auxotrophy_2$Compound == "Thr" & 
                                     Auxotrophy_2$Prototrophy == 0 & Auxotrophy_2$phylum == "Actinobacteriota")]) / numbauxogenomes  *100
nauxoThrBac <- nrow(Auxotrophy_2[c(Auxotrophy_2$Compound == "Thr" & 
                                     Auxotrophy_2$Prototrophy == 0 & Auxotrophy_2$phylum == "Bacteroidota")]) / numbauxogenomes  *100
nauxoThrBdel <- nrow(Auxotrophy_2[c(Auxotrophy_2$Compound == "Thr" & 
                                      Auxotrophy_2$Prototrophy == 0 & Auxotrophy_2$phylum == "Bdellovibrionota")]) / numbauxogenomes  *100
nauxoThrCamp <- nrow(Auxotrophy_2[c(Auxotrophy_2$Compound == "Thr" & 
                                      Auxotrophy_2$Prototrophy == 0 & Auxotrophy_2$phylum == "Campylobacterota")]) / numbauxogenomes  *100
nauxoThrCyan <- nrow(Auxotrophy_2[c(Auxotrophy_2$Compound == "Thr" & 
                                      Auxotrophy_2$Prototrophy == 0 & Auxotrophy_2$phylum == "Cyanobacteria")]) / numbauxogenomes  *100
nauxoThrDesulf <- nrow(Auxotrophy_2[c(Auxotrophy_2$Compound == "Thr" & 
                                        Auxotrophy_2$Prototrophy == 0 & Auxotrophy_2$phylum == "Desulfobacterota_A")]) / numbauxogenomes  *100
nauxoThrElusi <- nrow(Auxotrophy_2[c(Auxotrophy_2$Compound == "Thr" & 
                                       Auxotrophy_2$Prototrophy == 0 & Auxotrophy_2$phylum == "Elusimicrobiota")]) / numbauxogenomes  *100
nauxoThrEremi <- nrow(Auxotrophy_2[c(Auxotrophy_2$Compound == "Thr" & 
                                       Auxotrophy_2$Prototrophy == 0 & Auxotrophy_2$phylum == "Eremiobacterota")]) / numbauxogenomes  *100
nauxoThrEury <- nrow(Auxotrophy_2[c(Auxotrophy_2$Compound == "Thr" & 
                                      Auxotrophy_2$Prototrophy == 0 & Auxotrophy_2$phylum == "Euryarchaeota")]) / numbauxogenomes  *100
nauxoThrFibro <- nrow(Auxotrophy_2[c(Auxotrophy_2$Compound == "Thr" & 
                                       Auxotrophy_2$Prototrophy == 0 & Auxotrophy_2$phylum == "Fibrobacterota")]) / numbauxogenomes  *100
nauxoThrFirm <- nrow(Auxotrophy_2[c(Auxotrophy_2$Compound == "Thr" & 
                                      Auxotrophy_2$Prototrophy == 0 & Auxotrophy_2$phylum == "Firmicutes")]) / numbauxogenomes  *100
nauxoThrFirmA <- nrow(Auxotrophy_2[c(Auxotrophy_2$Compound == "Thr" & 
                                       Auxotrophy_2$Prototrophy == 0 & Auxotrophy_2$phylum == "Firmicutes_A")]) / numbauxogenomes  *100
nauxoThrFirmB <- nrow(Auxotrophy_2[c(Auxotrophy_2$Compound == "Thr" & 
                                       Auxotrophy_2$Prototrophy == 0 & Auxotrophy_2$phylum == "Firmicutes_B")]) / numbauxogenomes  *100
nauxoThrFirmG <- nrow(Auxotrophy_2[c(Auxotrophy_2$Compound == "Thr" & 
                                       Auxotrophy_2$Prototrophy == 0 & Auxotrophy_2$phylum == "Firmicutes_G")]) / numbauxogenomes  *100
nauxoThrFirmI <- nrow(Auxotrophy_2[c(Auxotrophy_2$Compound == "Thr" & 
                                       Auxotrophy_2$Prototrophy == 0 & Auxotrophy_2$phylum == "Firmicutes_I")]) / numbauxogenomes  *100
nauxoThrFuso <- nrow(Auxotrophy_2[c(Auxotrophy_2$Compound == "Thr" & 
                                      Auxotrophy_2$Prototrophy == 0 & Auxotrophy_2$phylum == "Fusobacteriota")]) / numbauxogenomes  *100
nauxoThrHalo <- nrow(Auxotrophy_2[c(Auxotrophy_2$Compound == "Thr" & 
                                      Auxotrophy_2$Prototrophy == 0 & Auxotrophy_2$phylum == "Halobacterota")]) / numbauxogenomes  *100
nauxoThrMyxo <- nrow(Auxotrophy_2[c(Auxotrophy_2$Compound == "Thr" & 
                                      Auxotrophy_2$Prototrophy == 0 & Auxotrophy_2$phylum == "Myxococcota")]) / numbauxogenomes  *100
nauxoThrPates <- nrow(Auxotrophy_2[c(Auxotrophy_2$Compound == "Thr" & 
                                       Auxotrophy_2$Prototrophy == 0 & Auxotrophy_2$phylum == "Patescibacteria")]) / numbauxogenomes  *100
nauxoThrProteo <- nrow(Auxotrophy_2[c(Auxotrophy_2$Compound == "Thr" & 
                                        Auxotrophy_2$Prototrophy == 0 & Auxotrophy_2$phylum == "Proteobacteria")]) / numbauxogenomes  *100
nauxoThrSpiro <- nrow(Auxotrophy_2[c(Auxotrophy_2$Compound == "Thr" & 
                                       Auxotrophy_2$Prototrophy == 0 & Auxotrophy_2$phylum == "Spirochaetota")]) / numbauxogenomes  *100
nauxoThrSyner <- nrow(Auxotrophy_2[c(Auxotrophy_2$Compound == "Thr" & 
                                       Auxotrophy_2$Prototrophy == 0 & Auxotrophy_2$phylum == "Synergistota")]) / numbauxogenomes  *100
nauxoThrTherm <- nrow(Auxotrophy_2[c(Auxotrophy_2$Compound == "Thr" & 
                                       Auxotrophy_2$Prototrophy == 0 & Auxotrophy_2$phylum == "Thermoplasmatota")]) / numbauxogenomes  *100
nauxoThrVerru <- nrow(Auxotrophy_2[c(Auxotrophy_2$Compound == "Thr" & 
                                       Auxotrophy_2$Prototrophy == 0 & Auxotrophy_2$phylum == "Verrucomicrobiota")]) / numbauxogenomes  *100
nauxoThrNa <- nrow(Auxotrophy_2[c(Auxotrophy_2$Compound == "Thr" & 
                                    Auxotrophy_2$Prototrophy == 0 & Auxotrophy_2$phylum == "NA")]) / numbauxogenomes  *100

#dataframe for auxotrophies Threonine:
auxoThr <- data.frame(nauxoThrAct, nauxoThrBac, nauxoThrBdel, nauxoThrCamp, nauxoThrCyan, nauxoThrDesulf, nauxoThrElusi, nauxoThrEremi,
                      nauxoThrEury, nauxoThrFibro, nauxoThrFirm, nauxoThrFirmA, nauxoThrFirmB, nauxoThrFirmG, nauxoThrFirmI, nauxoThrFuso,
                      nauxoThrHalo, nauxoThrMyxo, nauxoThrNa, nauxoThrPates, nauxoThrProteo, nauxoThrSpiro, nauxoThrSyner, nauxoThrTherm,
                      nauxoThrVerru) 
auxoThrphylum <-  t(auxoThr)
auxoThrphyla <- data.frame(auxoThrphylum)
colnames(auxoThrphyla) <- "Auxo"

#Glutamine (Gln)
nauxoGlnAct <- nrow(Auxotrophy_2[c(Auxotrophy_2$Compound == "Gln" & 
                                     Auxotrophy_2$Prototrophy == 0 & Auxotrophy_2$phylum == "Actinobacteriota")]) / numbauxogenomes  *100
nauxoGlnBac <- nrow(Auxotrophy_2[c(Auxotrophy_2$Compound == "Gln" & 
                                     Auxotrophy_2$Prototrophy == 0 & Auxotrophy_2$phylum == "Bacteroidota")]) / numbauxogenomes  *100
nauxoGlnBdel <- nrow(Auxotrophy_2[c(Auxotrophy_2$Compound == "Gln" & 
                                      Auxotrophy_2$Prototrophy == 0 & Auxotrophy_2$phylum == "Bdellovibrionota")]) / numbauxogenomes  *100
nauxoGlnCamp <- nrow(Auxotrophy_2[c(Auxotrophy_2$Compound == "Gln" & 
                                      Auxotrophy_2$Prototrophy == 0 & Auxotrophy_2$phylum == "Campylobacterota")]) / numbauxogenomes  *100
nauxoGlnCyan <- nrow(Auxotrophy_2[c(Auxotrophy_2$Compound == "Gln" & 
                                      Auxotrophy_2$Prototrophy == 0 & Auxotrophy_2$phylum == "Cyanobacteria")]) / numbauxogenomes  *100
nauxoGlnDesulf <- nrow(Auxotrophy_2[c(Auxotrophy_2$Compound == "Gln" & 
                                        Auxotrophy_2$Prototrophy == 0 & Auxotrophy_2$phylum == "Desulfobacterota_A")]) / numbauxogenomes  *100
nauxoGlnElusi <- nrow(Auxotrophy_2[c(Auxotrophy_2$Compound == "Gln" & 
                                       Auxotrophy_2$Prototrophy == 0 & Auxotrophy_2$phylum == "Elusimicrobiota")]) / numbauxogenomes  *100
nauxoGlnEremi <- nrow(Auxotrophy_2[c(Auxotrophy_2$Compound == "Gln" & 
                                       Auxotrophy_2$Prototrophy == 0 & Auxotrophy_2$phylum == "Eremiobacterota")]) / numbauxogenomes  *100
nauxoGlnEury <- nrow(Auxotrophy_2[c(Auxotrophy_2$Compound == "Gln" & 
                                      Auxotrophy_2$Prototrophy == 0 & Auxotrophy_2$phylum == "Euryarchaeota")]) / numbauxogenomes  *100
nauxoGlnFibro <- nrow(Auxotrophy_2[c(Auxotrophy_2$Compound == "Gln" & 
                                       Auxotrophy_2$Prototrophy == 0 & Auxotrophy_2$phylum == "Fibrobacterota")]) / numbauxogenomes  *100
nauxoGlnFirm <- nrow(Auxotrophy_2[c(Auxotrophy_2$Compound == "Gln" & 
                                      Auxotrophy_2$Prototrophy == 0 & Auxotrophy_2$phylum == "Firmicutes")]) / numbauxogenomes  *100
nauxoGlnFirmA <- nrow(Auxotrophy_2[c(Auxotrophy_2$Compound == "Gln" & 
                                       Auxotrophy_2$Prototrophy == 0 & Auxotrophy_2$phylum == "Firmicutes_A")]) / numbauxogenomes  *100
nauxoGlnFirmB <- nrow(Auxotrophy_2[c(Auxotrophy_2$Compound == "Gln" & 
                                       Auxotrophy_2$Prototrophy == 0 & Auxotrophy_2$phylum == "Firmicutes_B")]) / numbauxogenomes  *100
nauxoGlnFirmG <- nrow(Auxotrophy_2[c(Auxotrophy_2$Compound == "Gln" & 
                                       Auxotrophy_2$Prototrophy == 0 & Auxotrophy_2$phylum == "Firmicutes_G")]) / numbauxogenomes  *100
nauxoGlnFirmI <- nrow(Auxotrophy_2[c(Auxotrophy_2$Compound == "Gln" & 
                                       Auxotrophy_2$Prototrophy == 0 & Auxotrophy_2$phylum == "Firmicutes_I")]) / numbauxogenomes  *100
nauxoGlnFuso <- nrow(Auxotrophy_2[c(Auxotrophy_2$Compound == "Gln" & 
                                      Auxotrophy_2$Prototrophy == 0 & Auxotrophy_2$phylum == "Fusobacteriota")]) / numbauxogenomes  *100
nauxoGlnHalo <- nrow(Auxotrophy_2[c(Auxotrophy_2$Compound == "Gln" & 
                                      Auxotrophy_2$Prototrophy == 0 & Auxotrophy_2$phylum == "Halobacterota")]) / numbauxogenomes  *100
nauxoGlnMyxo <- nrow(Auxotrophy_2[c(Auxotrophy_2$Compound == "Gln" & 
                                      Auxotrophy_2$Prototrophy == 0 & Auxotrophy_2$phylum == "Myxococcota")]) / numbauxogenomes  *100
nauxoGlnPates <- nrow(Auxotrophy_2[c(Auxotrophy_2$Compound == "Gln" & 
                                       Auxotrophy_2$Prototrophy == 0 & Auxotrophy_2$phylum == "Patescibacteria")]) / numbauxogenomes  *100
nauxoGlnProteo <- nrow(Auxotrophy_2[c(Auxotrophy_2$Compound == "Gln" & 
                                        Auxotrophy_2$Prototrophy == 0 & Auxotrophy_2$phylum == "Proteobacteria")]) / numbauxogenomes  *100
nauxoGlnSpiro <- nrow(Auxotrophy_2[c(Auxotrophy_2$Compound == "Gln" & 
                                       Auxotrophy_2$Prototrophy == 0 & Auxotrophy_2$phylum == "Spirochaetota")]) / numbauxogenomes  *100
nauxoGlnSyner <- nrow(Auxotrophy_2[c(Auxotrophy_2$Compound == "Gln" & 
                                       Auxotrophy_2$Prototrophy == 0 & Auxotrophy_2$phylum == "Synergistota")]) / numbauxogenomes  *100
nauxoGlnTherm <- nrow(Auxotrophy_2[c(Auxotrophy_2$Compound == "Gln" & 
                                       Auxotrophy_2$Prototrophy == 0 & Auxotrophy_2$phylum == "Thermoplasmatota")]) / numbauxogenomes  *100
nauxoGlnVerru <- nrow(Auxotrophy_2[c(Auxotrophy_2$Compound == "Gln" & 
                                       Auxotrophy_2$Prototrophy == 0 & Auxotrophy_2$phylum == "Verrucomicrobiota")]) / numbauxogenomes  *100
nauxoGlnNa <- nrow(Auxotrophy_2[c(Auxotrophy_2$Compound == "Gln" & 
                                    Auxotrophy_2$Prototrophy == 0 & Auxotrophy_2$phylum == "NA")]) / numbauxogenomes  *100

#dataframe for auxotrophies Glutamine:
auxoGln <- data.frame(nauxoGlnAct, nauxoGlnBac, nauxoGlnBdel, nauxoGlnCamp, nauxoGlnCyan, nauxoGlnDesulf, nauxoGlnElusi, nauxoGlnEremi,
                      nauxoGlnEury, nauxoGlnFibro, nauxoGlnFirm, nauxoGlnFirmA, nauxoGlnFirmB, nauxoGlnFirmG, nauxoGlnFirmI, nauxoGlnFuso,
                      nauxoGlnHalo, nauxoGlnMyxo, nauxoGlnNa, nauxoGlnPates, nauxoGlnProteo, nauxoGlnSpiro, nauxoGlnSyner, nauxoGlnTherm,
                      nauxoGlnVerru) 
auxoGlnphylum <-  t(auxoGln)
auxoGlnphyla <- data.frame(auxoGlnphylum)
colnames(auxoGlnphyla) <- "Auxo"

#Glutamate
nauxoGluAct <- nrow(Auxotrophy_2[c(Auxotrophy_2$Compound == "Glu" & 
                                     Auxotrophy_2$Prototrophy == 0 & Auxotrophy_2$phylum == "Actinobacteriota")]) / numbauxogenomes  *100
nauxoGluBac <- nrow(Auxotrophy_2[c(Auxotrophy_2$Compound == "Glu" & 
                                     Auxotrophy_2$Prototrophy == 0 & Auxotrophy_2$phylum == "Bacteroidota")]) / numbauxogenomes  *100
nauxoGluBdel <- nrow(Auxotrophy_2[c(Auxotrophy_2$Compound == "Glu" & 
                                      Auxotrophy_2$Prototrophy == 0 & Auxotrophy_2$phylum == "Bdellovibrionota")]) / numbauxogenomes  *100
nauxoGluCamp <- nrow(Auxotrophy_2[c(Auxotrophy_2$Compound == "Glu" & 
                                      Auxotrophy_2$Prototrophy == 0 & Auxotrophy_2$phylum == "Campylobacterota")]) / numbauxogenomes  *100
nauxoGluCyan <- nrow(Auxotrophy_2[c(Auxotrophy_2$Compound == "Glu" & 
                                      Auxotrophy_2$Prototrophy == 0 & Auxotrophy_2$phylum == "Cyanobacteria")]) / numbauxogenomes  *100
nauxoGluDesulf <- nrow(Auxotrophy_2[c(Auxotrophy_2$Compound == "Glu" & 
                                        Auxotrophy_2$Prototrophy == 0 & Auxotrophy_2$phylum == "Desulfobacterota_A")]) / numbauxogenomes  *100
nauxoGluElusi <- nrow(Auxotrophy_2[c(Auxotrophy_2$Compound == "Glu" & 
                                       Auxotrophy_2$Prototrophy == 0 & Auxotrophy_2$phylum == "Elusimicrobiota")]) / numbauxogenomes  *100
nauxoGluEremi <- nrow(Auxotrophy_2[c(Auxotrophy_2$Compound == "Glu" & 
                                       Auxotrophy_2$Prototrophy == 0 & Auxotrophy_2$phylum == "Eremiobacterota")]) / numbauxogenomes  *100
nauxoGluEury <- nrow(Auxotrophy_2[c(Auxotrophy_2$Compound == "Glu" & 
                                      Auxotrophy_2$Prototrophy == 0 & Auxotrophy_2$phylum == "Euryarchaeota")]) / numbauxogenomes  *100
nauxoGluFibro <- nrow(Auxotrophy_2[c(Auxotrophy_2$Compound == "Glu" & 
                                       Auxotrophy_2$Prototrophy == 0 & Auxotrophy_2$phylum == "Fibrobacterota")]) / numbauxogenomes  *100
nauxoGluFirm <- nrow(Auxotrophy_2[c(Auxotrophy_2$Compound == "Glu" & 
                                      Auxotrophy_2$Prototrophy == 0 & Auxotrophy_2$phylum == "Firmicutes")]) / numbauxogenomes  *100
nauxoGluFirmA <- nrow(Auxotrophy_2[c(Auxotrophy_2$Compound == "Glu" & 
                                       Auxotrophy_2$Prototrophy == 0 & Auxotrophy_2$phylum == "Firmicutes_A")]) / numbauxogenomes  *100
nauxoGluFirmB <- nrow(Auxotrophy_2[c(Auxotrophy_2$Compound == "Glu" & 
                                       Auxotrophy_2$Prototrophy == 0 & Auxotrophy_2$phylum == "Firmicutes_B")]) / numbauxogenomes  *100
nauxoGluFirmG <- nrow(Auxotrophy_2[c(Auxotrophy_2$Compound == "Glu" & 
                                       Auxotrophy_2$Prototrophy == 0 & Auxotrophy_2$phylum == "Firmicutes_G")]) / numbauxogenomes  *100
nauxoGluFirmI <- nrow(Auxotrophy_2[c(Auxotrophy_2$Compound == "Glu" & 
                                       Auxotrophy_2$Prototrophy == 0 & Auxotrophy_2$phylum == "Firmicutes_I")]) / numbauxogenomes  *100
nauxoGluFuso <- nrow(Auxotrophy_2[c(Auxotrophy_2$Compound == "Glu" & 
                                      Auxotrophy_2$Prototrophy == 0 & Auxotrophy_2$phylum == "Fusobacteriota")]) / numbauxogenomes  *100
nauxoGluHalo <- nrow(Auxotrophy_2[c(Auxotrophy_2$Compound == "Glu" & 
                                      Auxotrophy_2$Prototrophy == 0 & Auxotrophy_2$phylum == "Halobacterota")]) / numbauxogenomes  *100
nauxoGluMyxo <- nrow(Auxotrophy_2[c(Auxotrophy_2$Compound == "Glu" & 
                                      Auxotrophy_2$Prototrophy == 0 & Auxotrophy_2$phylum == "Myxococcota")]) / numbauxogenomes  *100
nauxoGluPates <- nrow(Auxotrophy_2[c(Auxotrophy_2$Compound == "Glu" & 
                                       Auxotrophy_2$Prototrophy == 0 & Auxotrophy_2$phylum == "Patescibacteria")]) / numbauxogenomes  *100
nauxoGluProteo <- nrow(Auxotrophy_2[c(Auxotrophy_2$Compound == "Glu" & 
                                        Auxotrophy_2$Prototrophy == 0 & Auxotrophy_2$phylum == "Proteobacteria")]) / numbauxogenomes  *100
nauxoGluSpiro <- nrow(Auxotrophy_2[c(Auxotrophy_2$Compound == "Glu" & 
                                       Auxotrophy_2$Prototrophy == 0 & Auxotrophy_2$phylum == "Spirochaetota")]) / numbauxogenomes  *100
nauxoGluSyner <- nrow(Auxotrophy_2[c(Auxotrophy_2$Compound == "Glu" & 
                                       Auxotrophy_2$Prototrophy == 0 & Auxotrophy_2$phylum == "Synergistota")]) / numbauxogenomes  *100
nauxoGluTherm <- nrow(Auxotrophy_2[c(Auxotrophy_2$Compound == "Glu" & 
                                       Auxotrophy_2$Prototrophy == 0 & Auxotrophy_2$phylum == "Thermoplasmatota")]) / numbauxogenomes  *100
nauxoGluVerru <- nrow(Auxotrophy_2[c(Auxotrophy_2$Compound == "Glu" & 
                                       Auxotrophy_2$Prototrophy == 0 & Auxotrophy_2$phylum == "Verrucomicrobiota")]) / numbauxogenomes  *100
nauxoGluNa <- nrow(Auxotrophy_2[c(Auxotrophy_2$Compound == "Glu" & 
                                    Auxotrophy_2$Prototrophy == 0 & Auxotrophy_2$phylum == "NA")]) / numbauxogenomes  *100

#dataframe for auxotrophies Glutamine:
auxoGlu <- data.frame(nauxoGluAct, nauxoGluBac, nauxoGluBdel, nauxoGluCamp, nauxoGluCyan, nauxoGluDesulf, nauxoGluElusi, nauxoGluEremi,
                      nauxoGluEury, nauxoGluFibro, nauxoGluFirm, nauxoGluFirmA, nauxoGluFirmB, nauxoGluFirmG, nauxoGluFirmI, nauxoGluFuso,
                      nauxoGluHalo, nauxoGluMyxo, nauxoGluNa, nauxoGluPates, nauxoGluProteo, nauxoGluSpiro, nauxoGluSyner, nauxoGluTherm,
                      nauxoGluVerru) 
auxoGluphylum <-  t(auxoGlu)
auxoGluphyla <- data.frame(auxoGluphylum)
colnames(auxoGluphyla) <- "Auxo"
#Glycine
nauxoGlyAct <- nrow(Auxotrophy_2[c(Auxotrophy_2$Compound == "Gly" & 
                                     Auxotrophy_2$Prototrophy == 0 & Auxotrophy_2$phylum == "Actinobacteriota")]) / numbauxogenomes *100
nauxoGlyBac <- nrow(Auxotrophy_2[c(Auxotrophy_2$Compound == "Gly" & 
                                     Auxotrophy_2$Prototrophy == 0 & Auxotrophy_2$phylum == "Bacteroidota")]) / numbauxogenomes *100
nauxoGlyBdel <- nrow(Auxotrophy_2[c(Auxotrophy_2$Compound == "Gly" & 
                                      Auxotrophy_2$Prototrophy == 0 & Auxotrophy_2$phylum == "Bdellovibrionota")]) / numbauxogenomes *100
nauxoGlyCamp <- nrow(Auxotrophy_2[c(Auxotrophy_2$Compound == "Gly" & 
                                      Auxotrophy_2$Prototrophy == 0 & Auxotrophy_2$phylum == "Campylobacterota")]) / numbauxogenomes *100
nauxoGlyCyan <- nrow(Auxotrophy_2[c(Auxotrophy_2$Compound == "Gly" & 
                                      Auxotrophy_2$Prototrophy == 0 & Auxotrophy_2$phylum == "Cyanobacteria")]) / numbauxogenomes *100
nauxoGlyDesulf <- nrow(Auxotrophy_2[c(Auxotrophy_2$Compound == "Gly" & 
                                        Auxotrophy_2$Prototrophy == 0 & Auxotrophy_2$phylum == "Desulfobacterota_A")]) / numbauxogenomes *100
nauxoGlyElusi <- nrow(Auxotrophy_2[c(Auxotrophy_2$Compound == "Gly" & 
                                       Auxotrophy_2$Prototrophy == 0 & Auxotrophy_2$phylum == "Elusimicrobiota")]) / numbauxogenomes *100
nauxoGlyEremi <- nrow(Auxotrophy_2[c(Auxotrophy_2$Compound == "Gly" & 
                                       Auxotrophy_2$Prototrophy == 0 & Auxotrophy_2$phylum == "Eremiobacterota")]) / numbauxogenomes *100
nauxoGlyEury <- nrow(Auxotrophy_2[c(Auxotrophy_2$Compound == "Gly" & 
                                      Auxotrophy_2$Prototrophy == 0 & Auxotrophy_2$phylum == "Euryarchaeota")]) / numbauxogenomes *100
nauxoGlyFibro <- nrow(Auxotrophy_2[c(Auxotrophy_2$Compound == "Gly" & 
                                       Auxotrophy_2$Prototrophy == 0 & Auxotrophy_2$phylum == "Fibrobacterota")]) / numbauxogenomes *100
nauxoGlyFirm <- nrow(Auxotrophy_2[c(Auxotrophy_2$Compound == "Gly" & 
                                      Auxotrophy_2$Prototrophy == 0 & Auxotrophy_2$phylum == "Firmicutes")]) / numbauxogenomes *100
nauxoGlyFirmA <- nrow(Auxotrophy_2[c(Auxotrophy_2$Compound == "Gly" & 
                                       Auxotrophy_2$Prototrophy == 0 & Auxotrophy_2$phylum == "Firmicutes_A")]) / numbauxogenomes *100
nauxoGlyFirmB <- nrow(Auxotrophy_2[c(Auxotrophy_2$Compound == "Gly" & 
                                       Auxotrophy_2$Prototrophy == 0 & Auxotrophy_2$phylum == "Firmicutes_B")]) / numbauxogenomes *100
nauxoGlyFirmG <- nrow(Auxotrophy_2[c(Auxotrophy_2$Compound == "Gly" & 
                                       Auxotrophy_2$Prototrophy == 0 & Auxotrophy_2$phylum == "Firmicutes_G")]) / numbauxogenomes *100
nauxoGlyFirmI <- nrow(Auxotrophy_2[c(Auxotrophy_2$Compound == "Gly" & 
                                       Auxotrophy_2$Prototrophy == 0 & Auxotrophy_2$phylum == "Firmicutes_I")]) / numbauxogenomes *100
nauxoGlyFuso <- nrow(Auxotrophy_2[c(Auxotrophy_2$Compound == "Gly" & 
                                      Auxotrophy_2$Prototrophy == 0 & Auxotrophy_2$phylum == "Fusobacteriota")]) / numbauxogenomes *100
nauxoGlyHalo <- nrow(Auxotrophy_2[c(Auxotrophy_2$Compound == "Gly" & 
                                      Auxotrophy_2$Prototrophy == 0 & Auxotrophy_2$phylum == "Halobacterota")]) / numbauxogenomes *100
nauxoGlyMyxo <- nrow(Auxotrophy_2[c(Auxotrophy_2$Compound == "Gly" & 
                                      Auxotrophy_2$Prototrophy == 0 & Auxotrophy_2$phylum == "Myxococcota")]) / numbauxogenomes *100
nauxoGlyPates <- nrow(Auxotrophy_2[c(Auxotrophy_2$Compound == "Gly" & 
                                       Auxotrophy_2$Prototrophy == 0 & Auxotrophy_2$phylum == "Patescibacteria")]) / numbauxogenomes *100
nauxoGlyProteo <- nrow(Auxotrophy_2[c(Auxotrophy_2$Compound == "Gly" & 
                                        Auxotrophy_2$Prototrophy == 0 & Auxotrophy_2$phylum == "Proteobacteria")]) / numbauxogenomes *100
nauxoGlySpiro <- nrow(Auxotrophy_2[c(Auxotrophy_2$Compound == "Gly" & 
                                       Auxotrophy_2$Prototrophy == 0 & Auxotrophy_2$phylum == "Spirochaetota")]) / numbauxogenomes *100
nauxoGlySyner <- nrow(Auxotrophy_2[c(Auxotrophy_2$Compound == "Gly" & 
                                       Auxotrophy_2$Prototrophy == 0 & Auxotrophy_2$phylum == "Synergistota")]) / numbauxogenomes *100
nauxoGlyTherm <- nrow(Auxotrophy_2[c(Auxotrophy_2$Compound == "Gly" & 
                                       Auxotrophy_2$Prototrophy == 0 & Auxotrophy_2$phylum == "Thermoplasmatota")]) / numbauxogenomes *100
nauxoGlyVerru <- nrow(Auxotrophy_2[c(Auxotrophy_2$Compound == "Gly" & 
                                       Auxotrophy_2$Prototrophy == 0 & Auxotrophy_2$phylum == "Verrucomicrobiota")]) / numbauxogenomes *100
nauxoGlyNa <- nrow(Auxotrophy_2[c(Auxotrophy_2$Compound == "Gly" & 
                                    Auxotrophy_2$Prototrophy == 0 & Auxotrophy_2$phylum == "NA")]) / numbauxogenomes *100

#dataframe for auxotrophies Glytamine:
auxoGly <- data.frame(nauxoGlyAct, nauxoGlyBac, nauxoGlyBdel, nauxoGlyCamp, nauxoGlyCyan, nauxoGlyDesulf, nauxoGlyElusi, nauxoGlyEremi,
                      nauxoGlyEury, nauxoGlyFibro, nauxoGlyFirm, nauxoGlyFirmA, nauxoGlyFirmB, nauxoGlyFirmG, nauxoGlyFirmI, nauxoGlyFuso,
                      nauxoGlyHalo, nauxoGlyMyxo, nauxoGlyNa, nauxoGlyPates, nauxoGlyProteo, nauxoGlySpiro, nauxoGlySyner, nauxoGlyTherm,
                      nauxoGlyVerru) 
auxoGlyphylum <-  t(auxoGly)
auxoGlyphyla <- data.frame(auxoGlyphylum)
colnames(auxoGlyphyla) <- "Auxo"
#Serine
nauxoSerAct <- nrow(Auxotrophy_2[c(Auxotrophy_2$Compound == "Ser" & 
                                     Auxotrophy_2$Prototrophy == 0 & Auxotrophy_2$phylum == "Actinobacteriota")]) / numbauxogenomes  *100
nauxoSerBac <- nrow(Auxotrophy_2[c(Auxotrophy_2$Compound == "Ser" & 
                                     Auxotrophy_2$Prototrophy == 0 & Auxotrophy_2$phylum == "Bacteroidota")]) / numbauxogenomes  *100
nauxoSerBdel <- nrow(Auxotrophy_2[c(Auxotrophy_2$Compound == "Ser" & 
                                      Auxotrophy_2$Prototrophy == 0 & Auxotrophy_2$phylum == "Bdellovibrionota")]) / numbauxogenomes  *100
nauxoSerCamp <- nrow(Auxotrophy_2[c(Auxotrophy_2$Compound == "Ser" & 
                                      Auxotrophy_2$Prototrophy == 0 & Auxotrophy_2$phylum == "Campylobacterota")]) / numbauxogenomes  *100
nauxoSerCyan <- nrow(Auxotrophy_2[c(Auxotrophy_2$Compound == "Ser" & 
                                      Auxotrophy_2$Prototrophy == 0 & Auxotrophy_2$phylum == "Cyanobacteria")]) / numbauxogenomes  *100
nauxoSerDesulf <- nrow(Auxotrophy_2[c(Auxotrophy_2$Compound == "Ser" & 
                                        Auxotrophy_2$Prototrophy == 0 & Auxotrophy_2$phylum == "Desulfobacterota_A")]) / numbauxogenomes  *100
nauxoSerElusi <- nrow(Auxotrophy_2[c(Auxotrophy_2$Compound == "Ser" & 
                                       Auxotrophy_2$Prototrophy == 0 & Auxotrophy_2$phylum == "Elusimicrobiota")]) / numbauxogenomes  *100
nauxoSerEremi <- nrow(Auxotrophy_2[c(Auxotrophy_2$Compound == "Ser" & 
                                       Auxotrophy_2$Prototrophy == 0 & Auxotrophy_2$phylum == "Eremiobacterota")]) / numbauxogenomes  *100
nauxoSerEury <- nrow(Auxotrophy_2[c(Auxotrophy_2$Compound == "Ser" & 
                                      Auxotrophy_2$Prototrophy == 0 & Auxotrophy_2$phylum == "Euryarchaeota")]) / numbauxogenomes  *100
nauxoSerFibro <- nrow(Auxotrophy_2[c(Auxotrophy_2$Compound == "Ser" & 
                                       Auxotrophy_2$Prototrophy == 0 & Auxotrophy_2$phylum == "Fibrobacterota")]) / numbauxogenomes  *100
nauxoSerFirm <- nrow(Auxotrophy_2[c(Auxotrophy_2$Compound == "Ser" & 
                                      Auxotrophy_2$Prototrophy == 0 & Auxotrophy_2$phylum == "Firmicutes")]) / numbauxogenomes  *100
nauxoSerFirmA <- nrow(Auxotrophy_2[c(Auxotrophy_2$Compound == "Ser" & 
                                       Auxotrophy_2$Prototrophy == 0 & Auxotrophy_2$phylum == "Firmicutes_A")]) / numbauxogenomes  *100
nauxoSerFirmB <- nrow(Auxotrophy_2[c(Auxotrophy_2$Compound == "Ser" & 
                                       Auxotrophy_2$Prototrophy == 0 & Auxotrophy_2$phylum == "Firmicutes_B")]) / numbauxogenomes  *100
nauxoSerFirmG <- nrow(Auxotrophy_2[c(Auxotrophy_2$Compound == "Ser" & 
                                       Auxotrophy_2$Prototrophy == 0 & Auxotrophy_2$phylum == "Firmicutes_G")]) / numbauxogenomes  *100
nauxoSerFirmI <- nrow(Auxotrophy_2[c(Auxotrophy_2$Compound == "Ser" & 
                                       Auxotrophy_2$Prototrophy == 0 & Auxotrophy_2$phylum == "Firmicutes_I")]) / numbauxogenomes  *100
nauxoSerFuso <- nrow(Auxotrophy_2[c(Auxotrophy_2$Compound == "Ser" & 
                                      Auxotrophy_2$Prototrophy == 0 & Auxotrophy_2$phylum == "Fusobacteriota")]) / numbauxogenomes  *100
nauxoSerHalo <- nrow(Auxotrophy_2[c(Auxotrophy_2$Compound == "Ser" & 
                                      Auxotrophy_2$Prototrophy == 0 & Auxotrophy_2$phylum == "Halobacterota")]) / numbauxogenomes  *100
nauxoSerMyxo <- nrow(Auxotrophy_2[c(Auxotrophy_2$Compound == "Ser" & 
                                      Auxotrophy_2$Prototrophy == 0 & Auxotrophy_2$phylum == "Myxococcota")]) / numbauxogenomes  *100
nauxoSerPates <- nrow(Auxotrophy_2[c(Auxotrophy_2$Compound == "Ser" & 
                                       Auxotrophy_2$Prototrophy == 0 & Auxotrophy_2$phylum == "Patescibacteria")]) / numbauxogenomes  *100
nauxoSerProteo <- nrow(Auxotrophy_2[c(Auxotrophy_2$Compound == "Ser" & 
                                        Auxotrophy_2$Prototrophy == 0 & Auxotrophy_2$phylum == "Proteobacteria")]) / numbauxogenomes  *100
nauxoSerSpiro <- nrow(Auxotrophy_2[c(Auxotrophy_2$Compound == "Ser" & 
                                       Auxotrophy_2$Prototrophy == 0 & Auxotrophy_2$phylum == "Spirochaetota")]) / numbauxogenomes  *100
nauxoSerSyner <- nrow(Auxotrophy_2[c(Auxotrophy_2$Compound == "Ser" & 
                                       Auxotrophy_2$Prototrophy == 0 & Auxotrophy_2$phylum == "Synergistota")]) / numbauxogenomes  *100
nauxoSerTherm <- nrow(Auxotrophy_2[c(Auxotrophy_2$Compound == "Ser" & 
                                       Auxotrophy_2$Prototrophy == 0 & Auxotrophy_2$phylum == "Thermoplasmatota")]) / numbauxogenomes  *100
nauxoSerVerru <- nrow(Auxotrophy_2[c(Auxotrophy_2$Compound == "Ser" & 
                                       Auxotrophy_2$Prototrophy == 0 & Auxotrophy_2$phylum == "Verrucomicrobiota")]) / numbauxogenomes  *100
nauxoSerNa <- nrow(Auxotrophy_2[c(Auxotrophy_2$Compound == "Ser" & 
                                    Auxotrophy_2$Prototrophy == 0 & Auxotrophy_2$phylum == "NA")]) / numbauxogenomes  *100

#dataframe for auxotrophies Sertamine:
auxoSer <- data.frame(nauxoSerAct, nauxoSerBac, nauxoSerBdel, nauxoSerCamp, nauxoSerCyan, nauxoSerDesulf, nauxoSerElusi, nauxoSerEremi,
                      nauxoSerEury, nauxoSerFibro, nauxoSerFirm, nauxoSerFirmA, nauxoSerFirmB, nauxoSerFirmG, nauxoSerFirmI, nauxoSerFuso,
                      nauxoSerHalo, nauxoSerMyxo, nauxoSerNa, nauxoSerPates, nauxoSerProteo, nauxoSerSpiro, nauxoSerSyner, nauxoSerTherm,
                      nauxoSerVerru) 
auxoSerphylum <-  t(auxoSer)
auxoSerphyla <- data.frame(auxoSerphylum)
colnames(auxoSerphyla) <- "Auxo"
#Cysteine
nauxoCysAct <- nrow(Auxotrophy_2[c(Auxotrophy_2$Compound == "Cys" & 
                                     Auxotrophy_2$Prototrophy == 0 & Auxotrophy_2$phylum == "Actinobacteriota")]) / numbauxogenomes  *100
nauxoCysBac <- nrow(Auxotrophy_2[c(Auxotrophy_2$Compound == "Cys" & 
                                     Auxotrophy_2$Prototrophy == 0 & Auxotrophy_2$phylum == "Bacteroidota")]) / numbauxogenomes  *100
nauxoCysBdel <- nrow(Auxotrophy_2[c(Auxotrophy_2$Compound == "Cys" & 
                                      Auxotrophy_2$Prototrophy == 0 & Auxotrophy_2$phylum == "Bdellovibrionota")]) / numbauxogenomes  *100
nauxoCysCamp <- nrow(Auxotrophy_2[c(Auxotrophy_2$Compound == "Cys" & 
                                      Auxotrophy_2$Prototrophy == 0 & Auxotrophy_2$phylum == "Campylobacterota")]) / numbauxogenomes  *100
nauxoCysCyan <- nrow(Auxotrophy_2[c(Auxotrophy_2$Compound == "Cys" & 
                                      Auxotrophy_2$Prototrophy == 0 & Auxotrophy_2$phylum == "Cyanobacteria")]) / numbauxogenomes  *100
nauxoCysDesulf <- nrow(Auxotrophy_2[c(Auxotrophy_2$Compound == "Cys" & 
                                        Auxotrophy_2$Prototrophy == 0 & Auxotrophy_2$phylum == "Desulfobacterota_A")]) / numbauxogenomes  *100
nauxoCysElusi <- nrow(Auxotrophy_2[c(Auxotrophy_2$Compound == "Cys" & 
                                       Auxotrophy_2$Prototrophy == 0 & Auxotrophy_2$phylum == "Elusimicrobiota")]) / numbauxogenomes  *100
nauxoCysEremi <- nrow(Auxotrophy_2[c(Auxotrophy_2$Compound == "Cys" & 
                                       Auxotrophy_2$Prototrophy == 0 & Auxotrophy_2$phylum == "Eremiobacterota")]) / numbauxogenomes  *100
nauxoCysEury <- nrow(Auxotrophy_2[c(Auxotrophy_2$Compound == "Cys" & 
                                      Auxotrophy_2$Prototrophy == 0 & Auxotrophy_2$phylum == "Euryarchaeota")]) / numbauxogenomes  *100
nauxoCysFibro <- nrow(Auxotrophy_2[c(Auxotrophy_2$Compound == "Cys" & 
                                       Auxotrophy_2$Prototrophy == 0 & Auxotrophy_2$phylum == "Fibrobacterota")]) / numbauxogenomes  *100
nauxoCysFirm <- nrow(Auxotrophy_2[c(Auxotrophy_2$Compound == "Cys" & 
                                      Auxotrophy_2$Prototrophy == 0 & Auxotrophy_2$phylum == "Firmicutes")]) / numbauxogenomes  *100
nauxoCysFirmA <- nrow(Auxotrophy_2[c(Auxotrophy_2$Compound == "Cys" & 
                                       Auxotrophy_2$Prototrophy == 0 & Auxotrophy_2$phylum == "Firmicutes_A")]) / numbauxogenomes  *100
nauxoCysFirmB <- nrow(Auxotrophy_2[c(Auxotrophy_2$Compound == "Cys" & 
                                       Auxotrophy_2$Prototrophy == 0 & Auxotrophy_2$phylum == "Firmicutes_B")]) / numbauxogenomes  *100
nauxoCysFirmG <- nrow(Auxotrophy_2[c(Auxotrophy_2$Compound == "Cys" & 
                                       Auxotrophy_2$Prototrophy == 0 & Auxotrophy_2$phylum == "Firmicutes_G")]) / numbauxogenomes  *100
nauxoCysFirmI <- nrow(Auxotrophy_2[c(Auxotrophy_2$Compound == "Cys" & 
                                       Auxotrophy_2$Prototrophy == 0 & Auxotrophy_2$phylum == "Firmicutes_I")]) / numbauxogenomes  *100
nauxoCysFuso <- nrow(Auxotrophy_2[c(Auxotrophy_2$Compound == "Cys" & 
                                      Auxotrophy_2$Prototrophy == 0 & Auxotrophy_2$phylum == "Fusobacteriota")]) / numbauxogenomes  *100
nauxoCysHalo <- nrow(Auxotrophy_2[c(Auxotrophy_2$Compound == "Cys" & 
                                      Auxotrophy_2$Prototrophy == 0 & Auxotrophy_2$phylum == "Halobacterota")]) / numbauxogenomes  *100
nauxoCysMyxo <- nrow(Auxotrophy_2[c(Auxotrophy_2$Compound == "Cys" & 
                                      Auxotrophy_2$Prototrophy == 0 & Auxotrophy_2$phylum == "Myxococcota")]) / numbauxogenomes  *100
nauxoCysPates <- nrow(Auxotrophy_2[c(Auxotrophy_2$Compound == "Cys" & 
                                       Auxotrophy_2$Prototrophy == 0 & Auxotrophy_2$phylum == "Patescibacteria")]) / numbauxogenomes  *100
nauxoCysProteo <- nrow(Auxotrophy_2[c(Auxotrophy_2$Compound == "Cys" & 
                                        Auxotrophy_2$Prototrophy == 0 & Auxotrophy_2$phylum == "Proteobacteria")]) / numbauxogenomes  *100
nauxoCysSpiro <- nrow(Auxotrophy_2[c(Auxotrophy_2$Compound == "Cys" & 
                                       Auxotrophy_2$Prototrophy == 0 & Auxotrophy_2$phylum == "Spirochaetota")]) / numbauxogenomes  *100
nauxoCysSyner <- nrow(Auxotrophy_2[c(Auxotrophy_2$Compound == "Cys" & 
                                       Auxotrophy_2$Prototrophy == 0 & Auxotrophy_2$phylum == "Synergistota")]) / numbauxogenomes  *100
nauxoCysTherm <- nrow(Auxotrophy_2[c(Auxotrophy_2$Compound == "Cys" & 
                                       Auxotrophy_2$Prototrophy == 0 & Auxotrophy_2$phylum == "Thermoplasmatota")]) / numbauxogenomes  *100
nauxoCysVerru <- nrow(Auxotrophy_2[c(Auxotrophy_2$Compound == "Cys" & 
                                       Auxotrophy_2$Prototrophy == 0 & Auxotrophy_2$phylum == "Verrucomicrobiota")]) / numbauxogenomes  *100
nauxoCysNa <- nrow(Auxotrophy_2[c(Auxotrophy_2$Compound == "Cys" & 
                                    Auxotrophy_2$Prototrophy == 0 & Auxotrophy_2$phylum == "NA")]) / numbauxogenomes  *100

#dataframe for auxotrophies Cystamine:
auxoCys <- data.frame(nauxoCysAct, nauxoCysBac, nauxoCysBdel, nauxoCysCamp, nauxoCysCyan, nauxoCysDesulf, nauxoCysElusi, nauxoCysEremi,
                      nauxoCysEury, nauxoCysFibro, nauxoCysFirm, nauxoCysFirmA, nauxoCysFirmB, nauxoCysFirmG, nauxoCysFirmI, nauxoCysFuso,
                      nauxoCysHalo, nauxoCysMyxo, nauxoCysNa, nauxoCysPates, nauxoCysProteo, nauxoCysSpiro, nauxoCysSyner, nauxoCysTherm,
                      nauxoCysVerru) 
auxoCysphylum <-  t(auxoCys)
auxoCysphyla <- data.frame(auxoCysphylum)
colnames(auxoCysphyla) <- "Auxo"
#Asparagine
nauxoAsnAct <- nrow(Auxotrophy_2[c(Auxotrophy_2$Compound == "Asn" & 
                                     Auxotrophy_2$Prototrophy == 0 & Auxotrophy_2$phylum == "Actinobacteriota")]) / numbauxogenomes  *100
nauxoAsnBac <- nrow(Auxotrophy_2[c(Auxotrophy_2$Compound == "Asn" & 
                                     Auxotrophy_2$Prototrophy == 0 & Auxotrophy_2$phylum == "Bacteroidota")]) / numbauxogenomes  *100
nauxoAsnBdel <- nrow(Auxotrophy_2[c(Auxotrophy_2$Compound == "Asn" & 
                                      Auxotrophy_2$Prototrophy == 0 & Auxotrophy_2$phylum == "Bdellovibrionota")]) / numbauxogenomes  *100
nauxoAsnCamp <- nrow(Auxotrophy_2[c(Auxotrophy_2$Compound == "Asn" & 
                                      Auxotrophy_2$Prototrophy == 0 & Auxotrophy_2$phylum == "Campylobacterota")]) / numbauxogenomes  *100
nauxoAsnCyan <- nrow(Auxotrophy_2[c(Auxotrophy_2$Compound == "Asn" & 
                                      Auxotrophy_2$Prototrophy == 0 & Auxotrophy_2$phylum == "Cyanobacteria")]) / numbauxogenomes  *100
nauxoAsnDesulf <- nrow(Auxotrophy_2[c(Auxotrophy_2$Compound == "Asn" & 
                                        Auxotrophy_2$Prototrophy == 0 & Auxotrophy_2$phylum == "Desulfobacterota_A")]) / numbauxogenomes  *100
nauxoAsnElusi <- nrow(Auxotrophy_2[c(Auxotrophy_2$Compound == "Asn" & 
                                       Auxotrophy_2$Prototrophy == 0 & Auxotrophy_2$phylum == "Elusimicrobiota")]) / numbauxogenomes  *100
nauxoAsnEremi <- nrow(Auxotrophy_2[c(Auxotrophy_2$Compound == "Asn" & 
                                       Auxotrophy_2$Prototrophy == 0 & Auxotrophy_2$phylum == "Eremiobacterota")]) / numbauxogenomes  *100
nauxoAsnEury <- nrow(Auxotrophy_2[c(Auxotrophy_2$Compound == "Asn" & 
                                      Auxotrophy_2$Prototrophy == 0 & Auxotrophy_2$phylum == "Euryarchaeota")]) / numbauxogenomes  *100
nauxoAsnFibro <- nrow(Auxotrophy_2[c(Auxotrophy_2$Compound == "Asn" & 
                                       Auxotrophy_2$Prototrophy == 0 & Auxotrophy_2$phylum == "Fibrobacterota")]) / numbauxogenomes  *100
nauxoAsnFirm <- nrow(Auxotrophy_2[c(Auxotrophy_2$Compound == "Asn" & 
                                      Auxotrophy_2$Prototrophy == 0 & Auxotrophy_2$phylum == "Firmicutes")]) / numbauxogenomes  *100
nauxoAsnFirmA <- nrow(Auxotrophy_2[c(Auxotrophy_2$Compound == "Asn" & 
                                       Auxotrophy_2$Prototrophy == 0 & Auxotrophy_2$phylum == "Firmicutes_A")]) / numbauxogenomes  *100
nauxoAsnFirmB <- nrow(Auxotrophy_2[c(Auxotrophy_2$Compound == "Asn" & 
                                       Auxotrophy_2$Prototrophy == 0 & Auxotrophy_2$phylum == "Firmicutes_B")]) / numbauxogenomes  *100
nauxoAsnFirmG <- nrow(Auxotrophy_2[c(Auxotrophy_2$Compound == "Asn" & 
                                       Auxotrophy_2$Prototrophy == 0 & Auxotrophy_2$phylum == "Firmicutes_G")]) / numbauxogenomes  *100
nauxoAsnFirmI <- nrow(Auxotrophy_2[c(Auxotrophy_2$Compound == "Asn" & 
                                       Auxotrophy_2$Prototrophy == 0 & Auxotrophy_2$phylum == "Firmicutes_I")]) / numbauxogenomes  *100
nauxoAsnFuso <- nrow(Auxotrophy_2[c(Auxotrophy_2$Compound == "Asn" & 
                                      Auxotrophy_2$Prototrophy == 0 & Auxotrophy_2$phylum == "Fusobacteriota")]) / numbauxogenomes  *100
nauxoAsnHalo <- nrow(Auxotrophy_2[c(Auxotrophy_2$Compound == "Asn" & 
                                      Auxotrophy_2$Prototrophy == 0 & Auxotrophy_2$phylum == "Halobacterota")]) / numbauxogenomes  *100
nauxoAsnMyxo <- nrow(Auxotrophy_2[c(Auxotrophy_2$Compound == "Asn" & 
                                      Auxotrophy_2$Prototrophy == 0 & Auxotrophy_2$phylum == "Myxococcota")]) / numbauxogenomes  *100
nauxoAsnPates <- nrow(Auxotrophy_2[c(Auxotrophy_2$Compound == "Asn" & 
                                       Auxotrophy_2$Prototrophy == 0 & Auxotrophy_2$phylum == "Patescibacteria")]) / numbauxogenomes  *100
nauxoAsnProteo <- nrow(Auxotrophy_2[c(Auxotrophy_2$Compound == "Asn" & 
                                        Auxotrophy_2$Prototrophy == 0 & Auxotrophy_2$phylum == "Proteobacteria")]) / numbauxogenomes  *100
nauxoAsnSpiro <- nrow(Auxotrophy_2[c(Auxotrophy_2$Compound == "Asn" & 
                                       Auxotrophy_2$Prototrophy == 0 & Auxotrophy_2$phylum == "Spirochaetota")]) / numbauxogenomes  *100
nauxoAsnSyner <- nrow(Auxotrophy_2[c(Auxotrophy_2$Compound == "Asn" & 
                                       Auxotrophy_2$Prototrophy == 0 & Auxotrophy_2$phylum == "Synergistota")]) / numbauxogenomes  *100
nauxoAsnTherm <- nrow(Auxotrophy_2[c(Auxotrophy_2$Compound == "Asn" & 
                                       Auxotrophy_2$Prototrophy == 0 & Auxotrophy_2$phylum == "Thermoplasmatota")]) / numbauxogenomes  *100
nauxoAsnVerru <- nrow(Auxotrophy_2[c(Auxotrophy_2$Compound == "Asn" & 
                                       Auxotrophy_2$Prototrophy == 0 & Auxotrophy_2$phylum == "Verrucomicrobiota")]) / numbauxogenomes  *100
nauxoAsnNa <- nrow(Auxotrophy_2[c(Auxotrophy_2$Compound == "Asn" & 
                                    Auxotrophy_2$Prototrophy == 0 & Auxotrophy_2$phylum == "NA")]) / numbauxogenomes  *100

#dataframe for auxotrophies Asparagine:
auxoAsn <- data.frame(nauxoAsnAct, nauxoAsnBac, nauxoAsnBdel, nauxoAsnCamp, nauxoAsnCyan, nauxoAsnDesulf, nauxoAsnElusi, nauxoAsnEremi,
                      nauxoAsnEury, nauxoAsnFibro, nauxoAsnFirm, nauxoAsnFirmA, nauxoAsnFirmB, nauxoAsnFirmG, nauxoAsnFirmI, nauxoAsnFuso,
                      nauxoAsnHalo, nauxoAsnMyxo, nauxoAsnNa, nauxoAsnPates, nauxoAsnProteo, nauxoAsnSpiro, nauxoAsnSyner, nauxoAsnTherm,
                      nauxoAsnVerru) 
auxoAsnphylum <-  t(auxoAsn)
auxoAsnphyla <- data.frame(auxoAsnphylum)
colnames(auxoAsnphyla) <- "Auxo"
#Aspartate
nauxoAspAct <- nrow(Auxotrophy_2[c(Auxotrophy_2$Compound == "Asp" & 
                                     Auxotrophy_2$Prototrophy == 0 & Auxotrophy_2$phylum == "Actinobacteriota")]) / numbauxogenomes  *100
nauxoAspBac <- nrow(Auxotrophy_2[c(Auxotrophy_2$Compound == "Asp" & 
                                     Auxotrophy_2$Prototrophy == 0 & Auxotrophy_2$phylum == "Bacteroidota")]) / numbauxogenomes  *100
nauxoAspBdel <- nrow(Auxotrophy_2[c(Auxotrophy_2$Compound == "Asp" & 
                                      Auxotrophy_2$Prototrophy == 0 & Auxotrophy_2$phylum == "Bdellovibrionota")]) / numbauxogenomes  *100
nauxoAspCamp <- nrow(Auxotrophy_2[c(Auxotrophy_2$Compound == "Asp" & 
                                      Auxotrophy_2$Prototrophy == 0 & Auxotrophy_2$phylum == "Campylobacterota")]) / numbauxogenomes  *100
nauxoAspCyan <- nrow(Auxotrophy_2[c(Auxotrophy_2$Compound == "Asp" & 
                                      Auxotrophy_2$Prototrophy == 0 & Auxotrophy_2$phylum == "Cyanobacteria")]) / numbauxogenomes  *100
nauxoAspDesulf <- nrow(Auxotrophy_2[c(Auxotrophy_2$Compound == "Asp" & 
                                        Auxotrophy_2$Prototrophy == 0 & Auxotrophy_2$phylum == "Desulfobacterota_A")]) / numbauxogenomes  *100
nauxoAspElusi <- nrow(Auxotrophy_2[c(Auxotrophy_2$Compound == "Asp" & 
                                       Auxotrophy_2$Prototrophy == 0 & Auxotrophy_2$phylum == "Elusimicrobiota")]) / numbauxogenomes  *100
nauxoAspEremi <- nrow(Auxotrophy_2[c(Auxotrophy_2$Compound == "Asp" & 
                                       Auxotrophy_2$Prototrophy == 0 & Auxotrophy_2$phylum == "Eremiobacterota")]) / numbauxogenomes  *100
nauxoAspEury <- nrow(Auxotrophy_2[c(Auxotrophy_2$Compound == "Asp" & 
                                      Auxotrophy_2$Prototrophy == 0 & Auxotrophy_2$phylum == "Euryarchaeota")]) / numbauxogenomes  *100
nauxoAspFibro <- nrow(Auxotrophy_2[c(Auxotrophy_2$Compound == "Asp" & 
                                       Auxotrophy_2$Prototrophy == 0 & Auxotrophy_2$phylum == "Fibrobacterota")]) / numbauxogenomes  *100
nauxoAspFirm <- nrow(Auxotrophy_2[c(Auxotrophy_2$Compound == "Asp" & 
                                      Auxotrophy_2$Prototrophy == 0 & Auxotrophy_2$phylum == "Firmicutes")]) / numbauxogenomes  *100
nauxoAspFirmA <- nrow(Auxotrophy_2[c(Auxotrophy_2$Compound == "Asp" & 
                                       Auxotrophy_2$Prototrophy == 0 & Auxotrophy_2$phylum == "Firmicutes_A")]) / numbauxogenomes  *100
nauxoAspFirmB <- nrow(Auxotrophy_2[c(Auxotrophy_2$Compound == "Asp" & 
                                       Auxotrophy_2$Prototrophy == 0 & Auxotrophy_2$phylum == "Firmicutes_B")]) / numbauxogenomes  *100
nauxoAspFirmG <- nrow(Auxotrophy_2[c(Auxotrophy_2$Compound == "Asp" & 
                                       Auxotrophy_2$Prototrophy == 0 & Auxotrophy_2$phylum == "Firmicutes_G")]) / numbauxogenomes  *100
nauxoAspFirmI <- nrow(Auxotrophy_2[c(Auxotrophy_2$Compound == "Asp" & 
                                       Auxotrophy_2$Prototrophy == 0 & Auxotrophy_2$phylum == "Firmicutes_I")]) / numbauxogenomes  *100
nauxoAspFuso <- nrow(Auxotrophy_2[c(Auxotrophy_2$Compound == "Asp" & 
                                      Auxotrophy_2$Prototrophy == 0 & Auxotrophy_2$phylum == "Fusobacteriota")]) / numbauxogenomes  *100
nauxoAspHalo <- nrow(Auxotrophy_2[c(Auxotrophy_2$Compound == "Asp" & 
                                      Auxotrophy_2$Prototrophy == 0 & Auxotrophy_2$phylum == "Halobacterota")]) / numbauxogenomes  *100
nauxoAspMyxo <- nrow(Auxotrophy_2[c(Auxotrophy_2$Compound == "Asp" & 
                                      Auxotrophy_2$Prototrophy == 0 & Auxotrophy_2$phylum == "Myxococcota")]) / numbauxogenomes  *100
nauxoAspPates <- nrow(Auxotrophy_2[c(Auxotrophy_2$Compound == "Asp" & 
                                       Auxotrophy_2$Prototrophy == 0 & Auxotrophy_2$phylum == "Patescibacteria")]) / numbauxogenomes  *100
nauxoAspProteo <- nrow(Auxotrophy_2[c(Auxotrophy_2$Compound == "Asp" & 
                                        Auxotrophy_2$Prototrophy == 0 & Auxotrophy_2$phylum == "Proteobacteria")]) / numbauxogenomes  *100
nauxoAspSpiro <- nrow(Auxotrophy_2[c(Auxotrophy_2$Compound == "Asp" & 
                                       Auxotrophy_2$Prototrophy == 0 & Auxotrophy_2$phylum == "Spirochaetota")]) / numbauxogenomes  *100
nauxoAspSyner <- nrow(Auxotrophy_2[c(Auxotrophy_2$Compound == "Asp" & 
                                       Auxotrophy_2$Prototrophy == 0 & Auxotrophy_2$phylum == "Synergistota")]) / numbauxogenomes  *100
nauxoAspTherm <- nrow(Auxotrophy_2[c(Auxotrophy_2$Compound == "Asp" & 
                                       Auxotrophy_2$Prototrophy == 0 & Auxotrophy_2$phylum == "Thermoplasmatota")]) / numbauxogenomes  *100
nauxoAspVerru <- nrow(Auxotrophy_2[c(Auxotrophy_2$Compound == "Asp" & 
                                       Auxotrophy_2$Prototrophy == 0 & Auxotrophy_2$phylum == "Verrucomicrobiota")]) / numbauxogenomes  *100
nauxoAspNa <- nrow(Auxotrophy_2[c(Auxotrophy_2$Compound == "Asp" & 
                                    Auxotrophy_2$Prototrophy == 0 & Auxotrophy_2$phylum == "NA")]) / numbauxogenomes  *100

#dataframe for auxotrophies Asptamine:
auxoAsp <- data.frame(nauxoAspAct, nauxoAspBac, nauxoAspBdel, nauxoAspCamp, nauxoAspCyan, nauxoAspDesulf, nauxoAspElusi, nauxoAspEremi,
                      nauxoAspEury, nauxoAspFibro, nauxoAspFirm, nauxoAspFirmA, nauxoAspFirmB, nauxoAspFirmG, nauxoAspFirmI, nauxoAspFuso,
                      nauxoAspHalo, nauxoAspMyxo, nauxoAspNa, nauxoAspPates, nauxoAspProteo, nauxoAspSpiro, nauxoAspSyner, nauxoAspTherm,
                      nauxoAspVerru) 
auxoAspphylum <-  t(auxoAsp)
auxoAspphyla <- data.frame(auxoAspphylum)
colnames(auxoAspphyla) <- "Auxo"
#Chorismate
nauxoChorAct <- nrow(Auxotrophy_2[c(Auxotrophy_2$Compound == "Chor" & 
                                      Auxotrophy_2$Prototrophy == 0 & Auxotrophy_2$phylum == "Actinobacteriota")]) / numbauxogenomes  *100
nauxoChorBac <- nrow(Auxotrophy_2[c(Auxotrophy_2$Compound == "Chor" & 
                                      Auxotrophy_2$Prototrophy == 0 & Auxotrophy_2$phylum == "Bacteroidota")]) / numbauxogenomes  *100
nauxoChorBdel <- nrow(Auxotrophy_2[c(Auxotrophy_2$Compound == "Chor" & 
                                       Auxotrophy_2$Prototrophy == 0 & Auxotrophy_2$phylum == "Bdellovibrionota")]) / numbauxogenomes  *100
nauxoChorCamp <- nrow(Auxotrophy_2[c(Auxotrophy_2$Compound == "Chor" & 
                                       Auxotrophy_2$Prototrophy == 0 & Auxotrophy_2$phylum == "Campylobacterota")]) / numbauxogenomes  *100
nauxoChorCyan <- nrow(Auxotrophy_2[c(Auxotrophy_2$Compound == "Chor" & 
                                       Auxotrophy_2$Prototrophy == 0 & Auxotrophy_2$phylum == "Cyanobacteria")]) / numbauxogenomes  *100
nauxoChorDesulf <- nrow(Auxotrophy_2[c(Auxotrophy_2$Compound == "Chor" & 
                                         Auxotrophy_2$Prototrophy == 0 & Auxotrophy_2$phylum == "Desulfobacterota_A")]) / numbauxogenomes  *100
nauxoChorElusi <- nrow(Auxotrophy_2[c(Auxotrophy_2$Compound == "Chor" & 
                                        Auxotrophy_2$Prototrophy == 0 & Auxotrophy_2$phylum == "Elusimicrobiota")]) / numbauxogenomes  *100
nauxoChorEremi <- nrow(Auxotrophy_2[c(Auxotrophy_2$Compound == "Chor" & 
                                        Auxotrophy_2$Prototrophy == 0 & Auxotrophy_2$phylum == "Eremiobacterota")]) / numbauxogenomes  *100
nauxoChorEury <- nrow(Auxotrophy_2[c(Auxotrophy_2$Compound == "Chor" & 
                                       Auxotrophy_2$Prototrophy == 0 & Auxotrophy_2$phylum == "Euryarchaeota")]) / numbauxogenomes  *100
nauxoChorFibro <- nrow(Auxotrophy_2[c(Auxotrophy_2$Compound == "Chor" & 
                                        Auxotrophy_2$Prototrophy == 0 & Auxotrophy_2$phylum == "Fibrobacterota")]) / numbauxogenomes  *100
nauxoChorFirm <- nrow(Auxotrophy_2[c(Auxotrophy_2$Compound == "Chor" & 
                                       Auxotrophy_2$Prototrophy == 0 & Auxotrophy_2$phylum == "Firmicutes")]) / numbauxogenomes  *100
nauxoChorFirmA <- nrow(Auxotrophy_2[c(Auxotrophy_2$Compound == "Chor" & 
                                        Auxotrophy_2$Prototrophy == 0 & Auxotrophy_2$phylum == "Firmicutes_A")]) / numbauxogenomes  *100
nauxoChorFirmB <- nrow(Auxotrophy_2[c(Auxotrophy_2$Compound == "Chor" & 
                                        Auxotrophy_2$Prototrophy == 0 & Auxotrophy_2$phylum == "Firmicutes_B")]) / numbauxogenomes  *100
nauxoChorFirmG <- nrow(Auxotrophy_2[c(Auxotrophy_2$Compound == "Chor" & 
                                        Auxotrophy_2$Prototrophy == 0 & Auxotrophy_2$phylum == "Firmicutes_G")]) / numbauxogenomes  *100
nauxoChorFirmI <- nrow(Auxotrophy_2[c(Auxotrophy_2$Compound == "Chor" & 
                                        Auxotrophy_2$Prototrophy == 0 & Auxotrophy_2$phylum == "Firmicutes_I")]) / numbauxogenomes  *100
nauxoChorFuso <- nrow(Auxotrophy_2[c(Auxotrophy_2$Compound == "Chor" & 
                                       Auxotrophy_2$Prototrophy == 0 & Auxotrophy_2$phylum == "Fusobacteriota")]) / numbauxogenomes  *100
nauxoChorHalo <- nrow(Auxotrophy_2[c(Auxotrophy_2$Compound == "Chor" & 
                                       Auxotrophy_2$Prototrophy == 0 & Auxotrophy_2$phylum == "Halobacterota")]) / numbauxogenomes  *100
nauxoChorMyxo <- nrow(Auxotrophy_2[c(Auxotrophy_2$Compound == "Chor" & 
                                       Auxotrophy_2$Prototrophy == 0 & Auxotrophy_2$phylum == "Myxococcota")]) / numbauxogenomes  *100
nauxoChorPates <- nrow(Auxotrophy_2[c(Auxotrophy_2$Compound == "Chor" & 
                                        Auxotrophy_2$Prototrophy == 0 & Auxotrophy_2$phylum == "Patescibacteria")]) / numbauxogenomes  *100
nauxoChorProteo <- nrow(Auxotrophy_2[c(Auxotrophy_2$Compound == "Chor" & 
                                         Auxotrophy_2$Prototrophy == 0 & Auxotrophy_2$phylum == "Proteobacteria")]) / numbauxogenomes  *100
nauxoChorSpiro <- nrow(Auxotrophy_2[c(Auxotrophy_2$Compound == "Chor" & 
                                        Auxotrophy_2$Prototrophy == 0 & Auxotrophy_2$phylum == "Spirochaetota")]) / numbauxogenomes  *100
nauxoChorSyner <- nrow(Auxotrophy_2[c(Auxotrophy_2$Compound == "Chor" & 
                                        Auxotrophy_2$Prototrophy == 0 & Auxotrophy_2$phylum == "Synergistota")]) / numbauxogenomes  *100
nauxoChorTherm <- nrow(Auxotrophy_2[c(Auxotrophy_2$Compound == "Chor" & 
                                        Auxotrophy_2$Prototrophy == 0 & Auxotrophy_2$phylum == "Thermoplasmatota")]) / numbauxogenomes  *100
nauxoChorVerru <- nrow(Auxotrophy_2[c(Auxotrophy_2$Compound == "Chor" & 
                                        Auxotrophy_2$Prototrophy == 0 & Auxotrophy_2$phylum == "Verrucomicrobiota")]) / numbauxogenomes  *100
nauxoChorNa <- nrow(Auxotrophy_2[c(Auxotrophy_2$Compound == "Chor" & 
                                     Auxotrophy_2$Prototrophy == 0 & Auxotrophy_2$phylum == "NA")]) / numbauxogenomes  *100

#dataframe for auxotrophies Chortamine:
auxoChor <- data.frame(nauxoChorAct, nauxoChorBac, nauxoChorBdel, nauxoChorCamp, nauxoChorCyan, nauxoChorDesulf, nauxoChorElusi, nauxoChorEremi,
                       nauxoChorEury, nauxoChorFibro, nauxoChorFirm, nauxoChorFirmA, nauxoChorFirmB, nauxoChorFirmG, nauxoChorFirmI, nauxoChorFuso,
                       nauxoChorHalo, nauxoChorMyxo, nauxoChorNa, nauxoChorPates, nauxoChorProteo, nauxoChorSpiro, nauxoChorSyner, nauxoChorTherm,
                       nauxoChorVerru) 
auxoChorphylum <-  t(auxoChor)
auxoChorphyla <- data.frame(auxoChorphylum)
colnames(auxoChorphyla) <- "Auxo"

#alle auxophyla mergen 
install.packages("dplyr")
library(dplyr)
All_Auxos <- bind_rows(auxoAlaphyla, auxoValphyla, auxoMetphyla, auxoLeuphyla, 
                       auxoIlephyla, auxoProphyla, auxoTrpphyla, auxoPhephyla,
                       auxoLysphyla, auxoArgphyla,  auxoHisphyla,  auxoTyrphyla,
                       auxoThrphyla, auxoGlnphyla, auxoGluphyla, auxoGlyphyla,
                       auxoSerphyla, auxoCysphyla, auxoAsnphyla, auxoAspphyla,
                       auxoChorphyla)

All_Auxos
#add a cloumn with phylum (ist immer dasselbe, wahrscheinlich einmal einen Vektor mit allen Phyla erstellen und hinzufügen und der müsste sich 
#dann selber verdoppeln nach unten hin)
#merge irrelevant phyla to other
#others include: Bdellovibrionota, Eluismicrobiota, Eremiobacterota, Halobacterota, Myxococcota, NA, Patescibacteria, Thermoplasmatota
phylum <- c("Actinobacteriota", "Bacteroidota", "Other", "Campylobacterota",
            "Other", "Other", "Other", "Other",
            "Euryarchaeota", "Other", "Firmicutes", "Firmicutes", "Firmicutes",
            "Firmicutes", "Firmicutes", "Fusobacteriota", "Other", "Other", 
            "Other", "Other", "Proteobacteria", "Other", "Other",
            "Other", "Other")
All_Auxos$Phylum <- phylum

#phylum <- c("Actinobacteriota", "Bacteroidota", "Bdellovibrionota", "Campylobacterota",
#            "Cyanobacteria", "Desulfobacterota_A", "Elusimicrobiota", "Eremiobacterota",
#            "Euryarchaeota", "Fibrobacterota", "Firmicutes", "Firmicutes_A", "Firmicutes_B",
#            "Firmicutes_G", "Firmicutes_I", "Fusobacteriota", "Halobacterota", "Myxococcota", 
#           "NA", "Patescibacteria", "Proteobacteria", "Spirochaetota", "Synergistota",
#            "Thermoplasmatota", "Verrucomicrobiota")

#add a cloumn with the amino acids
Amino_acids <- rep(c("Ala", "Val", "Met", "Leu", "Ile",
                     "Pro", "Trp", "Phe", "Lys", "Arg",
                     "His", "Tyr", "Thr", "Gln", "Glu",
                     "Gly", "Ser", "Cys", "Asn", "Asp", 
                     "Chor"), each = 25)
All_Auxos$Aminoacids <- Amino_acids
View(All_Auxos)
#delete archaes
All_Auxos <- All_Auxos[!(All_Auxos$Phylum == "Euryarchaeota" | All_Auxos$Phylum == "Thermoplasmatota"), ]


All_Auxos
View(All_Auxos)
#visualization
brewer_palette <- c("#d73027", "#fc8d59", "#fee090", "#ffffbf", "#e0f3f8", "#91bfdb", "#4575b4")
library(ggplot2)
r2 <- ggplot(data=All_Auxos, aes(x=Aminoacids, y=Auxo, fill=Phylum)) +
  geom_bar(stat="identity") + 
  ylab("Auxotrophies [%]") +
  coord_cartesian(ylim = c(1,100)) +
  xlab("Amino acids") +
  theme(panel.background = element_rect(fill="white", colour= "black")) +
  scale_fill_manual(values = brewer_palette)+
  ggtitle("Abundancies of amino acid auxotrophies in HRGM genomes")
r2

r3 <- ggplot(data=All_Auxos, aes(x=Aminoacids, y=Auxo, fill=Phylum)) +
  geom_bar(stat="identity") + 
  ylab("Auxotrophien [%]") +
  coord_cartesian(ylim = c(1,100)) +
  xlab("Aminosäuren") +
  theme(panel.background = element_rect(fill="white", colour= "black")) +
  scale_fill_manual(values = brewer_palette)+
  ggtitle("Häufigkeiten von Aminosäure-Auxotrophien")
r3




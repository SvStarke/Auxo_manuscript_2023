###########################     SCFA production   ##############################
library(MicrobiomeGS2)
library(ggplot2)
library(data.table)
exchange <- get_exchanges(models)
View(exchange)
fwrite(exchange, file = "exchange.csv")

SCF <- exchange$name == "Butyrate" | exchange$name == "Propionate" | exchange$name == "Acetate" |
  exchange$name == "DL-Lactate" | exchange$name == "Formate" | exchange$name == "Fumarate" | 
  exchange$name == "Succinate"


SCFAs <- exchange[SCF, ]
is.data.table(SCFAs)
View(SCFAs)
remove(SCFA)

#get growth rates
m_gr <- lapply(models, FUN = get_growth)
head(m_gr)
m_growth <- data.table(Genome = names(m_gr),
                       Growth = unlist(m_gr))
fwrite(m_growth, file = "m_grwoth.csv")
#merge the files
auxo_growth <- merge(Auxotrophy_2, m_growth, by.x = "Genomes",
                     by.y = "Genome", allow.cartesian=TRUE)
fwrite(Auxotrophy_2, file="Auxotrophy_2.csv")

auxo_Prod1 <- merge(auxo_growth, SCFAs, by.x = "Genomes",
                    by.y = "model", allow.cartesian=TRUE)
is.data.table(auxo_Prod1)
View(auxo_Prod1)

#delete archaes
auxo_Prod1<- auxo_Prod1[!(auxo_Prod1$phylum == "Euryarchaeota" | auxo_Prod1$phylum == "Thermoplasmatota"), ]
View(auxo_Prod1)

#put all Firmicutes phyla to one Firmicutes phylum
auxo_Prod1$phylum[auxo_Prod1$phylum == "Firmicutes_A"] <- "Firmicutes"
auxo_Prod1$phylum[auxo_Prod1$phylum == "Firmicutes_B"] <- "Firmicutes"
auxo_Prod1$phylum[auxo_Prod1$phylum == "Firmicutes_G"] <- "Firmicutes"
auxo_Prod1$phylum[auxo_Prod1$phylum == "Firmicutes_I"] <- "Firmicutes"

#create a column with productionrates
auxo_Prod1$prod <- auxo_Prod1$flux/auxo_Prod1$Growth

############################       SCFA     ####################################
#first decide which SCFA is picked for the analysis
SCFA <- "Butyrate"
SCFA <- "Acetate"
SCFA <- "Propionate"
SCFA <- "DL-Lactate"
SCFA <- "Formate"
SCFA <- "Fumarate"
SCFA <- "Succinate"

#only if one SCFA was run before
remove(S2)
remove(Ala,Arg,Asn,Asp,Chor,Cys,Gln,Glu,Gly,His,Ile,Leu,Lys,Met,Phe,Pro,Ser,Thr,Trp,Tyr,Val)
remove(pvalueArg,pvalueAsn,pvalueChor,pvalueCys,
         pvalueGln,pvalueGly,pvalueHis,pvalueIle,pvalueLeu,
         pvalueLys,pvalueMet,pvaluePhe,pvaluePro,pvalueSer, pvalueThr,
         pvalueTrp,pvalueTyr,pvalueVal)
remove(pvalueArg,testAsn,testChor,testCys,
  testGln,testGly,testHis,testIle,testLeu,
  testLys,testMet,testPhe,testPro,testSer, testThr,
  testTrp,testTyr,testVal)

#analyze the production of SCFA
S2 <- auxo_Prod1[grepl(SCFA, name) & Prototrophy <=1]
S2$prod[S2$prod < 0] <- 0
S2$Production = ifelse(S2$flux > 0, "Yes", "No")
S2$new <- paste(S2$Compound, S2$Prototrophy)
View(S2)

########################    statistical analysis  ##############################

testAla <- wilcox.test(S2$prod[S2$new == "Ala 1"], S2$prod[S2$new == "Ala 0"],
                       paired = FALSE, conf.level = 0.95)
pvalueAla <- testAla$p.value

testArg <- wilcox.test(S2$prod[S2$new == "Arg 1"], S2$prod[S2$new == "Arg 0"],
                       paired = FALSE, conf.level = 0.95)
pvalueArg <- testArg$p.value
testAsn <- wilcox.test(S2$prod[S2$new == "Asn 1"], S2$prod[S2$new == "Asn 0"],
                       paired = FALSE, conf.level = 0.95)
pvalueAsn <- testAsn$p.value

testAsp <- wilcox.test(S2$prod[S2$new == "Asp 1"], S2$prod[S2$new == "Asp 0"],
                       paired = FALSE, conf.level = 0.95)
pvalueAsp <- testAsp$p.value
testChor <- wilcox.test(S2$prod[S2$new == "Chor 1"], S2$prod[S2$new == "Chor 0"],
                        paired = FALSE, conf.level = 0.95)
pvalueChor <- testChor$p.value

testCys <- wilcox.test(S2$prod[S2$new == "Cys 1"], S2$prod[S2$new == "Cys 0"],
                       paired = FALSE, conf.level = 0.95)
pvalueCys <- testCys$p.value
testGln <- wilcox.test(S2$prod[S2$new == "Gln 1"], S2$prod[S2$new == "Gln 0"],
                       paired = FALSE, conf.level = 0.95)
pvalueGln <- testGln$p.value

testGlu <- wilcox.test(S2$prod[S2$new == "Glu 1"], S2$prod[S2$new == "Glu 0"],
                       paired = FALSE, conf.level = 0.95)
pvalueGlu <- testGlu$p.value
testGly <- wilcox.test(S2$prod[S2$new == "Gly 1"], S2$prod[S2$new == "Gly 0"],
                       paired = FALSE, conf.level = 0.95)
pvalueGly <- testGly$p.value

testHis <- wilcox.test(S2$prod[S2$new == "His 1"], S2$prod[S2$new == "His 0"],
                       paired = FALSE, conf.level = 0.95)
pvalueHis <- testHis$p.value
testIle <- wilcox.test(S2$prod[S2$new == "Ile 1"], S2$prod[S2$new == "Ile 0"],
                       paired = FALSE, conf.level = 0.95)
pvalueIle <- testIle$p.value

testLeu <- wilcox.test(S2$prod[S2$new == "Leu 1"], S2$prod[S2$new == "Leu 0"],
                       paired = FALSE, conf.level = 0.95)
pvalueLeu <- testLeu$p.value
testLys <- wilcox.test(S2$prod[S2$new == "Lys 1"], S2$prod[S2$new == "Lys 0"],
                       paired = FALSE, conf.level = 0.95)
pvalueLys <- testLys$p.value

testMet <- wilcox.test(S2$prod[S2$new == "Met 1"], S2$prod[S2$new == "Met 0"],
                       paired = FALSE, conf.level = 0.95)
pvalueMet <- testMet$p.value
testPhe <- wilcox.test(S2$prod[S2$new == "Phe 1"], S2$prod[S2$new == "Phe 0"],
                       paired = FALSE, conf.level = 0.95)
pvaluePhe <- testPhe$p.value

testPro <- wilcox.test(S2$prod[S2$new == "Pro 1"], S2$prod[S2$new == "Pro 0"],
                       paired = FALSE, conf.level = 0.95)
pvaluePro <- testPro$p.value
testSer <- wilcox.test(S2$prod[S2$new == "Ser 1"], S2$prod[S2$new == "Ser 0"],
                       paired = FALSE, conf.level = 0.95)
pvalueSer <- testSer$p.value

testThr <- wilcox.test(S2$prod[S2$new == "Thr 1"], S2$prod[S2$new == "Thr 0"],
                       paired = FALSE, conf.level = 0.95)
pvalueThr <- testThr$p.value
testTrp <- wilcox.test(S2$prod[S2$new == "Trp 1"], S2$prod[S2$new == "Trp 0"],
                       paired = FALSE, conf.level = 0.95)
pvalueTrp <- testTrp$p.value

testTyr <- wilcox.test(S2$prod[S2$new == "Tyr 1"], S2$prod[S2$new == "Tyr 0"],
                       paired = FALSE, conf.level = 0.95)
pvalueTyr <- testTyr$p.value
testVal <- wilcox.test(S2$prod[S2$new == "Val 1"], S2$prod[S2$new == "Val 0"],
                       paired = FALSE, conf.level = 0.95)
pvalueVal <- testVal$p.value

###############  preparation of data for   fisher test #########################
#Ala 
Ala1 <- nrow(S2[c(S2$Prototrophy == 0 & S2$Production == "Yes" & S2$Compound == "Ala")])
Ala2<- nrow(S2[c(S2$Prototrophy == 0 & S2$Production == "No" & S2$Compound == "Ala")])
Ala3 <- nrow(S2[c(S2$Prototrophy == 1 & S2$Production == "Yes"& S2$Compound == "Ala")])
Ala4 <- nrow(S2[c(S2$Prototrophy == 1 & S2$Production == "No" & S2$Compound == "Ala")])
Ala <- data.frame(Ala1,Ala2,Ala3,Ala4)
colnames(Ala) <- c("a","b","c","d")

#Arg
Arg1 <- nrow(S2[c(S2$Prototrophy == 0 & S2$Production == "Yes" & S2$Compound == "Arg")])
Arg2<- nrow(S2[c(S2$Prototrophy == 0 & S2$Production == "No" & S2$Compound == "Arg")])
Arg3 <- nrow(S2[c(S2$Prototrophy == 1 & S2$Production == "Yes"& S2$Compound == "Arg")])
Arg4 <- nrow(S2[c(S2$Prototrophy == 1 & S2$Production == "No" & S2$Compound == "Arg")])
Arg <- data.frame(Arg1,Arg2,Arg3,Arg4)
colnames(Arg) <- c("a","b","c","d")
#Asn
Asn1 <- nrow(S2[c(S2$Prototrophy == 0 & S2$Production == "Yes" & S2$Compound == "Asn")])
Asn2<- nrow(S2[c(S2$Prototrophy == 0 & S2$Production == "No" & S2$Compound == "Asn")])
Asn3 <- nrow(S2[c(S2$Prototrophy == 1 & S2$Production == "Yes"& S2$Compound == "Asn")])
Asn4 <- nrow(S2[c(S2$Prototrophy == 1 & S2$Production == "No" & S2$Compound == "Asn")])
Asn <- data.frame(Asn1,Asn2,Asn3,Asn4)
colnames(Asn) <- c("a","b","c","d")
#Asp
Asp1 <- nrow(S2[c(S2$Prototrophy == 0 & S2$Production == "Yes" & S2$Compound == "Asp")])
Asp2<- nrow(S2[c(S2$Prototrophy == 0 & S2$Production == "No" & S2$Compound == "Asp")])
Asp3 <- nrow(S2[c(S2$Prototrophy == 1 & S2$Production == "Yes"& S2$Compound == "Asp")])
Asp4 <- nrow(S2[c(S2$Prototrophy == 1 & S2$Production == "No" & S2$Compound == "Asp")])
Asp <- data.frame(Asp1,Asp2,Asp3,Asp4)
colnames(Asp) <- c("a","b","c","d")
#Chor
Chor1 <- nrow(S2[c(S2$Prototrophy == 0 & S2$Production == "Yes" & S2$Compound == "Chor")])
Chor2<- nrow(S2[c(S2$Prototrophy == 0 & S2$Production == "No" & S2$Compound == "Chor")])
Chor3 <- nrow(S2[c(S2$Prototrophy == 1 & S2$Production == "Yes"& S2$Compound == "Chor")])
Chor4 <- nrow(S2[c(S2$Prototrophy == 1 & S2$Production == "No" & S2$Compound == "Chor")])
Chor <- data.frame(Chor1,Chor2,Chor3,Chor4)
colnames(Chor) <- c("a","b","c","d")
#Cys
Cys1 <- nrow(S2[c(S2$Prototrophy == 0 & S2$Production == "Yes" & S2$Compound == "Cys")])
Cys2<- nrow(S2[c(S2$Prototrophy == 0 & S2$Production == "No" & S2$Compound == "Cys")])
Cys3 <- nrow(S2[c(S2$Prototrophy == 1 & S2$Production == "Yes"& S2$Compound == "Cys")])
Cys4 <- nrow(S2[c(S2$Prototrophy == 1 & S2$Production == "No" & S2$Compound == "Cys")])
Cys <- data.frame(Cys1,Cys2,Cys3,Cys4)
colnames(Cys) <- c("a","b","c","d")
#Gln
Gln1 <- nrow(S2[c(S2$Prototrophy == 0 & S2$Production == "Yes" & S2$Compound == "Gln")])
Gln2<- nrow(S2[c(S2$Prototrophy == 0 & S2$Production == "No" & S2$Compound == "Gln")])
Gln3 <- nrow(S2[c(S2$Prototrophy == 1 & S2$Production == "Yes"& S2$Compound == "Gln")])
Gln4 <- nrow(S2[c(S2$Prototrophy == 1 & S2$Production == "No" & S2$Compound == "Gln")])
Gln <- data.frame(Gln1,Gln2,Gln3,Gln4)
colnames(Gln) <- c("a","b","c","d")
#Glu
Glu1 <- nrow(S2[c(S2$Prototrophy == 0 & S2$Production == "Yes" & S2$Compound == "Glu")])
Glu2<- nrow(S2[c(S2$Prototrophy == 0 & S2$Production == "No" & S2$Compound == "Glu")])
Glu3 <- nrow(S2[c(S2$Prototrophy == 1 & S2$Production == "Yes"& S2$Compound == "Glu")])
Glu4 <- nrow(S2[c(S2$Prototrophy == 1 & S2$Production == "No" & S2$Compound == "Glu")])
Glu <- data.frame(Glu1,Glu2,Glu3,Glu4)
colnames(Glu) <- c("a","b","c","d")
#Gly
Gly1 <- nrow(S2[c(S2$Prototrophy == 0 & S2$Production == "Yes" & S2$Compound == "Gly")])
Gly2<- nrow(S2[c(S2$Prototrophy == 0 & S2$Production == "No" & S2$Compound == "Gly")])
Gly3 <- nrow(S2[c(S2$Prototrophy == 1 & S2$Production == "Yes"& S2$Compound == "Gly")])
Gly4 <- nrow(S2[c(S2$Prototrophy == 1 & S2$Production == "No" & S2$Compound == "Gly")])
Gly <- data.frame(Gly1,Gly2,Gly3,Gly4)
colnames(Gly) <- c("a","b","c","d")
#His
His1 <- nrow(S2[c(S2$Prototrophy == 0 & S2$Production == "Yes" & S2$Compound == "His")])
His2<- nrow(S2[c(S2$Prototrophy == 0 & S2$Production == "No" & S2$Compound == "His")])
His3 <- nrow(S2[c(S2$Prototrophy == 1 & S2$Production == "Yes"& S2$Compound == "His")])
His4 <- nrow(S2[c(S2$Prototrophy == 1 & S2$Production == "No" & S2$Compound == "His")])
His <- data.frame(His1,His2,His3,His4)
colnames(His) <- c("a","b","c","d")
#Ile
Ile1 <- nrow(S2[c(S2$Prototrophy == 0 & S2$Production == "Yes" & S2$Compound == "Ile")])
Ile2<- nrow(S2[c(S2$Prototrophy == 0 & S2$Production == "No" & S2$Compound == "Ile")])
Ile3 <- nrow(S2[c(S2$Prototrophy == 1 & S2$Production == "Yes"& S2$Compound == "Ile")])
Ile4 <- nrow(S2[c(S2$Prototrophy == 1 & S2$Production == "No" & S2$Compound == "Ile")])
Ile <- data.frame(Ile1,Ile2,Ile3,Ile4)
colnames(Ile) <- c("a","b","c","d")
#Leu
Leu1 <- nrow(S2[c(S2$Prototrophy == 0 & S2$Production == "Yes" & S2$Compound == "Leu")])
Leu2<- nrow(S2[c(S2$Prototrophy == 0 & S2$Production == "No" & S2$Compound == "Leu")])
Leu3 <- nrow(S2[c(S2$Prototrophy == 1 & S2$Production == "Yes"& S2$Compound == "Leu")])
Leu4 <- nrow(S2[c(S2$Prototrophy == 1 & S2$Production == "No" & S2$Compound == "Leu")])
Leu <- data.frame(Leu1,Leu2,Leu3,Leu4)
colnames(Leu) <- c("a","b","c","d")
#Lys
Lys1 <- nrow(S2[c(S2$Prototrophy == 0 & S2$Production == "Yes" & S2$Compound == "Lys")])
Lys2<- nrow(S2[c(S2$Prototrophy == 0 & S2$Production == "No" & S2$Compound == "Lys")])
Lys3 <- nrow(S2[c(S2$Prototrophy == 1 & S2$Production == "Yes"& S2$Compound == "Lys")])
Lys4 <- nrow(S2[c(S2$Prototrophy == 1 & S2$Production == "No" & S2$Compound == "Lys")])
Lys <- data.frame(Lys1,Lys2,Lys3,Lys4)
colnames(Lys) <- c("a","b","c","d")
#Met
Met1 <- nrow(S2[c(S2$Prototrophy == 0 & S2$Production == "Yes" & S2$Compound == "Met")])
Met2<- nrow(S2[c(S2$Prototrophy == 0 & S2$Production == "No" & S2$Compound == "Met")])
Met3 <- nrow(S2[c(S2$Prototrophy == 1 & S2$Production == "Yes"& S2$Compound == "Met")])
Met4 <- nrow(S2[c(S2$Prototrophy == 1 & S2$Production == "No" & S2$Compound == "Met")])
Met <- data.frame(Met1,Met2,Met3,Met4)
colnames(Met) <- c("a","b","c","d")
#Phe
Phe1 <- nrow(S2[c(S2$Prototrophy == 0 & S2$Production == "Yes" & S2$Compound == "Phe")])
Phe2<- nrow(S2[c(S2$Prototrophy == 0 & S2$Production == "No" & S2$Compound == "Phe")])
Phe3 <- nrow(S2[c(S2$Prototrophy == 1 & S2$Production == "Yes"& S2$Compound == "Phe")])
Phe4 <- nrow(S2[c(S2$Prototrophy == 1 & S2$Production == "No" & S2$Compound == "Phe")])
Phe <- data.frame(Phe1,Phe2,Phe3,Phe4)
colnames(Phe) <- c("a","b","c","d")
#Pro
Pro1 <- nrow(S2[c(S2$Prototrophy == 0 & S2$Production == "Yes" & S2$Compound == "Pro")])
Pro2<- nrow(S2[c(S2$Prototrophy == 0 & S2$Production == "No" & S2$Compound == "Pro")])
Pro3 <- nrow(S2[c(S2$Prototrophy == 1 & S2$Production == "Yes"& S2$Compound == "Pro")])
Pro4 <- nrow(S2[c(S2$Prototrophy == 1 & S2$Production == "No" & S2$Compound == "Pro")])
Pro <- data.frame(Pro1,Pro2,Pro3,Pro4)
colnames(Pro) <- c("a","b","c","d")
#Ser
Ser1 <- nrow(S2[c(S2$Prototrophy == 0 & S2$Production == "Yes" & S2$Compound == "Ser")])
Ser2<- nrow(S2[c(S2$Prototrophy == 0 & S2$Production == "No" & S2$Compound == "Ser")])
Ser3 <- nrow(S2[c(S2$Prototrophy == 1 & S2$Production == "Yes"& S2$Compound == "Ser")])
Ser4 <- nrow(S2[c(S2$Prototrophy == 1 & S2$Production == "No" & S2$Compound == "Ser")])
Ser <- data.frame(Ser1,Ser2,Ser3,Ser4)
colnames(Ser) <- c("a","b","c","d")
#Thr
Thr1 <- nrow(S2[c(S2$Prototrophy == 0 & S2$Production == "Yes" & S2$Compound == "Thr")])
Thr2<- nrow(S2[c(S2$Prototrophy == 0 & S2$Production == "No" & S2$Compound == "Thr")])
Thr3 <- nrow(S2[c(S2$Prototrophy == 1 & S2$Production == "Yes"& S2$Compound == "Thr")])
Thr4 <- nrow(S2[c(S2$Prototrophy == 1 & S2$Production == "No" & S2$Compound == "Thr")])
Thr <- data.frame(Thr1,Thr2,Thr3,Thr4)
colnames(Thr) <- c("a","b","c","d")
#Trp
Trp1 <- nrow(S2[c(S2$Prototrophy == 0 & S2$Production == "Yes" & S2$Compound == "Trp")])
Trp2<- nrow(S2[c(S2$Prototrophy == 0 & S2$Production == "No" & S2$Compound == "Trp")])
Trp3 <- nrow(S2[c(S2$Prototrophy == 1 & S2$Production == "Yes"& S2$Compound == "Trp")])
Trp4 <- nrow(S2[c(S2$Prototrophy == 1 & S2$Production == "No" & S2$Compound == "Trp")])
Trp <- data.frame(Trp1,Trp2,Trp3,Trp4)
colnames(Trp) <- c("a","b","c","d")
#Tyr
Tyr1 <- nrow(S2[c(S2$Prototrophy == 0 & S2$Production == "Yes" & S2$Compound == "Tyr")])
Tyr2<- nrow(S2[c(S2$Prototrophy == 0 & S2$Production == "No" & S2$Compound == "Tyr")])
Tyr3 <- nrow(S2[c(S2$Prototrophy == 1 & S2$Production == "Yes"& S2$Compound == "Tyr")])
Tyr4 <- nrow(S2[c(S2$Prototrophy == 1 & S2$Production == "No" & S2$Compound == "Tyr")])
Tyr <- data.frame(Tyr1,Tyr2,Tyr3,Tyr4)
colnames(Tyr) <- c("a","b","c","d")
#Val
Val1 <- nrow(S2[c(S2$Prototrophy == 0 & S2$Production == "Yes" & S2$Compound == "Val")])
Val2<- nrow(S2[c(S2$Prototrophy == 0 & S2$Production == "No" & S2$Compound == "Val")])
Val3 <- nrow(S2[c(S2$Prototrophy == 1 & S2$Production == "Yes"& S2$Compound == "Val")])
Val4 <- nrow(S2[c(S2$Prototrophy == 1 & S2$Production == "No" & S2$Compound == "Val")])
Val <- data.frame(Val1,Val2,Val3,Val4)
colnames(Val) <- c("a","b","c","d")

#bind number of rows of every AA in one datatable
df1 <- rbind(Ala,Arg,Asn,Asp,Chor,Cys,Gln,Glu,Gly,His,Ile,Leu,Lys,Met,Phe,Pro,Ser,Thr,Trp,Tyr,Val)


#########################    fisher test   #####################################
#depending on the SCFA picked for the beginning you have to choose the next function

#butyrate
but <- t(apply(df1,1, function(x) {x1 <- fisher.test(matrix(x, nr=2))
c(x1$p.value, x1$conf.int, x1$estimate)}))
colnames(but) <- c("p-value", "CI left", "CI right", "odds_ratio")
butyrate <- as.data.frame(but)
butyrate$AA <- c("Ala", "Arg", "Asn", "Asp", "Chor", "Cys", "Gln", "Glu", "Gly", "His", "Ile", "Leu","Lys","Met",
            "Phe", "Pro", "Ser", "Thr", "Trp", "Tyr", "Val")
butyrate$SCFA <- SCFA
#making the amino acid with only no auxotrophy to NA
butyrate$Wilcox_pvalue2 <- c(NA,pvalueArg,pvalueAsn,NA,pvalueChor,pvalueCys,
                         pvalueGln,NA,pvalueGly,pvalueHis,pvalueIle,pvalueLeu,
                         pvalueLys,pvalueMet,pvaluePhe,pvaluePro,pvalueSer, pvalueThr,
                         pvalueTrp,pvalueTyr,pvalueVal)

Wilcox_pvalue <- c(NA,pvalueArg,pvalueAsn,NA,pvalueChor,pvalueCys,
  pvalueGln,NA,pvalueGly,pvalueHis,pvalueIle,pvalueLeu,
  pvalueLys,pvalueMet,pvaluePhe,pvaluePro,pvalueSer, pvalueThr,
  pvalueTrp,pvalueTyr,pvalueVal)

pvalues <- ifelse(Wilcox_pvalue < 0.05, "*","")
butyrate$Wilcox_pvalue <- pvalues

#acetate
#bind number of rows of every AA in one datatable
ace <- t(apply(df1,1, function(x) {x1 <- fisher.test(matrix(x, nr=2))
c(x1$p.value, x1$conf.int, x1$estimate)}))
colnames(ace) <- c("p-value", "CI left", "CI right", "odds_ratio")
acetate <- as.data.frame(ace)
acetate$AA <- c("Ala", "Arg", "Asn", "Asp", "Chor", "Cys", "Gln", "Glu", "Gly", "His", "Ile", "Leu","Lys","Met",
             "Phe", "Pro", "Ser", "Thr", "Trp", "Tyr", "Val")
acetate$SCFA <- SCFA
#making the amino acid with only no auxotrophy to NA
acetate$Wilcox_pvalue2 <- c(NA,pvalueArg,pvalueAsn,NA,pvalueChor,pvalueCys,
                         pvalueGln,NA,pvalueGly,pvalueHis,pvalueIle,pvalueLeu,
                         pvalueLys,pvalueMet,pvaluePhe,pvaluePro,pvalueSer, pvalueThr,
                         pvalueTrp,pvalueTyr,pvalueVal)

Wilcox_pvalue <- c(NA,pvalueArg,pvalueAsn,NA,pvalueChor,pvalueCys,
                   pvalueGln,NA,pvalueGly,pvalueHis,pvalueIle,pvalueLeu,
                   pvalueLys,pvalueMet,pvaluePhe,pvaluePro,pvalueSer, pvalueThr,
                   pvalueTrp,pvalueTyr,pvalueVal)

pvalues <- ifelse(Wilcox_pvalue < 0.05, "*","")
acetate$Wilcox_pvalue <- pvalues

#propionate
#bind number of rows of every AA in one datatable
df1 <- rbind(Ala,Arg,Asn,Asp,Chor,Cys,Gln,Glu,Gly,His,Ile,Leu,Lys,Met,Phe,Pro,Ser,Thr,Trp,Tyr,Val)

pro<- t(apply(df1,1, function(x) {x1 <- fisher.test(matrix(x, nr=2))
c(x1$p.value, x1$conf.int, x1$estimate)}))
colnames(pro) <- c("p-value", "CI left", "CI right", "odds_ratio")
propionate <- as.data.frame(pro)
propionate$AA <- c("Ala", "Arg", "Asn", "Asp", "Chor", "Cys", "Gln", "Glu", "Gly", "His", "Ile", "Leu","Lys","Met",
             "Phe", "Pro", "Ser", "Thr", "Trp", "Tyr", "Val")
propionate$SCFA <- SCFA
#making the amino acid with only no auxotrophy to NA
propionate$Wilcox_pvalue2 <- c(NA,pvalueArg,pvalueAsn,NA,pvalueChor,pvalueCys,
                         pvalueGln,NA,pvalueGly,pvalueHis,pvalueIle,pvalueLeu,
                         pvalueLys,pvalueMet,pvaluePhe,pvaluePro,pvalueSer, pvalueThr,
                         pvalueTrp,pvalueTyr,pvalueVal)

Wilcox_pvalue <- c(NA,pvalueArg,pvalueAsn,NA,pvalueChor,pvalueCys,
                   pvalueGln,NA,pvalueGly,pvalueHis,pvalueIle,pvalueLeu,
                   pvalueLys,pvalueMet,pvaluePhe,pvaluePro,pvalueSer, pvalueThr,
                   pvalueTrp,pvalueTyr,pvalueVal)

pvalues <- ifelse(Wilcox_pvalue < 0.05, "*","")
propionate$Wilcox_pvalue <- pvalues


#Lactate
#bind number of rows of every AA in one datatable
lac<- t(apply(df1,1, function(x) {x1 <- fisher.test(matrix(x, nr=2))
c(x1$p.value, x1$conf.int, x1$estimate, x1$Compound)}))
colnames(lac) <- c("p-value", "CI left", "CI right", "odds_ratio")
lactate <- as.data.frame(lac)
lactate$AA <- c("Ala", "Arg", "Asn", "Asp", "Chor", "Cys", "Gln", "Glu", "Gly", "His", "Ile", "Leu","Lys","Met",
                   "Phe", "Pro", "Ser", "Thr", "Trp", "Tyr", "Val")
lactate$SCFA <- SCFA
#making the amino acid with only no auxotrophy to NA
lactate$Wilcox_pvalue2 <- c(NA,pvalueArg,pvalueAsn,NA,pvalueChor,pvalueCys,
                               pvalueGln,NA,pvalueGly,pvalueHis,pvalueIle,pvalueLeu,
                               pvalueLys,pvalueMet,pvaluePhe,pvaluePro,pvalueSer, pvalueThr,
                               pvalueTrp,pvalueTyr,pvalueVal)

Wilcox_pvalue <- c(NA,pvalueArg,pvalueAsn,NA,pvalueChor,pvalueCys,
                   pvalueGln,NA,pvalueGly,pvalueHis,pvalueIle,pvalueLeu,
                   pvalueLys,pvalueMet,pvaluePhe,pvaluePro,pvalueSer, pvalueThr,
                   pvalueTrp,pvalueTyr,pvalueVal)

pvalues <- ifelse(Wilcox_pvalue < 0.05, "*","")
lactate$Wilcox_pvalue <- pvalues

#Formate
#bind number of rows of every AA in one datatable
form <- t(apply(df1,1, function(x) {x1 <- fisher.test(matrix(x, nr=2))
c(x1$p.value, x1$conf.int, x1$estimate)}))
colnames(form) <- c("p-value", "CI left", "CI right", "odds_ratio")
formate <- as.data.frame(form)
formate$AA <- c("Ala", "Arg", "Asn", "Asp", "Chor", "Cys", "Gln", "Glu", "Gly", "His", "Ile", "Leu","Lys","Met",
                "Phe", "Pro", "Ser", "Thr", "Trp", "Tyr", "Val")
formate$SCFA <- SCFA
#making the amino acid with only no auxotrophy to NA
formate$Wilcox_pvalue2 <- c(NA,pvalueArg,pvalueAsn,NA,pvalueChor,pvalueCys,
                            pvalueGln,NA,pvalueGly,pvalueHis,pvalueIle,pvalueLeu,
                            pvalueLys,pvalueMet,pvaluePhe,pvaluePro,pvalueSer, pvalueThr,
                            pvalueTrp,pvalueTyr,pvalueVal)

Wilcox_pvalue <- c(NA,pvalueArg,pvalueAsn,NA,pvalueChor,pvalueCys,
                   pvalueGln,NA,pvalueGly,pvalueHis,pvalueIle,pvalueLeu,
                   pvalueLys,pvalueMet,pvaluePhe,pvaluePro,pvalueSer, pvalueThr,
                   pvalueTrp,pvalueTyr,pvalueVal)

pvalues <- ifelse(Wilcox_pvalue < 0.05, "*","")
formate$Wilcox_pvalue <- pvalues

#Fumarate
#bind number of rows of every AA in one datatable
fum <- t(apply(df1,1, function(x) {x1 <- fisher.test(matrix(x, nr=2))
c(x1$p.value, x1$conf.int, x1$estimate)}))
colnames(fum) <- c("p-value", "CI left", "CI right", "odds_ratio")
fumarate <- as.data.frame(fum)
fumarate$AA <- c("Ala", "Arg", "Asn", "Asp", "Chor", "Cys", "Gln", "Glu", "Gly", "His", "Ile", "Leu","Lys","Met",
                "Phe", "Pro", "Ser", "Thr", "Trp", "Tyr", "Val")
fumarate$SCFA <- SCFA
#making the amino acid with only no auxotrophy to NA
fumarate$Wilcox_pvalue2 <- c(NA,pvalueArg,pvalueAsn,NA,pvalueChor,pvalueCys,
                            pvalueGln,NA,pvalueGly,pvalueHis,pvalueIle,pvalueLeu,
                            pvalueLys,pvalueMet,pvaluePhe,pvaluePro,pvalueSer, pvalueThr,
                            pvalueTrp,pvalueTyr,pvalueVal)

Wilcox_pvalue <- c(NA,pvalueArg,pvalueAsn,NA,pvalueChor,pvalueCys,
                   pvalueGln,NA,pvalueGly,pvalueHis,pvalueIle,pvalueLeu,
                   pvalueLys,pvalueMet,pvaluePhe,pvaluePro,pvalueSer, pvalueThr,
                   pvalueTrp,pvalueTyr,pvalueVal)

pvalues <- ifelse(Wilcox_pvalue < 0.05, "*","")
fumarate$Wilcox_pvalue <- pvalues

#Succinate
#bind number of rows of every AA in one datatable
succ <- t(apply(df1,1, function(x) {x1 <- fisher.test(matrix(x, nr=2))
c(x1$p.value, x1$conf.int, x1$estimate)}))
colnames(succ) <- c("p-value", "CI left", "CI right", "odds_ratio")
succinate <- as.data.frame(succ)
succinate$AA <- c("Ala", "Arg", "Asn", "Asp", "Chor", "Cys", "Gln", "Glu", "Gly", "His", "Ile", "Leu","Lys","Met",
                 "Phe", "Pro", "Ser", "Thr", "Trp", "Tyr", "Val")
succinate$SCFA <- SCFA
#making the amino acid with only no auxotrophy to NA
succinate$Wilcox_pvalue2 <- c(NA,pvalueArg,pvalueAsn,NA,pvalueChor,pvalueCys,
                             pvalueGln,NA,pvalueGly,pvalueHis,pvalueIle,pvalueLeu,
                             pvalueLys,pvalueMet,pvaluePhe,pvaluePro,pvalueSer, pvalueThr,
                             pvalueTrp,pvalueTyr,pvalueVal)

Wilcox_pvalue <- c(NA,pvalueArg,pvalueAsn,NA,pvalueChor,pvalueCys,
                   pvalueGln,NA,pvalueGly,pvalueHis,pvalueIle,pvalueLeu,
                   pvalueLys,pvalueMet,pvaluePhe,pvaluePro,pvalueSer, pvalueThr,
                   pvalueTrp,pvalueTyr,pvalueVal)

pvalues <- ifelse(Wilcox_pvalue < 0.05, "*","")
succinate$Wilcox_pvalue <- pvalues

#bind rows
stats <- rbind(butyrate,acetate, propionate, lactate, formate, fumarate, succinate)

###########################    visualization    ################################
library(MetBrewer)
library(ggplot2)

ggplot(stats, aes(AA, SCFA, fill = odds_ratio, label = Wilcox_pvalue)) +
  geom_tile() +
  geom_text() +
  scale_fill_gradientn(colours = met.brewer("Derain"))+
  labs(caption = " * pvalues <0.05 (Wilcox Test)") +
  ylab("Organic acid") +
  theme_minimal()

################## attempts to make a loop/function ############################
#preparation of data for fisher test
t <- aggregate(S2, list(S2$Compound, S2$Prototrophy, S2$Production), FUN=NROW)
View(t)
tu <- data.table(t$Group.1, t$Group.2, t$Group.3, t$Genomes)
colnames(tu) <- c("Compound", "Prototrophy","Production","nrows")
Compound <- c("Ala","Ala", "Asp", "Asp","Glu","Glu")
Prototrophy <- c(0,0,0,0,0,0)
Production <- c("No", "Yes", "No", "Yes","No", "Yes")
nrows <- c(0,0,0,0,0,0)
v <- data.table(Compound, Prototrophy, Production, nrows)
v1 <- rbind(tu,v)
v1 <- v1[order(v1$Compound),]
library(dplyr)
v1 %>% arrange(Compound)
a <- v1[grepl(0, Prototrophy) & Production == "Yes"]
a <- v1$nrows[v1$Prototrophy == 0 & v1$Production == "Yes"]
b <- v1$nrows[v1$Prototrophy == 0 & v1$Production == "No"]
c <- v1$nrows[v1$Prototrophy == 1 & v1$Production == "Yes"]
d <- v1$nrows[v1$Prototrophy == 1 & v1$Production == "No"]
df1 <- data.frame(a,b,c,d)

#control point
rowSums(df1)

################ create other plots (boxplot, violinplot) ######################

##########################     visualization    ################################
# ########################## auxotrophic microbiota #############################
# ggplot(S2[Prototrophy == 0],aes(Compound, prod)) +
#   geom_point(stat = "identity", position = position_fill())
# ggplot(S2[Prototrophy == 0],aes(Compound, prod)) +
#   geom_boxplot()
# 
# 
# ########################## prototrophic microbiota #############################
# ggplot(S2[Prototrophy == 1],aes(Compound, prod)) +
#   geom_point(stat = "identity", position = position_dodge2())
# ggplot(S2[Prototrophy == 1],aes(Compound, prod)) +
#   geom_boxplot()
# 
# 
# ###############################   both   #######################################
# #overview comparison auxotrophy and prototrophy
# ggplot(S2, aes(x=Prototrophy, y= prod, fill = Compound)) +
#   geom_bar(stat="identity", color ="black", position = position_dodge(width =0.9)) 



# S2$Butyrate = ifelse(S2$new== "Ala 1"& S2$Genomes == "HRGM_Genome_1856", NA,
#                      ifelse(S2$new == "Arg 1" & S2$Genomes == "HRGM_Genome_1856", pvalueArg, 
#                             ifelse(S2$new == "Asn 1" & S2$Genomes == "HRGM_Genome_1856", pvalueAsn, 
#                                    ifelse(S2$new == "Asp 1" & S2$Genomes == "HRGM_Genome_1856", NA, 
#                                           ifelse(S2$new == "Chor 1" & S2$Genomes == "HRGM_Genome_1856", pvalueChor, 
#                                                  ifelse(S2$new == "Cys 1" & S2$Genomes == "HRGM_Genome_1856", pvalueCys, 
#                                                         ifelse(S2$new == "Gln 1" & S2$Genomes == "HRGM_Genome_1856", pvalueGln, 
#                                                                ifelse(S2$new == "Glu 1" & S2$Genomes == "HRGM_Genome_1856", NA, 
#                                                                       ifelse(S2$new == "Gly 1" & S2$Genomes == "HRGM_Genome_1856", pvalueGly, 
#                                                                              ifelse(S2$new == "His 1" & S2$Genomes == "HRGM_Genome_1856", pvalueHis, 
#                                                                                     ifelse(S2$new == "Ile 1" & S2$Genomes == "HRGM_Genome_1856", pvalueIle, 
#                                                                                            ifelse(S2$new == "Leu 1" & S2$Genomes == "HRGM_Genome_1856", pvalueLeu,
#                                                                                                   ifelse(S2$new == "Lys 1" & S2$Genomes == "HRGM_Genome_1856", pvalueLys, 
#                                                                                                          ifelse(S2$new == "Met 1" & S2$Genomes == "HRGM_Genome_1856", pvalueMet, 
#                                                                                                                 ifelse(S2$new == "Phe 1" & S2$Genomes == "HRGM_Genome_1856", pvaluePhe, 
#                                                                                                                        ifelse(S2$new == "Pro 1" & S2$Genomes == "HRGM_Genome_1856", pvaluePro, 
#                                                                                                                               ifelse(S2$new == "Ser 1" & S2$Genomes == "HRGM_Genome_1856", pvalueSer, 
#                                                                                                                                      ifelse(S2$new == "Thr 1" & S2$Genomes == "HRGM_Genome_1856", pvalueThr, 
#                                                                                                                                             ifelse(S2$new == "Trp 1" & S2$Genomes == "HRGM_Genome_1856", pvalueTrp, 
#                                                                                                                                                    ifelse(S2$new == "Tyr 1" & S2$Genomes == "HRGM_Genome_1856", pvalueTyr, 
#                                                                                                                                                           ifelse(S2$new== "Val 1"& S2$Genomes == "HRGM_Genome_1856", pvalueVal, NA)))))))))))))))))))))
# S2$acetate = ifelse(S2$new== "Ala 1"& S2$Genomes == "HRGM_Genome_0063", NA,
#                     ifelse(S2$new == "Arg 1" & S2$Genomes == "HRGM_Genome_0063", pvalueArg, 
#                            ifelse(S2$new == "Asn 1" & S2$Genomes == "HRGM_Genome_0063", pvalueAsn, 
#                                   ifelse(S2$new == "Asp 1" & S2$Genomes == "HRGM_Genome_0063", NA, 
#                                          ifelse(S2$new == "Chor 1" & S2$Genomes == "HRGM_Genome_0063", pvalueChor, 
#                                                 ifelse(S2$new == "Cys 1" & S2$Genomes == "HRGM_Genome_0063", pvalueCys, 
#                                                        ifelse(S2$new == "Gln 1" & S2$Genomes == "HRGM_Genome_0063", pvalueGln, 
#                                                               ifelse(S2$new == "Glu 1" & S2$Genomes == "HRGM_Genome_0063", NA, 
#                                                                      ifelse(S2$new == "Gly 1" & S2$Genomes == "HRGM_Genome_0063", pvalueGly, 
#                                                                             ifelse(S2$new == "His 1" & S2$Genomes == "HRGM_Genome_0063", pvalueHis, 
#                                                                                    ifelse(S2$new == "Ile 1" & S2$Genomes == "HRGM_Genome_0063", pvalueIle, 
#                                                                                           ifelse(S2$new == "Leu 1" & S2$Genomes == "HRGM_Genome_0063", pvalueLeu,
#                                                                                                  ifelse(S2$new == "Lys 1" & S2$Genomes == "HRGM_Genome_0063", pvalueLys, 
#                                                                                                         ifelse(S2$new == "Met 1" & S2$Genomes == "HRGM_Genome_0063", pvalueMet, 
#                                                                                                                ifelse(S2$new == "Phe 1" & S2$Genomes == "HRGM_Genome_0063", pvaluePhe, 
#                                                                                                                       ifelse(S2$new == "Pro 1" & S2$Genomes == "HRGM_Genome_0063", pvaluePro, 
#                                                                                                                              ifelse(S2$new == "Ser 1" & S2$Genomes == "HRGM_Genome_0063", pvalueSer, 
#                                                                                                                                     ifelse(S2$new == "Thr 1" & S2$Genomes == "HRGM_Genome_0063", pvalueThr, 
#                                                                                                                                            ifelse(S2$new == "Trp 1" & S2$Genomes == "HRGM_Genome_0063", pvalueTrp, 
#                                                                                                                                                   ifelse(S2$new == "Tyr 1" & S2$Genomes == "HRGM_Genome_0063", pvalueTyr, 
#                                                                                                                                                          ifelse(S2$new== "Val 1"& S2$Genomes == "HRGM_Genome_0063", pvalueVal, NA)))))))))))))))))))))
# S2$propionate = ifelse(S2$new== "Ala 1"& S2$Genomes == "	 HRGM_Genome_1670", NA,
#                        ifelse(S2$new == "Arg 1" & S2$Genomes == "HRGM_Genome_1670", pvalueArg, 
#                               ifelse(S2$new == "Asn 1" & S2$Genomes == "HRGM_Genome_1670", pvalueAsn, 
#                                      ifelse(S2$new == "Asp 1" & S2$Genomes == "HRGM_Genome_1670", NA, 
#                                             ifelse(S2$new == "Chor 1" & S2$Genomes == "HRGM_Genome_1670", pvalueChor, 
#                                                    ifelse(S2$new == "Cys 1" & S2$Genomes == "HRGM_Genome_1670", pvalueCys, 
#                                                           ifelse(S2$new == "Gln 1" & S2$Genomes == "HRGM_Genome_1670", pvalueGln, 
#                                                                  ifelse(S2$new == "Glu 1" & S2$Genomes == "HRGM_Genome_1670", NA, 
#                                                                         ifelse(S2$new == "Gly 1" & S2$Genomes == "HRGM_Genome_1670", pvalueGly, 
#                                                                                ifelse(S2$new == "His 1" & S2$Genomes == "HRGM_Genome_1670", pvalueHis, 
#                                                                                       ifelse(S2$new == "Ile 1" & S2$Genomes == "HRGM_Genome_1670", pvalueIle, 
#                                                                                              ifelse(S2$new == "Leu 1" & S2$Genomes == "HRGM_Genome_1670", pvalueLeu,
#                                                                                                     ifelse(S2$new == "Lys 1" & S2$Genomes == "HRGM_Genome_1670", pvalueLys, 
#                                                                                                            ifelse(S2$new == "Met 1" & S2$Genomes == "HRGM_Genome_1670", pvalueMet, 
#                                                                                                                   ifelse(S2$new == "Phe 1" & S2$Genomes == "HRGM_Genome_1670", pvaluePhe, 
#                                                                                                                          ifelse(S2$new == "Pro 1" & S2$Genomes == "HRGM_Genome_1670", pvaluePro, 
#                                                                                                                                 ifelse(S2$new == "Ser 1" & S2$Genomes == "HRGM_Genome_1670", pvalueSer, 
#                                                                                                                                        ifelse(S2$new == "Thr 1" & S2$Genomes == "HRGM_Genome_1670", pvalueThr, 
#                                                                                                                                               ifelse(S2$new == "Trp 1" & S2$Genomes == "HRGM_Genome_1670", pvalueTrp, 
#                                                                                                                                                      ifelse(S2$new == "Tyr 1" & S2$Genomes == "HRGM_Genome_1670", pvalueTyr, 
#                                                                                                                                                             ifelse(S2$new== "Val 1"& S2$Genomes == "HRGM_Genome_1670", pvalueVal, NA)))))))))))))))))))))
# S2$lactate = ifelse(S2$new== "Ala 1"& S2$Genomes == "	 HRGM_Genome_0001", NA,
#                     ifelse(S2$new == "Arg 1" & S2$Genomes == "HRGM_Genome_0001", pvalueArg, 
#                            ifelse(S2$new == "Asn 1" & S2$Genomes == "HRGM_Genome_0001", pvalueAsn, 
#                                   ifelse(S2$new == "Asp 1" & S2$Genomes == "HRGM_Genome_0001", NA, 
#                                          ifelse(S2$new == "Chor 1" & S2$Genomes == "HRGM_Genome_0001", pvalueChor, 
#                                                 ifelse(S2$new == "Cys 1" & S2$Genomes == "HRGM_Genome_0001", pvalueCys, 
#                                                        ifelse(S2$new == "Gln 1" & S2$Genomes == "HRGM_Genome_0001", pvalueGln, 
#                                                               ifelse(S2$new == "Glu 1" & S2$Genomes == "HRGM_Genome_0001", NA, 
#                                                                      ifelse(S2$new == "Gly 1" & S2$Genomes == "HRGM_Genome_0001", pvalueGly, 
#                                                                             ifelse(S2$new == "His 1" & S2$Genomes == "HRGM_Genome_0001", pvalueHis, 
#                                                                                    ifelse(S2$new == "Ile 1" & S2$Genomes == "HRGM_Genome_0001", pvalueIle, 
#                                                                                           ifelse(S2$new == "Leu 1" & S2$Genomes == "HRGM_Genome_0001", pvalueLeu,
#                                                                                                  ifelse(S2$new == "Lys 1" & S2$Genomes == "HRGM_Genome_0001", pvalueLys, 
#                                                                                                         ifelse(S2$new == "Met 1" & S2$Genomes == "HRGM_Genome_0001", pvalueMet, 
#                                                                                                                ifelse(S2$new == "Phe 1" & S2$Genomes == "HRGM_Genome_0001", pvaluePhe, 
#                                                                                                                       ifelse(S2$new == "Pro 1" & S2$Genomes == "HRGM_Genome_0001", pvaluePro, 
#                                                                                                                              ifelse(S2$new == "Ser 1" & S2$Genomes == "0001", pvalueSer, 
#                                                                                                                                     ifelse(S2$new == "Thr 1" & S2$Genomes == "HRGM_Genome_0001", pvalueThr, 
#                                                                                                                                            ifelse(S2$new == "Trp 1" & S2$Genomes == "HRGM_Genome_0001", pvalueTrp, 
#                                                                                                                                                   ifelse(S2$new == "Tyr 1" & S2$Genomes == "HRGM_Genome_0001", pvalueTyr, 
#                                                                                                                                                          ifelse(S2$new== "Val 1"& S2$Genomes == "HRGM_Genome_0001", pvalueVal, NA)))))))))))))))))))))
# 
# ###########################   visualization ####################################
# #butyrate
# ggplot(S2, aes(x=new, y= prod, label = ifelse(Butyrate < 0.05, "*",""))) +
#   geom_violin(outlier.shape = NA) +
#   theme(axis.text.x = element_text (angle=90)) +
#   geom_text(vjust = -20) +
#   coord_cartesian(ylim = c(0,1)) +
#   ylab("exchange rates [mmol/gDW]\n") +
#   xlab("\nAmino acid Auxotrophy/Prototrophy") +
#   labs(caption = "\n\n0 – Auxotrophy, 1 - Prototrophy\n\n*p<0.05 (Wilcox-Test)") +
#   ggtitle("Comparison of butyrate productionrates by auxotrophic and prototrophic microbiota")
# 
# 
# #acetate
# ggplot(S2, aes(x=new, y= prod, label = ifelse(acetate < 0.05, "*",""))) +
#   geom_boxplot(outlier.shape = NA) +
#   theme(axis.text.x = element_text (angle=90)) +
#   geom_text(vjust = -5.8) +
#   coord_cartesian(ylim = c(0,50)) +
#   ylab("production rates [mmol/gDW]\n") +
#   xlab("\nAmino acid Auxotrophy/Prototrophy") +
#   labs(caption = "\n\n0 – Auxotrophy, 1 - Prototrophy\n\n*p<0.05 (Wilcox Test)") +
#   ggtitle("Comparison of acetate production rates by auxotrophic and prototrophic microbiota")
# 
# 
# #propionate
# ggplot(S2, aes(x=new, y= prod, label = ifelse(propionate < 0.05, "*",""))) +
#   geom_point (outlier.shape = NA) +
#   theme(axis.text.x = element_text (angle=90)) +
#   geom_text(vjust = -4, colour = "blue") +
#   coord_cartesian(ylim = c(0,50))+
#   ylab("production rates [mmol/gDW]\n") +
#   xlab("\nAmino acid Auxotrophy/Prototrophy") +
#   labs(caption = "\n\n0 – Auxotrophy, 1 - Prototrophy\n\n*p<0.05 (Wilcox Test)") +
#   ggtitle("Comparison of propionate production rates by auxotrophic and prototrophic microbiota")
# 
# 
# 
# 

#scatterplot for genome completeness and number of auxotrophies per genome
library(MicrobiomeGS2)
library(stringr)
install.packages("Hmisc")                 #######produces correlation matrices
install.packages("PerformanceAnalytics")  #######will be used to produce a figure that shows a lot!
install.packages("car")  

############################### create table ###################################


# fetch/ read models
models <- fetch_model_collection("/mnt/nuuk/2021/HRGM/models/")

models <- readRDS("")

# Anwendung (Hier Auxotrophies vorhersagen)
model.auxo <- lapply(models, FUN = predict_auxotrohies)


Metadata <- fread("/mnt/nuuk/2021/HRGM/REPR_Genomes_metadata.tsv")

#analyzing Auxotrophies
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
Auxotrophy_11 <- as.data.table(Auxotrophy)
View(Auxotrophy_11)
library(data.table)

Auxotrophy_20 <- merge(Auxotrophy_11, Metadata, by.x = "Genomes",
                       by.y = "HRGM name")
head(Auxotrophy_20)
View(Auxotrophy_20)

#count number of rows
Auxotrophy_11$count <- rowSums(Auxotrophy_11 == 0)
Auxotrophy_12 <- merge(Auxotrophy_11, Metadata, by.x = "Genomes",
                       by.y = "HRGM name")
names(Auxotrophy_12) [names(Auxotrophy_12) == "Completeness (%)"] <- "Completeness"


###########################  Correlation test "Spearman"  ######################

cor.test(Auxotrophy_12$Completeness, Auxotrophy_12$count, method = "spearman", exact = FALSE) 

############################ visualization #####################################

ggplot(Auxotrophy_12, aes(Completeness, count)) +
  geom_point()+
  geom_smooth() +
  xlab("Completeness [%]") +
  ylab("Number of Auxotrophies") +
  ggtitle("Correlation of Number of Auxotrophies and Completeness ") +
  theme_minimal() +
  xlim(50,100) +
  ylim(0,21)


install.packages("ggpubr")
library(ggpubr)
ggscatter(Auxotrophy_12, x="Completeness", y="count",
          add = "reg.line", conf.int =TRUE,
          cor.coef = TRUE, cor.method = "spearman",
          title="Spearman Correlation",
          xlab="Completeness [%]",
          ylab = "Number of Auxotrophies",
          ggtheme = theme_linedraw())


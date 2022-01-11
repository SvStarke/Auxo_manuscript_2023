### scatterplot for genome completeness and number of auxotrophies per genome ##

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


### scatterplot for genome completeness and number of auxotrophies per genome ##

###########################  Correlation test "Spearman"  ######################

cor.test(Auxotrophy_12$Completeness, Auxotrophy_12$count, method = "spearman", exact = FALSE) 

############################ visualization #####################################

c <- ggplot(Auxotrophy_12, aes(Completeness, count)) +
  geom_point()+
  geom_smooth() +
  xlab("Completeness [%]") +
  ylab("Number of Auxotrophies")  +
  theme_minimal() +
  xlim(50,100) +
  ylim(0,21)+
  theme(axis.title.x = element_text(size=20, colour = "black", margin = margin(10,0,0,0)),
        axis.title.y = element_text(size = 20, colour = "black", margin = margin(0,10,0,0)),
        axis.text.x = element_text(size=18, colour="black"),
        axis.text.y = element_text(size=18, colour="black"))
c

ggsave("/Users/Svenja/Desktop/Completeness_NumbAuxos.pdf", plot = c,
              width = 6, height = 5)





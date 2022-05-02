#######################  number of auxotrophies and genome length ##############


########## visualization ############
length_count <- ggplot(Auxotrophy_12, aes(x =`Genome Length (bp)`, y= count)) +
  geom_point(size = 0.1) +
  geom_smooth() +
  xlab("Genome Length (bp)") +
  ylab("Number of auxotrophies") +
  theme_minimal()+
  theme(axis.title.x = element_text(color = "black", size = 10, face = "bold", margin = margin(t = 10, r = 0, b= 0, l = 0))) +
  theme(axis.title.y = element_text(colour = "black", size = 10, face = "bold", margin = margin(t =0, r = 10, b= 0, l = 0))) +
  theme(axis.text.x = element_text(colour = "black", size = 9)) +
  theme(axis.text.y = element_text(colour = "black", size = 9)) +
  
  ylim(0,21)+
  coord_cartesian()
length_count 

ggsave("output/plots/Genome_length_numb_auxos.pdf", plot = length_count,
       width = 6, height = 4)

################   Abundancies of auxotrophies in the gut   ####################
data <- fread("/Users/Svenja/Downloads/FOCUS_HRGM_abundancies_2")

#filter for BL
data <- data[focus.call == "BL"]
View(data)

#filter
Auxotrophy_2[,c(4:17)] <- NULL

sub <- unique(data$subject)
relAA <- unique(Auxotrophy_2$Compound)
Auxotrophy_2
p <- list()
k <- 1

for (subi in sub) {
  print(subi) 
    for (AAi in relAA) {
      x <- data[subject == subi]
      y <- Auxotrophy_2[Compound == AAi]
      z <- merge(x,y, by.x = "model", by.y = "Genomes")
      t <- z[Prototrophy == 0]
      p[[k]] <- t
      k <- k +1
  }
}

u <- rbindlist(p) 
sumfreq <- aggregate(u$freq, by=list(subject=u$subject, AA=u$Compound), FUN=sum)
sumfreqAA <- aggregate(sumfreq$x, by=list(Aminoacid = sumfreq$AA), FUN=median)

#####################        visualization        ##############################
a <- ggplot( sumfreqAA, aes(Aminoacid, x)) +
  geom_bar(stat="identity", fill = "#7E0018DD",colour = "black", width= 1) +
  ylab("Abundance of Auxotrophy") +
  xlab("Amino Acids") +
  theme(axis.line = element_line(size=0.2, colour = "black")) +
  theme(panel.background = element_rect(fill="white", colour= "white")) +
  theme(axis.title.y = element_text(colour = "black", size = 16, face = "bold", margin = margin(0,10,0,0)))+
  theme(axis.title.x = element_text(colour = "black", size = 16, face = "bold", margin = margin(10,0,0,0))) +
  theme(axis.text.x = element_text(size=14, colour = "black", hjust = 1,  angle = 45, margin = margin(10,0,0,0))) +
  theme(axis.text.y = element_text(size = 14, colour = "black")) +
  theme(plot.margin= margin(0.5,0.5,0.5,0.5, "cm")) +
  coord_cartesian()

a

ggsave("/Users/Svenja/Desktop/barplot_frequence_gut.pdf", plot = a,
       width = 6, height = 4)


#################### boxplot with summary of statistical anyslsis ##############
#standard error
 r <- ggplot(sumfreq, aes(x=AA, y=x)) + geom_bar(stat="summary", fun="mean", fill = "#fdae6b",colour = "black", width= 1) + 
  geom_errorbar(stat="summary", fun.data="mean_se") +
  ylab("Abundance of Auxotrophy")+
  xlab("Amino Acids") +
  theme(axis.line = element_line(size=0.2, colour = "black")) +
  theme(panel.background = element_rect(fill="white", colour= "white")) +
  theme(axis.title.y = element_text(colour = "black", size = 16, face = "bold", margin = margin(0,10,0,0)))+
  theme(axis.title.x = element_text(colour = "black", size = 16, face = "bold", margin = margin(10,0,0,0))) +
  theme(axis.text.x = element_text(size=14, colour = "black", hjust = 1,  angle = 45, margin = margin(10,0,0,0))) +
  theme(axis.text.y = element_text(size = 14, colour = "black")) +
  theme(plot.margin= margin(0.5,0.5,0.5,0.5, "cm"))
  
ggsave("/Users/Svenja/Desktop/Standarderror_frequency_gut.pdf", plot = r,
       width = 6, height = 5)

#boxplot
b <- ggplot(sumfreq, aes(AA, x)) +
  geom_boxplot( outlier.shape = NA, width = 0.8) +
  ylab("Abundance of Auxotrophy")+
  xlab("Amino Acids") +
  theme(axis.line = element_line(size=0.2, colour = "black")) +
  theme(panel.background = element_rect(fill="white", colour= "white")) +
  theme(axis.title.y = element_text(colour = "black", size = 16, face = "bold", margin = margin(0,10,0,0)))+
  theme(axis.title.x = element_text(colour = "black", size = 16, face = "bold", margin = margin(10,0,0,0))) +
  theme(axis.text.x = element_text(size=14, colour = "black", hjust = 1,  angle = 45, margin = margin(10,0,0,0))) +
  theme(axis.text.y = element_text(size = 14, colour = "black")) +
  theme(plot.margin= margin(0.5,0.5,0.5,0.5, "cm"))

ggsave("/Users/Svenja/Desktop/boxplot_frequency_gut.pdf", plot = b,
       width = 6, height = 4)


##add standard deviation
c <- b + stat_summary(fun = mean, colour = "red", geom = "point") +
  stat_summary(fun.min = function(x) mean(x)- sd(x),
               fun.max = function(x) mean(x) + sd(x),
               geom = "errorbar",
               colour = "red",
               width = .3)

ggsave("/Users/Svenja/Desktop/boxplot_standarddeviation_frequency_gut.pdf", plot = c,
       width = 6, height = 4)



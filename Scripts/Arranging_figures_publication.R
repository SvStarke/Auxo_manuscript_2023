######################### arranging the figures for publication ################
######figure 1 #############
#done with Flowchartdesigner 3

############# figure 2##########
library(gridExtra)
library(cowplot)
fi <- grid.arrange(arrangeGrob(pt,o, ncol=2),
             nrow=2,
             abun)
fi2 <- as_ggplot(fi) +
  draw_plot_label(label= c("A", "B","C"), size =12,
                  x = c(0,0.5,0), y = c(1, 1, 0.5))
fi2

ggsave("output/plots/figure2.pdf", plot = fi2,
       width = 10, height = 8)
#################figure 3 #################
#completeness of the pathways
fi3 <- ggarrange(tr,hi,ch,se,le,il1,il2,il3,va,il4,il5,
                                 labels = c("A","B","C", "D","E"),
                                 hjust = c(-0.1,0.5,0.5,0.5,-0.1),
                                 ncol=4, nrow= 3, common.legend = TRUE, legend = "bottom",
                                     heights = c(1,1,1,1), widths= c(1,1,1,1))
fi3

ggsave("output/plots/figure3.pdf", plot = fi3,
       width = 9, height = 8)


################# figure 4 ###########
fi4 <- ggarrange(b, linear_health_BMI_age1,all2, ncol=3,
                   nrow=1, heights = c(1,1,1), widths= c(3,1,1.4),
                 labels = "AUTO", hjust = c(-0.5,1,-3))
fi4
ggsave("output/plots/figure4.pdf", plot = fi4,
       width = 9, height = 5)


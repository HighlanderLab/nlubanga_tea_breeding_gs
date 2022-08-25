# Merge results from 25% and 100% parent replacement scenarios. 
load("PlotProgress_finalPlots/Data25.RData")
load("PlotProgress_finalPlots/Data100.RData")

# Merge all line plots
a.25 <- a.25+theme(legend.position = "none")
b.25 <- b.25+theme(legend.position = "none")
c.25 <- c.25+theme(legend.position = "none")
summary25 <- ggarrange(a.25,b.25,c.25, nrow = 1, ncol=3, common.legend = T)
summary25 <- annotate_figure(summary25, top = text_grob("25% parents replaced", color = "black", face = "bold", size = 14))
a <- a+theme(legend.position = "none")
b <- b+theme(legend.position = "none")
c <- c+theme(legend.position = "none")
summary100 <- ggarrange(a,b,c, nrow = 1, ncol=3, common.legend = T)
summary100 <- annotate_figure(summary100, top = text_grob("100% parents replaced", color = "black", face = "bold", size = 14))

print(all <- ggarrange(summary25,summary100, nrow = 2, ncol=1, common.legend = TRUE, legend="bottom"))
ggsave(filename = paste("PlotProgress_finalPlots/Summary_all.png"), plot = all, width = 12, height = 8, scale = 1.1)

# Merge genetic gain plots
a.25 <- a.25 + ggtitle("25% parents replaced")
a <- a + ggtitle("100% parents replaced") + scale_y_continuous("",limits = c(4000,10500))
(gain <- ggarrange(a.25,a, nrow = 1, ncol=2, common.legend = T, legend="bottom"))
ggsave(filename = paste("PlotProgress_finalPlots/Summary_1gain.png"), plot = gain , width = 12, height = 6 , scale = 1.1)

# Merge genetic gain boxplots
a1.25 <- a1.25 + ggtitle("25% parents replaced")
a1 <- a1 + ggtitle("100% parents replaced")
print(gainbox <- ggarrange(a1.25,a1, nrow = 1, ncol=2, common.legend = F))
ggsave(filename = paste("PlotProgress_finalPlots/Summary_1gainbox.png"), plot = gainbox, width = 12, height = 6, scale = 1.1)

# Merge genetic variance plots
b.25 <- b.25 + ggtitle("25% parents replaced")
b <- b + ggtitle("100% parents replaced") + scale_y_continuous("",limits=c(0,140000))
(variance <- ggarrange(b.25,b, nrow = 1, ncol=2, common.legend = T, legend="bottom"))
ggsave(filename = paste("PlotProgress_finalPlots/Summary_2variance.png"), plot = variance , width = 12, height = 6 , scale = 1.1)

# Merge genetic accuracy plots
c.25 <- c.25 + ggtitle("25% parents replaced")
c <- c + ggtitle("100% parents replaced") + scale_y_continuous("",limits=c(0,1))
(accuracy <- ggarrange(c.25,c, nrow = 1, ncol=2, common.legend = T, legend="bottom"))
ggsave(filename = paste("PlotProgress_finalPlots/Summary_3accuracy.png"), plot = accuracy , width = 12, height = 6 , scale = 1.1)





library(ggplot2)
library(reshape2)
data=read.csv("MAG_pontential_modified.txt",header = T,sep='\t')
df=melt(data)
df$MAG=factor(
  df$MAG,ordered=T,
  levels = c(
    "Cadiz_LYX104","Cadiz_LYX52","Cadiz_LYX97","Cadiz_LYX36",
    "Cadiz_LYX64","Cadiz_LYX80","Cadiz_LYX1","Cadiz_LYX88",
    "Cadiz_LYX94","Cadiz_LYX81","Cadiz_LYX72","Cadiz_LYX26",
    "Cadiz_LYX31","Cadiz_LYX47","Cadiz_LYX43","Cadiz_LYX7",
    "Cadiz_LYX42","Cadiz_LYX74","Cadiz_LYX86","Cadiz_LYX46",
    "Cadiz_LYX58","Cadiz_LYX62","Cadiz_LYX49","Cadiz_LYX30",
    "Cadiz_LYX59","Cadiz_LYX40","Cadiz_LYX20","Cadiz_LYX105",
    "Cadiz_LYX10","Cadiz_LYX28","Cadiz_LYX76","Cadiz_LYX79",
    "Cadiz_LYX87","Cadiz_LYX93","Cadiz_LYX68","Cadiz_LYX17",
    "Cadiz_LYX6","Cadiz_LYX22","Cadiz_LYX39","Cadiz_LYX65",
    "Cadiz_LYX96","Cadiz_LYX67","Cadiz_LYX33"
  )
  )

ggplot(data = df, mapping = aes(x = variable, y = MAG, colour = value)) +
  geom_point()+
  theme_bw()+
  scale_color_gradientn(values = seq(0,1,0.5),
                        colours = c('white','gray','black'))+
  xlab("")+ylab("")+
  theme(
    legend.position = "none",
    axis.text.x = element_text(vjust = 0.5, hjust = 0.5, angle = 90),
    axis.text.y = element_text(vjust = 0, hjust = 0, angle = 0)
    )+
  theme(panel.grid =element_blank()) +  
  theme(axis.text = element_blank()) +   
  theme(axis.ticks = element_blank()) +   
  theme(panel.border = element_blank())
ggsave("Fig S2.pdf",units="in",width = 7.5,dpi=300)

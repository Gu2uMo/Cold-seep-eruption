library(ggplot2)
library(reshape2)
data=read.table("Fig 2 Bacteria_genus.txt",header = T)
df=melt(data)
colnames(df)=c("Genus","fac","Sample","relative abundance")
(p_b=ggplot(data = df, mapping = aes(x = Sample, y = Genus, colour = fac, size = `relative abundance`)) +
  geom_point()+
  theme_classic()+
  scale_x_discrete(labels=c("Origin","L8","L15","L30","L8II","H15","H30"))+
  scale_y_discrete(limits=c("Sulfurovum","Desulfosarcina","uncultured_Desulfosarcinaceae",
                            "Desulfuromonas","Izimaplasma","Proteiniclasticum",
                            "Dethiosulfatibacter","Fusibacter","Anaerovirgula","JTB215",
                            "uncultured_Desulfitobacteriales","SG8-4","uncultured_Rhizobiaceae",
                            "Marinobacter","Cupriavidus","Thiobacillus","Methylophaga",
                            "Halomonas","Kangiella","Pseudomonas","uncultured_Pseudomonadaceae",
                            "Others"))+
  theme(
        axis.title=element_text(colour='black', size=10),
        axis.text=element_text(colour='black',size=10),
        legend.text=element_text(size=10),
        legend.title=element_text(size=10),
        legend.position = "left"
        )+
  scale_radius()+scale_color_discrete(guide="none")
)

data=read.table("Fig 2 Archaea_genus.txt",header = T)
df=melt(data)
colnames(df)=c("Genus","Sample","relative abundance")
(p_a=ggplot(data = df, mapping = aes(x = Sample, y = Genus, colour = Genus,size = `relative abundance`)) +
  geom_point()+
  theme_classic()+
  scale_x_discrete(labels=c("Origin","L8","L15","L30","H8","H15","H30"))+
  scale_y_discrete(position = "right")+  
  theme(
        axis.title=element_text(colour='black', size=10),
        axis.text=element_text(colour='black',size=10),
        legend.text=element_text(size=10),
        legend.title=element_text(size=10)
  )+
  scale_radius()+scale_color_discrete(guide="none")
)

library(cowplot)
plot_grid(p_b, p_a, labels = c("Bacteria", "Archaea"))
ggsave(file="Fig 2.pdf",dpi=300,width = 15, units="in",limitsize = FALSE)




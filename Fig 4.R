library(ggrepel)
library(vegan)
library(ggplot2)
library(scales)
library(ade4)
library(tidyverse)

setwd("D:/workspace/MyNote/Experiment_record/by_project/2020-03-Project_Cadiz/source_data/06_MAG")
data=read.csv('abundant_MAG.txt',row.names = 1,sep='\t',header = T)

KO=read.csv('../07_TPM/07_MAG_module_TPM.txt',row.names = 1,sep='\t',header = T)
KO <- KO %>% filter(apply(KO,1,max) > 1000)

####PCoA
pcoa_dist=vegdist(t(data),method='euclidean')
pcoa<- dudi.pco(pcoa_dist, scan = FALSE,nf=3)
pcoa_eig <- (pcoa$eig)[1:2] / sum(pcoa$eig)

#提取样本点坐标（前两轴）
sample_site <- data.frame(data.frame({pcoa$li})[1:2])
sample_site$names <- rownames(sample_site)
names(sample_site)[1:2] <- c('PCoA1', 'PCoA2')

vec.pcoa<-envfit(sample_site[,c(1:2)],t(KO), perm=10000)
vec.pcoa.df<-as.data.frame(
  vec.pcoa$vectors$arrows*sqrt(vec.pcoa$vectors$r)
)
vec.pcoa.df$species<-rownames(vec.pcoa.df)

m1=as.data.frame(vec.pcoa$vectors$arrows)
m2=as.data.frame(vec.pcoa$vectors$pvals)
m_total=cbind(m1,m2)
m_total$species<-rownames(vec.pcoa.df)
cor_df=as.data.frame(m_total[m_total$`vec.pcoa$vectors$pvals`<0.01,])


module=read.table('../07_TPM/KEGG_module.tsv',row.names = 3,header = F,sep='\t')
cor_df$`KEGG_Module`=module[rownames(cor_df),]$V2
cor_df$`KEGG_Module2`=module[rownames(cor_df),]$V3
cor_df$`KEGG_Module3`=module[rownames(cor_df),]$V4

sample_site$PCoA1s=rescale_max(sample_site$PCoA1,to=c(-0.75,0.75))
sample_site$PCoA2s=rescale_max(sample_site$PCoA2,to=c(-0.75,0.75))
sample_site$sample=c("Origin","8MPa","15MPa","30MPa","8MPaII","15MPaII","30MPaII")
cor_df$PCoA1s=cor_df$PCoA1
cor_df$PCoA2s=cor_df$PCoA2

ggplot(sample_site, aes(PCoA1s, PCoA2s)) +
  theme_classic()+#去掉背景框
  geom_point(size = 5)+  #可在这里修改点的透明度、大小
  theme(panel.grid = element_line(color = 'gray', linetype = 2, size = 0.1), 
        panel.background = element_rect(color = 'black', fill = 'transparent'), 
        legend.title=element_blank()
  )+
  labs(x = paste('PCoA1: ', round(100 * pcoa_eig[1], 2), '%'), 
       y = paste('PCoA2: ', round(100 * pcoa_eig[2], 2), '%'))+
  geom_label_repel(data=sample_site,aes(label=sample))+
  geom_label_repel(data=cor_df,aes(label=species),label.size = 1)+
  geom_segment(data=cor_df,aes(x=0,xend=PCoA1,
                               y=0,yend=PCoA2,colour=KEGG_Module),
               arrow = arrow(length = unit(0.5, "cm")),
               linetype=1)
  

ggsave("PCoA_MAG_ab_moduleII.tiff",dpi=300)

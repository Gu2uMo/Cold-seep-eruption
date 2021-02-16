library(readxl)
library(grid)
library(reshape2)
library(ggplot2)
library("lubridate")
library("scales")
library("magrittr")
library("tidyr")
# cell numbers

data=read_excel("Table S1.concentration and cell number.xlsx",sheet="cell number",na="NA")
data$bac=data$`Bacteria (cell number/ ml)`/1000000
data$bsd=data$Bacteria_standard_deviation/1000000
data$arc=data$`Archaea  (cell number/ ml)`/50000
data$asd=data$Archaea_standard_deviation/50000

data$srb=data$`SRB_absolute_abundance (cell number/ ml)`/1000000
data$srb_sd=data$SRB_standard_deviation/1000000

data$ANME=data$`ANME_absolute_abundance  (cell number/ ml)`/50000
data$anme_sd=data$ANME_standard_deviation/50000

data$Met=data$`Methanolobus_absolute_abundance  (cell number/ ml)`/50000
data$Met_sd=data$Methanolobus_standard_deviation/50000

(p1=ggplot()+
    #bac
    geom_point(data = data,aes(x = factor(Sample),y = bac,group=1,colour="Bacteria"),size=1.5)+
    geom_errorbar(data = data,aes(x = factor(Sample),y = bac,ymin=bac-bsd, ymax=bac+bsd,colour="Bacteria"), width=.1) +
    #arc
    geom_point(data = data,aes(x = factor(Sample),y = arc ,group=1,colour="Archaea"),size=1.5)+
    geom_errorbar(data = data,aes(x = factor(Sample),y = arc,ymin=arc-asd, ymax=arc+asd,colour="Archaea"), width=.1) +
    #SRB
    geom_point(data = data,aes(x = factor(Sample),y =srb,group=1,colour="SRB"),size=1.5)+
    geom_errorbar(data = data,aes(x = factor(Sample),y = srb,ymin=srb-srb_sd, ymax=srb+srb_sd,colour="SRB"), width=.1) +
    #ANME
    geom_point(data = data,aes(x = factor(Sample),y =ANME,group=1,colour="ANME"),size=1.5)+
    geom_errorbar(data = data,aes(x = factor(Sample),y = ANME,ymin=ANME-anme_sd, ymax=ANME+anme_sd,colour="ANME"), width=.1)+
    #Mehanolobus
    geom_point(data = data,aes(x = factor(Sample),y =Met,group=1,colour="Methanolobus"),size=1.5)+
    geom_errorbar(data = data,aes(x = factor(Sample),y = Met,ymin=Met-Met_sd, ymax=Met+Met_sd,colour="Methanolobus"), width=.1)+
    #axis
    scale_y_continuous(limits=c(0,3),expand=c(0,0),
                       breaks=pretty_breaks(5),                           
                       sec.axis = sec_axis( ~rescale(.,c(0,1.5)),
                                            name = "Archaea 10x5 cell numbers / ml")
    )+
    theme_bw()+
    theme(plot.title = element_text(hjust = 0.5))+
    theme(axis.title.x = element_text(size = 10),axis.title.y = element_text(size = 10))+
    xlab("sample")+
    ylab("Bacteria 10x6 cell numbers / ml")+
    labs(title="cell number changes among samples")+
    scale_x_discrete(limits=c("Origin","L8","L15","L30","L8II", "H15", "H30"))
)

## acitivity
data=read_excel("Table S1.concentration and cell number.xlsx",sheet="key compounds concentration",na="NA")
######sulfide
sulfide=data[,c(1,2,3)]
colnames(sulfide)=c("days","Sample","value")
(p_sulfide=ggplot(sulfide, aes(x=days, y=value,group=sulfide$Sample))+ 
    geom_point()+
    geom_smooth(method = "lm",aes(group=sulfide$Sample))+
    theme_bw()+xlab('day')+ylab('rate (μmol/day)')+
    geom_vline(aes(xintercept=60), colour="#990000", linetype="dashed")+
    geom_vline(aes(xintercept=120), colour="#990000", linetype="dashed")+
    geom_vline(aes(xintercept=240), colour="#990000", linetype="dashed")+
    geom_vline(aes(xintercept=180), colour="#990000", linetype="dashed")+
    geom_vline(aes(xintercept=300), colour="#990000", linetype="dashed")+
    geom_vline(aes(xintercept=360), colour="#990000", linetype="dashed")+
    theme_bw()+
    annotate("text", x=30 , y=200,label='L8') +
    annotate("text", x=90 , y=200,label='L15') +
    annotate("text", x=150 , y=200,label='L30') +
    annotate("text", x=210 , y=200,label='L8II') +
    annotate("text", x=270 , y=200,label='H15') +
    annotate("text", x=330 , y=200,label='H30') +
    annotate("text", x=30 , y=155,label='μ=0.993') +
    annotate("text", x=90 , y=155,label='μ=2.297') +
    annotate("text", x=150 , y=155,label='μ=-12.264') +
    annotate("text", x=210 , y=155,label='μ=1.591') +
    annotate("text", x=270 , y=155,label='μ=8.208') +
    annotate("text", x=330 , y=155,label='μ=3.536') +
    labs(title='Production of Sulfide')+
    theme(plot.title = element_text(hjust = 0.5))+
    scale_x_continuous(breaks=seq(0, 360, 60))  
)


#######acetate
acetate=data[,c(1,2,4)]
colnames(acetate)=c("days","Sample","value")

(p_acetate=ggplot(acetate, aes(x=days, y=value,group=acetate$Sample))+ 
    geom_point()+
    geom_smooth(method = "lm",aes(group=acetate$Sample))+
    theme_bw()+xlab('day')+ylab('concentration (μM)')+
    geom_vline(aes(xintercept=60), colour="#990000", linetype="dashed")+
    geom_vline(aes(xintercept=120), colour="#990000", linetype="dashed")+
    geom_vline(aes(xintercept=240), colour="#990000", linetype="dashed")+
    geom_vline(aes(xintercept=180), colour="#990000", linetype="dashed")+
    geom_vline(aes(xintercept=300), colour="#990000", linetype="dashed")+
    geom_vline(aes(xintercept=360), colour="#990000", linetype="dashed")+
    theme_bw()+
    annotate("text", x=30 , y=7,label='L8') +
    annotate("text", x=90 , y=7,label='L15') +
    annotate("text", x=150 , y=7,label='L30') +
    annotate("text", x=210 , y=7,label='L8II') +
    annotate("text", x=270 , y=7,label='H15') +
    annotate("text", x=330 , y=7,label='H30') +
    annotate("text", x=30 , y=6,label='μ=2.559') +
    annotate("text", x=90 , y=6,label='μ=1.390') +
    annotate("text", x=150 , y=6,label='μ=1.680') +
    annotate("text", x=210 , y=6,label='μ=2.430') +
    annotate("text", x=270 , y=6,label='μ=1.229') +
    annotate("text", x=330 , y=6,label='μ=0.697') +
    labs(title='Concentration of acetate')+
    theme(plot.title = element_text(hjust = 0.5))+
    scale_x_continuous(breaks=seq(0, 360, 60))  
)





## merge all plot
library(cowplot)

plot_grid(p_sulfide, p_acetate,p1,
          labels = c("A", "B","C"),nrow=3)
ggsave(file="Fig 1.pdf",dpi=300,width = 7.3,units="in",limitsize = FALSE)

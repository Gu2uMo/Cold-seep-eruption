## 导入数据
library(tidyverse)
data <- read.table("Fig 5 MAG_module_TPM.txt",header = TRUE,row.names = 1,sep = "\t")
group <- read.table("Fig 5 Group.txt",header = FALSE,sep = "\t")
data <- data/10000
data <- data %>% filter(apply(data,1,max) > 0.05)
data <- t(data)
data1 <- data.frame(data,group$V2)
colnames(data1) <- c(colnames(data),"Group")
data1$Group <- as.factor(data1$Group)


diff <- data1 %>% 
  select_if(is.numeric) %>%
  map_df(~ broom::tidy(t.test(. ~ Group,data = data1)), .id = 'var')

diff$p.value <- p.adjust(diff$p.value,"none")
diff <- diff %>% filter(p.value < 0.1)


## 绘图数据构建
## 左侧条形图
abun.bar <- data1[,c(diff$var,"Group")] %>% 
  gather(variable,value,-Group) %>% 
  group_by(variable,Group) %>% 
  summarise(Mean = mean(value))

## 右侧散点图
diff.mean <- diff[,c("var","estimate","conf.low","conf.high","p.value")]
diff.mean$Group <- c(ifelse(diff.mean$estimate >0,levels(data1$Group)[1],
                            levels(data1$Group)[2]))
diff.mean <- diff.mean[order(diff.mean$estimate,decreasing = TRUE),]

## 左侧条形图
library(ggplot2)
cbbPalette <- c("#E69F00", "#56B4E9")
abun.bar$variable <- factor(abun.bar$variable,levels = rev(diff.mean$var))
p1 <- ggplot(abun.bar,aes(variable,Mean,fill = Group)) +
  scale_x_discrete(limits = levels(diff.mean$var)) +
  coord_flip() +
  xlab("") +
  ylab("Mean proportion (%)") +
  theme(panel.background = element_rect(fill = 'transparent'),
        panel.grid = element_blank(),
        axis.ticks.length = unit(0.4,"lines"), 
        axis.ticks = element_line(color='black'),
        axis.line = element_line(colour = "black"),
        axis.title.x=element_text(colour='black', size=12,face = "bold"),
        axis.text=element_text(colour='black',size=10,face = "bold"),
        legend.title=element_blank(),
        legend.text=element_text(size=12,face = "bold",colour = "black",
                                 margin = margin(r = 20)),
        legend.position = c(-1,-0.1),
        legend.direction = "horizontal",
        legend.key.width = unit(0.8,"cm"),
        legend.key.height = unit(0.5,"cm"))


for (i in 1:(nrow(diff.mean) - 1)) 
  p1 <- p1 + annotate('rect', xmin = i+0.5, xmax = i+1.5, ymin = -Inf, ymax = Inf, 
                      fill = ifelse(i %% 2 == 0, 'white', 'gray95'))

p1 <- p1 + 
  geom_bar(stat = "identity",position = "dodge",width = 0.7,colour = "black") +
  scale_fill_manual(values=cbbPalette)


## 右侧散点图
diff.mean$var <- factor(diff.mean$var,levels = levels(abun.bar$variable))
diff.mean$p.value <- signif(diff.mean$p.value,3)
diff.mean$p.value <- as.character(diff.mean$p.value)
p2 <- ggplot(diff.mean,aes(var,estimate,fill = Group)) +
  theme(panel.background = element_rect(fill = 'transparent'),
        panel.grid = element_blank(),
        axis.ticks.length = unit(0.4,"lines"), 
        axis.ticks = element_line(color='black'),
        axis.line = element_line(colour = "black"),
        axis.title.x=element_text(colour='black', size=12,face = "bold"),
        axis.text=element_text(colour='black',size=10,face = "bold"),
        axis.text.y = element_blank(),
        legend.position = "none",
        axis.line.y = element_blank(),
        axis.ticks.y = element_blank(),
        plot.title = element_text(size = 15,face = "bold",colour = "black",hjust = 0.5)) +
  scale_x_discrete(limits = levels(diff.mean$var)) +
  coord_flip() +
  xlab("") +
  ylab("Difference in mean proportions (%)") +
  labs(title="95% confidence intervals") 

for (i in 1:(nrow(diff.mean) - 1)) 
  p2 <- p2 + annotate('rect', xmin = i+0.5, xmax = i+1.5, ymin = -Inf, ymax = Inf, 
                      fill = ifelse(i %% 2 == 0, 'white', 'gray95'))

p2 <- p2 +
  geom_errorbar(aes(ymin = conf.low, ymax = conf.high), 
                position = position_dodge(0.8), width = 0.5, size = 0.5) +
  geom_point(shape = 21,size = 3) +
  scale_fill_manual(values=cbbPalette) +
  geom_hline(aes(yintercept = 0), linetype = 'dashed', color = 'black')


p3 <- ggplot(diff.mean,aes(var,estimate,fill = Group)) +
  geom_text(aes(y = 0,x = var),label = diff.mean$p.value,
            hjust = 0,fontface = "bold",inherit.aes = FALSE,size = 3) +
  geom_text(aes(x = nrow(diff.mean)/2 +0.5,y = 0.85),label = "P-value (corrected)",
            srt = 90,fontface = "bold",size = 5) +
  coord_flip() +
  ylim(c(0,1)) +
  theme(panel.background = element_blank(),
        panel.grid = element_blank(),
        axis.line = element_blank(),
        axis.ticks = element_blank(),
        axis.text = element_blank(),
        axis.title = element_blank())

## 图像拼接
library(patchwork)
(p <- p1 + p2 + p3 + plot_layout(widths = c(4,6,2)))

## 保存图像
ggsave("stamp_new_20201225.pdf",p,width = 20,units="in",dpi=300)






################################################
################################################
### 作者：果子
### 更新时间：2020-11-20
### 微信公众号:果子学生信
### 私人微信：guotosky
### 个人博客: https://codingsoeasy.com/
### 个人邮箱：hello_guozi@126.com

### 本节任务：画火山图，订制火山图
rm(list = ls())
##用ggplot2
library(ggplot2)
library(ggrepel)
library(dplyr)
###3load(file = "output/diffgeneAD.Rdata")
### 无论是什么名称，都改为data
load(file = "output/allDiffAD.Rdata")
data <- allDiff
### 如果没有gene这一列就需要添加
data$gene <- rownames(data)
## 仔细观察data数据
## 如果是你自己的数据，至少有三列
## logFC，P.Value，gene
cut_off_pvalue = 0.0000001  #统计显著性
cut_off_logFC = 0.58           #差异倍数值

# 根据阈值参数，上调基因设置为‘up’，下调基因设置为‘Down’，无差异设置为‘Stable’，并保存到change列中
data$change = ifelse(data$P.Value < cut_off_pvalue & abs(data$logFC) >= cut_off_logFC, 
                        ifelse(data$logFC> cut_off_logFC ,'Up','Down'),
                        'Stable')
head(data)
p <- ggplot(
  # 数据、映射、颜色
  data, aes(x = logFC, y = -log10(P.Value), colour=change)) +
  geom_point(alpha=0.4, size=3.5) +
  scale_color_manual(values=c("#00CC99", "#d2dae2","#ff4757"))+
  # 辅助线
  geom_vline(xintercept=c(-0.58,0.58),lty=4,col="black",lwd=0.8) +
  geom_hline(yintercept = -log10(cut_off_pvalue),lty=4,col="black",lwd=0.8) +
  # 坐标轴
  labs(x="log2(fold change)",
       y="-log10 (p-value)")+
  theme_bw()+
  # 图例
  theme(plot.title = element_text(hjust = 0.5), 
        legend.position="right", 
        legend.title = element_blank())
p







###########################################
ggplot(data=data, aes(x=logFC, y =-log10(P.Value))) +
  ## 三个部分分别画点
  geom_point(data=subset(data,abs(data$logFC) <= 1),aes(size=abs(logFC)),color="black",alpha=0.1) +
  geom_point(data=subset(data,data$P.Value<0.05 & data$logFC > 0.58),aes(size=abs(logFC)),color="red",alpha=0.2) +
  geom_point(data=subset(data,data$P.Value<0.05 & data$logFC < -0.58),aes(size=abs(logFC)),color="green",alpha=0.2) +
  ## 画线
  geom_hline(yintercept = -log10(0.05),lty=4,lwd=0.6,alpha=0.8)+
  geom_vline(xintercept = c(0.58,-0.58),lty=4,lwd=0.6,alpha=0.8)+
  ## 主题
  theme_bw()+
  theme(panel.border = element_blank(),
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),   
        axis.line = element_line(colour = "black"))+
  labs(x="log2 (fold change)",y="-log10 (q-value)")+
  theme(plot.title = element_text(hjust = 0.5))+
  theme(legend.position='none')+
  ## 标签
  geom_text_repel(data=subset(data, abs(logFC) >1.5), aes(label=gene),col="black",alpha = 0.8)



library(ggplot2)
library(ggrepel)
data <- allDiff
data$gene <- rownames(data)
significant <- as.factor(data$adj.P.Val<0.05 & abs(data$logFC) > 0.58)
p <- ggplot(data=data, aes(x=logFC, y =-log10(adj.P.Val),color=significant)) +
  geom_point(alpha=0.8, size=1.2,col="black")+
  geom_point(data=subset(data, logFC > 0.58),alpha=0.8, size=1.2,col="red")+
  geom_point(data=subset(data, logFC < -0.58),alpha=0.6, size=1.2,col="blue")+
  labs(x="log2 (fold change)",y="-log10 (adj.P.Val)")+
  theme(plot.title = element_text(hjust = 0.4))+
  geom_hline(yintercept = -log10(0.05),lty=4,lwd=0.6,alpha=0.8)+
  geom_vline(xintercept = c(0.58,-0.58),lty=4,lwd=0.6,alpha=0.8)+
  theme_bw()+
  theme(panel.border = element_blank(),
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),   
        axis.line = element_line(colour = "black")) +
  geom_point(data=subset(data, abs(logFC) >= 3),alpha=0.8, size=3,col="green")+
  geom_text_repel(data=subset(data, abs(logFC) > 3), 
                  aes(label=gene),col="black",alpha = 0.8)
p + coord_cartesian(xlim = c(-1.5,1.5))
###########################################################
##换一种风格
library(ggplot2)
library(ggrepel)
data <- allDiff
data$gene <- rownames(data)
data$significant <- as.factor(data$adj.P.Val<0.05 & abs(data$logFC) > 1)
data$gene <- rownames(data)
ggplot(data=data, aes(x=logFC, y =-log10(adj.P.Val),color=significant)) +
  geom_point(alpha=0.8, size=1.2,col="black")+
  geom_point(data=subset(data, logFC > 0.58),alpha=0.8, size=1.2,col="red")+
  geom_point(data=subset(data, logFC < -0.58),alpha=0.6, size=1.2,col="blue")+
  labs(x="log2 (fold change)",y="-log10 (adj.P.Val)")+
  theme(plot.title = element_text(hjust = 0.4))+
  geom_hline(yintercept = -log10(0.05),lty=4,lwd=0.6,alpha=0.8)+
  geom_vline(xintercept = c(0.58,-0.58),lty=4,lwd=0.6,alpha=0.8)+
  theme_bw()+
  theme(panel.border = element_blank(),
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),   
        axis.line = element_line(colour = "black")) +
  geom_point(data=subset(data, abs(logFC) >= 3),alpha=0.8, size=3,col="green")+
  geom_text_repel(data=subset(data, abs(logFC) > 3), 
                  aes(label=gene),col="black",alpha = 0.8)
## 加载R包
library(export)
## 导成PPT可编辑的格式
graph2ppt(file="output/volcano.pptx")
graph2pdf(file="output/volcano.pdf")

## GEO教程长期更新的链接是这个:
## https://codingsoeasy.com/archives/geo
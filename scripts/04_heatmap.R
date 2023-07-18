################################################
################################################
### 作者：果子
### 更新时间：2020-11-20
### 微信公众号:果子学生信
### 私人微信：guotosky
### 个人博客: https://codingsoeasy.com/
### 个人邮箱：hello_guozi@126.com

### 本节任务：画热图，以及掌握热图的用法
###############
####heatmap热图
### 热图就是表达量数据，行是差异基因，列是样本
### 想一想该如何提取

### 用行名提取数据
rm(list = ls())
## 1.加载表达数据
load(file = "output/exprSet_rmdupAD.Rdata")
## 2.加载差异列表
load(file = "output/allDiffAD.Rdata")
library(dplyr)
library(tibble)
diffgene <- allDiff %>% 
  rownames_to_column() %>% 
  filter(adj.P.Val < 0.05) %>% 
  filter(abs(logFC) >0.58) %>% 
  column_to_rownames()

## 3.用名称提取部分数据用作热图绘制
######heatdata <- exprSet[rownames(diffgene),]
heatdata<-exprSet[c("XK",
                     "PPP3CB",
                     "PPP3CA",
                     "TSPO",
                     "YWHAE",
                     "ATP1B1",
                     "FZD9",
                     "ITPKB",
                     "GPER1",
                     "SNX10",
                     "GRINA",
                     "ITPR1",
                     "ATP2B3",
                     "CCKBR",
                     "PTH1R",
                     "ATP2A2",
                     "NELL2",
                     "TMEM178A",
                     "NOTCH1",
                     "ATF4"),]
## 4.制作一个分组信息用于注释
######group <- c(rep("con",3),rep("treat",3)) 
group <- c(rep("con",81),rep("treat",80),"con","con","con","con","treat","treat","treat","treat","treat","treat","treat","treat","treat","treat","con","treat","treat","treat","treat","treat","treat","treat","con","con","con","con","con","con","con","con","con","con","con","con") 
annotation_col <- data.frame(group)
rownames(annotation_col) <- colnames(heatdata)

## 加载热图的R包
library(pheatmap)
## 颜色包 viridisLite
library(viridisLite)
## 直接作图
pheatmap(heatdata)

### 经过修饰的图
pheatmap(heatdata, #热图的数据
         cluster_rows = TRUE,#行聚类
         cluster_cols = F,#列聚类，可以看出样本之间的区分度
         annotation_col =annotation_col, #标注样本分类
         annotation_legend=T, # 显示注释
         show_rownames = T,# 显示行名
         scale = "row", #以行来标准化，这个功能很不错
         color = viridis(10, alpha = 1, begin = 0.5, end = 1, direction = 1),#调色
         cellwidth = 3, # 格子宽度
         cellheight = 20,# 格子高度
         fontsize = 10)# 字体大小
install.packages("RColorBrewer")
#加载RColorBrewer这个R包
library(RColorBrewer)
coul <- colorRampPalette(brewer.pal(10, "Spectral"))(50)
pheatmap(heatdata, #热图的数据
         cluster_rows = F,#行聚类
         cluster_cols = F,#列聚类，可以看出样本之间的区分r度
         annotation_col =annotation_col, #标注样本分类
         annotation_legend=T, # 显示注释
         show_rownames = T,
         show_colnames = F, 
         # 显示行名
         scale = "row", #以行来标准化，这个功能很不错
         col = coul,
         cellwidth = 2, # 格子宽度
         cellheight = 20,# 格子高度
         fontsize = 13)# 字体大小


############################################################
### 可以重新筛选，阈值设大一点
diffgene <- allDiff %>% 
  rownames_to_column() %>% 
  filter(adj.P.Val < 0.05) %>% 
  filter(abs(logFC) >3) %>% 
  column_to_rownames()

## 用名称提取部分数据用作热图绘制
heatdata <- exprSet[rownames(diffgene),]
## 制作一个分组信息用于注释
group <- c(rep("con",3),rep("treat",3)) 
annotation_col <- data.frame(group)
rownames(annotation_col) <- colnames(heatdata)

#如果注释出界，可以通过调整格子比例和字体修正
pheatmap(heatdata, #热图的数据
         cluster_rows = TRUE,#行聚类
         cluster_cols = TRUE,#列聚类，可以看出样本之间的区分度
         annotation_col =annotation_col, #标注样本分类
         annotation_legend=TRUE, # 显示注释
         show_rownames = T,# 显示行名
         scale = "row", #以行来标准化，这个功能很不错
         color = viridis(10, alpha = 1, begin = 0.5, end = 1, direction = 1),#调色
         #filename = "heatmap_F.pdf",#是否保存
         cellwidth = 30, # 格子宽度
         cellheight = 12,# 格子高度
         fontsize = 10 # 字体大小
         )

################################################
### 导出图
## 加载R包
library(export)
## 导成PPT可编辑的格式
graph2ppt(file="output/heatmap_modified.pptx")
graph2pdf(file="output/heatmap_modified.pdf")

##################################################################
## one more thing!!
## 热图的强悍教程
## https://jokergoo.github.io/ComplexHeatmap-reference/book/
## GEO教程长期更新的链接是这个:
## https://codingsoeasy.com/archives/geo

################################################
################################################
### 作者：果子
### 更新时间：2020-11-20
### 微信公众号:果子学生信
### 私人微信：guotosky
### 个人博客: https://codingsoeasy.com/
### 个人邮箱：hello_guozi@126.com

### 本节任务：使用limma来做芯片的差异分析
#############
#加载limma包，用于校正和比较差异
rm(list = ls())
library(limma)
### 加载数据，注意解决报错
exprSet<-exprSet[,-c(1:36)]
save(exprSet,file = "output/exprSet_rmdup_yanzheng.Rdata")
load(file = "output/exprSet_rmdup_GSE132903.Rdata")

### 1.创建分组
### 这一步根据样本来就行，原则就是: 跟样本匹配，取决于样本的排序
group <- c(rep("con",81),rep("AD",80),"con","con","con","con","AD","AD","AD","AD","AD","AD","AD","AD","AD","AD","AD","AD","AD","AD","AD","AD","AD","AD","con","con","con","con","con","con","con","con","con","con","con","con") 
##group <- c("con","con","treat","con","treat","treat") 
### 分组变成向量，并且限定leves的顺序
### levels里面，把对照组放在前面
group <- c(rep("treat",56),rep("con",44))
group <-factor(group,levels = c("con","AD"))

### 主成分分析PCA：提前预测结果
### 行是样本列是基因
#which(!na.omit(exprSet))
#oopsmat<-t(exprSet)
#which(apply(oopsmat, 2, var)==0)
#pcadata = oopsmat[ , which(apply(oopsmat, 2, var) != 0)] 
res.pca <- prcomp(t(exprSet), scale = TRUE)
##res.pca <- PCA(t(exprSet), graph = FALSE)

library(factoextra)       
fviz_pca_ind(res.pca,col.ind = group) ##,col.ind = group,label="none"
fviz_pca_ind(res.pca,col.ind = group,label="none",
             addEllipses=TRUE, ellipse.level=0.95,
             palette = c("#999999", "#E69F00", "#56B4E9"))








### 构建比较矩阵
design <- model.matrix(~factor(group))
### 比较矩阵命名
colnames(design) <- levels(group)
design

### 2.线性模型拟合
fit <- lmFit(exprSet,design)
### 3.贝叶斯检验
fit2 <- eBayes(fit)
### 4.输出差异分析结果,其中coef的数目不能操过design的列数
### 此处的2代表的是design中第二列和第一列的比较
allDiff=topTable(fit2,adjust='fdr',coef=2,number=Inf) 
### 这个数据很重要需要保存一下
save(allDiff,file = "output/allDiffAD.Rdata")

###################################################################################
### 定义差异基因：差异倍数2倍，矫正后的p值小于0.05

library(dplyr)
diffgene <- allDiff %>% 
  filter(P.Value < 0.05) %>% 
  filter(abs(logFC) >0.58)

### 如果出现行名丢失的情况，需要先把行名变成列，处理完毕后再把列变成行名
### 这个工作是由tibble这个包里面的rownames_to_column()和column_to_rownames()完成的
library(tibble)
diffgene <- allDiff %>% 
  rownames_to_column() %>% 
  filter(adj.P.Val < 0.05) %>% 
  filter(abs(logFC) >1) %>% 
  column_to_rownames()

### 可选方案:使用subset直接获取,&是and的意思
diffgene <- subset(allDiff,abs(logFC) >1 & adj.P.Val < 0.05)
test <- allDiff[allDiff$adj.P.Val < 0.05 & abs(allDiff$logFC)>1,]
### 该数据也需要保存，此处一次性保存两个数据，如果是多个，一次写入变量名称即可。
save(diffgene,group,file = "output/diffgene_yanzheng.Rdata")
### 到此差异基因的分析就结束了
####################################################################################
####################################################################################
## 作图环节
## 1.把现在数据调整成可以作图的格式
### 这个技能是data wrangling部分重点掌握的技能
### 复习一下流程：输入数据是表达量，经过三步
### 1.探针ID转换，2.行列转置，3，添加分组信息。最终获得的是数据框

### 行列转置
exprSet <- as.data.frame(t(exprSet))
### 添加分组信息
dd <- cbind(group=group,exprSet)
### 截取部分展示,这就是清洁数据
test = dd[,1:10]

## 2.作图展示
library(ggplot2)
ggplot(data = dd,aes(x=group,y=PPP3CA,fill=group))+
  geom_boxplot()+
  geom_point()+
  theme_bw()

## 3.steal plot
my_comparisons <- list(
  c("treat", "con")
)
library(ggpubr)
ggboxplot(
  dd, x = "group", y = "PPP3CA",
  color = "group", palette = c("#00AFBB", "#E7B800"),
  add = "jitter"
)+
  stat_compare_means(comparisons = my_comparisons, method = "t.test")

## 改写成函数
diffplot <- function(gene){
  my_comparisons <- list(
    c("treat", "con")
  )
  library(ggpubr)
  ggboxplot(
    dd, x = "group", y = gene,
    color = "group", palette = c("#00AFBB", "#E7B800"),
    add = "jitter"
  )+
    stat_compare_means(comparisons = my_comparisons, method = "t.test")
}

diffplot("ATP2A2")
diffplot("ITPR1")

## 4.多个基因作图查看
## 先把基因提取出来
##############genelist <- rownames(diffgene)[1:9]
diffgene$gene<-rownames(diffgene)
aa<-diffgene[diffgene$gene=="ATF4",]
 bb<-diffgene[diffgene$gene=="ITPKB",]
 cc<-diffgene[diffgene$gene=="GRINA",]
 zz<-diffgene[diffgene$gene=="PPP3CA",]
 ee<-rbind(aa,bb,cc)
 genelist<-rownames(ee)
## 再提取表达量，使用名称选取行
data <- dd[,c("group",genelist)]
## 用pivot_longer调整数据，数据变长，增加的是行
library(tidyr)
data <- data %>% 
  pivot_longer(cols=-1,
               names_to= "gene",
               values_to = "expression")
## 多基因作图
## 作图
library(ggpubr)
ggplot(data = data,aes(x=gene,y=expression,fill=group))+
  geom_boxplot()+
  geom_jitter()+
  theme_bw(base_size = 13)+
  stat_compare_means(aes(group=group), label = "p.signif", method = "t.test")
## 尝试更清晰的展示
ggplot(data = data,aes(x=group,y=expression,fill=group))+
  geom_boxplot()+
  geom_jitter()+
  theme_bw()+
  facet_grid(.~gene)+
  stat_compare_means(comparisons = my_comparisons, label = "p.signif", method = "t.test")

## 图片导出
library(export)
## 导成PPT可编辑的格式
graph2ppt(file="output/diffgenboxplot.pptx")
## 其他自己想要的格式
graph2pdf(file="output/diffgenboxplot.pdf")
graph2tif(file="output/diffgenboxplot.tif")
## 导成AI可以编辑的状态
graph2eps(file="output/diffgenboxplot.eps")

## 问题，如果差异基因少怎么办？
## 此处休息15分钟
####################################################
## 推荐阅读
## GEO芯片分析中的大坑，差异基因完全相反！
## https://mp.weixin.qq.com/s/rPiYACna0tkwgX0xGzNQgQ
## GEO芯片如果超过了两组，也可以一次搞定差异分析
## https://mp.weixin.qq.com/s/YcTbLpMQ3gnCC3hXsWVmnw
## GEO芯片中配对样本如何做差异分析
## https://mp.weixin.qq.com/s/tz0CsJDumvnzY8WLqvoT-A
## 因子(factor)就像贤内助，让你始终分清主次，拨开云雾。
## https://mp.weixin.qq.com/s/d_DjXdxJydapb6_2KLwekg

## GEO的样本名称太多而且排序不规则，你们都是手动分组的么？
## https://mp.weixin.qq.com/s/m3hc4slyV9arw70kKPOClg
## GEO教程长期更新的链接是这个:
## https://codingsoeasy.com/archives/geo
####design
design<-model.matrix(~0+factor(group))
colnames(design)<-levels(factor(group))
rownames(design)<-colnames(exprSet)

################################################
################################################
### 作者：果子
### 更新时间：2020-11-20
### 微信公众号:果子学生信
### 私人微信：guotosky
### 个人博客: https://codingsoeasy.com/
### 个人邮箱：hello_guozi@126.com

### 本节任务：数据预处理，探针ID转换，探针去重

################################################
## 自动log化，解释
## 注意此时数据的格式是矩阵
ex <- exprSet
qx <- as.numeric(quantile(ex, c(0., 0.25, 0.5, 0.75, 0.99, 1.0), na.rm=T))
LogC <- (qx[5] > 100) ||
  (qx[6]-qx[1] > 50 && qx[2] > 0) ||
  (qx[2] > 0 && qx[2] < 1 && qx[4] > 1 && qx[4] < 2)

## 开始判断
if (LogC) { 
  ex[which(ex <= 0)] <- NaN
  ## 取log2
  exprSet <- log2(ex)
  print("log2 transform finished")
  }else{
    print("log2 transform not needed")
  }

library(limma) 
boxplot(exprSet,outline=FALSE, notch=T, las=2)
### 该函数默认使用quntile 矫正差异 
exprSet=normalizeBetweenArrays(exprSet)
boxplot(exprSet,outline=FALSE, notch=T, las=2)
## 这步把矩阵转换为数据框很重要
exprSet <- as.data.frame(exprSet)

#################################################################
## 探针基因名转换
##platformMap 中有常见的平台个R注释包的对应关系，这是我整理的。
## 读取，这都是我们已经讲过的
platformMap <- data.table::fread("D:/资料/double/07_GEO_best_practice/resource/platformMap.txt",data.table = F)

## 平台的名称如何知道?
index <- "10558"
## 数据储存在bioc_package这一列中
paste0(platformMap$bioc_package[grep(index,platformMap$gpl)],".db")

## 安装R包,可以直接安装，这里用了判断
if(!requireNamespace("illuminaHumanv4.db")){
  options(BioC_mirror="https://mirrors.ustc.edu.cn/bioc/")
  BiocManager::install("illuminaHumanv4.db",update = F,ask = F,force = TRUE)
} 

######不知道怎能么办
GPL23040_anno<-data.table::fread("D:/资料/double/07_GEO_best_practice/data/GSE122063_family.soft",skip ="ID")
library(dplyr)
library(tidyr)


probe2symbol_df <- GPL23040_anno %>% 
  
  select(ID,mrna_assignment) %>% 
  
  filter(mrna_assignment != "---") %>% 
  
  separate(mrna_assignment,c("drop","drop1","symbol"),sep="//") %>%
  select(ID,symbol)

write.table (probe2symbol_df, file ="D:/资料/double/07_GEO_best_practice/data/probe2symbol_df.txt",  row.names =F, col.names =TRUE) 
probe2symbol_df <- probe2symbol_df %>% 
  separate(symbol,c("symbol","drop"),sep=",") %>%
  select(ID,symbol)
colnames(probe2symbol_df)[1] = 'probe_id'

## 加载R包
library(illuminaHumanv4.db)
## 获取探针和基因的对应关系：这是探针注释的关键步骤
probe2symbol_df <- toTable(get("illuminaHumanv4SYMBOL"))
## 探针有多少个？
length(unique(probe2symbol_df$probe_id))
## 这么多行中，基因名称有重复的么？
length(unique(probe2symbol_df$symbol))

####################################################################
### 探针转换以及去重，获得最终的表达矩阵
### 拆分体会,要学会分步探索排查，#号的使用
library(dplyr)
library(tibble)
exprSet <- exprSet %>% 
  ## 行名转列名,因为只有变成数据框的列,才可以用inner_join
   rownames_to_column("probe_id") %>% 
  ## 合并探针的信息
  inner_join(probe2symbol_df,by="probe_id") %>% 
  ## 去掉多余信息
  dplyr::select(-probe_id) %>%  
  ## 重新排列
  dplyr::select(symbol,everything()) %>%  
  ## rowMeans求出行的平均数(这边的.代表上面传入的数据)
  ## .[,-1]表示去掉出入数据的第一列，然后求行的平均值
  mutate(rowMean =rowMeans(.[,-1])) %>% 
  ## 把表达量的平均值按从大到小排序
  arrange(desc(rowMean)) %>% 
  ## 去重，symbol留下第一个
  distinct(symbol,.keep_all = T) %>% 
  ## 反向选择去除rowMean这一列
  dplyr::select(-rowMean) %>% 
  ## 列名转行名
  column_to_rownames("symbol")
save(exprSet,file = "D:/资料/double/07_GEO_best_practice/output/exprSet_rmdup_GSE132903.Rdata")
write.csv(exprSet, file="D:/资料/double/07_GEO_best_practice/output/GSE132903.csv")  
###有基因symbol的处理
GPL23040_anno<-data.table::fread("D:/资料/double/07_GEO_best_practice/data/GSE122063_family.soft",skip ="ID")
library(dplyr)
library(tidyr)
"ID_REF"<- GPL23040_anno$ID
"symbol"<-GPL23040_anno$GENE_SYMBOL
"ID_REF"<-as.character(ID_REF)
probe2symbol_df<-data.frame(ID_REF,symbol)
exprSet$ID_REF<-rownames(exprSet)

exprSet<-inner_join(probe2symbol_df,exprSet,by="ID_REF")

exprSet<-dplyr::select(exprSet,-ID_REF)

#exprSet<-dplyr::select(exprSet,symbol,everything())
rowMean =rowMeans(exprSet[,-1])

exprSet$rowMean<-rowMean

exprSet<-arrange(exprSet,dplyr::desc(rowMean))

exprSet<-distinct(exprSet,symbol,.keep_all = T)

exprSet[1:10,1:10]

exprSet<-dplyr::select(exprSet,-rowMean)

rownames(exprSet)=exprSet[,1]  #取出第一列

exprSet=exprSet[,-1] 
### 保存数据
save(exprSet,file = "output/exprSet_rmdup_GSE132903.Rdata")

################################################################################
## 补充阅读部分:很重要
## 十分重要，解决剩下的20% 问题
## platformap中没有注释的R包怎么办！！！！
## 探针对应的信息可以从平台文件获取
## https://mp.weixin.qq.com/s/nWbMO4mULgN__nPjooRDlg
## https://mp.weixin.qq.com/s/CSHdvRK6xoNJU91tpper_w
## https://mp.weixin.qq.com/s/DlioHHXQd-W-96tXLWrQvA
## NM_，NR_开头的识别号如何转换成基因名称
## https://mp.weixin.qq.com/s/FdCcliMCYj4Yb4grzIQMaA
## 非编码序列如何转换
## https://mp.weixin.qq.com/s/X8rUnEasKy3Dk-EoUAvC2A
## 如何让基因名称在多个数据库间随意转换？
## https://mp.weixin.qq.com/s/wsiceQmNVveoggiqeDSlmQ

### GZ02_批次矫正视频课程-包括GEO和TCGA  
### https://weidian.com/item.html?itemID=3528879676
### GZ03_GEO芯片的探针ID转换，包括mRNA，lncRNA和环状RNA 
### https://weidian.com/item.html?itemID=3584919356

## GEO教程长期更新的链接是这个:
## https://codingsoeasy.com/archives/geo

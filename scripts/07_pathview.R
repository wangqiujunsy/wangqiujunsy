################################################
################################################
### 作者：果子
### 更新时间：2020-11-20
### 微信公众号:果子学生信
### 私人微信：guotosky
### 个人博客: https://codingsoeasy.com/
### 个人邮箱：hello_guozi@126.com

### 本节任务：使用pathview 这个R包可视化KEGG的结果
################################################
rm(list = ls())
load(file = "output/EGG.Rdata")
KEGG_df <- as.data.frame(EGG)

## 安装这个R包
if(!requireNamespace("pathview")){
  options(BioC_mirror="https://mirrors.ustc.edu.cn/bioc/")
  BiocManager::install("pathview",update = F,ask = F)
} 

### 这个分析需要什么啊？？
#############################################################################
### 很重要的部分，制作geneList
### 什么是geneList

library(clusterProfiler)
load(file = "output/diffgene.Rdata")
gene <- rownames(diffgene)
## 基因名称转换，从symbol 到ENTREZID
gene = bitr(gene, fromType="SYMBOL", toType="ENTREZID", OrgDb="org.Hs.eg.db")
## 去重
gene <- dplyr::distinct(gene,SYMBOL,.keep_all=TRUE)
## 新的数据框
gene_df <- data.frame(logFC=diffgene$logFC,SYMBOL = rownames(diffgene))
## 根据基因名称合并数据
gene_df <- merge(gene_df,gene,by="SYMBOL")

## geneList 三部曲
## 1.获取基因logFC
geneList <- gene_df$logFC
## 2.命名
names(geneList) = gene_df$ENTREZID
## 3.排序很重要
geneList = sort(geneList, decreasing = TRUE)

head(geneList)
############################################################################
### pathview可视化
library(pathview)
pathway.id = "hsa04110"
pv.out <- pathview(gene.data  = geneList,
                     pathway.id = pathway.id,
                     species    = "hsa",
                     limit      = list(gene=max(abs(geneList)), cpd=1))
## 改变倍数的颜色
pv.out <- pathview(gene.data  = geneList,
                   pathway.id = pathway.id,
                   species    = "hsa"
                   #limit      = list(gene=max(abs(geneList)), cpd=1)
                   )
## 改变构图
pv.out <- pathview(gene.data  = geneList,
                   pathway.id = pathway.id,
                   species    = "hsa",
                   #limit      = list(gene=max(abs(geneList)), cpd=1),
                   #改变构图
                   kegg.native=F)

#########################################################################
### pathview批量画图
### 新建文件夹，
dir.create("output/pathview_out")
### 设置工作目录到想要的地方
getwd()
### 然后循环绘图
for (pathway.id in KEGG_df$ID ){
  pathview(gene.data  = geneList,
           pathway.id = pathway.id,
           species    = "hsa"
  )
}

### 结束后切换回原来的工作目录
### 应该是07_GEO_best_practice
getwd()

###################################################################
## GEO教程长期更新的链接是这个:
## https://codingsoeasy.com/archives/geo
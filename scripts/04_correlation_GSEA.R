################################################
################################################
### 作者：果子
### 更新时间：2020-7-2
### 微信公众号:果子学生信
### 私人微信：guotosky
### 个人博客: https://codingsoeasy.com/

### 单基因相关性GSEA
### 注释任意基因，编码，非编码，microRNA，甲基化位点
rm(list = ls())
load(file = "exprSet_arrange.Rdata")

y <- as.numeric(exprSet[1,])
rownames <- rownames(exprSet)
cor_data_df <- data.frame(rownames)
### 批量相关性分析
for (i in 1:length(rownames)){
  print(i)
  test <- cor.test(as.numeric(exprSet[i,]),y,method="spearman")
  cor_data_df[i,2] <- test$estimate
  cor_data_df[i,3] <- test$p.value
}
names(cor_data_df) <- c("symbol","correlation","pvalue")

## geneList 三部曲
## 1.获取基因logFC
geneList <- cor_data_df$correlation
## 2.命名
names(geneList) = cor_data_df$symbol
## 3.排序很重要
geneList = sort(geneList, decreasing = TRUE)


### 基因集
library(msigdbr)
dd <- msigdbr(species = "Homo sapiens")
hallmarks <- dd %>% 
  filter(gs_cat == "H") %>% 
  select(gs_name,gene_symbol)

### GSEA分析
library(clusterProfiler)
y <- GSEA(geneList,TERM2GENE =hallmarks)

### 看整体分布
dotplot(y,showCategory=12,split=".sign")+facet_grid(~.sign)

### 又是神器！基于单基因批量相关性分析的GSEA
### https://mp.weixin.qq.com/s/sZJPW8OWaLNBiXXrs7UYFw
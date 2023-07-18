################################################
################################################
### 作者：果子
### 更新时间：2020-7-2
### 微信公众号:果子学生信
### 私人微信：guotosky
### 个人博客: https://codingsoeasy.com/

### 更详细的GSEA，leading edge analysis
### xPierGSEA
rm(list = ls())

### 基于pi的GSEA分析
options("repos"=c(CRAN="https://mirrors.tuna.tsinghua.edu.cn/CRAN/"))
options(BioC_mirror="https://mirrors.tuna.tsinghua.edu.cn/bioconductor")
if(!require("RCircos")) install.packages("RCircos",update = F,ask = F)
if(!require("ggnetwork")) install.packages("ggnetwork",update = F,ask = F)
if(!require("tibble")) install.packages("tibble",update = F,ask = F)
if(!require("gridExtra")) install.packages("gridExtra",update = F,ask = F)
if(!require("Pi")) BiocManager::install("Pi",update = F,ask = F)

install.packages("https://cran.r-project.org/src/contrib/Archive/XGR/XGR_1.1.7.tar.gz",repos=NULL,type="source")


 ### 数据是什么？
load(file = "D:/资料/double/07_GEO_best_practice/output/allDiffAD.Rdata")

library(dplyr)
library(tibble)
data <- allDiff %>%
  rownames_to_column("gene") %>% 
  select(gene,logFC,adj.P.Val) %>%
  mutate(
    offset = logFC,
    offset_abs = abs(logFC),
  ) %>% 
  mutate(rank = rank(-offset), priority = offset) %>%
  select(priority, rank, gene) %>% 
  column_to_rownames("gene")


# GSEA using MsigdbH
####################
### GSEA 分析
library(Pi)
index <- c("GWAS2EF", "GWAS_LD", "IlluminaHumanHT",
           "IlluminaOmniExpress", "ig.DO",
           "ig.EF", "ig.GOBP", "ig.GOCC", "ig.GOMF", "ig.HPCM", "ig.HPMA",
           "ig.HPMI", "ig.HPPA",
           "ig.MP", "org.Hs.eg", "org.Hs.egDGIdb", "org.Hs.egDO", "org.Hs.egGOBP",
           "org.Hs.egGOCC", "org.Hs.egGOMF", "org.Hs.egHPCM", "org.Hs.egHPMA",
           "org.Hs.egHPMI",
           "org.Hs.egHPPA", "org.Hs.egMP", "org.Hs.egMsigdbC1",
           "org.Hs.egMsigdbC2BIOCARTA",
           "org.Hs.egMsigdbC2CGP", "org.Hs.egMsigdbC2CPall",
           "org.Hs.egMsigdbC2CP",
           "org.Hs.egMsigdbC2KEGG", "org.Hs.egMsigdbC2REACTOME",
           "org.Hs.egMsigdbC3MIR",
           "org.Hs.egMsigdbC3TFT", "org.Hs.egMsigdbC4CGN", "org.Hs.egMsigdbC4CM",
           "org.Hs.egMsigdbC5BP", "org.Hs.egMsigdbC5CC", "org.Hs.egMsigdbC5MF",
           "org.Hs.egMsigdbC6", "org.Hs.egMsigdbC7", "org.Hs.egMsigdbH",
           "org.Hs.egPS",
           "org.Hs.egSF", "org.Hs.egPfam", "org.Hs.string", "org.Hs.PCommons_DN",
           "org.Hs.PCommons_UN")

for (i in index) {
  print(i)
  geneset = xRDataLoader(RData=i)
  assign(i,geneset)
  save(list = i,file = paste0("ontology_Rdata/",i,".Rdata"))
}
eGSEA<-
  Pi::xPierGSEA(
    data,
    ontology = "MsigdbH",
    size.range = c(20, 5000),
    nperm = 10000,
    fast = F,
    RData.location = paste0(getwd(),"/ontology_Rdata"))


##eGSEA <-
  Pi::xPierGSEA(
    data,
    ontology = "MsigdbH",
    size.range = c(20, 5000),
    nperm = 20000,
    fast = F,
    RData.location = "http://galahad.well.ox.ac.uk/bigdata")

####################
### 使用xGSEAdotplot来画图
### leading edge 分析
Pi::xGSEAdotplot(
  eGSEA,
  top = 1,
  leading = T,
  leading.edge.only = F,
  max.overlaps =100,
  colormap = "white-orange-yellow",
  zlim = c(-1, 3),
  peak.color = 'black',
  clab = 'LFC\n(ipiNivo - pembro)',
  signature = FALSE
) + theme(
  plot.title = element_text(hjust = 0.5, size = 8),
  plot.subtitle = element_text(hjust = 0.5, size = 6)
)

### 批量作图画
ls_gp <- xGSEAdotplot(eGSEA, top=1:4, signature=F,colormap = "gbr",
                      subtitle ="both",leading = T)
library(gridExtra)
grid.arrange(grobs=ls_gp, ncol=2)

####################
### 画个火山图
df_summary <- eGSEA$df_summary %>%
  mutate(label = gsub('HALLMARK_', '', setID)) %>%
  mutate(flag = ifelse(adjp < 0.01, 'Y', 'N'))

### 火山图
ggplot(df_summary, aes(x = nes, y = -log10(adjp))) +
  geom_point(alpha = 0.6, shape = 16) +
  xlab("NES") +
  ylab(expression(-log[10]("FDR"))) +
  theme_bw() + 
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank()) +
  ggrepel::geom_text_repel(
    data = subset(df_summary, flag == 'Y'),
    aes(label = label),
    size = 2,
    show.legend = F,
    segment.alpha = 0.5,
    segment.color = "grey50",
    segment.size = 0.2,
    arrow = arrow(length = unit(0.01, 'npc'))
  )

### 使用火山图呈现GSEA富集分析的结果
### https://mp.weixin.qq.com/s/jolWmKoLic5m_M5F1E8K2g
### 使用举例
### https://www.nature.com/articles/s41591-019-0734-6
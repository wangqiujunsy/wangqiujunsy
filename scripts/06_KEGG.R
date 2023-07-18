################################################
################################################
### 作者：果子
### 更新时间：2020-11-20
### 微信公众号:果子学生信
### 私人微信：guotosky
### 个人博客: https://codingsoeasy.com/
### 个人邮箱：hello_guozi@126.com

### 本节任务：KEGG分析
rm(list = ls())
library(clusterProfiler)
load(file = "output/diffgeneAD.Rdata")
### 这个分析需要什么数据？
### 获得基因列表
df<-data.frame(diffgene[c("XK","PPP3CB","PPP3CA","TSPO","YWHAE","ATP1B1","FZD9","ITPKB",
                          "GPER1","SNX10","GRINA","ITPR1","ATP2B3","CCKBR","PTH1R","ATP2A2",
                          "NELL2","TMEM178A","NOTCH1","ATF4"),])
gene <- rownames(df)
#基因名称转换，返回的是数据框
gene = bitr(gene, 
            fromType="SYMBOL", 
            toType="ENTREZID", 
            OrgDb="org.Hs.eg.db")
head(gene)
##############################################################
#**KEGG分析**
##############################################################
library(KEGG.db)
### organism = 'hsa'
### http://rest.kegg.jp/list/organism
### https://mp.weixin.qq.com/s/PwrdQAkG3pTlwMB6Mj8wXQ
 library(createKEGGdb)
species <-c("mmu","hsa")
create_kegg_db(species)

EGG <- enrichKEGG(gene = gene$ENTREZID,
                 organism = 'hsa',
                 pvalueCutoff = 0.05,
                 use_internal_data =T)

## 画图
barplot(EGG)
dotplot(EGG)

######################################################################
### KEGG的富集分析比较特殊，他的背后是个网站
KEGG_df <- as.data.frame(EGG)
symboldata <- setReadable(EGG, OrgDb="org.Hs.eg.db", keyType = "ENTREZID")
symboldata  <- as.data.frame(symboldata )

browseKEGG(EGG, 'hsa04110')
save(EGG,file = "output/EGG.Rdata")

#############################################################################
## 1.横坐标使用富集分数,rich factor,参考附赠视频
## 2.富集分析原理，参考附赠视频

## GEO教程长期更新的链接是这个:
## https://codingsoeasy.com/archives/geo
###GO富集分析
library(AnnotationDbi)
library(org.Hs.eg.db)#基因注释包
library(clusterProfiler)#富集包
library(dplyr)
library(ggplot2)#画图包
df$SYMBOL<-rownames(df)

gene.df <- bitr(df$SYMBOL,fromType="SYMBOL",toType="ENTREZID", OrgDb = org.Hs.eg.db)
gene <- gene.df$ENTREZID                
                
                ego_ALL <- enrichGO(gene = gene,#我们上面定义了
                    OrgDb=org.Hs.eg.db,
                    keyType = "ENTREZID",
                    ont = "ALL",#富集的GO类型
                    pAdjustMethod = "BH",#这个不用管，一般都用的BH
                    minGSSize = 1,
                    pvalueCutoff = 0.05,#P值可以取0.05
                    qvalueCutoff = 0.05,
                    readable = TRUE)
ego_CC <- enrichGO(gene = gene,
                   OrgDb=org.Hs.eg.db,
                   keyType = "ENTREZID",
                   ont = "CC",
                   pAdjustMethod = "BH",
                   minGSSize = 1,
                   pvalueCutoff = 0.05,
                   qvalueCutoff = 0.05,
                   readable = TRUE)
ego_BP <- enrichGO(gene = gene,
                   OrgDb=org.Hs.eg.db,
                   keyType = "ENTREZID",
                   ont = "BP",
                   pAdjustMethod = "BH",
                   minGSSize = 1,
                   pvalueCutoff = 0.05,
                   qvalueCutoff = 0.05,
                   readable = TRUE)

ego_MF <- enrichGO(gene = gene,
                   OrgDb=org.Hs.eg.db,
                   keyType = "ENTREZID",
                   ont = "MF",
                   pAdjustMethod = "BH",
                   minGSSize = 1,
                   pvalueCutoff = 0.05,
                   qvalueCutoff = 0.05,
                   readable = TRUE)
ego_ALL <- as.data.frame(ego_ALL)
ego_result_BP <- as.data.frame(ego_BP)
ego_result_CC <- as.data.frame(ego_CC)
ego_result_MF <- as.data.frame(ego_MF)
display_number = c(5, 5, 5)#这三个数字分别代表选取的BP、CC、MF的通路条数，这个自己设置就行了
ego_result_BP <- as.data.frame(ego_BP)[1:display_number[1], ]
ego_result_CC <- as.data.frame(ego_CC)[1:display_number[2], ]
ego_result_MF <- as.data.frame(ego_MF)[1:display_number[3], ]

go_enrich_df <- data.frame(
  ID=c(ego_result_BP$ID, ego_result_CC$ID, ego_result_MF$ID),                         Description=c(ego_result_BP$Description,ego_result_CC$Description,ego_result_MF$Description),
  GeneNumber=c(ego_result_BP$Count, ego_result_CC$Count, ego_result_MF$Count),
  type=factor(c(rep("biological process", display_number[1]), 
                rep("cellular component", display_number[2]),
                rep("molecular function", display_number[3])), 
              levels=c("biological process", "cellular component","molecular function" )))
for(i in 1:nrow(go_enrich_df)){
  description_splite=strsplit(go_enrich_df$Description[i],split = " ")
  description_collapse=paste(description_splite[[1]][1:5],collapse = " ") #这里的5就是指5个单词的意思，可以自己更改
  go_enrich_df$Description[i]=description_collapse
  go_enrich_df$Description=gsub(pattern = "NA","",go_enrich_df$Description)
}
go_enrich_df$type_order=factor(rev(as.integer(rownames(go_enrich_df))),labels=rev(go_enrich_df$Description))#这一步是必须的，为了让柱子按顺序显示，不至于很乱
COLS <- c("#66C3A5", "#8DA1CB", "#FD8D62")#设定颜色
##竖着的图
library(stringr)
####ggplot(data=go_enrich_df, aes(x=type_order,y=GeneNumber, fill=type)) + #横纵轴取值
  geom_bar(stat="identity", width=0.8) + #柱状图的宽度，可以自己设置
  scale_fill_manual(values = COLS) + ###颜色
  coord_flip() + ##这一步是让柱状图横过来，不加的话柱状图是竖着的
  labs(title = "The Most Enriched GO Terms",
       x = "GO term", 
       y = "Gene numbers")+
  theme(axis.title.x = element_text(size = 24),
        legend.title = element_text(size=24),
        legend.text = element_text(size=18))+
  theme_bw(base_size = 13)+
  scale_x_discrete(labels = function(x) str_wrap(x, width = 30) )
#横着的图
go_enrich_df$type_order=factor(rev(as.integer(rownames(go_enrich_df))),labels=rev(go_enrich_df$Description))
COLS <- c("#66C3A5", "#8DA1CB", "#FD8D62")
ggplot(data=go_enrich_df, aes(x=type_order,y=GeneNumber, fill=type)) + 
  geom_bar(stat="identity", width=0.8) + 
  scale_fill_manual(values = COLS) + 
  theme_bw() + 
  xlab("GO term") + 
  ylab("Num of Genes") + 
  labs(title = "The Most Enriched GO Terms")+ 
  theme(axis.text.x=element_text(face = "bold", color="gray50",angle = 70,vjust = 1, hjust = 1,size=6 ))

###########KEGG

kk <- enrichKEGG(gene = gene,keyType = "kegg",organism= "human", qvalueCutoff = 0.05, pvalueCutoff=0.05)

#2、可视化
###柱状图
hh <- as.data.frame(kk)#自己记得保存结果哈！
rownames(hh) <- 1:nrow(hh)
hh$order=factor(rev(as.integer(rownames(hh))),labels = rev(hh$Description))
hh<-hh[1:10,]
ggplot(hh,aes(y=order,x=Count,fill=p.adjust))+
  geom_bar(stat = "identity",width=0.7)+####柱子宽度
  #coord_flip()+##颠倒横纵轴
  scale_fill_gradient(low ="#FD8D62",high = "#66C3A5")+#颜色自己可以换
  labs(title = "KEGG Pathways Enrichment",
       x = "Gene numbers", 
       y = "Pathways")+
  theme(axis.title.x = element_text(face = "bold",size = 24),
        axis.title.y = element_text(face = "bold",size = 24),
        legend.title = element_text(face = "bold",size = 24))+
  theme_bw(base_size = 13)


##气泡图

hh <- as.data.frame(kk)
rownames(hh) <- 1:nrow(hh)
hh$order=factor(rev(as.integer(rownames(hh))),labels = rev(hh$Description))
ggplot(hh,aes(y=order,x=Count))+
  geom_point(aes(size=Count,color=-1*p.adjust))+# 修改点的大小
  scale_color_gradient(low="green",high = "red")+
  labs(color=expression(p.adjust,size="Count"), 
       x="Gene Number",y="Pathways",title="KEGG Pathway Enrichment")+
  theme_bw()

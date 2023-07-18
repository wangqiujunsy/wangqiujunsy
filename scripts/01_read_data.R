################################################
################################################
### 作者：果子
### 更新时间：2021-04-10
### 微信公众号:果子学生信
### 私人微信：guotosky
### 个人博客: https://codingsoeasy.com/
### 个人邮箱：hello_guozi@126.com

### 本节目标：下载geo的数据

###########
#先解压GSE42872_series_matrix.txt.gz，注意解压到当前目录，再读入
#comment.char="!" 意思是！后面的内容不要读取，可以打开文件看一下?read.table

exprSet <- read.table("D:/资料/double/07_GEO_best_practice/data/GSE132903_series_matrix.txt",
                      comment.char="!",
                      stringsAsFactors=F,
                      header=T)

### 高频操作来了：第一列变成行名
### 分两步操作
rownames(exprSet) <- exprSet[,1]
exprSet <- exprSet[,-1]
##############
write.csv(cor_data_df, file="D:/资料/double/07_GEO_best_practice/output/AD/GSEA/GSEAATF4.csv")


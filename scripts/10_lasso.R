#x <- read.csv("expr.xls",header = T,row.names = 1,sep="\t")
#y <- read.csv("data.xls",header = T,row.names = 1,sep="\t")

#typeof(x)

load(file = "output/exprSet_rmdup_AD.Rdata")

ex<-t(exprSet)

ex<-as.data.frame(ex)


group<-c(rep("con",81),rep("treat",80),"con","con","con","con","treat","treat","treat","treat","treat","treat","treat","treat","treat","treat","con","treat","treat","treat","treat","treat","treat","treat","con","con","con","con","con","con","con","con","con","con","con","con") 
####group=factor(group,levels=c('con','treat'),labels=c('0',"1"))


ex$group<-group

library(glmnet)
ex$ID<-rownames(ex)
ex[1:6,19456:19461]
####提取列名
cols_remain2<- c("XK","PPP3CB",
                 "PPP3CA","TSPO","YWHAE",
                 "ATP1B1",	"FZD9","ITPKB","GPER1","SNX10","GRINA","ITPR1","ATP2B3"
                 ,"CCKBR","PTH1R","ATP2A2","NELL2","TMEM178A","NOTCH1","ATF4") #列名
x<- exprSet[ ,colnames(exprSet) %in% cols_remain2]
modle<-ex[,19460:19461]

meta['event'][meta['event'] == 'treat'] = 1 
meta['event'][meta['event'] == 'con'] = 0
meta<-modle
colnames(meta)<-c("event","ID")
y<-meta
#data <- data.matrix(data)
library(glmnet)
library(dplyr)
library(caret)
library(pROC)
library(ggplot2)
library(patchwork)
#加载数据
#data.matrix(x)
#生成随机数种子（如果想让结果可复现）：
#set.seed(233)
inA <- createDataPartition(y$event,
p = 0.8,#切割比例，即从总体中抽取的样本比例，一般常见分割还有8（训练集）:2（测试集）
times = 1,#拆分的数量
list = F)#不将数据以列表的方式返回
#训练集
#head(inA)

train <-x[inA,]
train<- data.matrix(train)
train_meta <- y[inA,]
#测试集
test <-x[-inA,]
test <- data.matrix(test)
test_meta <- y[-inA,]

test_meta$event<-as.numeric(test_meta$event)
train_meta$event<-as.numeric(train_meta$event)
#colnames(train_cl) <- colnames(y)
#rownames(train_cl) <-rownames(train_exp)

#train_cl

#prop.table(table(train_cl$group))
#prop.table(table(test_cl$group))
prop.table(table(train_meta$event))
prop.table(table(test_meta$event))
#构建模型和筛选最佳λ值：
cv_fit <- cv.glmnet(x = train,y = train_meta$event,nlambda = 1000,alpha = 1)

pdf('D:/资料/double/07_GEO_best_practice/output/AD/LASSO/image_1.pdf')
plot(cv_fit)
dev.off()

#提取两个λ值：
lambda.min <- cv_fit$lambda.min
lambda.1se <- cv_fit$lambda.1se

print("提取两个λ值")
model_lasso_min <- glmnet(train, train_meta$event , alpha = 1, lambda = lambda.min)
model_lasso_1se <- glmnet(train, train_meta$event , alpha = 1, lambda = lambda.1se)
model_lasso <- glmnet(train,train_meta$event, nlambda=1000, alpha=1)
#每一行代表一个所构建的模型；
#Df：自由度，代表了非零的线性模型拟合系数的个数；
#%Dev：由模型解释的残差的比例(评估模型好坏)，一般来说越接近1说明模型构建的越好；
#Lambda：λ值，lasso回归复杂度调整的程度是由参数λ来控制的，λ越大，对变量较多的线性模型的惩罚力度就越大，从而最终获得一个变量较少的模型。
pdf('D:/资料/double/07_GEO_best_practice/output/AD/LASSO/lambda.pdf')
plot(model_lasso, xvar="lambda", label=TRUE)
plot(model_lasso, xvar="norm", label=TRUE) #或xvar="dev"
dev.off()

#拎出模型使用的基因(存放在beta中)：
gene_min <- rownames(model_lasso_min$beta)[as.numeric(model_lasso_min$beta)!=0]
gene_1se <- rownames(model_lasso_1se$beta)[as.numeric(model_lasso_1se$beta)!=0]

print("重要的标志物")
print(gene_min)
print(gene_1se)
#predict用训练集构建的模型cv_fit，来预测测试集的数据（test_exp和test_cl）

lasso.prob <- predict(cv_fit, newx = train,s =c(lambda.min,lambda.1se) )

df <- as.data.frame(cbind(train_meta$event ,lasso.prob))

colnames(df) <- c("event","pre_min","pre_1se")

#使用箱线图比较min和1se的预测结果：
library(ggplot2)
df$event <- as.factor(df$event)

#lambda.min预测箱线图：
p <- ggplot(df,
aes(x = event,y = pre_min,group = event))+
geom_boxplot(df,
mapping = aes(color = event,fill = event,group = event),
alpha = 0.4,
show.legend = F)+
geom_jitter(aes(fill = event,color = event),
alpha = 0.3,width = 0.2,shape = 21,size = 1.5)+
scale_y_continuous(limits = c(0,1))+
theme_bw()
#lambda.1se预测箱线图：
p1 <- ggplot(df,
aes(x = event,y = pre_1se,group = event))+
geom_boxplot(df,
mapping = aes(color = event,fill = event,group = event),
alpha = 0.4,
show.legend = F)+
geom_jitter(aes(fill = event,color = event),
alpha = 0.3,width = 0.2,shape = 21,size = 1.5)+
scale_y_continuous(limits = c(0,1))+
theme_bw()


pdf('box_test.pdf')
p+p1
dev.off()
test_df <- df
head(test_df)

#predict用训练集构建的模型cv_fit，来预测训练集的数据（train_exp和train_cl）
lasso.prob <- predict(cv_fit, newx = train,s = c(lambda.min,lambda.1se) )

df <- as.data.frame(cbind(train_meta$event ,lasso.prob))

colnames(df) <- c("event","pre_min","pre_1se")

#使用箱线图比较min和1se的预测结果：
library(ggplot2)
df$event <- as.factor(df$event)

#lambda.min预测箱线图：
p <- ggplot(df,
aes(x = event,y = pre_min,group = event))+
geom_boxplot(df,
mapping = aes(color = event,fill = event,group = event),
alpha = 0.4,
show.legend = F)+
geom_jitter(aes(fill = event,color = event),
alpha = 0.3,width = 0.2,shape = 21,size = 1.5)+
scale_y_continuous(limits = c(0,1))+
theme_bw()

#lambda.1se预测箱线图：
p1 <- ggplot(df,
aes(x = event,y = pre_1se,group = event))+
geom_boxplot(df,
mapping = aes(color = event,fill = event,group = event),
alpha = 0.4,
show.legend = F)+
geom_jitter(aes(fill = event,color = event),
alpha = 0.3,width = 0.2,shape = 21,size = 1.5)+
scale_y_continuous(limits = c(0,1))+
theme_bw()

pdf('box_train.pdf')
p+p1
dev.off()

train_df <- df
head(train_exp)

#构建ROC对象+计算AUC（对pre_min和pre_1se批量处理）：
roc1 <- roc(event ~ pre_min+pre_1se, train_df,direction = "<")
#direction默认"auto"时，通过计算control和cases的中位数判断方向(默认偏向更高AUC值方向)
#自定义direction：当control预测值高于cases：选">";当control预测值低于cases：选"<"

pdf('D:/资料/double/07_GEO_best_practice/output/AD/LASSO/roc1.pdf')
plot.roc(roc1$pre_1se, #取子集（pre_1se）绘制对应模型曲线
print.auc=TRUE,#显示AUC值
auc.polygon=T, #是否显示AUC多边形区域
auc.polygon.col="#52c3ca")
#ci(roc1$pre_min)#查看95%CI范围

p <- plot.roc(event ~ pre_min, train_df,
ci = T,
percent=TRUE,
print.auc=TRUE)
#计算置信区间；
ci <- ci.se(p,specificities=seq(0, 100, 1))
#ROC曲线中显示置信区间：

plot(ci, type="shape", col="lightblue")


#绘制pre_min和pre_1se的ROC曲线组合图：
plot.roc(roc1$pre_min, main="Train ROC", col="#FF99CC90")
plot.roc(roc1$pre_1se, col="#99CC0090", add = T) #add:是否将ROC添加到现有plot上
#添加图例：
legend("bottomright", bty = "n",lwd=3,
legend = c("pre_min = 0.9522", "pre_1se = 0.8276"),
col = c("#FF99CC", "#99CC00"))

dev.off()

#测试集预测效果
roc2 <- roc(event ~ pre_min + pre_1se, test_df,direction = "<")

#绘制pre_min预测的ROC曲线：
pdf('D:/资料/double/07_GEO_best_practice/output/AD/LASSO/roc2.pdf')
plot.roc(roc2$pre_1se,
print.auc=TRUE,
auc.polygon=T,
auc.polygon.col="#52c3ca")
dev.off()
#显示置信区间：
ci(roc2$pre_min)#查看95%CI范围
p <- plot.roc(event ~ pre_min, test_df,
ci = T,
percent=TRUE,
print.auc=TRUE)
#计算置信区间；
ci <- ci.se(p,specificities=seq(0, 100, 1))
#ROC曲线中显示置信区间：
plot(ci, type="shape", col="#52c3ca")
#绘制pre_min和pre_1se的ROC曲线组合图：
plot.roc(roc2$pre_min, main="Test ROC", col="#52c3ca")
plot.roc(roc2$pre_1se, col="#fee882", add = T) #add:是否将ROC添加到现有plot上
#添加图例：
legend("bottomright", bty = "n",lwd=3,
legend = c("pre_min = 0.691", "pre_1se = 0.683"),
col = c("#52c3ca", "#52c3ca"))
dev.off()





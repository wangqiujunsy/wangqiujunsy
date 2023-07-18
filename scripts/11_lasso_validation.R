load(file = "D:/Rfile/double/07_GEO_best_practice/output//exprSet_rmdup_GSE5281.Rdata")
yanzheng<-exprSet
ex<-t(yanzheng)

ex<-as.data.frame(ex)


group<-c(rep("con",12),rep("treat",16)) 
ex$group<-group

library(glmnet)
ex$ID<-rownames(ex)
ex[1:6,23520:23523]

x<- ex[ ,c("XK","PPP3CB",
           "PPP3CA","TSPO","YWHAE","ATP1B1",	"FZD9","SNX10","GRINA","ITPR1","CCKBR","PTH1R","ATP2A2","NELL2",
           "NOTCH1","ATF4","TMEM178A","ATP2B3", "GPER1","ITPKB")]
#colnames(x)<-c("XK","PPP3CB",
"PPP3CA","TSPO","YWHAE","ATP1B1",	"FZD9","SNX10","GRINA","ITPR1","CCKBR","PTH1R","ATP2A2","NELL2",
"NOTCH1","ATF4","TMEM178A","ATP2B3", "GPER1","ITPKB")
modle<-ex[,23522:23523]
meta<-modle

colnames(meta)<-c("event","ID")
meta['event'][meta['event'] == 'treat'] = 1 
meta['event'][meta['event'] == 'con'] = 0
y<-meta
yanzheng_meta<-y

yanzheng_meta$event<-as.numeric(yanzheng_meta$event)

x<-as.matrix(x)

lasso.prob <- predict(cv_fit, newx = x,s = c(lambda.min,lambda.1se) )

df <- as.data.frame(cbind(yanzheng_meta$event ,lasso.prob))

colnames(df) <- c("event","pre_min","pre_1se")
yanzheng_df<-df

#验证集预测效果
roc3 <- roc(event ~ pre_min + pre_1se, yanzheng_df,direction = "auto")

#绘制pre_min预测的ROC曲线：
#pdf('D:/资料/double/07_GEO_best_practice/output/AD/LASSO/roc2.pdf')
plot.roc(roc3$pre_min,
         print.auc=TRUE,
         auc.polygon=T,
         auc.polygon.col="#52c3ca")
dev.off()
#显示置信区间：
ci(roc3$pre_min)#查看95%CI范围
p <- plot.roc(event ~ pre_min, yanzheng_df,
              ci = T,
              percent=TRUE,
              print.auc=TRUE)
#计算置信区间；
ci <- ci.se(p,specificities=seq(0, 100, 1))
#ROC曲线中显示置信区间：
plot(ci, type="shape", col="#52c3ca")

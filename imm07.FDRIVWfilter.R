pvalFilter=0.05      #pvalue过滤条件
fdrFilter=0.2        #FDR过滤条件
mrFile="table.MRresult.csv"        #孟德尔随机化分析的结果文件
pleFile="table.pleiotropy.csv"     #多效性的结果文件
setwd("E:\\K_课题设计\\MR课题\\免疫\\弥漫大B细胞淋巴瘤diffuse large B-cell lymphoma\\07.FDRIVWfilter")     #设置工作目录

#读取孟德尔随机化的结果文件
rt=read.csv(mrFile, header=T, sep=",", check.names=F)

#对pvalue进行矫正
ivwRT=rt[rt$method=="Inverse variance weighted",]
ivwPval=ivwRT[,"pval"]
fdr=p.adjust(ivwPval, method="fdr")
ivwRT=cbind(ivwRT, fdr)
View(ivwRT)

ivwRT=ivwRT[((ivwRT$pval<pvalFilter) & (ivwRT$fdr<fdrFilter)),]

#提取五种方法OR值方向一致的代谢物
ivw=data.frame()
for(metabolite in unique(ivwRT$exposure)){
  metData=rt[rt$exposure==metabolite,]
  if(sum(metData$or>1)==nrow(metData) | sum(metData$or<1)==nrow(metData)){
    ivw=rbind(ivw, ivwRT[ivwRT$exposure==metabolite,])
  }
}

#读取多效性的结果文件
pleRT=read.csv(pleFile, header=T, sep=",", check.names=F)
#剔除pvalue小于0.05的代谢物
pleRT=pleRT[pleRT$pval>0.05,]
metLists=as.vector(pleRT$exposure)
outTab=ivw[ivw$exposure %in% metLists,]
write.csv(outTab, file="IVW.filter.csv", row.names=F)


######生信自学网: https://www.biowolf.cn/
######课程链接1: https://shop119322454.taobao.com
######课程链接2: https://ke.biowolf.cn
######课程链接3: https://ke.biowolf.cn/mobile
######光俊老师邮箱: seqbio@foxmail.com
######光俊老师微信: eduBio


######Video source: https://ke.biowolf.cn
######生信自学网: https://www.biowolf.cn/
######微信公众号：biowolf_cn
######合作邮箱：biowolf@foxmail.com
######答疑微信: 18520221056

pvalFilter=0.05      #pvalue过滤条件(0.05   0.01   0.001)
mrFile="table.MRresult.csv"        #孟德尔随机化分析的结果文件
pleFile="table.pleiotropy.csv"     #多效性的结果文件
setwd("E:\\K_课题设计\\MR课题\\免疫\\弥漫大B细胞淋巴瘤diffuse large B-cell lymphoma\\07.FDRIVWfilter")     #设置工作目录

#读取孟德尔随机化的结果文件
rt=read.csv(mrFile, header=T, sep=",", check.names=F)

#根据pvalue对孟德尔随机化分析的结果进行过滤
ivwRT=rt[rt$method=="Inverse variance weighted",]
ivwRT=ivwRT[ivwRT$pval<pvalFilter,]

#提取五种方法OR值方向一致的免疫细胞
ivw=data.frame()
for(cell in unique(ivwRT$exposure)){
  cellData=rt[rt$exposure==cell,]
  if(sum(cellData$or>1)==nrow(cellData) | sum(cellData$or<1)==nrow(cellData)){
    ivw=rbind(ivw, ivwRT[ivwRT$exposure==cell,])
  }
}

#读取多效性的结果文件
pleRT=read.csv(pleFile, header=T, sep=",", check.names=F)
#剔除pvalue小于0.05的免疫细胞
pleRT=pleRT[pleRT$pval>0.05,]
cellLists=as.vector(pleRT$exposure)
outTab=ivw[ivw$exposure %in% cellLists,]
write.csv(outTab, file="immune-disease.IVWfilter.csv", row.names=F)


######Video source: https://ke.biowolf.cn
######生信自学网: https://www.biowolf.cn/
######微信公众号：biowolf_cn
######合作邮箱：biowolf@foxmail.com
######答疑微信: 18520221056



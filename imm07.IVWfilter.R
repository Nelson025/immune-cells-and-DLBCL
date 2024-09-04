mrFile="table.MRresult.csv"        #孟德尔随机化分析的结果文件
pleFile="table.pleiotropy.csv"     #多效性的结果文件
setwd("E:\\K_课题设计\\MR课题\\免疫\\弥漫大B细胞淋巴瘤diffuse large B-cell lymphoma\\07.IVWfilter")     #设置工作目录

#读取孟德尔随机化的结果文件
rt=read.csv(mrFile, header=T, sep=",", check.names=F)
#提取IVW和Wm方法pvalue<0.05的免疫细胞
ivw=data.frame()
for (immuneCell in unique(rt$exposure)) {
  immData <- rt[rt$exposure == immuneCell, ]
  if (sum(immData$or > 1) == nrow(immData) | sum(immData$or < 1) == nrow(immData)) {
    if (sum(immData$method == "Inverse variance weighted" & immData$pval < 0.05) > 0 & 
        sum(immData$method == "Weighted median" & immData$pval < 0.05) > 0) {
      ivw <- rbind(ivw, immData)
    }
  }
}

#读取多效性的结果文件
pleRT=read.csv(pleFile, header=T, sep=",", check.names=F)
#剔除pvalue小于0.05的免疫细胞
pleRT=pleRT[pleRT$pval>0.05,]
immuneLists=as.vector(pleRT$exposure)
outTab=ivw[ivw$exposure %in% immuneLists,]
write.csv(outTab, file="IVW.filter.csv", row.names=F)



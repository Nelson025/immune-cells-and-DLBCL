#install.packages("reshape2")
#install.packages("circlize")

#if (!require("BiocManager", quietly = TRUE))
#    install.packages("BiocManager")
#BiocManager::install("ComplexHeatmap")


#引用包
library(reshape2)
library(circlize)
library(ComplexHeatmap)

mrFile="table.MRresult.csv"      #孟德尔随机化分析的结果文件
setwd("E:\\K_课题设计\\MR课题\\免疫\\弥漫大B细胞淋巴瘤diffuse large B-cell lymphoma\\07.Circos")      #设置工作目录

#读取孟德尔随机化的结果文件
rt=read.csv(mrFile, header=T, sep=",", check.names=F)

#对数据进行过滤
#selRT=rt[rt$pval<0.05,]
#rt=rt[(rt$id.exposure %in% as.vector(selRT$id.exposure)),]
#rt$method[rt$method=="Inverse variance weighted"]="IVW"

# 筛选出method为"Inverse variance weighted"且pval<0.05的行
ivw_sel <- rt[rt$method == "Inverse variance weighted" & rt$pval < 0.05,]

# 筛选出method为"Weighted median"且pval<0.05的行
wm_sel <- rt[rt$method == "Weighted median" & rt$pval < 0.05,]

# 找到id.exposure在两个条件下都满足的行
common_id_exposure <- intersect(ivw_sel$id.exposure, wm_sel$id.exposure)

# 筛选原始数据框中id.exposure满足上述条件的行
rt <- rt[rt$id.exposure %in% common_id_exposure,]

# 替换method名称
rt$method[rt$method == "Inverse variance weighted"] <- "IVW"


#数据整理
mat=acast(rt, exposure~method, value.var="pval")

library(reshape2) # 加载reshape2包以使用acast函数
# 删除指定行名的数据
rows_to_remove <- c("HLA DR on plasmacytoid DC", "CD64 on CD14+ CD16+ monocyte", "CD39+ CD8br %T cell","CD8br NKT AC","B cell % CD3- lymphocyte","NK %CD3- lymphocyte")
mat <- mat[!rownames(mat) %in% rows_to_remove, ]

#设置图形颜色
col_color = colorRamp2(c(0, 0.05, 0.5, 1), c("red", "white", "skyblue", "blue"))

#绘制圈图
circos.par(gap.after=c(20))
pdf(file="circos1.pdf", width=8, height=8)
circos.heatmap(mat, col = col_color,
               dend.side = "inside",         #聚类图位于圈图内侧
               rownames.side = "outside",    #代谢物名称位于圈图外侧
               bg.border = "black")          #背景边框颜色

#添加孟德尔随机化分析方法的名称
cn = colnames(mat)     #孟德尔随机化分析方法的名称
n = length(cn)
circos.text(rep(CELL_META$cell.xlim[2], n) + 
              convert_x(1, "mm"), 2+(n:1)*0.92,     #x轴和y轴范围
            cn, cex = 0.4,      #展示的列名和列名字体大小
            adj = c(0, 0.5),
            facing = "inside")
circos.clear()

#添加图例
lgd = Legend(title="Pvalue", col_fun=col_color)
grid.draw(lgd)
dev.off()


######生信自学网: https://www.biowolf.cn/
######课程链接1: https://shop119322454.taobao.com
######课程链接2: https://ke.biowolf.cn
######课程链接3: https://ke.biowolf.cn/mobile
######光俊老师邮箱: seqbio@foxmail.com
######光俊老师微信: eduBio


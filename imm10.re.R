# 引用包
library(TwoSampleMR)

# 设置变量
exposureFile = "finngen_R10_C3_DLBCL_EXALLC"  # 暴露数据ID
ivwFile = "19IVW.filter.csv"  # IVW方法过滤的结果文件

# 设置工作目录
setwd("E:\\K_课题设计\\MR课题\\免疫\\弥漫大B细胞淋巴瘤diffuse large B-cell lymphoma\\10.revise")  

# 读取结局数据文件 (IVW方法过滤的结果文件)
rt = read.csv(ivwFile, header = TRUE, sep = ",", check.names = FALSE)
row.names(rt) = rt$id.exposure

####################################################1.提取暴露----
# 提取暴露数据
exposureData = read_exposure_data(filename = exposureFile, sep = "\t",
                                  snp_col = "rsids",
                                  beta_col = "beta",
                                  se_col = "sebeta",
                                  effect_allele_col = "alt",
                                  other_allele_col = "ref",
                                  pval_col = "pval",
                                  eaf_col = "af_alt")

# 添加样本信息
exposureData$ncase.exposure <- 1050
exposureData$ncontrol.exposure <- 314193
exposureData$samplesize.exposure <- 1050 + 314193

exposure <- "diffuse large B-cell lymphoma"
exposureData$exposure <- exposure

# 根据pvalue<5e-08对结果进行过滤
exposure_dat <- subset(exposureData, pval.exposure < 5e-08)

# 去除连锁不平衡的SNP
exposure_dat_clumped <- clump_data(exposure_dat, clump_kb = 10000, clump_r2 = 0.001)

# F值过滤
Ffilter = 10  # F值过滤条件
dat <- exposure_dat_clumped

# 计算F检验值
N = dat[1, "samplesize.exposure"]  # 获取样品的数目
dat = transform(dat, R2 = 2 * ((beta.exposure)^2) * eaf.exposure * (1 - eaf.exposure))  # 计算R2
dat = transform(dat, F = (N - 2) * R2 / (1 - R2))  # 计算F检验值

# 根据F值>10进行过滤, 删除弱工具变量
exposure_dat = dat[dat$F > Ffilter, ]
write.csv(exposure_dat, "exposure.F.csv", row.names = FALSE)

#########################################################2.读取结局，进行MR分析----
# 对结局数据进行循环(免疫细胞)
revPvalVec = c()
for (i in row.names(rt)) {
  # 提取结局数据
  outcome_dat <- extract_outcome_data(snps = exposure_dat$SNP, outcomes = i)
  
  # 将暴露数据和结局数据合并
  dat <- harmonise_data(exposure_dat, outcome_dat)
  
  # MR-PRESSO异常值检测（偏倚的SNP）
  # presso = run_mr_presso(dat)
  # write.csv(presso[[1]]$`MR-PRESSO results`$`Global Test`, file = paste0(i, ".table.MR-PRESSO_Global.csv"))
  # write.csv(presso[[1]]$`MR-PRESSO results`$`Outlier Test`, file = paste0(i, ".table.MR-PRESSO_Outlier.csv"))
  
  # 孟德尔随机化分析
  mrResult = mr(dat)
  mrTab = generate_odds_ratios(mrResult)
  # 输出孟德尔随机化分析的结果
  write.csv(mrTab, file = paste0(i, ".table.MRresult.csv"), row.names = FALSE)
  revPvalVec = c(revPvalVec, mrResult$pval[3])
  
  # 对IVW方法pvalue小于0.05的结果可视化
  if (!is.na(mrResult$pval[3]) && mrResult$pval[3] < 0.05) {
    # 输出用于孟德尔随机化的工具变量
    outTab = dat[dat$mr_keep == "TRUE", ]
    write.csv(outTab, file = paste0(i, ".table.SNP.csv"), row.names = FALSE)
    
    # 异质性检验
    heterTab = mr_heterogeneity(dat)
    write.csv(heterTab, file = paste0(i, ".table.heterogeneity.csv"), row.names = FALSE)
    
    # 多效性检验
    pleioTab = mr_pleiotropy_test(dat)
    write.csv(pleioTab, file = paste0(i, ".table.pleiotropy.csv"), row.names = FALSE)
    
    # 绘制散点图
    pdf(file = paste0(i, ".scatter_plot.pdf"), width = 7, height = 6.5)
    p1 = mr_scatter_plot(mrResult, dat)
    print(p1)
    dev.off()
    
    # 森林图
    res_single = mr_singlesnp(dat)  # 得到每个工具变量对结局的影响
    pdf(file = paste0(i, ".forest.pdf"), width = 6.5, height = 5)
    p2 = mr_forest_plot(res_single)
    print(p2)
    dev.off()
    
    # 漏斗图
    pdf(file = paste0(i, ".funnel_plot.pdf"), width = 6.5, height = 6)
    p3 = mr_funnel_plot(singlesnp_results = res_single)
    print(p3)
    dev.off()
    
    # 留一法敏感性分析
    pdf(file = paste0(i, ".leaveoneout.pdf"), width = 6.5, height = 5)
    p4 = mr_leaveoneout_plot(leaveoneout_results = mr_leaveoneout(dat))
    print(p4)
    dev.off()
  }
}

# 输出反向孟德尔随机化分析的结果文件
outTab = cbind(rt, revPvale = revPvalVec)
outTab = outTab[as.numeric(outTab$revPvale) > 0.05, ]
write.csv(outTab, file = "immune-disease.REVfilter.csv", row.names = FALSE)

library(openxlsx)
library(lipidr)

d <- as_lipidomics_experiment(read.csv("data_matrix.csv"))
non_parsed_molecules(d)
d <- add_sample_annotation(d, "data_clin.csv")

####案例1 两个分组
##将文件读入 R
datadir = system.file("extdata", package="lipidr")
filelist = list.files(datadir, "data.csv", full.names = TRUE) 
d = read_skyline(filelist)
print(d)
##添加示例注释
clinical_file = system.file("extdata", "clin.csv", package="lipidr")
d = add_sample_annotation(d, clinical_file)
colData(d)
##数据子集
d_subset = d[1:10, 1:10]
rowData(d_subset)
colData(d)

d_qc = d[, d$SampleType == "QC"]
rowData(d_qc)
colData(d_qc)

pc_lipids = rowData(d)$Class %in% c("PC", "PCO", "PCP")
d_pc = d[pc_lipids,]
rowData(d_pc)
colData(d_pc)

lipid_classes = rowData(d)$Class %in% c("Cer", "PC", "LPC")
groups = d$BileAcid != "DCA"
d = d[lipid_classes, groups]
#QC sample subset
d_qc = d[, d$group == "QC"]

##原始数据质量检查
plot_samples(d, type = "tic", log = TRUE)
plot_molecules(d_qc, "sd", measure = "Retention Time", log = FALSE)
plot_molecules(d_qc, "cv", measure = "Area")
plot_lipidclass(d, "boxplot")
##互动图
use_interactive_graphics()
use_interactive_graphics(interactive=FALSE)
##总结过渡
# d_summarized = summarize_transitions(d, method = "average")
d_summarized<-d
##归一化
#概率商归一化 (PQN),PQN 方法通过比较样品和参考光谱之间的商分布来确定每个样品的稀释系数，然后使用该稀释系数进行样品归一化。
d_normalized = normalize_pqn(d_summarized, measure = "Area", exclude = "blank", log = TRUE)
plot_samples(d_normalized, "boxplot")##通过指定exclude = "blank"，将自动检测空白运行并将其排除在标准化过程之外。
##内标归一化
#内标归一化校正样品之间的脂质类别特异性差异。使用相同脂质类别的相应内标对脂质类别进行标准化。如果没有找到相应的内标，则使用所有测量的内标的平均值
d_normalized_istd = normalize_istd(d_summarized, measure = "Area", exclude = "blank", log = TRUE)
##多变量分析,您可以使用 PCA 或 PCoA（经典 MDS）调查样本变异。
mvaresults = mva(d_normalized, measure="Area", method="PCA")
plot_mva(mvaresults, color_by="SampleType", components = c(1,2))#通过指定components参数可以绘制其他组件。例如components = c(2,3)绘制第二个和第三个分量。
##监督多变量分析
#可以执行受监督的多变量分析，例如 OPLS 和 OPLS-DA，以确定哪些脂质与感兴趣的组（y 变量）相关。在这个例子中，我们使用“饮食”作为分组，并在分数图中显示结果。
mvaresults = mva(d_normalized, method = "OPLS-DA", group_col = "SampleType", groups=c("RA", "SLE"))
plot_mva(mvaresults, color_by="SampleType")
#我们还可以绘制负载并显示有助于不同（饮食）组之间分离的重要脂质。
plot_mva_loadings(mvaresults, color_by="Class", top.n=10)
#或者，我们可以提取前N 个脂质及其注释。
top_lipids(mvaresults, top.n=10)
##差异分析
#工作流的这一步需要limma安装包。应使用归一化和对数转换数据。
de_results = de_analysis(
  data=d_normalized, 
  HighFat_water - NormalDiet_water,
  measure="Area"
)
head(de_results)
significant_molecules(de_results)
#我们可以使用火山图将差异分析结果可视化。在这种情况下，我们关闭标记，因为在不同条件之间有大量显着改变的脂质。
plot_results_volcano(de_results, show.labels = FALSE)
##复杂的实验设计
#对于可能处理因子调整或交互作用的更复杂的实验设计，建议使用该de_design函数。用户可以提供一个设计矩阵或一个公式来创建一个。
# Using formula
de_design(
  data=d_normalized,
  design = ~ group, 
  coef="groupHighFat_water", 
  measure="Area")

# Using design matrix
design = model.matrix(~ group, data=colData(d_normalized))
de_design(
  data=d_normalized,
  design = design,
  coef="groupHighFat_water",
  measure="Area")
##富集分析
#lipidr根据脂质类别、总链长和不饱和度自动生成脂质组。然后使用差异分析的结果，按倍数变化或 p 值对测量的脂质进行排序。
enrich_results = lsea(de_results, rank.by = "logFC")
significant_lipidsets(enrich_results)
#富集分析结果的可视化。突出显示了丰富的脂质类别和链不饱和度。
plot_enrichment(de_results, significant_lipidsets(enrich_results), annotation="class")
plot_enrichment(de_results, significant_lipidsets(enrich_results), annotation="unsat")
#lipidr还可以绘制每类脂质的倍数变化，显示总链长和不饱和度。图中的数字表示用相同的链特性测量的多种脂质。
plot_chain_distribution(de_results)



####案例2 多分组
use_interactive_graphics()
##加载数据
list_mw_studies(keyword = "lipid")
d = fetch_mw_study("ST001111")
#修改无法识别的名字
non_parsed_molecules(d)
non_parsed <- non_parsed_molecules(d)
new_names <- sub("^.* \\(", "Cer (", non_parsed)
d <- update_molecule_names(d, old = non_parsed, new = new_names)
non_parsed_molecules(d)
#查看临床数据
colData(d)
#接下来，我们告诉lipidr我们的数据集已经过规范化和记录。
d <- set_logged(d, "Area", TRUE)
d <- set_normalized(d, "Area", TRUE)
##质控
#我们查看每个样品的总离子浓度 (TIC) 和分布（箱线图）。
plot_samples(d, "tic")
plot_samples(d, "boxplot")#由于此数据集中没有 QC 样本或技术重复，我们无法评估分子的 %CV。此外，无需summarize_transitions在非靶标数据集中。
##主成分分析
mvaresults = mva(d, measure="Area", method="PCA")
plot_mva(mvaresults, color_by="SampleType", components = c(1,2))#由PC1和解释的低方差PC2（在图中累积显示为R2X）表明这些临床样本中的脂质谱变化很大。
keep_samples <- !colnames(d) %in% c("18", "42")
d <- d[, keep_samples]
##单变量分析两组比较
#为了简单分析，我们可以比较癌症与良性以及癌症与转移。从 PCA 图中，我们可以预期癌症和转移之间的差异很小。
two_group <- de_analysis(d, Cancer-Benign, Cancer-Metastasis)
plot_results_volcano(two_group)
##单变量分析多组比较
#我们可以使用de_design函数进行 ANOVA 式的多组比较，允许用户提供自定义设计矩阵。
#在这个例子中，我们将使用癌症分期作为我们的分组变量。在我们的数据集中，我们可以看到 I-IV 阶段的样本。使用方差分析式分析，我们将识别所有可能在癌症阶段不同的脂质分子。
multi_group <- de_design(d, ~ Stage)
#在这里，我们使用公式波浪号来定义我们的预测变量。~ Stage公式表明我们对与 相关的特征（脂质分子）感兴趣Stage。
significant_molecules(multi_group)
##因子分析
#在复杂的实验设计和临床样本中，我们可能需要纠正混杂变量。这只需通过将变量添加到 中的公式即可完成de_design。例如，下面，我们在校正种族效应时对癌症效应感兴趣。
factorial_de <- de_design(d, ~ Race + SampleType, coef = "SampleTypeCancer")
significant_molecules(factorial_de)
plot_results_volcano(factorial_de)
##多变量分析-正交多变量分析
#可以执行受监督的多变量分析，例如 OPLS 和 OPLS-DA，以确定哪些脂质与感兴趣的组（y 变量）相关。
mvaresults = mva(d, method = "OPLS-DA", group_col = "SampleType", groups=c("Benign", "Cancer"))
plot_mva(mvaresults, color_by="SampleType")
#我们还可以绘制负载并显示有助于不同（饮食）组之间分离的重要脂质。
plot_mva_loadings(mvaresults, color_by="Class", top.n=10)
#或者，我们可以提取前N 个脂质及其注释。
top_lipids(mvaresults, top.n=10)
##多变量分析-具有连续响应变量的监督多变量分析
#OPLS-DA 只能应用于两组比较设置。在某些情况下，我们可能对脂质分子感兴趣……
#在此示例中，我们将癌症阶段的格式设置为数字向量。
stage <- d$Stage
stage[stage == "I"] <- 1
stage[stage == "II"] <- 2
stage[stage == "III"] <- 3
stage[stage == "IV"] <- 4
stage <- as.numeric(stage)
stage
#我们可以看到stage包含缺失值。我们应该先过滤掉它们。
d_filtered <- d[, !is.na(stage)]
stage <- stage[!is.na(stage)]

mvaresults = mva(d_filtered, method = "OPLS", group_col = stage )
plot_mva(mvaresults)

use_interactive_graphics(FALSE)
plot_mva_loadings(mvaresults, color_by="Class", top.n=10)
##富集分析
enrich_results = lsea(two_group, rank.by = "logFC")
significant_lipidsets(enrich_results)
#富集分析结果的可视化。突出显示了丰富的脂质类别。
plot_enrichment(two_group, significant_lipidsets(enrich_results), annotation="class")
##或者，我们可以突出显示显着丰富的链长
plot_enrichment(two_group, significant_lipidsets(enrich_results), annotation="length")
##脂链分析
plot_trend(two_group)
##聚类
plot_heatmap(d, cluster_rows="none")




d<-read.csv("data_matrix.csv",row.names = 1)
for(ii in rownames(d)){
  result<-sapply(as.numeric(d[ii,]),function(i){log10(i)})
  d[ii,]<-result
}
d <- as_lipidomics_experiment(d)
non_parsed_molecules(d)
d <- add_sample_annotation(d, "data_clin.csv")
#####自用代码整理
#接下来，我们告诉lipidr我们的数据集已经过规范化和记录。
d <- set_logged(d, "Area", TRUE)
d <- set_normalized(d, "Area", TRUE)
##质控
#我们查看每个样品的总离子浓度 (TIC) 和分布（箱线图）。
plot_samples(d, "tic")
plot_samples(d, "boxplot")#由于此数据集中没有 QC 样本或技术重复，我们无法评估分子的 %CV。此外，无需summarize_transitions在非靶标数据集中。
##主成分分析
mvaresults = mva(d, measure="Area", method="PCA")
pdf("PCA.pdf",width = 8,height = 8)
plot_mva(mvaresults, color_by="SampleType", components = c(1,2))#由PC1和解释的低方差PC2（在图中累积显示为R2X）表明这些临床样本中的脂质谱变化很大。
dev.off()

##单变量分析两组比较
#为了简单分析，我们可以比较癌症与良性以及癌症与转移。从 PCA 图中，我们可以预期癌症和转移之间的差异很小。
two_group <- de_analysis(d, SLE-HC, RA-HC,SLE-RA)
pdf("Volcano plot.pdf",width = 12,height = 9)
plot_results_volcano(two_group)
dev.off()
# ##单变量分析多组比较
# #我们可以使用de_design函数进行 ANOVA 式的多组比较，允许用户提供自定义设计矩阵。
# #在这个例子中，我们将使用癌症分期作为我们的分组变量。在我们的数据集中，我们可以看到 I-IV 阶段的样本。使用方差分析式分析，我们将识别所有可能在癌症阶段不同的脂质分子。
# multi_group <- de_design(d, ~ SampleType)
# #在这里，我们使用公式波浪号来定义我们的预测变量。~ Stage公式表明我们对与 相关的特征（脂质分子）感兴趣Stage。
# significant_molecules(multi_group)
# ##因子分析
# #在复杂的实验设计和临床样本中，我们可能需要纠正混杂变量。这只需通过将变量添加到 中的公式即可完成de_design。例如，下面，我们在校正种族效应时对癌症效应感兴趣。
# factorial_de <- de_design(d, ~ Race + SampleType, coef = "SampleTypeCancer")
# significant_molecules(factorial_de)
# plot_results_volcano(factorial_de)
##多变量分析-正交多变量分析
#可以执行受监督的多变量分析，例如 OPLS 和 OPLS-DA，以确定哪些脂质与感兴趣的组（y 变量）相关。
mvaresults1 = mva(d, method = "OPLS-DA", group_col = "SampleType", groups=c("HC", "SLE"))
pdf("oplsda-score-hc-sle.pdf",width = 8,height = 8)
plot_mva(mvaresults1, color_by="SampleType")
dev.off()
#我们还可以绘制负载并显示有助于不同（饮食）组之间分离的重要脂质。
pdf("oplsda-loading-hc-sle.pdf",width = 8,height = 8)
plot_mva_loadings(mvaresults1, color_by="Class", top.n=10)
dev.off()
#或者，我们可以提取前N 个脂质及其注释。
top_lipids(mvaresults, top.n=10)

mvaresults2 = mva(d, method = "OPLS-DA", group_col = "SampleType", groups=c("HC", "RA"))
pdf("oplsda-score-hc-ra.pdf",width = 8,height = 8)
plot_mva(mvaresults2, color_by="SampleType")
dev.off()
pdf("oplsda-loading-hc-ra.pdf",width = 8,height = 8)
plot_mva_loadings(mvaresults2, color_by="Class", top.n=10)
dev.off()

mvaresults2 = mva(d, method = "OPLS-DA", group_col = "SampleType", groups=c("RA", "SLE"))
pdf("oplsda-score-ra-sle.pdf",width = 8,height = 8)
plot_mva(mvaresults2, color_by="SampleType")
dev.off()
pdf("oplsda-loading-ra-sle.pdf",width = 8,height = 8)
plot_mva_loadings(mvaresults2, color_by="Class", top.n=10)
dev.off()
##多变量分析-具有连续响应变量的监督多变量分析
#OPLS-DA 只能应用于两组比较设置。在某些情况下，我们可能对脂质分子感兴趣……
# #在此示例中，我们将癌症阶段的格式设置为数字向量。
# stage <- d$Stage
# stage[stage == "I"] <- 1
# stage[stage == "II"] <- 2
# stage[stage == "III"] <- 3
# stage[stage == "IV"] <- 4
# stage <- as.numeric(stage)
# stage
# #我们可以看到stage包含缺失值。我们应该先过滤掉它们。
# d_filtered <- d[, !is.na(stage)]
# stage <- stage[!is.na(stage)]
# 
# mvaresults = mva(d_filtered, method = "OPLS", group_col = stage )
# plot_mva(mvaresults)
# 
# use_interactive_graphics(FALSE)
# plot_mva_loadings(mvaresults, color_by="Class", top.n=10)
##富集分析
enrich_results = lsea(two_group, rank.by = "logFC")
significant_lipidsets(enrich_results)
#富集分析结果的可视化。突出显示了丰富的脂质类别。
pdf("class.pdf",width = 20,height = 10)
plot_enrichment(two_group, significant_lipidsets(enrich_results), annotation="class")
dev.off()
##或者，我们可以突出显示显着丰富的链长
pdf("chainlong.pdf",width = 20,height = 10)
plot_enrichment(two_group, significant_lipidsets(enrich_results), annotation="length")
dev.off()
##脂链分析
pdf("chain.pdf",width = 15,height = 10)
plot_trend(two_group)
dev.off()
##聚类
plot_heatmap(d, cluster_rows="none")


pdf("chainunsat.pdf",width = 20,height = 10)
plot_enrichment(two_group, significant_lipidsets(enrich_results), annotation="unsat")
dev.off()


pdf("distribution.pdf",width = 30,height = 30)
plot_chain_distribution(two_group)
dev.off()




data<-two_group[two_group$contrast=="SLE - HC",]
data$Molecule<-str_replace_all(data$Molecule,"-","/")
data$class1<-raw$`Class-I`[match(data$Molecule,raw$Molecule)]


data<-data[order(data$total_cl),]
data$total_cl<-as.factor(data$total_cl)
data$class1[data$P.Value>0.05]<-"not significant"

color<-jet(length(unique(data$class1)))
names(color)<-unique(data$class1)
color["not significant"]<-"grey"

p <- ggplot(data, aes(x=total_cl, y=logFC,color=class1)) + 
  geom_jitter(aes(color=class1,size=abs(logFC)),shape=16,width = 0.35,alpha=0.4) + ##geom_quasirandom
  scale_color_manual(values = color)+
  scale_size_continuous(
    name = 'logFC', range = c(4, 8),
    breaks = c(min(abs(data$logFC)), max(abs(data$logFC))),
    labels = c(min(abs(data$logFC)), round(max(abs(data$logFC)), 3)),
    guide = guide_legend(title.hjust = 0.5)
  )+
  theme_bw() + theme(panel.grid.major = element_blank(),
                     panel.grid.minor = element_blank(), legend.key=element_blank(),
                     axis.title = element_text(size = 20),
                     axis.text = element_text(size = 20),
                     # axis.ticks = element_blank(),
                     legend.title = element_text(size = 20),
                     legend.text = element_text(size = 20))+
  # legend.position = c(0.85,0.8))+
  # coord_cartesian(ylim = c(0.1,0.55))+
  labs(x = "total chain long", y = "logfc")
# theme(legend.position="none")
ggsave("jitterplot2.png",p, width = 24, height = 8,units = "in",dpi = 600)
ggsave("jitterplot2.pdf",p, width = 24, height = 8)

write.xlsx(enrich_results,"enrich_results.xlsx",rowNames=T)
write.xlsx(two_group,"result.xlsx",rowNames=T,overwrite = T)

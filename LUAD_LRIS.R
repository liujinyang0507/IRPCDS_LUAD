rm(list=ls())
#引用包
library(data.table)
library(dplyr)
library(tidyverse)
library(gridExtra)
library(psych)
library(pheatmap)
library(RColorBrewer)
library(randomcoloR)
library(ggsci)


geneSets_im <- geneSets
names(geneSets)

sets <- c("28_Immune_Cell_Infiltration_GeneSets","7_Steps_of_the_Tumor_Immune_Cycle","75_ICDs_Gene_pathway")
immune_associated_genesets <- geneSets_im[names(geneSets_im) %in% sets]
names(immune_associated_genesets)
save(immune_associated_genesets,file = "immune_associated_genesets.Rdata")


# 获取 NPG 颜色调色板的前 10 种颜色
mycolpal2 <- pal_npg()(10)

names(geneSets$`28_Immune_Cell_Infiltration_GeneSets`)

roc_1_3_5_col <- c(C1="#fdc68a",C2="#05a59e",C3="#813e98")
#roc_1_3_5_col <- c(C1="#abe2b3",C2="#40c3ec",C3="#eae8d7")
clust_col <- c(C1="#C4B4CE",C2="#46AC85",C3="#60599C")
risk_col_c <- c("Low-risk"="#05a59e","High-risk"="#813e98")
risk_col_l <- c("#05a59e","#813e98")
names(risk_col_l) <- c("Low-risk","High-risk")
mycolpanel <- c("#abe2b3","#40c3ec","#eae8d7","#fdc68a","#05a59e","#813e98")
colorbrewers <- c("BrBG","PiYG","PRGn","PuOr","RdBu","RdYlBu","RdYlGn")
tumor_normal_col <- c(normal="#fdc68a",tumor="#05a59e")
stage_col <- c("Stage I"="#abe2b3","Stage II"="#40c3ec",
               "Stage III"="#eae8d7","Stage IV"="#fdc68a")
annotation_cohorts <- c("TCGA-LUAD","CPTAC-LUAD","GPL570_GSE31210","GPL570_GSE50081","GPL570_GSE37745")

#加载数据
root_dir <- "验证数据集\\LUAD\\clean_data"
load(file=file.path(root_dir,"all_cohorts_raw_data_mRNA_os_clean.Rdata"))
load(file=file.path(root_dir,"all_cohorts_raw_data_mRNA_tumor_normal_expression_clean.Rdata"))
load(file=file.path(root_dir,"all_cohorts_all_genes_unicoxs_res_sig.Rdata"))
lapply(all_cohorts_raw_data_mRNA_os_clean,function(x){dim(x)})

names(all_cohorts_raw_data_mRNA_os_clean)

#绘制PCD基因圈图

library(ggplot2)
library(RColorBrewer)

pcd_genesets <- geneSets$`18PCD`
pcd_gene_size <- lapply(pcd_genesets, function(x){length(x)})
length(names(pcd_genesets))
uni <- data.frame(num=unlist(pcd_gene_size),pcd=names(pcd_genesets))
uni <- uni %>% mutate(Proportion=paste0(round(num/sum(num)*100,digits = 2),"%"))
uni$num <- log(uni$num)
uni <- uni[order(uni$num,decreasing = F),]
uni

pal <- brewer.pal(10, colorbrewers[1])
colors = c(pal[4],pal[7])
uni$id <- seq(nrow(uni):1)

#计算标签角度：
number <- nrow(uni) #条形图柱子的数量
angle <- 90 - 360*(uni$id - 0.5)/number #减去0.5保证标签在中心位置，如果在极左或极右，则为0与1；
#选择标签对齐方式(左对齐or右对齐，后期在图中添加标签时需要)
uni$hjust <- ifelse(angle < -90, 1, 0) #根据角度判断
#翻转角度，便于显示：
uni$angle <- ifelse(angle < -90, angle + 180, angle)
uni$label <- uni$pcd
head(uni)
uni$label <- paste0(uni$pcd,"\n","(",uni$Proportion,")")

head(uni)
p <- ggplot(uni, aes(x = as.factor(id), y = num)) +
  geom_bar(stat = "identity", aes(fill = num)) +
  scale_fill_gradient(low = colors[1], high = colors[2]) +
  geom_text(aes(x = id, y = num - 1.2, label = label, hjust = hjust), color = "black", size = 1) +
  theme_minimal() +
  theme(axis.text = element_blank(), axis.title = element_blank(), panel.grid = element_blank())  +
  coord_polar(start = 0) + theme(legend.position = NULL)
p

ggsave(filename = file.path("PCD_roles_in_cancers","pcd_gene_size.pdf"),p,height = 4,width = 6)
write.csv(uni,file.path("PCD_roles_in_cancers","pcd_gene_size.csv"))


####################PCD相关基因在肺癌发生发展中的作用###########################
#癌和癌旁PCD相关通路的比较

dir.create("PCD_roles_in_cancers")
pcd_genesets <- geneSets$`18PCD`
#pcd_ssGSEA <- list()
# for (cohort in names(all_cohorts_raw_data_mRNA_tumor_normal_expression_clean)){
#   #cohort <- names(all_cohorts_raw_data_mRNA_tumor_normal_expression_clean)[1]
#   sur_exp <- all_cohorts_raw_data_mRNA_tumor_normal_expression_clean[[cohort]]
#   ssGSEA.res <- ssGSEA(data = sur_exp,need_log = F,genesets = pcd_genesets,scale = T)
#   ssGSEA.res <- data.frame(t(ssGSEA.res))
#   ssGSEA.res[1:4,1:4]
#   normal.ssGSEA.res <- ssGSEA.res[grepl("normal",row.names(ssGSEA.res)),]
#   tumor.ssGSEA.res <- ssGSEA.res[grepl("tumor",row.names(ssGSEA.res)),]
#   pcd_ssGSEA[[cohort]] <- list(tumor=tumor.ssGSEA.res,normal=normal.ssGSEA.res)
# }
# save(pcd_ssGSEA,file="pcd_ssGSEA.Rdata")
load(file="pcd_ssGSEA.Rdata")
tumor_normal_cohorts <- c("TCGA-LUAD","CPTAC-LUAD","GPL570_GSE31210")

#分析每个PCD在癌和癌旁中的表达改变的情况
pathway_diff_matrix <- matrix(NA,nrow = 18,ncol = 0)
row.names(pathway_diff_matrix) <- names(pcd_genesets)
for (cohort in tumor_normal_cohorts){
  #cohort <- tumor_normal_cohorts[1]
  tumor.ssGSEA.res <- pcd_ssGSEA[[cohort]]$tumor
  normal.ssGSEA.res <- pcd_ssGSEA[[cohort]]$normal
  data <- rbind(normal.ssGSEA.res,tumor.ssGSEA.res)
  
  group_dat <- data.frame(group=sapply(row.names(data),function(x){strsplit(x,"_")[[1]][2]}))
  row.names(group_dat) <- row.names(data) 
  head(group_dat)
  diff.res <- two.group.test(data=data,group=group_dat,controlgroup="normal",
                             experimentalgroup="tumor",method="limma")
  pathway_fg <- diff.res[names(pcd_genesets),c("logFC","adj.P.Val"),drop=F]
  pathway_fg[,"alteration"]  <- ifelse(pathway_fg$adj.P.Val>0.05,"P>0.05",ifelse(pathway_fg$logFC>0,"up","down"))
  pathway_fg <- pathway_fg[names(pcd_genesets),"alteration",drop=F]
  colnames(pathway_fg) <- cohort
  pathway_diff_matrix <- cbind(pathway_diff_matrix,pathway_fg)
  data[,"type"] <- sapply(row.names(data),function(x){strsplit(x,"_")[[1]][2]})
  data[,"id"] <- row.names(data)
  data_long <- data %>% melt(id=c("id","type"),variable.names="cell",value.name  = "value")
  
  p <- multi_item_group_violin_boxplot_plot(data=data_long,
                                            featureVariable="variable",
                                            valueVariable="value",
                                            groupVariable="type",
                                            colors=tumor_normal_col,
                                            ylab="PCD pathway score",
                                            type="boxplot")
  p
  ggsave(filename = file.path("PCD_roles_in_cancers",paste0(cohort,"_tumor_normal_pcd_pathways_score_diff_boxplot.pdf")),
         height = 5,width = 8 )
}
data <- pathway_diff_matrix
# 将矩阵转换为数值型，以便进行热图绘制
data_numeric <- matrix(NA, nrow = nrow(data), ncol = ncol(data))
row.names(data_numeric) <- row.names(data)
colnames(data_numeric) <- colnames(data)
data_numeric[data == "up"] <- 1
data_numeric[data == "down"] <- -1
data_numeric[data == "P>0.05"] <- 0

library(pheatmap)
#color_mapping <- colorRampPalette(c("#BC102B", "white", "red"))(3)
color_mapping = c( "up" = "#BC102B","down"="#05a59e","P>0.05" = "#EFEFEF")
# 绘制热图
pdf(file=file.path("PCD_roles_in_cancers",paste0("Tumor_normal_pcd_pathways_expression_heatmap.pdf")),height = 4,width = 8)
pheatmap(data_numeric, 
         color = color_mapping,cluster_rows = T,cluster_cols = F, 
         legend = F, show_rownames = T,show_colnames = T,display_numbers = FALSE,
         cellwidth = 6,cellheight = 6,fontsize_row=6,fontsize_col=6)
dev.off()

###################PCD相关通路在病理分期上的差异################################
for (cohort in annotation_cohorts){
  #cohort <- annotation_cohorts[1]
  tumor.ssGSEA.res <- pcd_ssGSEA[[cohort]]$tumor
  row.names(tumor.ssGSEA.res) <- gsub("_tumor","",row.names(tumor.ssGSEA.res))
  tumor.ssGSEA.res[1:4,1:4]
  clinical <- read.csv(file.path(root_data_dir,paste0(cohort,"_clinical.csv")),row.names = 1)
  clinical[1:4,1:4]
  clinical <- clinical[,c("Stage","pT","pN","pM")]
  clinical_pcd <- merge(clinical,tumor.ssGSEA.res,by.x=0,by.y=0)
  
  for (i in c("Stage","pT","pN","pM")){
    #i="pT"
    clinical_pcd_tmp <- clinical_pcd[!is.na(clinical_pcd[i]),]
    clinical_pcd_tmp <- clinical_pcd_tmp[!clinical_pcd_tmp[i]=="Uknown",]
    clinical_pcd_tmp$Stage <- factor(clinical_pcd_tmp$Stage,levels = c("Stage I","Stage II","Stage III","Stage IV"))
    clinical_pcd_tmp$pT <- factor(clinical_pcd_tmp$pT,levels = c("T1","T2","T3","T4"))
    clinical_pcd_tmp$pN <- factor(clinical_pcd_tmp$pN,levels = c("N0","N1","N2","N3"))
    clinical_pcd_tmp$pM <- factor(clinical_pcd_tmp$pM,levels = c("M0","M1"))
    clinical_pcd_tmp <- clinical_pcd_tmp[order(clinical_pcd_tmp[,i]),]
    
    if ((length(clinical_pcd_tmp[,i])>30) & (length(unique(clinical_pcd_tmp[,i]))>1)){
      plot_list <- list()
      for (j in colnames(clinical_pcd_tmp)[6:ncol(clinical_pcd_tmp)]){
        #j <- colnames(clinical_pcd_tmp)[6:ncol(clinical_pcd_tmp)][1]
        p <- boxtools(data=clinical_pcd_tmp,
                      valueVariable=j,
                      groupVariable=i,
                      colors=sample(mycolpanel,length(unique(clinical_pcd[,i]))),
                      ylab=j,
                      plot_tyep="boxplot",
                      signif_label="p.format")
        plot_list[[j]] <- p
      }
      p <- grid.arrange(grobs = plot_list,ncol=3 )
      p
      ggsave(filename = file.path("PCD_roles_in_cancers","病理分期差异",
                                  paste0(cohort,"_",i,"_pcd_pathway_diff.pdf")),p,width=20,height = 20,dpi=300)
    }
  }
}

#单因素分析PCD相关通路在每个队列中是风险因子还是保护因子
unicox.pcd.all.cohorts <- matrix(NA,nrow = 18,ncol = 0)
row.names(unicox.pcd.all.cohorts) <- names(pcd_genesets)
for (cohort in annotation_cohorts){
  #cohort <- annotation_cohorts[1]
  tumor.ssGSEA.res <- pcd_ssGSEA[[cohort]]$tumor
  sur_exp <- all_cohorts_raw_data_mRNA_os_clean[[cohort]][]
  sur <- sur_exp[grep("tumor",row.names(sur_exp)),c(1:2)]
  pathway_score_sur <- merge(sur,tumor.ssGSEA.res,by.x = 0,by.y = 0)
  pathway_score_sur <- pathway_score_sur %>% column_to_rownames("Row.names")
  head(pathway_score_sur)
  res <- unicox(pathway_score_sur)
  result <- data.frame(res$res)
  row.names(result) <- result$Variable
  head(result)
  unicox.res <- result[names(pcd_genesets),c("HR","p.value"),drop=F]
  unicox.res[,"alteration"]  <- ifelse(unicox.res$p.value>=0.05,"p>0.05",ifelse(unicox.res$HR>0,"Risky","Protective"))
  unicox.res <- unicox.res[names(pcd_genesets),"alteration",drop=F]
  colnames(unicox.res) <- cohort
  head(unicox.res)
  unicox.pcd.all.cohorts <- cbind(unicox.pcd.all.cohorts,unicox.res)
}
head(unicox.pcd.all.cohorts)
data <- unicox.pcd.all.cohorts
#data <- t(unicox.pcd.all.cohorts)
# 将矩阵转换为数值型，以便进行热图绘制
data_numeric <- matrix(NA, nrow = nrow(data), ncol = ncol(data))
row.names(data_numeric) <- row.names(data)
colnames(data_numeric) <- colnames(data)
data_numeric[data == "Risky"] <- 1
data_numeric[data == "Protective"] <- -1
data_numeric[data == "p>0.05"] <- 0
head(data_numeric)

library(pheatmap)
#color_mapping <- colorRampPalette(c("#BC102B", "white", "red"))(3)
color_mapping = c( "Risky" = "#BC102B","Protective"="#05a59e","p>0.05" = "#EFEFEF")
# 绘制热图
pdf(file=file.path("PCD_roles_in_cancers",paste0("Pcd_pathways_unicox_heatmap_T.pdf")),height = 8,width = 4)
pheatmap(data_numeric, 
         color = color_mapping,cluster_rows = T,cluster_cols = F, 
         legend = F, show_rownames = T,show_colnames = T,display_numbers = FALSE,
         cellwidth = 6,cellheight = 6,fontsize_row=6,fontsize_col=6)
dev.off()

#分析pcd相关通路和10个癌症经典通路的相关性
#计算10个经典癌症信号通路活性
# ten_common_cancer_pathway_genesets <- geneSets$`10_Classic_Cancer_Signaling_Pathways`
# ten_common_cancer_pathway_ssGSEA <- list()
# for (cohort in names(all_cohorts_raw_data_mRNA_tumor_normal_expression_clean)){
#   #cohort <- names(all_cohorts_raw_data_mRNA_tumor_normal_expression_clean)[1]
#   sur_exp <- all_cohorts_raw_data_mRNA_tumor_normal_expression_clean[[cohort]]
#   ssGSEA.res <- ssGSEA(data = sur_exp,need_log = F,genesets = ten_common_cancer_pathway_genesets,scale = T)
#   ssGSEA.res <- data.frame(t(ssGSEA.res))
#   ssGSEA.res[1:4,1:4]
#   normal.ssGSEA.res <- ssGSEA.res[grepl("normal",row.names(ssGSEA.res)),]
#   tumor.ssGSEA.res <- ssGSEA.res[grepl("tumor",row.names(ssGSEA.res)),]
#   ten_common_cancer_pathway_ssGSEA[[cohort]] <- list(tumor=tumor.ssGSEA.res,normal=normal.ssGSEA.res)
# }
# save(ten_common_cancer_pathway_ssGSEA,file="ten_common_cancer_pathway_ssGSEA.Rdata")
load(file="ten_common_cancer_pathway_ssGSEA.Rdata")
library(corrplot)
library(corrgram)
library(randomcoloR)

#绘制PCD相关通路得分和10个癌症经典信号通路得分
cor.pcd.common.cancer.matrix <- matrix(NA,nrow = 18,ncol = 0)
row.names(cor.pcd.common.cancer.matrix) <- names(geneSets$`18PCD`)
pvalue.pcd.common.cancer.matrix <- cor.pcd.common.cancer.matrix
for (cohort in annotation_cohorts){
  #cohort <- annotation_cohorts[1]
  tumor.ssGSEA.res <- pcd_ssGSEA[[cohort]]$tumor
  tumor.common_cancer_pathway.res <- ten_common_cancer_pathway_ssGSEA[[cohort]]$tumor
  colnames(tumor.common_cancer_pathway.res) <- gsub("\\."," ",colnames(tumor.common_cancer_pathway.res))
  head(tumor.common_cancer_pathway.res)
  head(tumor.ssGSEA.res)
  #计算相关性
  pp <- corr.test(tumor.ssGSEA.res, tumor.common_cancer_pathway.res,method = "spearman", adjust = "fdr")
  cor <- pp$r # 获取相关系数矩阵
  pvalue <- pp$p # 获取p-value矩阵
  cor <- cor[names(geneSets$`18PCD`),names(geneSets$`10_Classic_Cancer_Signaling_Pathways`)]
  colnames(cor) <- paste0(cohort,"#",colnames(cor))
  cor.pcd.common.cancer.matrix <- cbind(cor.pcd.common.cancer.matrix,cor)
  
  pvalue <- pvalue[names(geneSets$`18PCD`),names(geneSets$`10_Classic_Cancer_Signaling_Pathways`)]
  colnames(pvalue) <- paste0(cohort,"#",colnames(pvalue))
  pvalue.pcd.common.cancer.matrix <- cbind(pvalue.pcd.common.cancer.matrix,pvalue)
}
dim(cor.pcd.common.cancer.matrix)
cor.pcd.common.cancer.matrix.sort <- matrix(NA, nrow = nrow(cor.pcd.common.cancer.matrix), ncol = 0)
row.names(cor.pcd.common.cancer.matrix.sort) <- row.names(cor.pcd.common.cancer.matrix.sort)
pathways <- unique(sapply(colnames(cor.pcd.common.cancer.matrix),function(x){strsplit(x,"#")[[1]][2]}))
length(pathways)
for (i in pathways){
  #i=pathways[1]
  tmp <- cor.pcd.common.cancer.matrix[,grepl(i,pathways)]
  cor.pcd.common.cancer.matrix.sort <- cbind(cor.pcd.common.cancer.matrix.sort,tmp)
}
pvalue.pcd.common.cancer.matrix.sort <- pvalue.pcd.common.cancer.matrix[,colnames(cor.pcd.common.cancer.matrix.sort)]
pvalue.pcd.common.cancer.matrix.sort <- apply(pvalue.pcd.common.cancer.matrix.sort,2,function(x){ifelse(x > 0.05,"",".")})

#pvalue.pcd.common.cancer.matrix.sort <- apply(pvalue.pcd.common.cancer.matrix.sort,2,function(x){ifelse(x >0.05,"",ifelse(x>0.01,"*",ifelse(x>0.001,"**","***")))})
head(pvalue.pcd.common.cancer.matrix.sort)
cohort_annotation <- sapply(colnames(cor.pcd.common.cancer.matrix.sort),function(x){strsplit(x,"#")[[1]][1]})
pathway_annotation <- sapply(colnames(cor.pcd.common.cancer.matrix.sort),function(x){strsplit(x,"#")[[1]][2]})
annotation_col <- data.frame(cohort=cohort_annotation,pathway=pathway_annotation)
row.names(annotation_col) = colnames(cor.pcd.common.cancer.matrix.sort)
color_cohort <- mycolpanel[1:length(unique(cohort_annotation))]
names(color_cohort) <- c(unique(annotation_col$cohort))
color_pathway <- distinctColorPalette(length(unique(pathway_annotation)))
names(color_pathway) <- c(unique(annotation_col$pathway))

pdf(file.path("PCD_roles_in_cancers",paste0("Pcd_pathways_10_common_cancer_pathway_cor_heatmap.pdf")),height = 10,width = 10)
pheatmap(cor.pcd.common.cancer.matrix.sort,show_colnames = F,show_rownames = T,
         scale = "none",cluster_rows=T,cluster_cols=F,
         clustering_method = 'complete',
         color = colorRampPalette(c('#21b6af','white','#eeba4d'))(100),
         #color = colorRampPalette(c('#21b6af','white','#eeba4d'))(100),
         annotation_col = annotation_col,
         annotation_colors = list(cohort=color_cohort,pathway=color_pathway),
         gaps_col = c(5,10,15,20,25,30,35,40,45),
         display_numbers=pvalue.pcd.common.cancer.matrix.sort,fontsize_number=1,
         cellwidth = 6,cellheight = 6,fontsize_row=6,fontsize_col=6)
dev.off()
  
##计算28个免疫细胞浸润得分  
# Immune_Cell_Infiltration_genesets <- geneSets$`28_Immune_Cell_Infiltration_GeneSets`
# Immune_Cell_Infiltration_ssGSEA <- list()
# for (cohort in names(all_cohorts_raw_data_mRNA_tumor_normal_expression_clean)){
#   #cohort <- names(all_cohorts_raw_data_mRNA_tumor_normal_expression_clean)[1]
#   sur_exp <- all_cohorts_raw_data_mRNA_tumor_normal_expression_clean[[cohort]]
#   ssGSEA.res <- ssGSEA(data = sur_exp,need_log = F,genesets = Immune_Cell_Infiltration_genesets,scale = T)
#   ssGSEA.res <- data.frame(t(ssGSEA.res))
#   ssGSEA.res[1:4,1:4]
#   normal.ssGSEA.res <- ssGSEA.res[grepl("normal",row.names(ssGSEA.res)),]
#   tumor.ssGSEA.res <- ssGSEA.res[grepl("tumor",row.names(ssGSEA.res)),]
#   Immune_Cell_Infiltration_ssGSEA[[cohort]] <- list(tumor=tumor.ssGSEA.res,normal=normal.ssGSEA.res)
# }
# save(Immune_Cell_Infiltration_ssGSEA,file="28_Immune_Cell_Infiltration_ssGSEA.Rdata")
load(file="28_Immune_Cell_Infiltration_ssGSEA.Rdata") 
#绘制PCD相关通路得分和28个免疫细胞浸润之间的相关性
cor.pcd.immune.infiltration.matrix <- matrix(NA,nrow = 18,ncol = 0)
row.names(cor.pcd.immune.infiltration.matrix) <- names(geneSets$`18PCD`)
pvalue.pcd.immune.infiltration.matrix <- cor.pcd.immune.infiltration.matrix
for (cohort in annotation_cohorts){
  #cohort <- annotation_cohorts[1]
  tumor.ssGSEA.res <- pcd_ssGSEA[[cohort]]$tumor
  immune_Cell_Infiltration_ssGSEA <- Immune_Cell_Infiltration_ssGSEA[[cohort]]$tumor
  #colnames(immune_Cell_Infiltration_ssGSEA) <- gsub("\\."," ",colnames(immune_Cell_Infiltration_ssGSEA))
  #colnames(immune_Cell_Infiltration_ssGSEA) <- gsub(" Adaptive| Innate","",colnames(immune_Cell_Infiltration_ssGSEA))
  head(immune_Cell_Infiltration_ssGSEA)
  head(tumor.ssGSEA.res)
  colnames(immune_Cell_Infiltration_ssGSEA)
  #col <- sapply(names(geneSets$`28_Immune_Cell_Infiltration_GeneSets`),function(x){strsplit(x,"#")[[1]][1]})
  #计算相关性
  pp <- corr.test(tumor.ssGSEA.res, immune_Cell_Infiltration_ssGSEA,method = "spearman", adjust = "fdr")
  cor <- pp$r # 获取相关系数矩阵
  head(cor)
  pvalue <- pp$p # 获取p-value矩阵
  #cor <- cor[names(geneSets$`18PCD`),col]
  colnames(cor) <- paste0(cohort,"#",colnames(cor))
  cor.pcd.immune.infiltration.matrix <- cbind(cor.pcd.immune.infiltration.matrix,cor)
  #pvalue <- pvalue[names(geneSets$`18PCD`),col]
  colnames(pvalue) <- paste0(cohort,"#",colnames(pvalue))
  pvalue.pcd.immune.infiltration.matrix <- cbind(pvalue.pcd.immune.infiltration.matrix,pvalue)
}
dim(cor.pcd.immune.infiltration.matrix)
head(cor.pcd.immune.infiltration.matrix)

#row.names(cor.pcd.immune.infiltration.matrix.sort) <- row.names(cor.pcd.immune.infiltration.matrix.sort)
pathways <- unique(sapply(colnames(cor.pcd.immune.infiltration.matrix),function(x){strsplit(x,"#")[[1]][2]}))
pathway_Innate <- pathways[grepl("Innate",pathways)]
pathway_Adaptive <- pathways[grepl("Adaptive",pathways)]
#分别手动绘制适应性免疫和获得性免疫热图
#a=pathway_Adaptive
a=pathway_Innate
cor.pcd.immune.infiltration.matrix.sort <- matrix(NA, nrow = nrow(cor.pcd.immune.infiltration.matrix), ncol = 0)
for (i in a){
  #i=pathways[1]
  tmp <- cor.pcd.immune.infiltration.matrix[,grepl(i,pathways)]
  cor.pcd.immune.infiltration.matrix.sort <- cbind(cor.pcd.immune.infiltration.matrix.sort,tmp)
}
pvalue.pcd.immune.infiltration.matrix.sort <- pvalue.pcd.immune.infiltration.matrix[,colnames(cor.pcd.immune.infiltration.matrix.sort)]
pvalue.pcd.immune.infiltration.matrix.sort <- apply(pvalue.pcd.immune.infiltration.matrix.sort,2,function(x){ifelse(x >0.05,"",".")})
#pvalue.pcd.immune.infiltration.matrix.sort <- apply(pvalue.pcd.immune.infiltration.matrix.sort,2,function(x){ifelse(x >0.05,"",ifelse(x>0.01,"*",ifelse(x>0.001,"**","***")))})
head(pvalue.pcd.immune.infiltration.matrix.sort)

cohort_annotation <- sapply(colnames(cor.pcd.immune.infiltration.matrix.sort),function(x){strsplit(x,"#")[[1]][1]})
pathway_annotation <- sapply(colnames(cor.pcd.immune.infiltration.matrix.sort),function(x){strsplit(x,"#")[[1]][2]})
pathway_annotation <- sapply(pathway_annotation,function(x){gsub("\\.Innate|\\.Adaptive","",x)})
annotation_col <- data.frame(cohort=cohort_annotation,pathway=pathway_annotation)
head(annotation_col)
row.names(annotation_col) = colnames(cor.pcd.immune.infiltration.matrix.sort)
color_cohort <- mycolpanel[1:length(unique(cohort_annotation))]
names(color_cohort) <- c(unique(annotation_col$cohort))
color_pathway <- distinctColorPalette(length(unique(pathway_annotation)))
names(color_pathway) <- c(unique(annotation_col$pathway))
length(unique(annotation_col$pathway))
head(cor.pcd.immune.infiltration.matrix.sort)
pdf(file.path("PCD_roles_in_cancers",paste0("Pcd_pathways_28_immune_cell_",a,"_infiltration_cor_heatmap_t.pdf")),height = 20,width = 20)
# colnames(cor.pcd.immune.infiltration.matrix.sort) <- sapply(colnames(cor.pcd.immune.infiltration.matrix.sort),function(x){
#   strsplit(x,"#")[[1]][2]
# })
pheatmap(cor.pcd.immune.infiltration.matrix.sort,show_colnames = F,show_rownames = T,
         scale = "none",cluster_rows=T,cluster_cols=F,
         clustering_method = 'complete',
         color = colorRampPalette(c('#21b6af','white','#813e98'))(100),
         #color = colorRampPalette(c('#21b6af','white','#eeba4d'))(100),
         annotation_col = annotation_col,
         annotation_colors = list(cohort=color_cohort,pathway=color_pathway),
         #gaps_col = seq(5,70,5),
         gaps_col = seq(5,60,5),
         display_numbers=pvalue.pcd.immune.infiltration.matrix.sort,fontsize_number=1,
         cellwidth = 6,cellheight = 6,fontsize_row=6,fontsize_col=6)
dev.off()

########################一致性聚类######################################
dir_consensusClust <- "consensus_clust"
dir.create(dir_consensusClust)

########TCGA-LUAD##############
dir_clust <- file.path(dir_consensusClust,"TCGA-LUAD_Clust")
dir.create(dir_clust)

#TCGA-LUAD所有样本存活表达矩阵
TCGA_LUAD_sursive_exp <- all_cohorts_raw_data_mRNA_os_clean$`TCGA-LUAD`
TCGA_LUAD_sursive_exp[1:4,1:4]

TCGA_LUAD_tumor_sursive_exp <- TCGA_LUAD_sursive_exp[grepl("tumor",row.names(TCGA_LUAD_sursive_exp)),]
TCGA_LUAD_tumor_sursive_exp[1:4,1:4]

#计算免疫浸润得分
ssGSEA.res <- ssGSEA(data = TCGA_LUAD_tumor_sursive_exp[,-c(1:2)],need_log = F,genesets = geneSets$`28_Immune_Cell_Infiltration_GeneSets`,scale = T)
ssGSEA.res <- data.frame(t(ssGSEA.res))
row.names(ssGSEA.res) <- gsub("_tumor","",row.names(ssGSEA.res) ) 
ssGSEA.res[1:4,1:4]

ssGSEA.res

data <- ssGSEA.res
#进行一致性聚类
consensusClust(ssGSEA.res,heatmap_color = c("#fdc68a","#05a59e"),max_k = 6,outdir=dir_clust)

#读取聚类结果
clust.res <- read.csv(file.path(dir_clust,"cluster3.csv"))
table(clust.res$Cluster)
#按照免疫得分高-中，低进行从新排序
clust.res[,"Cluster_rename"] <- ifelse(clust.res$Cluster=="C1","C2",
                                       ifelse(clust.res$Cluster=="C2","C1","C3"))
#clust.res[,"Cluster_rename"] <- ifelse(clust.res$Cluster=="C1","C2","C1")
table(clust.res$Cluster_rename)
clust.res$Cluster <- NULL
colnames(clust.res) <- c("Sample","Cluster")
clust.res <- clust.res[order(clust.res$Cluster),]
#clust.res$Cluster <- factor(clust.res$Cluster,levels = c("C1","C2"))
clust.res$Cluster <- factor(clust.res$Cluster,levels = c("C1","C2","C3"))
row.names(clust.res) <- clust.res$Sample
clust.res$Sample <- NULL
head(clust.res)
table(clust.res$Cluster)

#合并免疫得分和免疫分型
imssGSEA_clust <- merge(clust.res,ssGSEA.res,by.x=0,by.y=0)
imssGSEA_clust <- imssGSEA_clust %>% column_to_rownames("Row.names")
imssGSEA_clust <- imssGSEA_clust[order(imssGSEA_clust$Cluster),]
#imssGSEA_clust$Cluster <- factor(imssGSEA_clust$Cluster,levels = c("C1","C2"))
imssGSEA_clust$Cluster <- factor(imssGSEA_clust$Cluster,levels = c("C1","C2","C3"))
write.csv(imssGSEA_clust,file.path(dir_clust,"TCGA_LUAD_imssGSEA_clust.csv"))

imssGSEA_clust <- read.csv(file.path(dir_clust,"TCGA_LUAD_imssGSEA_clust.csv"),row.names = 1)
head(imssGSEA_clust)
table(imssGSEA_clust$Cluster)

#合并临床和免疫聚类分型
TCGA_LUAD_clinical_file <- "TCGA-LUAD_clinical.csv"
clinical <- read.csv(TCGA_LUAD_clinical_file,header = 1,row.names = 1)
head(clinical)
# clinical_tmp <- merge(clinical,TCGA_LUAD_estimateScore,by.x=0,by.y=0)
# row.names(clinical_tmp) <- clinical_tmp$Row.names
# clinical_tmp$Row.names <- NULL
# write.csv(clinical_tmp,"TCGA-LUAD_clinical.csv")
clinical_imclust <- merge(clinical,clust.res,by.x = 0, by.y = 0)
clinical_imclust <- clinical_imclust %>% column_to_rownames("Row.names")
clinical_imclust <- clinical_imclust[order(clinical_imclust$Cluster),]
clinical_imclust$Cluster <- factor(clinical_imclust$Cluster,levels = c("C1","C2","C3"))
clinical_imclust[1:4,1:4]
write.csv(clinical_imclust,file.path(dir_clust,"TCGA_LUAD_clinical_imclust.csv"))

head(clinical_imclust)
clinical_imclust <-  read.csv(file.path(dir_clust,"TCGA_LUAD_clinical_imclust.csv"),row.names = 1)
table(clinical_imclust$Cluster)
#免疫浸润比较分析
colnames(clinical_imclust)
p <- boxtools(clinical_imclust,valueVariable="ImmuneScore_estimate",groupVariable="Cluster"
              ,colors=clust_col,ylab="immuneScore",plot_tyep="violin",prefix = "TCGA-LUAD")

p
ggsave(filename = file.path(dir_clust,"TCGA_LUAD_immuneScore_clust.pdf"),p)

ssGSEA.res <- imssGSEA_clust[,-1] 
clust.res <- clust.res[,1,drop=F]
p <- pca.plot(ssGSEA.res,clust.res,color = clust_col)
p
ggsave(filename = file.path(dir_clust,"TCGA_LUAD_imssGSEA_clust_pca.pdf"),p)

table(imssGSEA_clust$Cluster)

###28免疫浸润比较小提琴图
imssGSEA_clust_tmp_long <- melt(imssGSEA_clust,id="Cluster",variable.names = "cell",value.name  = "value")
head(imssGSEA_clust_tmp_long)
imssGSEA_clust_tmp_long[,"variable"] <- gsub("\\.Adaptive|\\.Innate","",imssGSEA_clust_tmp_long$variable)
table(imssGSEA_clust_tmp_long$Cluster)

p <- multi_item_group_violin_boxplot_plot(data=imssGSEA_clust_tmp_long,
                                          featureVariable="variable",
                                          valueVariable="value",
                                          groupVariable="Cluster",
                                          colors=clust_col,
                                          ylab="Infiltration",
                                          type="boxplot")

p
ggsave(filename = file.path(dir_clust,"TCGA_LUAD_infiltration_difference_between_clust_boxplot.pdf"),p,height = 3,width = 5.5)

#绘制免疫浸润热图
library(pheatmap)
imssGSEA_clust_tmp <- imssGSEA_clust
imssGSEA_clust_tmp <- imssGSEA_clust_tmp[order(imssGSEA_clust_tmp$Cluster),]
colnames(clinical_imclust)
col_names <- c("Cluster","Gender","Age","Stage")
for (i in col_names){
  print(i)
  print(unique(clinical_imclust[,i]))
}

col_annotation <- clinical_imclust[,col_names]
color_Gender <- c(Female='#E0864A',Male='rosybrown',"Unknown"="White")
color_Age <- c("<60"="#7ac7e2",">=60"="#54beaa","Unknown"="White")
color_Stage <- c("Stage I"='cornsilk',"Stage II"='paleturquoise',"Stage III"='goldenrod',"Stage IV"='firebrick',"Unknown"='White')

infiltration <- imssGSEA_clust_tmp[,-1]
colnames(infiltration) <- gsub("\\."," ",colnames(infiltration) )
type_annotation <- sapply(colnames(infiltration),function(x){tail(strsplit(x," ")[[1]],1)})
row_annotation <- data.frame(type=type_annotation,row.names = colnames(infiltration))
col_type <- c(Adaptive = "#B4DAAB", Innate = "#85ceb7")
head(row_annotation)

#clust_col <- c(C1="#C4B4CE",C2="#46AC85")

infiltration <- t(infiltration)
infiltration[1:4,1:4]
graphics.off()
pdf(file.path(dir_clust,"TCGA_LUAD_infiltration_heatmap.pdf"),height = 8,width = 8)
pheatmap(infiltration,show_colnames = F,show_rownames = T,
         scale = "row",cluster_rows=T,cluster_cols=F,
         clustering_method = 'ward',
         color = colorRampPalette(c('#21b6af','white','#eeba4d'))(100),
         #color = colorRampPalette(c('#21b6af','white','#eeba4d'))(100),
         annotation_names_row = F,annotation_names_col = F,
         annotation_row = row_annotation,
         annotation_col = col_annotation,
         annotation_colors = list(Cluster=clust_col,type =col_type,Gender=color_Gender,
                                  Age=color_Age,Stage=color_Stage))
dev.off()

#############################WGCNA分析########################3

dir_wgcna <- "WGCNA"
dir.create(dir_wgcna)

pcd.genes <- unique(unlist(geneSets$`18PCD`))
length(pcd.genes)

########TCGA-LUAD#######
dir_wgcna_TCGA_LUAD <- file.path(dir_wgcna,"TCGA-LUAD")
dir.create(dir_wgcna_TCGA_LUAD)

#进行编码
clinical_code <- read.csv(file.path(root_dir,"TCGA-LUAD_clinical_code.csv"),row.names = 1)
head(clinical_code)

clinical_imclust <- read.csv(file.path(dir_clust,"TCGA_LUAD_clinical_imclust.csv"),row.names = 1)
head(clinical_imclust)
clinical_code <- cbind(clinical_imclust[row.names(clinical_code),"Cluster"],clinical_code)
colnames(clinical_code)[1] <- "Cluster"
clinical_code <- clinical_code[!is.na(clinical_code$Cluster),!colnames(clinical_code) %in% c("OS","OS.time")]
clinical_code$Cluster
clinical_code$Cluster <- as.numeric(gsub("C","",clinical_code$Cluster))
table(clinical_code$Cluster)
head(clinical_code)
sum(is.na(clinical_code))
str(clinical_code)
 
#表达文件
TCGA_LUAD_tumor_sursive_exp[1:4,1:4]
TCGA_LUAD_tumor_wgcna_input <- TCGA_LUAD_tumor_sursive_exp[,colnames(TCGA_LUAD_tumor_sursive_exp) %in% pcd.genes]
row.names(TCGA_LUAD_tumor_wgcna_input) <- gsub("_tumor","",row.names(TCGA_LUAD_tumor_wgcna_input))
TCGA_LUAD_tumor_wgcna_input[1:4,1:4]
dim(TCGA_LUAD_tumor_wgcna_input)

same_sample <- intersect(row.names(clinical_code),row.names(TCGA_LUAD_tumor_wgcna_input))
TCGA_LUAD_tumor_wgcna_input <- TCGA_LUAD_tumor_wgcna_input[same_sample,]
clinical_code <- clinical_code[same_sample,]
TCGA_LUAD_tumor_wgcna_input[1:4,1:4]
head(clinical_code)

#样品聚类
sampleTree <- preWGCNA1(expData=TCGA_LUAD_tumor_wgcna_input,outdir=dir_wgcna_TCGA_LUAD)

#查看样本聚类结果，设置截断高度的cutoff值
clust <- preWGCNA2(sampleTree,cutHeight=50,minSize=10)
clust_sample_size <- data.frame(table(clust))
clust_sample_size
keep_clust <- as.numeric(as.character(clust_sample_size[clust_sample_size$Freq>20,"clust"]))
keepSamples <- clust %in% keep_clust
#查看每个聚类中的样本数，将保留的簇传递给runWGCNA
runWGCNA(expData=TCGA_LUAD_tumor_wgcna_input,traitData=clinical_code,keepSamples=keepSamples,outdir=dir_wgcna_TCGA_LUAD)

#查看和免疫分型相关性最高得模块，并获得模块内得基因
best_models <- "9_blue_genes.txt"
WGCNA_genes <- read.table(file.path(dir_wgcna_TCGA_LUAD,best_models),sep="\t")
WGCNA_genes <- WGCNA_genes[,1]
WGCNA_genes

best_models <- "9_red_genes.txt"
WGCNA_genes2 <- read.table(file.path(dir_wgcna_TCGA_LUAD,best_models),sep="\t")
WGCNA_genes2 <- WGCNA_genes2[,1]
WGCNA_genes2

WGCNA_genes <- c(WGCNA_genes,WGCNA_genes2)
length(WGCNA_genes)

#进一步对WGCNA_genes进行univariable cox
unicox_input <- TCGA_LUAD_sursive_exp[grepl("tumor",row.names(TCGA_LUAD_sursive_exp)),c("OS","OS.time",WGCNA_genes)]
unicox_input[1:4,1:4]
dim(unicox_input)

unicox.res<- multiple.item.unicox.analysis(unicox_input,Binning=TRUE)
unicox.res_tmp <- unicox.res[!unicox.res$pvalue=="",]
unicox.res.sig <- unicox.res_tmp[as.numeric(unicox.res_tmp$pvalue) < 0.05,]
dim(unicox.res.sig)

dir_wgcna_TCGA_LUAD <- "WGCNA\\TCGA-LUAD"
write.csv(unicox.res.sig,file.path(dir_wgcna_TCGA_LUAD,"TCGA_LUAD_WGCNA_genes_unicox_sig.csv"))
wgcna.unicox.res.sig <- read.csv(file.path(dir_wgcna_TCGA_LUAD,"TCGA_LUAD_WGCNA_genes_unicox_sig.csv"))
wgcna.unicox.res.sig.genes <- wgcna.unicox.res.sig$item
wgcna.unicox.res.sig.genes

###进行101中机器学习算法建模
dir_ml_1 <- "101_ml"
dir.create(dir_ml_1)

#获取肿瘤样本的存活表达矩阵
all_cohorts_raw_data_mRNA_tumor_os_clean <- lapply(all_cohorts_raw_data_mRNA_os_clean, function(x){x[grepl("tumor",row.names(x)),]})
#获取感兴趣的特征的存活表达矩阵
list_train_vali_Data <- lapply(all_cohorts_raw_data_mRNA_tumor_os_clean, function(x){x[,c("OS","OS.time",wgcna.unicox.res.sig.genes)]})
lapply(list_train_vali_Data,function(x){x[1:4,1:4]})
lapply(list_train_vali_Data,function(x){dim(x)})
lapply(list_train_vali_Data,function(x){sum(is.na(x))})
sum(is.na(list_train_vali_Data[["TCGA-LUAD"]]$OS))

#模型构建
Machine_Learning_101_Algorithms_for_Prognostic_Models(list_train_vali_Data=list_train_vali_Data,
                                                      train_cohort_name="TCGA-LUAD",
                                                      barCol="#FDA481",
                                                      heatmapCol=c("#46AC85", "#FFFFFF","#C4B4CE" ),
                                                      outdir="101_ml")

#选择最好的模型进行生存和ROC曲线分析
#设置输出目录
dir_survival <- file.path("101_ml","survival")
dir.create(dir_survival )
ROC_1_3_5_colors <- c("#F05C3BFF","#5C88DAFF","#5CB85CFF")

risk_matrix_tmp <- read.table(file.path(dir_ml,"model.riskMatrix.txt"),header=T, sep="\t", row.names=1, check.names=F, stringsAsFactors=F)
head(risk_matrix_tmp)
best_ml_model <- "Lasso+plsRcox"
risk_matrix <- risk_matrix_tmp[,c("Cohort","OS","OS.time",best_ml_model)]
head(risk_matrix)

#结果保存
risk_matrix_group <- data.frame()
for (cohort in unique(risk_matrix$Cohort)){
  #cohort=unique(risk_matrix$Cohort)[2]
  risk_matrix_tmp <- risk_matrix[risk_matrix$Cohort==cohort,c("OS","OS.time",best_ml_model)]
  colnames(risk_matrix_tmp)[3] <- "RS"
  
  #调用绘制生存曲线的函数
  survival_res <- surviva_ml(risk_matrix_tmp,color=risk_col_c,title = cohort)
  risk_matrix_tmp <- cbind(risk_matrix_tmp,survival_res[["risk-group"]])
  head(risk_matrix_tmp)
  colnames(risk_matrix_tmp)[4] <- "RR"
  pdf(file=file.path(dir_survival,paste0(cohort,"_survival_",survival_res[["cutoff"]],"_cutoff.pdf")),height = 6,width = 6)
  print(survival_res[["plot"]],newpage = FALSE)
  dev.off()
  
  #保存生存风险分组
  risk_matrix_tmp <- risk_matrix_tmp %>% dplyr::mutate(cohort=rep(cohort,nrow(risk_matrix_tmp)))
  risk_matrix_group <- rbind(risk_matrix_group,risk_matrix_tmp) 
  
  #调用绘制ROC的函数
  ROC_Curves_for_Risk_Scores_at_1_3_5_Years(risk_matrix_tmp,futimeVariable="OS.time",
                                            fustatVariable="OS",
                                            riskscoreVariable="RS",
                                            colors=ROC_1_3_5_colors,
                                            outdir=dir_survival,prefix = cohort)
}
write.csv(risk_matrix_group,file.path("101_ml","model.riskGroup.csv"),quote = F)

#计算cindex
#选择最好的模型进行生存和ROC曲线分析
dir_survival <- file.path("101_ml","survival")
dir.create(dir_survival )
ROC_1_3_5_colors <- roc_1_3_5_col

risk_matrix_tmp <- read.table(file.path("101_ml","model.riskMatrix.txt"),header=T, sep="\t", row.names=1, check.names=F, stringsAsFactors=F)
head(risk_matrix_tmp)
best_ml_model <- "Lasso+plsRcox"
risk_matrix <- risk_matrix_tmp[,c("Cohort","OS","OS.time",best_ml_model)]
head(risk_matrix)

#结果保存
risk_matrix_group <- data.frame()
for (cohort in unique(risk_matrix$Cohort)){
  #cohort=unique(risk_matrix$Cohort)[2]
  risk_matrix_tmp <- risk_matrix[risk_matrix$Cohort==cohort,c("OS","OS.time",best_ml_model)]
  colnames(risk_matrix_tmp)[3] <- "RS"
  
  #调用绘制生存曲线的函数
  survival_res <- surviva_ml(risk_matrix_tmp,color=risk_col_c,title = cohort)
  risk_matrix_tmp <- cbind(risk_matrix_tmp,survival_res[["risk-group"]])
  head(risk_matrix_tmp)
  colnames(risk_matrix_tmp)[4] <- "RR"
  pdf(file=file.path(dir_survival,paste0(cohort,"_survival_",survival_res[["cutoff"]],"_cutoff.pdf")),height = 3,width = 3.5)
  print(survival_res[["plot"]],newpage = FALSE)
  dev.off()
  
  #保存生存风险分组
  risk_matrix_tmp <- risk_matrix_tmp %>% dplyr::mutate(cohort=rep(cohort,nrow(risk_matrix_tmp)))
  risk_matrix_group <- rbind(risk_matrix_group,risk_matrix_tmp) 
  
  #调用绘制ROC的函数
  ROC_Curves_for_Risk_Scores_at_1_3_5_Years(risk_matrix_tmp,futimeVariable="OS.time",
                                            fustatVariable="OS",
                                            riskscoreVariable="RS",
                                            colors=ROC_1_3_5_colors,
                                            outdir=dir_survival,prefix = cohort)
}
risk_matrix_group$RR <- gsub("high","High-risk",risk_matrix_group$RR)
risk_matrix_group$RR <- gsub("low","Low-risk",risk_matrix_group$RR)
risk_matrix_group <- risk_matrix_group[,c("cohort","OS","OS.time","RS","RR")]
row.names(risk_matrix_group) <- gsub("_tumor","",row.names(risk_matrix_group))
write.csv(risk_matrix_group,file.path("101_ml","model.riskGroup.csv"),quote = F)


#模型中的系数
dir_gene_role <- "model_gene_role"
dir.create(dir_gene_role)
model <- readRDS(file.path("101_ml","model.MLmodels.rds"))
model_tmp <- model$`Lasso+plsRcox`

dat <- read.csv(file.path("101_ml","Coefficients.csv"))
dat <- dat[order(dat$coefficient,decreasing = T),]
dat$gene <- factor(dat$gene,levels = unique(dat$gene))
dat$group <- ifelse(dat$coefficient>0,"high","low")
head(dat)
p <- ggplot(dat, aes(x=gene, y=coefficient)) +
  geom_segment( aes(x=gene, xend=gene, y=0, yend=coefficient), color="grey",linewidth = 2) +
  geom_point( aes(color = group), size=4) +
  scale_color_manual(values = c("#fdc68a","#05a59e"))+ 
  theme_light() + theme(panel.grid.major = element_blank(),
                        panel.grid.minor = element_blank(),
                        axis.line = element_line(color = "black",size = 0.4),
                        axis.text.y = element_text(color = "black",size = 6),
                        axis.text.x = element_text(color = "black",size = 6,angle = 45),
                        axis.title = element_text(color = "black",size = 6),
                        legend.position = "none",
                        plot.title = element_text(hjust = 0.5,size = 8))
p
ggsave(filename = file.path(dir_gene_role,"model_coefficient.pdf"),p,height = 3,width = 5)

#模型中基因单因素分析
genes_in_model <- as.character(dat$gene)
TCGA_LUAD_tumor_sur_exp <- all_cohorts_raw_data_mRNA_os_clean$`TCGA-LUAD`
TCGA_LUAD_tumor_sur_exp <- TCGA_LUAD_tumor_sur_exp[grepl("tumor",row.names(TCGA_LUAD_tumor_sur_exp)),]
TCGA_LUAD_tumor_sur_exp[1:4,1:6]
data <- TCGA_LUAD_tumor_sur_exp[,c("OS","OS.time",genes_in_model)]

library(ezcox)
#修改因子水平，将contrast放到前面，ref_level放到后面
zz <- ezcox(data,covariates = colnames(data)[-c(1,2)],
            time = "OS.time",status = "OS",return_models=TRUE)
mds <- get_models(zz)
pdf(file.path(dir_gene_role,paste0("TCGA-LUAD","_unicox.pdf")),height = 10,width = 8)
show_models(mds)
dev.off()

#模型中基因来源弦图
#整理PCD相关基因来源
pcd.genes <- geneSets$`18PCD`
df <- data.frame()
for (i in names(pcd.genes)){
  #i=names(pcd.genes)[1]
  tmp <- data.frame(pcd=rep(i,length(pcd.genes[[i]])),genes=pcd.genes[[i]])
  df <- rbind(df,tmp)
}
#获取模型中的基因
dat <- read.csv(file.path("101_ml","Coefficients.csv"))
dat <- dat[order(dat$coefficient,decreasing = T),]
order_gene <- factor(dat$gene,levels = unique(dat$gene))
model_in_genes <- dat$gene
pdata <- df[df$genes %in% model_in_genes,]
colnames(pdata) <- c("pathway","gene")
table(pdata$pathway)
head(pdata)
write.csv(pdata,file = file.path("model_gene_role","gene_source.csv"))
chord_plot(pdata=pdata,outdir = "model_gene_role","101_model_gene_pcd")
gc()
names(all_cohorts_raw_data_mRNA_tumor_normal_expression_clean)
#比较模型中的基因在癌和癌症中是否有差异
dir.create(file.path("model_gene_role","tumor_normal_exp_diff"))
for (cohort in c("TCGA−LUAD","CPTAC−LUAD","GPL570_GSE31210")){
  model_gene_exp <- all_cohorts_raw_data_mRNA_tumor_normal_expression_clean$GPL570_GSE31210
  model_gene_exp <- model_gene_exp[,model_in_genes]
  model_gene_exp[1:4,1:4]
  model_gene_exp[,"group"] <- ifelse(grepl("tumor",row.names(model_gene_exp)),"tumor","normal")
  model_gene_exp_long <- melt(model_gene_exp,ids="group",value.name = "value",variable.names="genes")
  model_gene_exp_long$variable <- factor(model_gene_exp_long$variable,levels = order_gene)
  head(model_gene_exp_long)
  p <- multi_item_group_violin_boxplot_plot(model_gene_exp_long,featureVariable = "variable",
                                            valueVariable = "value",groupVariable = "group",colors = tumor_normal_col,ylab = "gene expression",type = "boxplot")
  
  ggsave(file=file.path("model_gene_role","tumor_normal_exp_diff",paste0("GPL570_GSE31210","_all_expression_diff.pdf")),plot = p,width = 5.5,height = 3)
  # for (gene in model_in_genes){
  #   #gene=model_in_genes[1]
  #   p <- boxtools(model_gene_exp,valueVariable = gene,groupVariable = "group",
  #                 colors =tumor_normal_col,ylab = gene,plot_tyep = "slimboxplot",prefix = gene )
  #   ggsave(file=file.path("model_gene_role","tumor_normal_exp_diff",paste0(gene,"_expression_diff.pdf")),plot = p,width = 3,height = 4)
  # }
}


#比较模型中的基因CPTAC-LUAD蛋白表达的差异
length(order_gene)
protein_exp <- data.frame(fread(file.path("CPTAC_LUAD_tumor_normal_protein_impute_knn_log2.csv"),sep=",",check.names = F,stringsAsFactors = F),row.names=1)
protein_exp[1:4,1:4]
row.names(protein_exp)
order_gene <- intersect(order_gene,colnames(protein_exp))
order_gene <- factor(order_gene,levels = order_gene)
order_gene
model_gene_protein <- protein_exp[,c("PTTG1IP","GNAI3","KRT8","CDCP1","BTK","FUCA1","ATP6V1B2","SIRT2","RNF5","PEBP1")]
#model_gene_protein <- protein_exp[,order_gene]
head(model_gene_protein)
model_gene_protein[,"group"] <- ifelse(grepl("tumor",row.names(model_gene_protein)),"tumor","normal")
head(model_gene_protein)
model_gene_protein_long <- melt(model_gene_protein,ids="group",value.name = "value",variable.names="genes")
head(model_gene_protein_long)
#model_gene_protein_long$variable <- factor(model_gene_protein_long$variable,levels = order_gene)
head(model_gene_protein_long)
p <- multi_item_group_violin_boxplot_plot(model_gene_protein_long,featureVariable = "variable",
                                          valueVariable = "value",groupVariable = "group",colors = tumor_normal_col,ylab = "protein expression",type = "splitviolin")

p
ggsave(file=file.path("model_gene_role","tumor_normal_exp_diff",paste0("CPTAC-LUAD","_all_protein_expression_diff.pdf")),plot = p,width = 10,height = 6)

#CPTC-LUAD蛋白表达和预后的关系
CPTAC_LUAD_survival <- all_cohorts_raw_data_mRNA_os_clean$`CPTAC-LUAD`[,c(1,2)]
CPTAC_LUAD_survival <- CPTAC_LUAD_survival[grepl("tumor",row.names(CPTAC_LUAD_survival)),]
head(CPTAC_LUAD_survival)

proteins <- c("PTTG1IP","GNAI3","KRT8","CDCP1","BTK","FUCA1","ATP6V1B2","SIRT2","RNF5","PEBP1")
protein_exp_tmp <- protein_exp[grepl("tumor",row.names(protein_exp)),proteins]
head(protein_exp_tmp)
protein_exp_emp_survival <- merge(CPTAC_LUAD_survival,protein_exp_tmp,by.x = 0,by.y=0)
protein_exp_emp_survival <- protein_exp_emp_survival %>% tibble::column_to_rownames("Row.names")
head(protein_exp_emp_survival)

dir.create(file.path("model_gene_role","CPTAC-LUAD_protein_exp_survival"))
for (i in proteins){
  cut <- surv_cutpoint(protein_exp_emp_survival,'OS.time','OS',i,minprop = 0.1)
  cat <- surv_categorize(cut)
  head(cat)
  pdf(file=file.path("model_gene_role","CPTAC-LUAD_protein_exp_survival",paste0(i,"_protein_exp_survival.pdf")),height = 5,width = 5)
  res <- single_item_survival_analysis(cat,groupVariable=i,need_binary=F,colors=tumor_normal_col,ylab="survival probability (OS)",risk.table=F)
  print(res$plot,newpage = F)
  dev.off()
}


#CPTC-LUAD基因表达和预后的关系
CPTAC_LUAD_survival <- all_cohorts_raw_data_mRNA_os_clean$`CPTAC-LUAD`[,]
CPTAC_LUAD_survival <- CPTAC_LUAD_survival[grepl("tumor",row.names(CPTAC_LUAD_survival)),c("OS","OS.time","PTTG1IP","CDCP1","KRT8")]
CPTAC_LUAD_survival[1:4,1:4]

for (i in c("PTTG1IP","CDCP1","KRT8")){
  cut <- surv_cutpoint(CPTAC_LUAD_survival,'OS.time','OS',i,minprop = 0.1)
  cat <- surv_categorize(cut)
  head(cat)
  pdf(file=file.path("model_gene_role",paste0(i,"_mRNA_exp_survival.pdf")),height = 5,width = 5)
  res <- single_item_survival_analysis(cat,groupVariable=i,need_binary=F,colors=tumor_normal_col,ylab="survival probability (OS)",risk.table=F)
  print(res$plot,newpage = F)
  dev.off()
}

#预测免疫治疗疗效
dir_predict_im_res <- "predict_icd_response"
dir.create(dir_predict_im_res)
model_path <- file.path("101_ml","model.MLmodels.rds")
method <- "Lasso+plsRcox"
#调用函数
predict_immunecheckpoint_inhibitor_response(model=model_path,method=method,outdir=dir_predict_im_res,mycolpanel=mycolpanel)


##模型性能比较
risk_matrix_group <- read.csv(file.path("101_ml","model.riskGroup.csv"))
library(tibble)
dir_compare <- "compare_performance"
dir.create(dir_compare)
dir_compre_clinical <- file.path(dir_compare,"compre_clinical")
dir.create(dir_compre_clinical )

#合并临床编码文件
clinical_code <- data.frame()
for (cohort in unique(risk_matrix_group$cohort)){
  #cohort <- unique(risk_matrix_group$cohort)[1]
  clinical_code_tmp <-  read.csv(file.path(clinical_code_dir,paste0(cohort,"_clinical_code.csv")),row.names = 1)
  clinical_code_tmp <- clinical_code_tmp[,!colnames(clinical_code_tmp) %in% c("OS","OS.time")]
  clinical_code <- bind_rows(clinical_code,clinical_code_tmp)
}
head(clinical_code)
write.csv(clinical_code,file.path(dir_compre_clinical,"clinical_code.csv"))

#合并临床编码和风险分组文件
risk_matrix_group <- read.csv(file.path("101_ml","model.riskGroup.csv"),row.names=1)
head(risk_matrix_group)
risk_clinical <- merge(risk_matrix_group,clinical_code,by.x=0,by.y=0)
risk_clinical <- risk_clinical %>% tibble::column_to_rownames("Row.names")
dim(risk_clinical)
head(risk_clinical)
#进行临床特征性能比较
for (cohort in unique(risk_clinical$cohort)){
  #cohort <- unique(risk_clinical$cohort)[2]
  risk_clinical_tmp <- risk_clinical[risk_clinical$cohort==cohort,]
  risk_clinical_tmp <- risk_clinical_tmp[,!apply(risk_clinical_tmp,2,
                                               function(x){all(is.na(x))})]
  risk_clinical_tmp$RR <- NULL
  cindex_compared_with_clinical_barplot(risk_clinical=risk_clinical_tmp,
                                              index_risk_score=4,
                                              outdir=dir_compre_clinical,
                                              prefix=cohort,only_plot = T)
}

risk_matrix_group <- read.csv(file.path("101_ml","model.riskGroup.csv"),row.names=1)
unique(risk_matrix_group$cohort)
head(risk_matrix_group)
dir_compre_signature <- file.path("compare_performance","compred_with_signature")
dir.create(dir_compre_signature)

#模型和已发表的特征性能比较
cohorts <- unique(risk_matrix_group$cohort)
for (i in 1:length(cohorts)){
  #i=1
  pal <- brewer.pal(10, colorbrewers[i])
  cohort <- cohorts[i]
  #cohort <- unique(risk_matrix_group$cohort)[2]
  risk_matrix_group_tmp <-risk_matrix_group[risk_matrix_group$cohort==cohort,]
  risk_matrix_group_tmp$RR <- NULL
  risk_matrix_group_tmp$cohort <- NULL
  head(risk_matrix_group_tmp)
  
  #因为多个队列merge矩阵gene太少，读取原始survival_exp矩阵
  files <- list.files(path = clean_data_dir ,pattern = "*mRNA*")
  file <- files[grepl(cohort,files)]
  survival_exp_tmp <- data.frame(fread(file.path(clean_data_dir,file),sep=","),row.names=1)
  
  #survival_exp_tmp <- all_cohorts_raw_data_mRNA_os_clean[[cohort]]
  survival_exp_tmp <- survival_exp_tmp[grepl("tumor",row.names(survival_exp_tmp)),]
  row.names(survival_exp_tmp) <- gsub("_tumor","",row.names(survival_exp_tmp))
  survival_exp_tmp[1:4,1:4]
  
  cindex_compared_with_paper_signature(cohort_name = cohort,
                                       rise_score = risk_matrix_group_tmp,
                                       survival_exp = survival_exp_tmp,
                                       index_riskscore = 3,
                                       plot_type = "lbbt",
                                       #colors = c(pal[4],pal[7]),
                                       colors = mycolpanel[i],
                                       outdir = dir_compre_signature,
                                       sheet_name_in_baseline = "LUAD",
                                       del_model=F,
                                       pmids = keep_pmids,
                                       only_plot = T)
}

###########################风险得分和免疫微环境相关性分析########################
#准备TIDE输入文件
library(data.table)
library(tidyverse)
dir.create("TIDE")
for (cohort in names(all_cohorts_raw_data_mRNA_tumor_normal_expression_clean)){
  #cohort = names(all_cohorts_raw_data_mRNA_tumor_normal_expression_clean)[1]
  exp <- all_cohorts_raw_data_mRNA_tumor_normal_expression_clean[[cohort]]
  tumor_exp <- exp[grepl("tumor",row.names(exp)),]
  tumor_exp <- t(tumor_exp)
  mydat <- t(apply(tumor_exp, 1, function(x)x-(mean(x))))
  tmp <- data.frame(gene=row.names(mydat))
  mydat <- cbind(tmp,mydat)
  mydat[1:4,1:4]
  write.table(mydat,file = file.path("TIDE",paste0(cohort,'TIDE.input.txt')),sep = '\t',row.names = F) ##记得用制表符分隔
  
}

#基质得分和免疫得分差异分析
dir.create("TIME")
dir.create(file.path("TIME","Estimate"))
for (cohort in annotation_cohorts){
  #cohort <- annotation_cohorts[1]
  riskScore <- read.csv(file.path("101_ml","model.riskGroup.csv"),row.names = 1)
  rs <- riskScore[riskScore$cohort==cohort,]
  rs <- rs[,c("cohort","RS","RR")]
  clinical <- read.csv(file.path(root_data_dir,paste0(cohort,"_clinical.csv")),row.names=1)
  clinical_riskScore <- merge(clinical,rs,by.x=0,by.y = 0)
  write.csv(clinical_riskScore,file=paste0(cohort,"_clinical_riskscore.csv"))
  score_col <- c("StromalScore_estimate","ImmuneScore_estimate","TumorPurity_estimate")
  data <- clinical_riskScore[,c("RS","RR",score_col)]
  data <- data[order(data$RR,decreasing = T),]
  data$RR <- factor(data$RR,levels=c("Low-risk","High-risk"))
  colnames(data) <- c("RS","RR","StromalScore","ImmuneScore","TumourPurity")
  head(data)
  for (j in colnames(data)[-c(1:2)]){
    plot_list <- list()
    #j = colnames(data)[-c(1:2)][3]
    p <- boxtools(data=data,valueVariable=j,groupVariable="RR",
                  colors=risk_col_l,j,"violin",signif_label = "p.format",prefix=cohort)
    p
    plot_list[["diff"]] <- p
    p2 <- vector_vector_association_analysis(data=data,variable_x="RS",variable_y=j,
                                             groupvariable="RR",
                                             group_color = risk_col_l,
                                             point_border_color="#231815",
                                             point_fill_color="#cdcac9",
                                             line_color="#00aeef",
                                             line_margin_color="#b2e7fa",
                                             xlab="RS",
                                             ylab=j,margin=FALSE,
                                             prefix = cohort)
    p2
    plot_list[["cor"]] <- p2$plot
    p <- grid.arrange(grobs = plot_list,ncol=2 )
    ggsave(filename = file.path("TIME","Estimate",paste0(cohort,"_",j,"_diff_between_high_low_risk.pdf")),p,width=5.5,height = 2,dpi=300)
  }
}

#TIDE差异分析
dir.create(file.path("TIME","TIDE"))
annotation_cohorts <- c("TCGA-LUAD","CPTAC-LUAD","GPL570_GSE31210","GPL570_GSE50081","GPL570_GSE37745")
for (cohort in annotation_cohorts){
  #cohort <- annotation_cohorts[1]
  riskScore <- read.csv(file.path("101_ml","model.riskGroup.csv"),row.names = 1)
  rs <- riskScore[riskScore$cohort==cohort,]
  rs <- rs[,c("cohort","RS","RR")]
  head(rs)
  TIDE <- read.csv(file.path("E:\\流程搭建\\多组学聚类分型\\文档\\TIDE",paste0(cohort,"_TIDE.csv")),row.names=1)
  row.names(TIDE) <- gsub("_tumor","",row.names(TIDE))
  TIDE_riskScore <- merge(TIDE,rs,by.x=0,by.y = 0)
  #write.csv(TIDE_riskScore,file=paste0(cohort,"_TIDE_riskscore.csv"))
  data <- TIDE_riskScore[,c("Row.names","RR","RS","TIDE")]
  data <- data[order(data$RR,decreasing = T),]
  data$RR <- factor(data$RR,levels=c("Low-risk","High-risk"))
  head(data)
  plot_list <- list()
  p1 <- boxtools(data=data,valueVariable="TIDE",groupVariable="RR",
                colors=risk_col_l,"TIDE score","box_violin",signif_label = "p.format",prefix=cohort)
  p1
  plot_list[["diff"]] <- p1

  p2 <- vector_vector_association_analysis(data=data,variable_x="RS",variable_y="TIDE",
                                          groupvariable="RR",
                                          group_color = risk_col_l,
                                          point_border_color="#231815",
                                          point_fill_color="#cdcac9",
                                          line_color="#00aeef",
                                          line_margin_color="#b2e7fa",
                                          xlab="RS",
                                          ylab="TIDE",margin=FALSE,prefix = cohort)
  p2$plot
  plot_list[["cor"]] <- p2$plot
  p <- grid.arrange(grobs = plot_list,ncol=2 )
  ggsave(filename = file.path("TIME","TIDE",paste0(cohort,"_Risk_TIDE.result.pdf")),p,width=5.5,height = 2,dpi=300)
}

#TMB和HRD和风险得分差异关联分析
dir.create(file.path("TIME","TMB"))
for (cohort in c("TCGA-LUAD", "CPTAC-LUAD")){
  #cohort <- "TCGA-LUAD"
  riskScore <- read.csv(file.path("101_ml","model.riskGroup.csv"),row.names = 1)
  rs <- riskScore[riskScore$cohort==cohort,]
  rs <- rs[,c("cohort","RS","RR")]
  head(rs)
  clinical <- read.csv(file.path(root_data_dir,paste0(cohort,"_clinical.csv")),row.names=1)
  colnames(clinical)
  for (i in c("HRD_Score","Nonsynonymous_Mutations")){
   # i = "HRD_Score"
    plot_list <- list()
    if (sum(!is.na(clinical[i]))>30){
      clinical_tmp <- clinical[!is.na(clinical[i]),i,drop=F]
      head(clinical_tmp)
      data <- merge(clinical_tmp,rs,by.x=0,by.y = 0)
      data <- data[order(data$RR,decreasing = T),]
      data$RR <- factor(data$RR,levels=c("Low-risk","High-risk"))
      p <- boxtools(data=data,valueVariable=i,groupVariable="RR",
                    colors=risk_col_l,i,"box_violin",signif_label = "p.format",prefix = cohort)
      plot_list[["diff"]] <- p
      p2 <- vector_vector_association_analysis(data=data,variable_x="RS",variable_y=i,
                                               groupvariable="RR",
                                               group_color = risk_col_l,
                                               point_border_color="#231815",
                                               point_fill_color="#cdcac9",
                                               line_color="#00aeef",
                                               line_margin_color="#b2e7fa",
                                               xlab="RS",
                                               ylab=i,margin=FALSE,prefix = cohort)
      plot_list[["cor"]] <- p2$plot
      p <- grid.arrange(grobs = plot_list,ncol=2)
      ggsave(filename = file.path("TIME","TMB",paste0(cohort,"_",i,"_diff_between_high_low_risk.pdf")),p,width=5.5,height = 2,dpi=300)
    }
  }
}

type <- unlist(lapply(row.names(TCGA_raw_exp_fpkm),function(x){strsplit(x,"-")[[1]][4]}))
TCGA_raw_exp_fpkm_tumor <- TCGA_raw_exp_fpkm[grepl("01A",type),]
TCGA_raw_exp_fpkm_tumor[,"sample"] <- substr(row.names(TCGA_raw_exp_fpkm_tumor),1,12)
TCGA_raw_exp_fpkm_tumor <- TCGA_raw_exp_fpkm_tumor %>% dplyr::distinct(sample,.keep_all=T)
row.names(TCGA_raw_exp_fpkm_tumor) <- substr(row.names(TCGA_raw_exp_fpkm_tumor),1,12)
dim(TCGA_raw_exp_fpkm_tumor)
TCGA_raw_exp_fpkm_tumor[1:4,1:4]
length(unlist(geneSets$ICDs_GeneSet))

icd_genes_tmp <- unlist(geneSets$`75_ICDs_GeneSet`)
icd_genes <- intersect(icd_genes_tmp,colnames(TCGA_raw_exp_fpkm_tumor))
length(icd_genes)

col_annotation <- geneSets$`75_ICDs_Gene_pathway`
col_annotation <- col_annotation[icd_genes,,drop=F]
head(col_annotation)
col_annotation_color <- mycolpal2[1:length(type)]
names(col_annotation_color) <- unique(row_annotation$pathway)

icd_genes_exp <- TCGA_raw_exp_fpkm_tumor[,icd_genes]
icd_genes_exp[1:4,1:4]
cohort <- "TCGA-LUAD"
#cohort <- annotation_cohorts[1]
riskScore <- read.csv(file.path("101_ml","model.riskGroup.csv"),row.names = 1)
rs <- riskScore[riskScore$cohort==cohort,]
head(rs)

same_sample <- intersect(row.names(icd_genes_exp),row.names(rs))
rs <- rs[same_sample,"RS",drop=F]
head(rs)
icd_genes_exp <- icd_genes_exp[same_sample,]
icd_genes_exp[1:4,1:4]

cor.res <- cor_test_between_two_matrixes(matrix1=rs,
                                         matrix2=icd_genes_exp,
                                         matrix1_omic_name="Risk score",
                                         matrix2_omic_name="ICDs",outdir=NULL,outprefix=NULL,
                                                     cluster_row = TRUE,cluster_col = TRUE,plot=F)
head(cor.res)
cor <- cor.res %>% dplyr::select(c("ICDs","r")) %>% tibble::column_to_rownames("ICDs") %>% t()
head(cor)

p <- cor.res %>% dplyr::select(c("ICDs","p_signif")) %>% tibble::column_to_rownames("ICDs") %>% t()
head(p)

# 绘制热图
pdf(file=file.path("TIME",paste0("TCGA-LUAD_ICD_Risk_score_cor_heatmap2.pdf")),height = 4,width = 10)
pheatmap(cor, color = colorRampPalette(c('#21b6af','white','#813e98'))(100),,cluster_rows = F,cluster_cols = F, 
         legend = T, show_rownames = F,show_colnames = T,display_numbers = p,
         cellwidth = 8,cellheight = 8,fontsize_row=6,fontsize_col=6,annotation_col = row_annotation)
dev.off()

#高低风险组在免疫细胞浸润上的差异
##免疫浸润热图####
library(ComplexHeatmap)
library(dplyr)
library(circlize)
library(psych)
library(data.table)
library(ggpubr)
library(RColorBrewer)
library(patchwork)

clinical <- read.csv("TCGA-LUAD_clinical_riskscore.csv",row.names = 2)

infiltration_estimation_file <- "TCGA_infiltration_estimation_for_tcga(01样本).csv"
infilt_dat <- read.csv(infiltration_estimation_file,sep=",",row.names = 1)
colnames(infilt_dat) <- gsub("\\.\\.","\\.",colnames(infilt_dat))
infilt_dat[1:4,1:4]
#进一步筛选进行展示
#infilt_dat <- infilt_dat[,grepl("TIMER$|MCPCOUNTER$|EPIC$|ESTIMATE$",colnames(infilt_dat))] 

#合并ESTIMATE得分
TCGA_ESTIMATE <- clinical[,grepl("estimate",colnames(clinical))]
colnames(TCGA_ESTIMATE) <- gsub("estimate","ESTIMATE",colnames(TCGA_ESTIMATE))
infilt_dat <- merge(infilt_dat,TCGA_ESTIMATE,by.x = 0,by.y = 0)
infilt_dat <- infilt_dat %>% column_to_rownames("Row.names")
colnames(infilt_dat)

#获得共有的样本
sample <- intersect(row.names(infilt_dat),row.names(clinical))
clinical <- clinical[sample,,drop=F]
clinical <- clinical[order(clinical$RR,decreasing = F),]
clinical$RR

infilt_dat <- infilt_dat[row.names(clinical),]
group <-clinical[,"RR",drop=F]
head(infilt_dat)
dim(group)

dir.create(file.path("TIME","different_method_infiltration"))
#不同方法比较两组免疫浸润小提琴图
methods = unique(sapply(strsplit(colnames(infilt_dat),"_"),tail,1))
for (method in methods){
  #method <- methods[6]
  infilt_dat_tmp <- infilt_dat[,grepl(method,colnames(infilt_dat))]
  infilt_dat_tmp <- merge(group,infilt_dat_tmp,by.x = 0,by.y = 0)
  infilt_dat_tmp_long <- melt(infilt_dat_tmp,id.vars = c("Row.names","RR"),value.name = "value",variable.name = "cell.type")
  head(infilt_dat_tmp_long)
  infilt_dat_tmp_long$RR <- factor(infilt_dat_tmp_long$RR,levels = c("Low-risk","High-risk"))
  infilt_dat_tmp_long$cell.type <- gsub(method,"",infilt_dat_tmp_long$cell.type)
  infilt_dat_tmp_long$cell.type <- gsub("\\.$|_","",infilt_dat_tmp_long$cell.type)
  n_cell_type <- length(unique(infilt_dat_tmp_long$cell.type))
  width <- ifelse(n_cell_type>=30,12,ifelse(n_cell_type>=20,10,ifelse(n_cell_type>=10,8,6)))
  p <- multi_item_group_violin_boxplot_plot(data=infilt_dat_tmp_long,
                                            featureVariable="cell.type",
                                            valueVariable="value",
                                            groupVariable="RR",
                                            colors=risk_col_c,
                                            ylab="infiltration",
                                            type="boxplot")
  p
  ggsave(filename = file.path("TIME","different_method_infiltration",paste0("TCGA-LUAD","high_low_risk_group_infiltration_diff_boxplot_",method,".pdf")),
         height = 4,width = width )
}
###获取CIBERSORT方法的免疫浸润得分有差异的箱线图###
infilt_dat_CIBERSORT <- infilt_dat[,grepl("CIBERSORT$",colnames(infilt_dat))]


#查看矩阵中最大和最小的前n个值,自己定义,然后丢弃 对应的(行)样本
# min_values <- sort(infilt_dat, decreasing = FALSE)[1:20]
# min_values
# max_values <-  sort(infilt_dat, decreasing = TRUE)[1:20]
# max_values
# discard_values <- c(max_values)
# discard_values <- c(min_values,max_values)
# discard_samples_index <- c()
# for (i in discard_values){
#   #i=discard_values[2]
#   indices <- which(infilt_dat == i, arr.ind = TRUE)
#   discard_index <- indices[1,1]
#   discard_samples_index <- c(discard_samples_index,discard_index)
# }
# discard_samples_index <- unique(discard_samples_index)
# infilt_dat <- infilt_dat[,-c(discard_samples_index)]
# dim(infilt_dat)

sample <- intersect(colnames(infilt_dat),row.names(clinical))
clinical <- clinical[sample,,drop=F]
clinical <- clinical[order(clinical$RS,decreasing = F),]
infilt_dat <- infilt_dat[,row.names(clinical)]
group <-clinical[,"RR",drop=F]
head(infilt_dat)
dim(group)

t.test.res = as.data.frame(matrix(nrow = nrow(infilt_dat),ncol = length(unique(group$RR))))
rownames(t.test.res) = rownames(infilt_dat) ; colnames(t.test.res) = c("p.val",'adj.p')
RRs <- sort(unique(group$RR))
for ( i in rownames(infilt_dat)){
  tmp0 = group$RR ; names(tmp0) = colnames(infilt_dat)
  tmp = data.frame(gene = as.numeric(infilt_dat[i,]),
                   group = group$RR)
  if (length(unique(group$RR))>2){
    res = summary(aov(gene~group,data = tmp))
    t.test.res[i,"p.val"] = res[[1]][["Pr(>F)"]][1]
  }else{
    data_subtype1 <- infilt_dat[i,group$RR==RRs[1]]
    data_subtype2 <- infilt_dat[i,group$RR==RRs[2]]
    res = wilcox.test(data_subtype1,data_subtype2)
    t.test.res[i,"p.val"] = res$p.value
  }
}
t.test.res$adj.p = p.adjust(t.test.res$p.val,method = "fdr")
t.test.res$p <- ifelse(t.test.res$adj.p >0.05,"NS",ifelse(t.test.res$adj.p>0.01,"*",
                                                          ifelse(t.test.res$adj.p>0.001,"**","***")))
write.csv(t.test.res,file = file.path("TIME","aov_res.csv"))

# 行注释数据
methods = paste(gsub("^","_",unique(sapply(strsplit(row.names(infilt_dat),"_"),tail,1))),collapse ="|")
labels_row = paste(t.test.res$p,gsub(methods,"",row.names(infilt_dat)),sep=" ")
labels_row <- gsub("_"," ",labels_row)
annotation_row <- data.frame(row.names = row.names(infilt_dat),
                             methods = (sapply(strsplit(row.names(infilt_dat),"_"),tail,1)),
                             check.names = F)
annotation_row$methods <- factor(annotation_row$methods,levels = unique(annotation_row$methods))

methods_col <- distinctColorPalette(length(unique(annotation_row$methods)))
names(methods_col) <- unique(annotation_row$methods)
#methods_col <- c("TIMER"="#abe2b3","MCPCOUNTER"="#40c3ec","EPIC"="#eae8d7","ESTIMATE"="#fdc68a")
# # 添加临床注释颜色信息
annotation_col <- clinical[,c("RS","RR"),drop=F]
annotation_col$RR <- factor(annotation_col$RR,levels = c("Low-risk","High-risk"))
RS_color <- colorRampPalette(c("white","#84E8EF", "#09A6B2"))(100)
RR_color <- c("Low-risk"="#05a59e","High-risk"="#813e98")

annotation_colors <- list(methods=methods_col,RS=RS_color,RR=RR_color)

col_fun = colorRampPalette(c( "greenyellow","white","red"))(100)
#col_fun = colorRampPalette(c("#fdc68a","#F7F5F5","#09A6B2"))(100)
table(annotation_col$RR)
#设置行注释的gap列
gap_row_tmp <- as.data.frame(table(annotation_row$methods))[,2]
gap_row_tmp <- gap_row_tmp[1:(length(gap_row_tmp)-1)]
gap_row <- c(gap_row_tmp[1])
for (i in gap_row_tmp[-1]){
  gap_row <- c(gap_row,gap_row[length(gap_row)] + i)
}

#设置列注释的gap列
gap_col_tmp <- as.data.frame(table(annotation_col$RR))[,2]
gap_col_tmp <- gap_col_tmp[1:(length(gap_col_tmp)-1)]
gap_col <- c(gap_col_tmp[1])
for (i in gap_col_tmp[-1]){
  gap_col <- c(gap_col,gap_col[length(gap_col)] + i)
}


pdf(file.path("TIME",paste0("TCGA-LUAD","_infiltration.pdf")),height=15,width=15)
pheatmap(infilt_dat,show_colnames = F,show_rownames = T,
         scale = "row",cluster_rows=F,cluster_cols=F,
         clustering_method = 'complete',
         color = col_fun,
         annotation_col = annotation_col,
         annotation_row = annotation_row,
         annotation_colors = annotation_colors,
         fontsize_row=6,fontsize_col=6,
         labels_row=labels_row,gaps_col = gap_col,
         gaps_row = gap_row)
dev.off()


#计算高低风险组在抗肿瘤免疫步骤上的差异
# tumor_immune_cycle_pathway_genesets <- geneSets$`7_Steps_of_the_Tumor_Immune_Cycle`
# tumor_immune_cycle_pathway_ssGSEA <- list()
# for (cohort in names(all_cohorts_raw_data_mRNA_tumor_normal_expression_clean)){
#   #cohort <- names(all_cohorts_raw_data_mRNA_tumor_normal_expression_clean)[1]
#   sur_exp <- all_cohorts_raw_data_mRNA_tumor_normal_expression_clean[[cohort]]
#   ssGSEA.res <- ssGSEA(data = sur_exp,need_log = F,genesets = tumor_immune_cycle_pathway_genesets,scale = T)
#   ssGSEA.res <- data.frame(t(ssGSEA.res))
#   ssGSEA.res[1:4,1:4]
#   normal.ssGSEA.res <- ssGSEA.res[grepl("normal",row.names(ssGSEA.res)),]
#   tumor.ssGSEA.res <- ssGSEA.res[grepl("tumor",row.names(ssGSEA.res)),]
#   tumor_immune_cycle_pathway_ssGSEA[[cohort]] <- list(tumor=tumor.ssGSEA.res,normal=normal.ssGSEA.res)
# }
# save(tumor_immune_cycle_pathway_ssGSEA,file="tumor_immune_cycle_pathway_ssGSEA.Rdata")
load(file="tumor_immune_cycle_pathway_ssGSEA.Rdata")
for (cohort in annotation_cohorts){
  #cohort <- annotation_cohorts[1]
  tumor.ssGSEA.res <- tumor_immune_cycle_pathway_ssGSEA[[cohort]]$tumor
  row.names(tumor.ssGSEA.res) <- gsub("_tumor","",row.names(tumor.ssGSEA.res))
  tumor.ssGSEA.res[1:4,1:4]
  
  riskScore <- read.csv(file.path("101_ml","model.riskGroup.csv"),row.names = 1)
  rs <- riskScore[riskScore$cohort==cohort,"RR",drop=F]
  data <- merge(tumor.ssGSEA.res,rs,by.x=0,by.y=0)
  data <- data[order(data$RR,decreasing = T),]
  data$RR <- factor(data$RR,levels = c("Low-risk","High-risk"))
  
  data_long <- data %>% melt(id=c("RR","Row.names"),value.name = "value",variable.names="type")
  head(data_long)
  p <- multi_item_group_violin_boxplot_plot(data=data_long,
                                            featureVariable="variable",
                                            valueVariable="value",
                                            groupVariable="RR",
                                            colors=risk_col_c,
                                            ylab="Anti-cancer immunity",
                                            type="boxplot")
  p
  ggsave(filename = file.path("TIME","抗肿瘤免疫",paste0(cohort,"_anti_cancer_immunity_boxplot.pdf")),
         height = 5,width = 8 )
}

#计算免疫治疗得分
dir.create(file.path("TIME","CYT"))
for (cohort in annotation_cohorts){
  #cohort <- names(all_cohorts_raw_data_mRNA_os_clean)[1]
  sur_exp <- all_cohorts_raw_data_mRNA_os_clean[[cohort]]
  sur_exp <- sur_exp[grepl("tumor",row.names(sur_exp)),-c(1:2)]
  row.names(sur_exp) <- gsub("_tumor","",row.names(sur_exp))
  sur_exp[1:4,1:4]
  
  res <- immunotherap_signature_score(sur_exp)
  names(res)
  for (i in names(res)){
    plot_list <- list()
    #i=names(res)[1]
    result <- res[[i]]
    
    riskScore <- read.csv(file.path("101_ml","model.riskGroup.csv"),row.names = 1)
    rs <- riskScore[riskScore$cohort==cohort,c("RS","RR"),drop=F]
    data <- merge(result,rs,by.x=0,by.y=0)
    head(rs)
    p <- boxtools(data=data,valueVariable=i,groupVariable="RR",
                  colors=risk_col_l,i,"box_violin",signif_label = "p.format",prefix = cohort)
    p
    plot_list[["diff"]] <- p
    p2 <- vector_vector_association_analysis(data=data,variable_x="RS",variable_y=i,
                                             groupvariable="RR",
                                             group_color = risk_col_l,
                                             point_border_color="#231815",
                                             point_fill_color="#cdcac9",
                                             line_color="#00aeef",
                                             line_margin_color="#b2e7fa",
                                             xlab="RS",
                                             ylab=i,margin=FALSE,prefix = cohort)
    p2$plot
    plot_list[["cor"]] <- p2$plot
    p <- grid.arrange(grobs = plot_list,ncol=2)
    ggsave(filename = file.path("TIME","CYT",paste0(cohort,"_",i,"_diff_cor_between_high_low_risk.pdf"))
           ,p,width=5.5,height = 2,dpi=300)
  }
}


################基因组####NULL################基因组#############################
genome <- file.path("genome")
dir.create(genome,recursive = T)

library(dplyr)
library(ComplexHeatmap)
library(maftools)
library(reshape2)
library(data.table)
library(RColorBrewer)
library(circlize)
library(dplyr)

#CPTAC数据绘制TMB-H,TMB-L生存分析
clinical <- read.csv("CPTAC-LUAD_clinical_riskscore.csv",row.names = 2)
clinical <- clinical[order(clinical$RR,decreasing = T),]
clinical$RR <- factor(clinical$RR,levels = c("Low-risk","High-risk"))
colnames(clinical)
clinical$OS.time <- as.numeric(clinical$OS.time)
clinical$OS <- as.numeric(clinical$OS)
clinical <- clinical %>% filter(!is.na(OS)) %>% filter(!is.na(OS.time))
clinical_tmp <- clinical %>% filter(TMB_State %in% c("TMB-H","TMB-L"))
#绘制TMB-H,TMB-L两组生存分析
p <- single_item_survival_analysis(clinical_tmp,futimeVariable="OS.time",fustatVariable="OS",
                                   groupVariable="TMB_State",need_binary=F,colors=mycolpanel[1:2],ylab="survival probability (OS)",
                                   risk.table=F)

p$plot

clinical <- read.csv("TCGA-LUAD_clinical_riskscore.csv",row.names = 2)
clinical <- clinical[order(clinical$RR,decreasing = T),]
clinical$RR <- factor(clinical$RR,levels = c("Low-risk","High-risk"))
colnames(clinical)

clinical_tmp <- clinical %>% filter(TMB_State %in% c("TMB-H","TMB-L"))
#绘制TMB-H,TMB-L两组生存分析
p <- single_item_survival_analysis(clinical_tmp,futimeVariable="OS.time",fustatVariable="OS",
                              groupVariable="TMB_State",need_binary=F,colors=mycolpanel[1:2],ylab="survival probability (OS)",
                                          risk.table=F)

p$plot

table(clinical$HRD_Label)
clinical_tmp <- clinical %>% filter(HRD_Label %in% c("Negative","Positive"))
#绘制HRD-R,HRD-P两组生存分析
p <- single_item_survival_analysis(clinical_tmp,futimeVariable="OS.time",fustatVariable="OS",
                                   groupVariable="HRD_Label",need_binary=F,colors=mycolpanel[1:2],ylab="survival probability (OS)",
                                   risk.table=F)

p$plot



maf_file <- file.path(root_data_dir,"TCGA-LUAD.mutect.somatic.maf")
mutation_data <- read.data(maf_file,"snv")

mutation <- mutation_data$mutation
mutation_dcast <- mutation_data$mutation_dcast
mutation_dcast_type <- mutation_data$mutation_dcast_type
head(mutation)
mutation_dcast[1:4,1:4]
mutation_dcast_type[1:4,1:4]

same_sample <- intersect(row.names(clinical),mutation$Patient)
clinical <- clinical[same_sample,]
mutation <- mutation[mutation$Patient %in% same_sample,]
mutation_dcast <- mutation_dcast[,same_sample]
mutation_dcast_type <- mutation_dcast_type[,same_sample]

table(clinical$RR)
low_risk_patient <- row.names(clinical)[clinical$RR=="Low-risk"]
high_risk_patient <- row.names(clinical)[clinical$RR=="High-risk"]
table(mutation_dcast["PNLIP",low_risk_patient])
table(mutation_dcast["PNLIP",high_risk_patient])

#对感兴趣的基因进行可视化，这里默认所有的基因
mutation_tmp <- mutation  %>% dplyr::select(Patient,Hugo_Symbol) 
mut2 <- sort(table(mutation_tmp$Hugo_Symbol),decreasing = TRUE)
mut2 <- data.frame(mut2)
row.names(mut2) <- mut2$Var1
mut2$Var1 <- NULL
head(mut2)

genes <- row.names(mut2)
mut2_tmp <- mut2 %>% filter(row.names(mut2) %in% genes )
mutation_dcast_tmp <- mutation_dcast[genes,]
mutation_dcast_tmp[1:4,1:4]

######### 基因进行亚型间突变差异比较######
mutation_dcast_tmp <- apply(mutation_dcast_tmp, 2, function(x){ifelse(x==0, "WT", "Mut")} )
mutation_dcast_tmp <- t(mutation_dcast_tmp)
mutation_dcast_tmp[1:3,1:3]
identical(row.names(mutation_dcast_tmp),row.names(clinical))

stat_res <- data.frame(matrix(NA,nrow = ncol(mutation_dcast_tmp),ncol = length(unique(clinical$RR))))
row.names(stat_res) <- colnames(mutation_dcast_tmp)
colnames(stat_res) <- c(sort(unique(clinical$RR)))
stat_res$p.value <- stat_res[,1]
head(stat_res)
#i="PNLIP"
for (i in colnames(mutation_dcast_tmp)){
  mat <- data.frame(gene = mutation_dcast_tmp[,i],group = sort(clinical$RR),row.names = row.names(mutation_dcast_tmp))
  head(mat)
  table_mutation <- table(mat$gene,mat$group)
  table_mutation
  if (nrow(table_mutation)>1){
    res = fisher.test(table(mat$gene,mat$group))
    #res = chisq.test(table(mat$gene,mat$group))
    tmp = as.matrix.data.frame(table(mat$gene,mat$group))
    head(tmp)
    for (j in 1:(ncol(stat_res)-1)){
      stat_res[i,j] = tmp[1,j]/sum(tmp[,j])
    }
    stat_res[i,"p.value"] = res[["p.value"]]
  }
}


# 将前 ncol - 1 列转换为百分比
for (i in 1:(ncol(stat_res) - 1)) {
  stat_res[, i] <- paste0(round(as.numeric(as.character(stat_res[, i])) * 100, 0), "%")
}

# 对P值进行校正
# stat_res$adj.p <- p.adjust(stat_res$p.value,method = "fdr")

# 对pvalue列进行显著性替换
stat_res$p <- ifelse(stat_res$p.value >0.05,"NS",ifelse(stat_res$p.value>0.01,"*",
                                                        ifelse(stat_res$p.value>0.001,"**","***")))

annotation_gene <- row.names(stat_res)[grepl("\\*",stat_res$p)]
annotation_list <- geneSets$commom_pawhtways_genesets
annotation_col_name <- "commom_pawhtways_genesets"
new_data <- gene_annotation(data=stat_res,annotation_gene = annotation_gene,
                annotation_list = annotation_list,annotation_col_name = annotation_col_name)

annotation_list <- geneSets$`10_Classic_Cancer_Signaling_Pathways_Activated_Repressed`
annotation_col_name <- "Classic_Cancer_Signaling_Pathways"
new_data <- gene_annotation(data=new_data,annotation_gene = annotation_gene,
                            annotation_list = annotation_list,annotation_col_name = annotation_col_name)


annotation_list <- geneSets$msigdb_H
annotation_col_name <- "hallmarker_pawhtways_genesets"
new_data <- gene_annotation(data=new_data,annotation_gene = annotation_gene,
                            annotation_list = annotation_list,annotation_col_name = annotation_col_name)

annotation_list <- geneSets$`470_pancer_driver_genes`
annotation_col_name <- "470_pancer_driver_genes"
new_data <- gene_annotation(data=new_data,annotation_gene = annotation_gene,
                            annotation_list = annotation_list,annotation_col_name = annotation_col_name)

write.csv(new_data,file.path("genome","mutation_compare_across_subtypes_significent_pathway_frequency_static.csv"))

#查看突变有差异的基因是否包含模型中的基因
dat <- read.csv(file.path("101_ml","Coefficients.csv"))
model_gene <- dat$gene
intersect(annotation_gene,model_gene)
#查看模型中的基因突变频率
model_gene_mutation <- new_data[model_gene,]
model_gene_mutation 


#####突变有差异的基因进行生存分析
TCGA_survival_exp <- all_cohorts_raw_data_mRNA_os_clean$`TCGA-LUAD`
TCGA_tumor_survival <- TCGA_survival_exp[grepl("tumor",row.names(TCGA_survival_exp)),c(1,2)]
row.names(TCGA_tumor_survival) <- gsub("_tumor","",row.names(TCGA_tumor_survival))
head(TCGA_tumor_survival)

stat_res <- read.csv(file.path("genome","mutation_compare_across_subtypes_significent_pathway_frequency_static.csv"),row.names = 1)
head(stat_res)
diff_sig_genes <- row.names(stat_res)[grepl("\\*",stat_res$p)]
diff_sig_genes
mutation_dcast_tmp_sig <- mutation_dcast_tmp[,diff_sig_genes]
mutation_dcast_tmp_sig_sur <- merge(TCGA_tumor_survival,mutation_dcast_tmp_sig,by.x=0,by.y=0)
mutation_dcast_tmp_sig_sur <- mutation_dcast_tmp_sig_sur %>% column_to_rownames("Row.names")
head(mutation_dcast_tmp_sig_sur)
res <- multiple.item.survival.analysis(mutation_dcast_tmp_sig_sur,need_binning=FALSE)
res.sig <- res[res$P.Value<0.05,]

stat_res <- res.sig
row.names(stat_res) <- stat_res$item
annotation_gene <- row.names(stat_res)
#运行上面基因注释代码
annotation_list <- geneSets$commom_pawhtways_genesets
annotation_col_name <- "commom_pawhtways_genesets"
new_data <- gene_annotation(data=stat_res,annotation_gene = annotation_gene,
                            annotation_list = annotation_list,annotation_col_name = annotation_col_name)

annotation_list <- geneSets$`10_Classic_Cancer_Signaling_Pathways_Activated_Repressed`
annotation_col_name <- "Classic_Cancer_Signaling_Pathways"
new_data <- gene_annotation(data=new_data,annotation_gene = annotation_gene,
                            annotation_list = annotation_list,annotation_col_name = annotation_col_name)


annotation_list <- geneSets$msigdb_H
annotation_col_name <- "hallmarker_pawhtways_genesets"
new_data <- gene_annotation(data=new_data,annotation_gene = annotation_gene,
                            annotation_list = annotation_list,annotation_col_name = annotation_col_name)

annotation_list <- geneSets$`470_pancer_driver_genes`
annotation_col_name <- "470_pancer_driver_genes"
new_data <- gene_annotation(data=new_data,annotation_gene = annotation_gene,
                            annotation_list = annotation_list,annotation_col_name = annotation_col_name)

write.csv(new_data,file.path("genome","mutation_compare_across_subtypes_significent_pathway_prognosis.csv"))

#查看在两组中的突变频率和生存的关系
stat_res_tmp <- stat_res[res.sig$item,]
stat_res_tmp
for (gene in res.sig$item){
  #gene <- res.sig$item[1]
  sur_res <- single_item_survival_analysis(mutation_dcast_tmp_sig_sur,futimeVariable="OS.time",fustatVariable="OS",
                                            groupVariable=gene,need_binary=F,colors=c("#fdc68a","#05a59e"),ylab="survival probability (OS)",
                                            risk.table=F)
  pdf(file=file.path(genome,paste0(gene,"mutation_wt_survival.pdf")),height = 4,width = 4)
  print(sur_res$plot,newpage = F)
  dev.off()
}

#查看特定癌症中常见的基因突变情况
commom_mutation <- c("TP53","EGFR","HER2","ALK","MET","ROS1","RET","BRAF")
stat_res_tmp <- stat_res[commom_mutation,]
stat_res_tmp

#查看mutation_compare_across_subtypes_significent_pathway_frequency_static.csv 文件，
#根据突变差异和生存预后差异及注释的差异选择感兴趣的基因进行展示
clinical$RR
plot_genes <- c("TP53","MUC16","MUC17","COL11A1","FMN2","NRCAM","CTNND2","SMARCA4","MYCBP2","BSN","ATRX","SPTBN1","NCOR2","KCNA5","MMP2")
mutation_dcast_type_plot <- mutation_dcast_type[plot_genes,row.names(clinical)]
mutation_dcast_type_plot[1:4,1:4]

stat_res_plot <- stat_res[plot_genes,]
stat_res_plot

###画图设置

color_Gender <- c(Female='#AECDE1',Male='#F3C17B',"Unknown"="White")
color_Age <- c("<60"="#7ac7e2",">=60"="#54beaa","Unknown"="White")
color_Stage <- c("Stage I"="#abe2b3","Stage II"="#40c3ec","Stage III"="#eae8d7","Stage IV"="#fdc68a","Unknown"='White')
RR_color <- c("Low-risk"="#05a59e","High-risk"="#813e98")


# 定义颜色
col = c("Missense_Mutation" = "#366A9C", "Nonsense_Mutation" = "#BBDE93",
        "In_Frame_Del" = "#EE8632" ,"In_Frame_Ins" = "#D0342B",
        "Multi_Hit" = "black", "Splice_Site" = "#F3C17B", "Frame_Shift_Del" = "#AECDE1",
        "Frame_Shift_Ins" = "#ED9E9B",
        "Nonstop_Mutation" = "#339900",
        "Unknown" = "#eaeaea")
wid = 3
hei = 2

alter_fun <- list(
  background = function(x, y, w, h) {
    grid.rect(x, y, w-unit(wid, "pt"), h-unit(hei, "pt"),
              gp = gpar(fill = "white", col = "white"))
  },
  Missense_Mutation = function(x, y, w, h) {
    grid.rect(x, y, w-unit(wid, "pt"), h-unit(hei, "pt"),
              gp = gpar(fill = col["Missense_Mutation"], col = col["Missense_Mutation"]))
  },
  Nonsense_Mutation = function(x, y, w, h) {
    grid.rect(x, y, w-unit(wid, "pt"), h-unit(hei, "pt"),
              gp = gpar(fill = col["Nonsense_Mutation"], col = col["Nonsense_Mutation"]))
  },
  In_Frame_Del = function(x, y, w, h) {
    grid.rect(x, y, w-unit(wid, "pt"), h-unit(hei, "pt"),
              gp = gpar(fill = col["In_Frame_Del"], col = col["In_Frame_Del"]))
  },
  In_Frame_Ins = function(x, y, w, h) {
    grid.rect(x, y, w-unit(wid, "pt"), h-unit(hei, "pt"),
              gp = gpar(fill = col["In_Frame_Ins"], col = col["In_Frame_Ins"]))
  },
  Multi_Hit = function(x, y, w, h) {
    grid.rect(x, y, w-unit(wid, "pt"), h-unit(hei, "pt"),
              gp = gpar(fill = col["Multi_Hit"], col = col["Multi_Hit"]))
  },
  Splice_Site = function(x, y, w, h) {
    grid.rect(x, y, w-unit(wid, "pt"), h-unit(hei, "pt"),
              gp = gpar(fill = col["Splice_Site"], col = col["Splice_Site"]))
  },
  Frame_Shift_Del = function(x, y, w, h) {
    grid.rect(x, y, w-unit(wid, "pt"), h-unit(hei, "pt"),
              gp = gpar(fill = col["Frame_Shift_Del"], col = col["Frame_Shift_Del"]))
  },
  Frame_Shift_Ins = function(x, y, w, h) {
    grid.rect(x, y, w-unit(wid, "pt"), h-unit(hei, "pt"),
              gp = gpar(fill = col["Frame_Shift_Ins"], col = col["Frame_Shift_Ins"]))
  },
  Nonstop_Mutation = function(x, y, w, h) {
    grid.rect(x, y, w-unit(wid, "pt"), h-unit(hei, "pt"),
              gp = gpar(fill = col["Nonstop_Mutation"], col = col["Nonstop_Mutation"]))
  },
  Unknown = function(x, y, w, h) {
    grid.rect(x, y, w-unit(wid, "pt"), h-unit(hei, "pt"),
              gp = gpar(fill = col["Unknown"], col = col["Unknown"]))
  }
)

# 定义图例
heatmap_legend_param = list(title = "Alternations", 
                            at = c("Missense_Mutation", "Nonsense_Mutation", 
                                   "Splice_Site","Frame_Shift_Del",
                                   "Frame_Shift_Ins","In_Frame_Del",
                                   "In_Frame_Ins","Nonstop_Mutation",
                                   "Translation_Start_Site"
                            ), 
                            labels = c("Missense_Mutation", "Nonsense_Mutation", 
                                       "Splice_Site","Frame_Shift_Del",
                                       "Frame_Shift_Ins","In_Frame_Del",
                                       "In_Frame_Ins","Nonstop_Mutation",
                                       "Multi_Hit"
                            )
)

# 添加临床信息
columnanno <- HeatmapAnnotation(Gender=clinical$Gender,
                                Age=clinical$Age,
                                Stage=clinical$Stage,
                                Subtypes=clinical$RR,
                                RS=clinical$RS,
                                show_annotation_name = T)
sum(is.na(mutation_dcast[,row.names(clinical)]))

p1 <- oncoPrint(mutation_dcast_type_plot[,row.names(clinical)], alter_fun_is_vectorized = FALSE, 
                alter_fun = alter_fun, col = col, show_pct = T, pct_side = 'right',
                column_title = '',
                row_names_side = "left", show_column_names = F, remove_empty_columns = F,
                heatmap_legend_param = list(title = "Mutation"),border = T,
                row_names_gp = gpar(fontsize = 6), column_title_gp = gpar(fontsize = 0),
                column_split = factor(clinical$RR, levels = sort(unique(clinical$RR))),column_gap = unit(5, "mm"),
                top_annotation = HeatmapAnnotation(Age=clinical$Age,
                                                   Gender=clinical$Gender,
                                                   Stage=clinical$Stage,
                                                   Subtypes=clinical$RR,
                                                   annotation_name_side = "left", annotation_name_gp = gpar(fontsize = 6),
                                                   col = list( Subtypes = RR_color,Gender = color_Gender,
                                                               Age=color_Age,Stage = color_Stage) ) )

p1
####将每个亚型的突变频率添加到右侧注释,不知道为什么用循环不能成功，只显示最后一个。只能分开添加。
p1 <-  p1 + rowAnnotation(x = anno_text(stat_res_plot$Low.risk, just = "center", 
                                        location = unit(0.5, "npc"), show_name = TRUE), annotation_name_rot = 0)
p1 <-  p1 + rowAnnotation(x = anno_text(stat_res_plot$High.risk, just = "center", 
                                        location = unit(0.5, "npc"), show_name = TRUE), annotation_name_rot = 0)
p1 <-  p1 + rowAnnotation(x = anno_text(stat_res_plot$p, just = "center", 
                                        location = unit(0.5, "npc"), show_name = TRUE), annotation_name_rot = 0)
export::graph2pdf(p1,file =  file.path("genome","Mut_plot.pdf"),width = 14, height = 7)

####CNA####
###CNA下载：UCSC Xena,下载copy number (gene-level)，gistic2 thresholded  文件
#gene-level copy number estimates. GISTIC2 further thresholded the estimated values to -2,-1,0,1,2, 
#representing homozygous deletion, single copy deletion, diploid normal copy, 
#low-level copy number amplification, or high-level copy number amplification.
cna_file <- file.path(root_dir,"TCGA-LUAD_Gistic2_CopyNumber_Gistic2_all_thresholded.by_genes")
cna <- data.frame(fread(cna_file,sep = "\t",),row.names = 1)
colnames(cna) <- substr(colnames(cna),1,12)
colnames(cna) <- gsub("\\.","-",colnames(cna))
amp_mtx <- t(cna)
amp_mtx[1:4,1:4]

####临床中包含clust列，聚类结果列，下面确保cna中每个样本和每个样本聚类结果对应
cnv_clinical_samples <- intersect(row.names(amp_mtx),row.names(clinical))
amp_mtx <- amp_mtx[cnv_clinical_samples,]
cnv_clinical <- clinical[cnv_clinical_samples,]

amp_mtx <- as.matrix(amp_mtx)
amp_mtx <- apply(amp_mtx, 2, function(x) {ifelse(x == 2 | x == 1, "Amp", "non_Amp")})
amp_mtx[1:3,1:4]
dim(amp_mtx)
# 统计每个基因扩增的频率
# 统计每列中值等于"amp"的个数
amp_count <- colSums(amp_mtx == "Amp", na.rm = TRUE)
result_df <- data.frame(Gene = colnames(amp_mtx), Amp_Count = amp_count)
result_df <- arrange(result_df, desc(Amp_Count))
row.names(result_df) <- result_df$Gene
result_df$Gene <- NULL

# 选择top200的基因进行亚型间差异分析,chisq.test 检验。
genes <- row.names(result_df)
result_df <- result_df %>% filter(row.names(result_df ) %in% genes)
head(result_df)

amp_mtx_tmp <- amp_mtx[,genes ]
amp_mtx_tmp[1:3,1:4] 

#统计不同亚型之间基因扩增的差异
stat_res <- data.frame(matrix(nrow = ncol(amp_mtx_tmp),ncol = length(unique(cnv_clinical$RR))))
row.names(stat_res) <- colnames(amp_mtx_tmp)
colnames(stat_res) <- c(sort(unique(cnv_clinical$RR)))
stat_res$p.value <- stat_res[,1]
head(stat_res)
for (i in colnames(amp_mtx_tmp)){
  mat <- data.frame(gene = amp_mtx_tmp[,i],group = sort(cnv_clinical$RR),row.names = row.names(amp_mtx_tmp))
  head(mat)
  res = chisq.test(table(mat$gene,mat$group))
  tmp = as.matrix.data.frame(table(mat$gene,mat$group))
  head(tmp)
  for (j in 1:(ncol(stat_res)-1)){
    stat_res[i,j] = tmp[1,j]/sum(tmp[,j])
  }
  stat_res[i,"p.value"] = res[["p.value"]]
}

# 将前 ncol - 1 列转换为百分比
for (i in 1:(ncol(stat_res) - 1)) {
  stat_res[, i] <- paste0(round(as.numeric(as.character(stat_res[, i])) * 100, 0), "%")
}

# 对P值进行校正
# stat_res$adj.p <- p.adjust(stat_res$p.value,method = "fdr")

# 对pvalue列进行显著性替换
stat_res$p <- ifelse(stat_res$p.value >0.05,"NS",ifelse(stat_res$p.value>0.01,"*",
                                                      ifelse(stat_res$p.value>0.001,"**","***")))
stat_res <- merge(stat_res,result_df,by=0,all=F)
stat_res  <- stat_res %>% arrange(Amp_Count) 
row.names(stat_res) <- stat_res$Row.names
stat_res$Row.names = NULL
head(stat_res)  


annotation_gene <- row.names(stat_res)[grepl("\\*",stat_res$p)]
annotation_list <- geneSets$commom_pawhtways_genesets
annotation_col_name <- "commom_pawhtways_genesets"
new_data <- gene_annotation(data=stat_res,annotation_gene = annotation_gene,
                            annotation_list = annotation_list,annotation_col_name = annotation_col_name)

annotation_list <- geneSets$`10_Classic_Cancer_Signaling_Pathways_Activated_Repressed`
annotation_col_name <- "Classic_Cancer_Signaling_Pathways"
new_data <- gene_annotation(data=new_data,annotation_gene = annotation_gene,
                            annotation_list = annotation_list,annotation_col_name = annotation_col_name)


annotation_list <- geneSets$msigdb_H
annotation_col_name <- "hallmarker_pawhtways_genesets"
new_data <- gene_annotation(data=new_data,annotation_gene = annotation_gene,
                            annotation_list = annotation_list,annotation_col_name = annotation_col_name)
write.csv(new_data,file.path("genome","cna_compare_across_subtypes_using_chisq_test.csv"))

####选择可视化的基因，这里，可以自己根据top_200_cna_compare_across_subtypes_using_chisq_test.csv 文件中的基因自己定义，
#目前代码中选择adj.p <=0.05的基因，然后选择高频的基因进行组合，共显示5个。

plot_genes <- row.names(stat_res[which(stat_res$adj.p<=0.05),])
need_supplement_gene_num <- 5 - length(plot_genes)
if (need_supplement_gene_num > 0){
  diff_gene <- setdiff(row.names(stat_res)[1:20],plot_genes)
  diff_gene <- diff_gene[1:need_supplement_gene_num]
  plot_genes <- c(plot_genes,diff_gene)
} 
stat_res_plot <- stat_res[plot_genes,]
stat_res_plot$gene <- row.names(stat_res_plot)
stat_res_plot <- stat_res_plot %>% arrange(desc(Amp_Count))
row.names(stat_res_plot) <- stat_res_plot$gene
stat_res_plot$gene <- NULL
stat_res_plot <- stat_res_plot[1:5,]
head(stat_res_plot)

amp_mtx2 = t(amp_mtx[,row.names(stat_res_plot)])
head(amp_mtx2)
col2 = c( "Amp" = "#BC102B","non_Amp" = "white")
wid = 3
hei = 2

alter_fun2 <- list(
  background = function(x, y, w, h) {
    grid.rect(x, y, w-unit(wid, "pt"), h-unit(hei, "pt"),
              gp = gpar(fill = "white", col = "white"))
  },
  Amp = function(x, y, w, h) {
    grid.rect(x, y, w-unit(wid, "pt"), h-unit(hei, "pt"),
              gp = gpar(fill = col2["Amp"], col = col2["Amp"]))
  },
  non_Amp = function(x, y, w, h) {
    grid.rect(x, y, w-unit(wid, "pt"), h-unit(hei, "pt"),
              gp = gpar(fill = col2["non_Amp"], col = col2["non_Amp"]))
  }
)

# 提取 clust 列的唯一值并排序
clust_levels <- sort(unique(clinical$RR))

# 创建因子变量，指定因子水平的顺序
p2 <- oncoPrint(
  amp_mtx2[, row.names(clinical)],
  alter_fun_is_vectorized = FALSE,
  alter_fun = alter_fun2,
  col = col2,
  show_pct = TRUE,
  show_column_names = FALSE,
  na_col = "#eaeaea",
  border = TRUE,
  row_names_gp = gpar(fontsize = 9),
  heatmap_legend_param = list(title = "CNA"),
  row_names_side = "left",
  pct_side = 'right',
  top_annotation = NULL,
  column_split = factor(clinical$RR, levels = clust_levels),
  column_gap = unit(5, "mm")
)

####将每个亚型的扩增频率添加到右侧注释,不知道为什么用循环不能成功，只显示最后一个。只能分开添加。
p2 <-  p2 + rowAnnotation(x = anno_text(stat_res_plot$CS1, just = "center", 
                                        location = unit(0.5, "npc"), show_name = TRUE), annotation_name_rot = 0)
p2 <-  p2 + rowAnnotation(x = anno_text(stat_res_plot$CS2, just = "center", 
                                        location = unit(0.5, "npc"), show_name = TRUE), annotation_name_rot = 0)

#hp_list = p1 %v% p2
export::graph2pdf(p2,file =  file.path("genome","CNA_plot.pdf"),width = 14, height = 7)


#突变特征比较
TCGA_LUAD_mutation_sig <- read.csv(file.path(root_dir,"TCGA-LUAD_mutation_signature.csv"),row.names = 1)
head(TCGA_LUAD_mutation_sig)
risk <- clinical[,c("RS","RR")]
head(risk)
risk_mut_signature <- merge(risk,TCGA_LUAD_mutation_sig,by.x=0,by.y=0)
head(risk_mut_signature)
data <- risk_mut_signature
data$RR <- factor(data$RR,levels=c("Low-risk","High-risk"))
signatures <-  colnames(risk_mut_signature)[grepl("Signature",colnames(risk_mut_signature))]
signatures <-  c("Signature.4","Signature.5","Signature.13","Signature.20")
for (i in signatures){
  #i=signatures[2]
  plot_list <- list()
  p <- boxtools(data=data,valueVariable=i,groupVariable="RR",
                colors=risk_col_l,i,"box_violin",signif_label = "p.format",prefix = "")
  p
  plot_list[["diff"]] <- p
  p2 <- vector_vector_association_analysis(data=data,variable_x="RS",variable_y=i,
                                           groupvariable="RR",
                                           group_color = risk_col_l,
                                           point_border_color="#231815",
                                           point_fill_color="#cdcac9",
                                           line_color="#00aeef",
                                           line_margin_color="#b2e7fa",
                                           xlab="RS",
                                           ylab=i,margin=FALSE,prefix = "")
  plot_list[["cor"]] <- p2$plot
  p <- grid.arrange(grobs = plot_list,ncol=2)
  ggsave(filename = file.path(genome,paste0("TCGA-LUAD_Mutation_signate_diff",i,".pdf")),p,width=5.5,height = 2,dpi=300)
}


#TCGA-LUAD免疫浸润单因素cox回归分析
outdir <- "scRNA"
dir.create(outdir)
names(all_cohorts_raw_data_mRNA_tumor_os_clean)
TCGA_LUAD_sursive_exp <- all_cohorts_raw_data_mRNA_os_clean$`TCGA-LUAD`
TCGA_LUAD_sursive_exp <- TCGA_LUAD_sursive_exp[grepl("tumor",row.names(TCGA_LUAD_sursive_exp)),]
TCGA_LUAD_sursive <- TCGA_LUAD_sursive_exp[,c(1,2)]
row.names(TCGA_LUAD_sursive) <- gsub("_tumor","",row.names(TCGA_LUAD_sursive))
head(TCGA_LUAD_sursive)

infiltration_estimation_file <- "E:\\流程搭建\\多组学聚类分型\\文档\\TCGA所有样本新生成特征文件\\TCGA_infiltration_estimation_for_tcga(01样本).csv"
infilt_dat <- read.csv(infiltration_estimation_file,sep=",",row.names = 1)
colnames(infilt_dat) <- gsub("\\.\\.","\\.",colnames(infilt_dat))
infilt_dat[1:4,1:4]
infilt_dat_tmp <- infilt_dat[,grepl("CIBERSORT.ABS",colnames(infilt_dat))]
infilt_dat[1:4,1:4]
colnames(infilt_dat_tmp)

data <- merge(TCGA_LUAD_sursive,infilt_dat_tmp,by.x=0,by.y = 0)
library(ezcox)
#修改因子水平，将contrast放到前面，ref_level放到后面
#TCGA_LIHC_cox_imput_tmp$risk <- factor(TCGA_LIHC_cox_imput_tmp$risk,levels = c("Low-risk","High-risk"))
zz <- ezcox(data,covariates = colnames(data)[-c(1,2)],
            time = "OS.time",status = "OS",return_models=TRUE)
mds <- get_models(zz)
pdf(file.path(outdir,paste0("infiltration","_unicox.pdf")),height = 10,width = 10)
show_models(mds)
dev.off()


#风险平分在单细胞层面验证
library(qs)
library(Seurat)

dat <- read.csv("Coefficients.csv",row.names = 1)
head(dat)

seurat.data <- qread("seurat.data.006.qs")
expression <- seurat.data@assays$RNA$data
expression_tmp <- expression[row.names(dat),]
expression_tmp[,1:4]

risk_score_sc_mart <- expression_tmp * dat$coefficient
risk_score_sc <- data.frame(colSums(risk_score_sc_mart))
colnames(risk_score_sc) <- "risk_score"

metadata <- seurat.data@meta.data
metadata_risk_score <- cbind(metadata,risk_score_sc[row.names(metadata),,drop=F])

umap <- seurat.data@reductions$umap@cell.embeddings
head(umap)
merge_risk_meta <- merge(metadata_risk_score,umap,by.x = 0,by.y = 0)
head(merge_risk_meta)

col_fun = colorRampPalette(c( "#1eaed3","#73d4ea","#ef8a7f","#ea3b28"))(1000)
theme_with_legend <- theme(panel.grid.major = element_blank(),
                           panel.grid.minor = element_blank(),
                           axis.text.y = element_text(color = "black",size = 6),
                           axis.text.x = element_text(color = "black",size = 6,angle = 45,hjust = 1),
                           axis.title = element_text(color = "black",size = 6),
                           legend.position = "right",
                           plot.title = element_text(hjust = 0.5,size = 8),
                           legend.title = element_blank(),
                           legend.background = element_rect(fill = 'transparent'),
                           legend.text = element_text(size = 6),
                           legend.key.size = unit(0.5, "cm"))  + theme_bw(base_rect_size = 0.5)

p <- ggplot(data = merge_risk_meta,aes(UMAP_1,UMAP_2,colour = risk_score)) +
      geom_point(size=0.1) +
      scale_color_gradientn(colours = col_fun) +
  theme_bw() + theme_with_legend

p
ggsave(file=file.path(dir_standard_process,'risk_score.pdf'),p,height = 6,width = 7)

p <- ggplot(data = merge_risk_meta,aes(UMAP_1,UMAP_2,colour = celltype)) +
  geom_point(size=0.1) +
  scale_fill_npg() + 
  theme_bw() + theme_with_legend

p
ggsave(file=file.path(dir_standard_process,'umap_cell_annotation.pdf'),p,height = 6,width = 8)


p <- ggplot(data = merge_risk_meta,aes(celltype,risk_score,fill = celltype)) +
  geom_boxplot(outlier.size = 1) + 
  scale_fill_npg() + 
  theme_bw() + 
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        axis.text.y = element_text(color = "black",size = 6),
        axis.text.x = element_text(color = "black",size = 6,angle = 45,hjust = 1),
        axis.title = element_text(color = "black",size = 6),
        legend.position = "right",
        plot.title = element_text(hjust = 0.5,size = 8),
        legend.title = element_blank(),
        legend.background = element_rect(fill = 'transparent'),
        legend.text = element_text(size = 6),
        legend.key.size = unit(0.5, "cm"))

p
ggsave(file=file.path(dir_standard_process,'risk_score_boxplot.pdf'),p,height = 5,width = 8)

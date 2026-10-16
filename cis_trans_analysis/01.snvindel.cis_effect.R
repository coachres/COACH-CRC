library(dplyr)
library(caret)
library(tidyverse)
library(tidyr)
library(data.table)
library(stats)

## lm
linear.regression = function(Y,X,gene_list,formula,Omics_name,outliers_list=list()){
  out_a = 1
  data = NULL
  for(i in gene_list){
    try({
      print(i)
      tmp_data = Y[,c("Sample",i,"ESTIMATE TumorPurity")]
      colnames(tmp_data) = c("Sample","test_target_Y","ESTIMATE_Purity")
      tmp_data$test_target_Y = as.numeric(tmp_data$test_target_Y)
      test_target_X = subset(X,select=c('Sample',i))
      colnames(test_target_X) = c('Sample',"test_target_X")
      tmp_data = merge(test_target_X, tmp_data,by = 'Sample', all = TRUE)

      ## remove outliers
      outliers_flag = "NO"
      if (i %in% outliers_list) {
        print("removing outliers...")
        tmp_data_WT <- tmp_data[tmp_data[,"test_target_X"]==0,]
        tmp_data_Mut <- tmp_data[tmp_data[,"test_target_X"]==1,]
        tmp_data_WT = tmp_data_WT[which(!tmp_data_WT[,"test_target_Y"] %in% boxplot.stats(tmp_data_WT[,"test_target_Y"])$out),]
        tmp_data_Mut = tmp_data_Mut[which(!tmp_data_Mut[,"test_target_Y"] %in% boxplot.stats(tmp_data_Mut[,"test_target_Y"])$out),]
        tmp_data = rbind(tmp_data_WT,tmp_data_Mut)
        outliers_flag = "YES"
      }
      
      if(length(unique(tmp_data$test_target_X))!=1){
          results = lm(formula=formula,data=tmp_data,na.action="na.omit")
          cooksd <- cooks.distance(results)
          results = summary(results)
          B <- results$coefficients[2,1]
          q <- results$coefficients[2,4]
          q2 <- results$coefficients[3,4]
          R2 <- results$adj.r.squared ## updates to adj R square
          AdjR = sqrt(abs(R2))   ## use the absolute values to calculate the square root
          B_AdjR = B * AdjR
          signed_AdjR = (B/abs(B)) * AdjR
          fvalue <- results$fstatistic
          pvalue <- 1-pf(fvalue[1],fvalue[2],fvalue[3])
          # value = c(i,Omics_name,R2,pvalue,B,q,AdjR,B_AdjR,signed_AdjR,signed_log10_adj_p)
          value = c(i,Omics_name,R2,B,AdjR,B_AdjR,signed_AdjR,pvalue,q,q2,outliers_flag)
          # print(cooksd)
          ##print
          if(out_a==1){ 
            data <- value
          }else{
            data <- rbind(data,value)
          }
          out_a = out_a+1
      }
    },FALSE)
  }
  if(is.null(dim(data))){
    return(data)
  }else{
    colnames(data) = c("Gene","Omics_name","AdjR_squared","B","AdjR","B_times_AdjR","signed_AdjR","pvalue","q","purity_q","outliers_filter")
    data = as.data.frame(data)
    data$AdjR = as.numeric(data$AdjR)
    data$q = as.numeric(data$q)
    data$B = as.numeric(data$B)
    data$Adj_pvalue = p.adjust(data$pvalue,method="fdr")
    data$Adj_q = p.adjust(data$q,method="fdr")
    data$purity_q = as.numeric(data$purity_q)
    data$Adj_purity_q = p.adjust(data$purity_q,method="fdr")

    data <- mutate(data,
      signed_log10_adj_p = (B/abs(B)) * (-log10(Adj_pvalue+(10**-200))),
      signed_log10_adj_q = (B/abs(B)) * (-log10(Adj_q+(10**-200))),
      signed_log10_adj_purity_q = (B/abs(B)) * (-log10(Adj_purity_q+(10**-200))),
      signed_log10_q = (B/abs(B)) * (-log10(q+(10**-200)))
    )
    return(data)
  }
}

linear.regression_phos = function(Y,X,gene_list,formula,Omics_name,outliers_list=list()){
  out_a = 1
  data = NULL
  for(i in gene_list){
    try({
      split_string <- strsplit(i, "_")  
      gene_prefix <- split_string[[1]][1]
      print(i)
      tmp_data = Y[,c("Sample",i,"ESTIMATE TumorPurity")]
      colnames(tmp_data) = c("Sample","test_target_Y","ESTIMATE_Purity")
      tmp_data$test_target_Y = as.numeric(tmp_data$test_target_Y)
      set.seed(123)
      noise <- rnorm(length(tmp_data$test_target_Y), mean = 0, sd = 0.5)
      tmp_data$noise = noise
      test_target_X = subset(X,select=c('Sample',gene_prefix))
      colnames(test_target_X) = c('Sample',"test_target_X")
      tmp_data = merge(test_target_X, tmp_data,by = 'Sample', all.y = TRUE)

      outliers_flag = "NO"
      if (i %in% outliers_list) {
        print("removing outliers...")
        tmp_data_WT <- tmp_data[tmp_data[,"test_target_X"]==0,]
        tmp_data_Mut <- tmp_data[tmp_data[,"test_target_X"]==1,]
        tmp_data_WT = tmp_data_WT[which(!tmp_data_WT[,"test_target_Y"] %in% boxplot.stats(tmp_data_WT[,"test_target_Y"])$out),]
        tmp_data_Mut = tmp_data_Mut[which(!tmp_data_Mut[,"test_target_Y"] %in% boxplot.stats(tmp_data_Mut[,"test_target_Y"])$out),]
        tmp_data = rbind(tmp_data_WT,tmp_data_Mut)
        outliers_flag = "YES"
      }

      if(length(unique(tmp_data$test_target_X))!=1){
          results = lm(formula=formula,data=tmp_data,na.action="na.omit")
          results = summary(results)
          B <- results$coefficients[2,1]
          q <- results$coefficients[2,4]
          q2 <- results$coefficients[3,4]
          R2 <- results$adj.r.squared ## updates to adj R square
          AdjR = sqrt(abs(R2))   ## use the absolute values to calculate the square root
          B_AdjR = B * AdjR
          signed_AdjR = (B/abs(B)) * AdjR
          # signed_log10_adj_p = (B/abs(B)) * (-log10(q))
          fvalue <- results$fstatistic
          pvalue <- 1-pf(fvalue[1],fvalue[2],fvalue[3])
          value = c(i,Omics_name,R2,B,AdjR,B_AdjR,signed_AdjR,pvalue,q,q2,outliers_flag)
          ##print
          if(out_a==1){ 
            data <- value
          }else{
            data <- rbind(data,value)
          }
          out_a = out_a+1
      }
    },FALSE)
  }
  if(is.null(dim(data))){
    return(data)
  }else{
    colnames(data) = c("Gene","Omics_name","AdjR_squared","B","AdjR","B_times_AdjR","signed_AdjR","pvalue","q","purity_q","outliers_filter")
    data = as.data.frame(data)
    data$AdjR = as.numeric(data$AdjR)
    data$q = as.numeric(data$q)
    data$B = as.numeric(data$B)
    data$Adj_pvalue = p.adjust(data$pvalue,method="fdr")
    data$Adj_q = p.adjust(data$q,method="fdr")
    data$purity_q = as.numeric(data$purity_q)
    data$Adj_purity_q = p.adjust(data$purity_q,method="fdr")
    data <- mutate(data,
      signed_log10_adj_p = (B/abs(B)) * (-log10(Adj_pvalue+(10**-200))),
      signed_log10_adj_q = (B/abs(B)) * (-log10(Adj_q+(10**-200))),
      signed_log10_adj_purity_q = (B/abs(B)) * (-log10(Adj_purity_q+(10**-200))),
      signed_log10_q = (B/abs(B)) * (-log10(q+(10**-200)))
    )
    return(data)
  }
}

## clinic info
crc1005_clinical_genetic_profile = read.table("clinical_genetic_profile.txt",sep = "\t",header = T,check.names = F)
is_hyper <- crc1005_clinical_genetic_profile$Hypermutated == "hyper"  
df_filtered <- crc1005_clinical_genetic_profile[!is_hyper, ]  
nonhyper_sample <- df_filtered$Sample

# RNA vs SNV_indel
# RNA expression
RNA_expression = read.table("FPKM.matrix.trans.log2_plus1.tsv",sep = "\t",header = T,check.names = F)
RNA_expression = RNA_expression[!is.na(RNA_expression$Sample),]
RNA_expression <- mutate_at(RNA_expression, colnames(RNA_expression), ~replace(., is.na(.), 0))

## SNV indel mutation
non_silent_snvindel = read.table("gene_count_nonsilent.trans.tsv",sep = "\t",header = T,check.names = F)

## select samples
non_silent_snvindel = non_silent_snvindel[non_silent_snvindel$Sample %in% RNA_expression$Sample,]
non_silent_snvindel = non_silent_snvindel[non_silent_snvindel$Sample %in% nonhyper_sample,]
RNA_expression = RNA_expression[RNA_expression$Sample %in% non_silent_snvindel$Sample,]
RNA_expression = RNA_expression[RNA_expression$Sample %in% nonhyper_sample,]

# head(RNA_expression)
# head(non_silent_snvindel)
## gene list
gene_lines <- readLines("snvindel_39_driver.gene.list")
gene_list <- strsplit(gene_lines, "\n")

## outliers list
outliers_lines <- readLines("outliers.RNA.gene.list")
outliers_list <- strsplit(outliers_lines, "\n")

formula_estimate_purity = test_target_Y ~ test_target_X + ESTIMATE_Purity
data = linear.regression(RNA_expression,non_silent_snvindel,gene_list,formula_estimate_purity,"SNV_RNA_estimate_purity",outliers_list)
write.table(data, file ="SNV_RNA_lm_estimate_purity.result.tsv", quote=FALSE, sep = "\t",row.names = F, col.names = TRUE)

#####
#### Protein vs SNV_indel 
protein = read.table("protein.trans.log2.tsv", sep = "\t", header = T, check.names = F)
protein = protein[!is.na(protein$Sample),]

## SNV indel mutation
non_silent_snvindel = read.table("gene_count.trans.tsv",sep = "\t",header = T,check.names = F)

## select samples
non_silent_snvindel = non_silent_snvindel[non_silent_snvindel$Sample %in% protein$Sample,]
non_silent_snvindel = non_silent_snvindel[non_silent_snvindel$Sample %in% nonhyper_sample,]
protein = protein[protein$Sample %in% non_silent_snvindel$Sample,]
protein = protein[protein$Sample %in% nonhyper_sample,]

## gene list
gene_lines <- readLines("snvindel_driver.gene.list")
gene_list <- strsplit(gene_lines, "\n")

## outliers list
outliers_lines <- readLines("outliers.protein.gene.list")
outliers_list <- strsplit(outliers_lines, "\n")

formula_estimate_purity = test_target_Y ~ test_target_X + ESTIMATE_Purity
data = linear.regression(protein,non_silent_snvindel,gene_list,formula_estimate_purity,"SNV_protein_protein_estimate_purity",outliers_list)
write.table(data, file ="SNV_protein_lm_estimate_purity.result.tsv", quote=FALSE, sep = "\t",row.names = F, col.names = TRUE)

########
## Phospho sites vs SNV_indel
print("phosphorylation analysis")
phosphate = read.table("phosSite.trans.log2.tsv",sep= "\t", header = T, check.names = F)

# ## SNV indel mutation
non_silent_snvindel = read.table("gene_count.trans.tsv",sep = "\t",header = T,check.names = F)

# ## select samples
non_silent_snvindel = non_silent_snvindel[non_silent_snvindel$Sample %in% phosphate$Sample,]
non_silent_snvindel = non_silent_snvindel[non_silent_snvindel$Sample %in% nonhyper_sample,]
phosphate = phosphate[phosphate$Sample %in% non_silent_snvindel$Sample,]
phosphate = phosphate[phosphate$Sample %in% nonhyper_sample,]

## gene list
gene_lines <- readLines("snvindel_39_driver.gene_phosphate.list")
gene_list <- strsplit(gene_lines, "\n")

## outliers list
outliers_lines <- readLines("outliers.phosphorylation.gene.list")
outliers_list <- strsplit(outliers_lines, "\n")

formula_estimate_purity = test_target_Y ~ test_target_X + ESTIMATE_Purity
data = linear.regression_phos(phosphate,non_silent_snvindel,gene_list,formula_estimate_purity,"SNV_phosphate_protein_estimate_purity",outliers_list)
write.table(data, file ="SNV_phosphate_protein_lm_estimate_purity.result.tsv", quote=FALSE, sep = "\t",row.names = F, col.names = TRUE)

########
## methyl
print("methylation analysis")
methyl = read.table("methyaltion_level.tsv",sep= "\t", header = T, check.names = F)
# ## SNV indel mutation
non_silent_snvindel = read.table("gene_count.trans.tsv",sep = "\t",header = T,check.names = F)

## select samples
non_silent_snvindel = non_silent_snvindel[non_silent_snvindel$Sample %in% methyl$Sample,]
non_silent_snvindel = non_silent_snvindel[non_silent_snvindel$Sample %in% nonhyper_sample,]
methyl = methyl[methyl$Sample %in% non_silent_snvindel$Sample,]
methyl = methyl[methyl$Sample %in% nonhyper_sample,]

## gene list
gene_lines <- readLines("snvindel_driver.gene.list")
gene_list <- strsplit(gene_lines, "\n")

## outliers list
outliers_lines <- readLines("outliers.methyl.gene.list")
outliers_list <- strsplit(outliers_lines, "\n")

## regression
formula_estimate_purity = test_target_Y ~ test_target_X + ESTIMATE_Purity
data = linear.regression(methyl,non_silent_snvindel,gene_list,formula_estimate_purity,"SNV_methyl_estimate_purity",outliers_lines)
write.table(data, file ="SNV_methyl_lm_estimate_purity.result.tsv", quote=FALSE, sep = "\t",row.names = F, col.names = TRUE)


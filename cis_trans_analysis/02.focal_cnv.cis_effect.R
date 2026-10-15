library(dplyr)

## lm
linear.regression = function(Y,X,gene_list,formula,Omics_name,cnv,outliers_list){
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
      tmp_data = merge(test_target_X, tmp_data,by = 'Sample', all.y = TRUE)

      cnv_gene = cnv[,c("Sample",i)]
      colnames(cnv_gene) = c('Sample',"mutation")
      tmp_data = merge(cnv_gene,tmp_data,by = 'Sample', all.y = TRUE)

      outliers_flag = "NO"
      sample_count <- nrow(tmp_data)
      sample_rmoutliers_count <- "NA"
      if (i %in% outliers_list) {
        print("removing outliers...")
        tmp_data_WT <- tmp_data[tmp_data[,"mutation"]==0,]
        tmp_data_Mut <- tmp_data[tmp_data[,"mutation"]==1,]
        tmp_data_WT = tmp_data_WT[which(!tmp_data_WT[,"test_target_Y"] %in% boxplot.stats(tmp_data_WT[,"test_target_Y"])$out),]
        tmp_data_WT = tmp_data_WT[which(!tmp_data_WT[,"test_target_X"] %in% boxplot.stats(tmp_data_WT[,"test_target_X"])$out),]
        tmp_data_Mut = tmp_data_Mut[which(!tmp_data_Mut[,"test_target_Y"] %in% boxplot.stats(tmp_data_Mut[,"test_target_Y"])$out),]
        tmp_data_Mut = tmp_data_Mut[which(!tmp_data_Mut[,"test_target_X"] %in% boxplot.stats(tmp_data_Mut[,"test_target_X"])$out),]
        tmp_data = rbind(tmp_data_WT,tmp_data_Mut)
        outliers_flag = "YES"
        sample_rmoutliers_count <- nrow(tmp_data)
      }
        number_of_zero <- nrow(tmp_data[tmp_data$test_target_Y == 0,])

      if(length(unique(tmp_data$test_target_X))!=1){
          results = lm(formula=formula,data=tmp_data,na.action="na.omit")
          results = summary(results)
          B <- results$coefficients[2,1]
          q <- results$coefficients[2,4]
          q2 <- results$coefficients[3,4]
          R2 <- results$r.squared ## updates to adj R square
          R = sqrt(abs(R2)) 
          signed_R = (B/abs(B)) * R 
          AdjR2 <- results$adj.r.squared ## updates to adj R square
          AdjR = sqrt(abs(AdjR2))   ## use the absolute values to calculate the square root
          B_AdjR = B * AdjR
          signed_AdjR = (B/abs(B)) * AdjR
          fvalue <- results$fstatistic
          pvalue <- 1-pf(fvalue[1],fvalue[2],fvalue[3])
          value = c(i,Omics_name,R2,AdjR2,B,R,AdjR,signed_R,signed_AdjR,q,outliers_flag,sample_count,sample_rmoutliers_count,number_of_zero)
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
    colnames(data) = c("Gene","Omics_name","R_squared","AdjR_squared","B","R","AdjR","signed_R","signed_AdjR","q","outliers_filter","sample_num","sample_num_rmoutliers","Target_Y_0_number")
    data = as.data.frame(data)
    data$AdjR = as.numeric(data$AdjR)
    data$q = as.numeric(data$q)
    data$B = as.numeric(data$B)
    # print(data$B)
    data$Adj_q = p.adjust(data$q,method="fdr")
    data <- mutate(data,
      signed_log10_adj_q = (B/abs(B)) * (-log10(Adj_q+(10**-500))),
    )
    return(data)
  }
}

linear.regression_phos = function(Y,X,gene_list,formula,Omics_name,cnv,outliers_list){
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
      test_target_X = subset(X,select=c('Sample',gene_prefix))
      colnames(test_target_X) = c('Sample',"test_target_X")
      tmp_data = merge(test_target_X, tmp_data,by = 'Sample', all.y = TRUE)

      cnv_gene = cnv[,c("Sample",gene_prefix)]
      colnames(cnv_gene) = c('Sample',"mutation")
      tmp_data = merge(cnv_gene,tmp_data,by = 'Sample', all.y = TRUE)

      outliers_flag = "NO"
      sample_count <- nrow(tmp_data)
      sample_rmoutliers_count <- "NA"

      if (i %in% outliers_list) {
        print("removing outliers...")
        tmp_data_WT <- tmp_data[tmp_data[,"mutation"]==0,]
        tmp_data_Mut <- tmp_data[tmp_data[,"mutation"]==1,]
        tmp_data_WT = tmp_data_WT[which(!tmp_data_WT[,"test_target_Y"] %in% boxplot.stats(tmp_data_WT[,"test_target_Y"])$out),]
        tmp_data_WT = tmp_data_WT[which(!tmp_data_WT[,"test_target_X"] %in% boxplot.stats(tmp_data_WT[,"test_target_X"])$out),]
        tmp_data_Mut = tmp_data_Mut[which(!tmp_data_Mut[,"test_target_Y"] %in% boxplot.stats(tmp_data_Mut[,"test_target_Y"])$out),]
        tmp_data_Mut = tmp_data_Mut[which(!tmp_data_Mut[,"test_target_X"] %in% boxplot.stats(tmp_data_Mut[,"test_target_X"])$out),]
        tmp_data = rbind(tmp_data_WT,tmp_data_Mut)
        outliers_flag = "YES"
        sample_rmoutliers_count <- nrow(tmp_data)
        number_of_zero <- nrow(tmp_data[tmp_data$test_target_Y==0,])
      }


      if(length(unique(tmp_data$test_target_X))!=1){
          results = lm(formula=formula,data=tmp_data,na.action="na.omit")
          results = summary(results)
          B <- results$coefficients[2,1]
          q <- results$coefficients[2,4]
          q2 <- results$coefficients[3,4]
          R2 <- results$r.squared
          R = sqrt(abs(R2))
          signed_R = (B/abs(B)) * R

          AdjR2 <- results$adj.r.squared
          AdjR = sqrt(abs(AdjR2))
          B_AdjR = B * AdjR
          signed_AdjR = (B/abs(B)) * AdjR
          fvalue <- results$fstatistic
          pvalue <- 1-pf(fvalue[1],fvalue[2],fvalue[3])
          #value = c(i,Omics_name,R2,AdjR2,B,R,AdjR,B_AdjR,signed_R,signed_AdjR,pvalue,q,q2,outliers_flag,sample_count,sample_rmoutliers_count)
          value = c(i,Omics_name,R2,AdjR2,B,R,AdjR,signed_R,signed_AdjR,q,outliers_flag,sample_count,sample_rmoutliers_count,number_of_zero)
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
    colnames(data) = c("Gene","Omics_name","R_squared","AdjR_squared","B","R","AdjR","signed_R","signed_AdjR","q","outliers_filter","sample_num","sample_num_rmoutliers","Target_Y_0_number")
    data = as.data.frame(data)
    data$AdjR = as.numeric(data$AdjR)
    data$q = as.numeric(data$q)
    data$B = as.numeric(data$B)
    data$Adj_q = p.adjust(data$q,method="fdr")
    data <- mutate(data,
      signed_log10_adj_q = (B/abs(B)) * (-log10(Adj_q+(10**-500))),
    )
    return(data)
  }
}



## cnv mutation states
cnv = read.table("focal_CNV.gene.matrix.tsv",sep = "\t",header = T,check.names = F)

## RNA vs SNV_indel
# RNA expression
print("RNA ciseffect")
RNA_expression = read.table("PKM.matrix.trans.log2_plus1.tsv",sep = "\t",header = T,check.names = F)
RNA_expression = RNA_expression[!is.na(RNA_expression$Sample),]
RNA_expression <- mutate_at(RNA_expression, colnames(RNA_expression), ~replace(., is.na(.), 0))

## focal cnv mutation
focal_cnv = read.table("cnv.selected.trans.rsBygene",sep = "\t",header = T,check.names = F)

# ## select samples
focal_cnv =focal_cnv[focal_cnv$Sample %in% RNA_expression$Sample,]
RNA_expression = RNA_expression[RNA_expression$Sample %in% focal_cnv$Sample,]

# ## gene list
gene_lines <- readLines("input_gene.RNA.list")
gene_list <- strsplit(gene_lines, "\n")


## outliers list
outliers_lines <- readLines("outliers.RNA.gene.list")
outliers_list <- strsplit(outliers_lines, "\n")

formula_estimate_purity = test_target_Y ~ test_target_X + ESTIMATE_Purity
data = linear.regression(RNA_expression,focal_cnv,gene_list,formula_estimate_purity,"purity_logratio_cnv_RNA_estimate_purity",cnv,outliers_list)
write.table(data, file ="estimate_purity.result.tsv", quote=FALSE, sep = "\t",row.names = F, col.names = TRUE)


# ########
# #### Protein vs SNV_indel
print("protein ciseffect") 
protein = read.table("protein.trans.log2.tsv", sep = "\t", header = T, check.names = F)
protein = protein[!is.na(protein$Sample),]

## focal cnv mutation
focal_cnv = read.table("cnv.selected.trans.rsBygene",sep = "\t",header = T,check.names = F)

# ## select samples
focal_cnv = focal_cnv[focal_cnv$Sample %in% protein$Sample,]
protein = protein[protein$Sample %in% focal_cnv$Sample,]

# ## gene list
gene_lines <- readLines("input_gene.protein.list")
gene_list <- strsplit(gene_lines, "\n")

## outliers list
outliers_lines <- readLines("outliers.protein.gene.list")
outliers_list <- strsplit(outliers_lines, "\n")

formula_estimate_purity = test_target_Y ~ test_target_X + ESTIMATE_Purity
data = linear.regression(protein,focal_cnv,gene_list,formula_estimate_purity,"purity_logratio_cnv_protein_estimate_purity",cnv,outliers_list)
write.table(data, file ="estimate_purity.result.tsv", quote=FALSE, sep = "\t",row.names = F, col.names = TRUE)

########
# ## Phospho sites vs SNV_indel
print("phosphorylation ciseffect")
phosphate = read.table("phosSite.20230907.trans.log2.tsv",sep= "\t", header = T, check.names = F)
# # ## SNV indel mutation
focal_cnv = read.table("cnv.selected.trans.rsBygene",sep = "\t",header = T,check.names = F)

# # ## select samples
focal_cnv = focal_cnv[focal_cnv$Sample %in% phosphate$Sample,]
phosphate = phosphate[phosphate$Sample %in% focal_cnv$Sample,]

# ## gene list
gene_lines <- readLines("input_gene.phosSite.list")
gene_list <- strsplit(gene_lines, "\n")

## outliers list
outliers_lines <- readLines("outliers.phosphorylation.gene.list")
outliers_list <- strsplit(outliers_lines, "\n")

formula_estimate_purity = test_target_Y ~ test_target_X + ESTIMATE_Purity
data = linear.regression_phos(phosphate,focal_cnv,gene_list,formula_estimate_purity,"purity_logratio_cnv_phosphate_protein_estimate_purity",cnv,outliers_list)
write.table(data, file ="02.purity_logratio_cnv_phosphate_protein_lm_estimate_purity.result.tsv", quote=FALSE, sep = "\t",row.names = F, col.names = TRUE)

q()
########
## methyl
print("methyl ciseffect")
methyl = read.table("methyaltion_level.tsv",sep= "\t", header = T, check.names = F)
# ## mutation
focal_cnv = read.table("cnv.selected.trans.rsBygene",sep = "\t",header = T,check.names = F)

## select samples
focal_cnv = focal_cnv[focal_cnv$Sample %in% methyl$Sample,]
methyl = methyl[methyl$Sample %in% focal_cnv$Sample,]

## gene list
gene_lines <- readLines("methyaltion_level_cnv.fomatted.gene.list")
gene_list <- strsplit(gene_lines, "\n")

## outliers list
outliers_lines <- readLines("outliers.methyl.gene.list")
outliers_list <- strsplit(outliers_lines, "\n")

## regression
formula_estimate_purity = test_target_Y ~ test_target_X + ESTIMATE_Purity
data = linear.regression(methyl,focal_cnv,gene_list,formula_estimate_purity,"purity_logratio_cnv_methyl_estimate_purity",cnv,outliers_list)
write.table(data, file ="purity_logratio_cnv_methyl_lm_estimate_purity.result.tsv", quote=FALSE, sep = "\t",row.names = F, col.names = TRUE)


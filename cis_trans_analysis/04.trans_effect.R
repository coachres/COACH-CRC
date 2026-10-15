library(dplyr)
library(tidyr)
library(data.table)
library(foreach)  
library(doParallel)


## lm for trans effects
linear.regression = function(Y,X,gene_list,test_target_gene_list,formula,Omics_name){
  out_a = 1
  data = NULL
  for(gene in gene_list){
    print(gene)
    for(i in test_target_gene_list){
      if(i!=gene){
        try({
          print(i)
          ## test target Y
          test_target_Y_gene = i
          tmp_data = Y[,c("Sample",test_target_Y_gene,"ESTIMATE TumorPurity")]
          print("1")
          colnames(tmp_data) = c("Sample","test_target_Y","ESTIMATE_Purity")
          tmp_data$test_target_Y = as.numeric(tmp_data$test_target_Y)

          ## test target X
          print("3")
          test_target_X_gene = gene
          test_target_X = subset(X,select=c('Sample',test_target_X_gene))
          print("2")
          colnames(test_target_X) = c('Sample',"test_target_X")

          ## set noise
          set.seed(123)
          noise <- rnorm(length(tmp_data$test_target_Y), mean = 0, sd = 0.5)
          tmp_data$noise = noise
          
          ## merge test data table
          tmp_data = merge(test_target_X, tmp_data, by = 'Sample', all.y = TRUE)

          ## lm calculation
          if(length(unique(tmp_data$test_target_X))!=1){
              results = lm(formula=formula,data=tmp_data,na.action="na.omit")
              results = summary(results)
              B <- results$coefficients[2,1]
              q <- results$coefficients[2,4]
              R2 <- results$adj.r.squared ## updates to adj R square
              AdjR = sqrt(abs(R2))
              B_AdjR = B * AdjR
              signed_AdjR = (B/abs(B)) * AdjR
              fvalue <- results$fstatistic
              pvalue <- 1-pf(fvalue[1],fvalue[2],fvalue[3])
              value = c(test_target_X_gene,test_target_Y_gene,Omics_name,R2,pvalue,B,q,AdjR,B_AdjR,signed_AdjR)
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
    }
  }
  if(is.null(dim(data))){
    return(data)
  }else{
    colnames(data) = c("Gene","trans_target_gene","Omics_name","AdjR_squared","pvalue","B","q","AdjR","B_times_AdjR","signed_AdjR")
    data = as.data.frame(data)
    data$AdjR = as.numeric(data$AdjR)
    data$q = as.numeric(data$q)
    data$B = as.numeric(data$B)
    data$Adj_pvalue = p.adjust(data$pvalue,method="fdr")
    data$Adj_q = p.adjust(data$q,method="fdr")
    data <- mutate(data,
      signed_log10_adj_p = (B/abs(B)) * (-log10(Adj_pvalue+(10**-200))),
      signed_log10_adj_q = (B/abs(B)) * (-log10(Adj_q+(10**-200))),
      signed_log10_q = (B/abs(B)) * (-log10(q+(10**-200)))
    )
    return(data)
  }
}

linear.regression_par = function(Y,X,gene,test_target_gene_list,formula,Omics_name){
  out_a = 1
  data = NULL
  for(i in test_target_gene_list){
    if(i!=gene){
      try({
        ## test target Y
        test_target_Y_gene = i
        tmp_data = Y[,c("Sample",test_target_Y_gene,"ESTIMATE TumorPurity")]
        colnames(tmp_data) = c("Sample","test_target_Y","ESTIMATE_Purity")
        tmp_data$test_target_Y = as.numeric(tmp_data$test_target_Y)

        ## test target X
        test_target_X_gene = gene
        test_target_X = subset(X,select=c('Sample',test_target_X_gene))
        colnames(test_target_X) = c('Sample',"test_target_X")

        ## set noise
        set.seed(123)
        noise <- rnorm(length(tmp_data$test_target_Y), mean = 0, sd = 0.5)
        tmp_data$noise = noise
        
        ## merge test data table
        tmp_data = merge(test_target_X, tmp_data, by = 'Sample', all.y = TRUE)

        ## lm calculation
        if(length(unique(tmp_data$test_target_X))!=1){
            results = lm(formula=formula,data=tmp_data,na.action="na.omit")
            results = summary(results)
            B <- results$coefficients[2,1]
            q <- results$coefficients[2,4]
            R2 <- results$r.squared
            AdjR = sqrt(R2)
            B_AdjR = B * AdjR
            signed_AdjR = (B/abs(B)) * AdjR
            signed_log10_adj_p = (B/abs(B)) * (-log10(q))
            fvalue <- results$fstatistic
            pvalue <- 1-pf(fvalue[1],fvalue[2],fvalue[3])
            value = c(test_target_X_gene,test_target_Y_gene,Omics_name,R2,pvalue,B,q,AdjR,B_AdjR,signed_AdjR,signed_log10_adj_p)
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
  }
  if(is.null(dim(data))){
    return(data)
  }else{
    colnames(data) = c("Gene","trans_target_gene","Omics_name","AdjR_squared","pvalue","B","q","AdjR","B_times_AdjR","signed_AdjR","signed_log10_adj_p")
    data = as.data.frame(data)
    data$AdjR = as.numeric(data$AdjR)
    data$q = as.numeric(data$q)
    return(data)
  }
}

## for trans effect phosphorylation
linear.regression_phos = function(Y,X,gene_list,phosphate_sites_list,formula,Omics_name){
  out_a = 1
  data = NULL
  for(gene in gene_list){
    print(gene)
    for(i in phosphate_sites_list){
      split_string <- strsplit(i, "_")  
      gene_prefix <- split_string[[1]][1]
      ## trans effect 
      print(gene_prefix)
      if(gene_prefix!=gene){
        try({
          ## test target Y
          test_target_Y_phosphate = i
          tmp_data = Y[,c("Sample",test_target_Y_phosphate,"ESTIMATE TumorPurity","Purity")]
          colnames(tmp_data) = c("Sample","test_target_Y","ESTIMATE_Purity","Purity")
          tmp_data$test_target_Y = as.numeric(tmp_data$test_target_Y)

          ## test target X
          test_target_X_gene = gene
          test_target_X = subset(X,select=c('Sample',test_target_X_gene))
          colnames(test_target_X) = c('Sample',"test_target_X")

          ## set noise
          set.seed(123)
          noise <- rnorm(length(tmp_data$test_target_Y), mean = 0, sd = 0.5)
          tmp_data$noise = noise

          ## merge test data table
          tmp_data = merge(test_target_X, tmp_data,by = 'Sample', all.y = TRUE)
          #lm calculating
          if(length(unique(tmp_data$test_target_X))!=1){
              results = lm(formula=formula,data=tmp_data,na.action="na.omit")
              results = summary(results)
              B <- results$coefficients[2,1]
              q <- results$coefficients[2,4]
              # R2 <- results$r.squared
              # AdjR = sqrt(R2)
              R2 <- results$adj.r.squared ## updates to adj R square
              AdjR = sqrt(abs(R2)) 
              B_AdjR = B * AdjR
              signed_AdjR = (B/abs(B)) * AdjR
              # signed_log10_adj_p = (B/abs(B)) * (-log10(q))
              fvalue <- results$fstatistic
              pvalue <- 1-pf(fvalue[1],fvalue[2],fvalue[3])
              value = c(test_target_X_gene,test_target_Y_phosphate,Omics_name,R2,pvalue,B,q,AdjR,B_AdjR,signed_AdjR)
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
    }
  }
  if(is.null(dim(data))){
    return(data)
  }else{
    colnames(data) = c("Gene","trans_target_phosphate_sites","Omics_name","AdjR_squared","pvalue","B","q","AdjR","B_times_AdjR","signed_AdjR")
    data = as.data.frame(data)
    data$AdjR = as.numeric(data$AdjR)
    data$q = as.numeric(data$q)
    data$B = as.numeric(data$B)
    data$Adj_pvalue = p.adjust(data$pvalue,method="fdr")
    data$Adj_q = p.adjust(data$q,method="fdr")
    data <- mutate(data,
      signed_log10_adj_p = (B/abs(B)) * (-log10(Adj_pvalue+(10**-200))),
      signed_log10_adj_q = (B/abs(B)) * (-log10(Adj_q+(10**-200))),
      signed_log10_q = (B/abs(B)) * (-log10(q+(10**-200)))
    )
    return(data)
  }
}

linear.regression_phos_par = function(Y,X,gene,phosphate_sites_list,formula,Omics_name){
  out_a = 1
  data = NULL
  for(i in phosphate_sites_list){
    split_string <- strsplit(i, "_")  
    gene_prefix <- split_string[[1]][1]
    ## trans effect 
    if(gene_prefix!=gene){
      try({
        ## test target Y
        test_target_Y_phosphate = i
        tmp_data = Y[,c("Sample",test_target_Y_phosphate,"ESTIMATE TumorPurity","Purity")]
        colnames(tmp_data) = c("Sample","test_target_Y","ESTIMATE_Purity","Purity")
        tmp_data$test_target_Y = as.numeric(tmp_data$test_target_Y)

        ## test target X
        test_target_X_gene = gene
        test_target_X = subset(X,select=c('Sample',test_target_X_gene))
        colnames(test_target_X) = c('Sample',"test_target_X")

        ## set noise
        set.seed(123)
        noise <- rnorm(length(tmp_data$test_target_Y), mean = 0, sd = 0.5)
        tmp_data$noise = noise

        ## merge test data table
        tmp_data = merge(test_target_X, tmp_data,by = 'Sample', all.y = TRUE)
        
        #lm calculating
        if(length(unique(tmp_data$test_target_X))!=1){
            results = lm(formula=formula,data=tmp_data,na.action="na.omit")
            results = summary(results)
            B <- results$coefficients[2,1]
            q <- results$coefficients[2,4]
            R2 <- results$r.squared
            AdjR = sqrt(R2)
            B_AdjR = B * AdjR
            signed_AdjR = (B/abs(B)) * AdjR
            signed_log10_adj_p = (B/abs(B)) * (-log10(q))
            fvalue <- results$fstatistic
            pvalue <- 1-pf(fvalue[1],fvalue[2],fvalue[3])
            value = c(test_target_X_gene,test_target_Y_phosphate,Omics_name,R2,pvalue,B,q,AdjR,B_AdjR,signed_AdjR,signed_log10_adj_p)
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
  }
  if(is.null(dim(data))){
    return(data)
  }else{
    colnames(data) = c("Gene","trans_target_phosphate_sites","Omics_name","AdjR_squared","pvalue","B","q","AdjR","B_times_AdjR","signed_AdjR","signed_log10_adj_p")
    data = as.data.frame(data)
    data$AdjR = as.numeric(data$AdjR)
    data$q = as.numeric(data$q)
    return(data)
  }
}

### parallel
args <- commandArgs(trailingOnly = TRUE) 

thread_index <- which(args == "--threads")  

if (!is.na(thread_index)) {  
  threads <- as.numeric(args[thread_index + 1])  
} else {  
  threads <- 1 
}

input_dir_index <- which(args == "--input-dir")  
input_dir <- args[input_dir_index + 1]

mut_input_index <- which(args == "--mut-input")  
mut_input <- args[mut_input_index + 1]
mut_input <- paste(input_dir, mut_input, sep="/")

gene_list_index <- which(args == "--gene")  
input_gene <- args[gene_list_index + 1]

output_prefix <- paste("SNVindel",input_gene,sep="_")

if ("nonhyper" %in% args) {
  clinic <- paste(input_dir,"clinical_genetic_profile.txt",sep="/")
  crc1005_clinical_genetic_profile = read.table(clinic,sep = "\t",header = T,check.names = F)
  is_hyper <- crc1005_clinical_genetic_profile$Hypermutated == "hyper"  
  df_filtered <- crc1005_clinical_genetic_profile[!is_hyper, ]  
  nonhyper_sample <- df_filtered$Sample
}

# RNA vs SNV_indel
## RNA expression
if ("RNA" %in% args) {  
  print("calculating correlation versus RNA data")
  RNA_input <- paste(input_dir, "FPKM.matrix.trans.log2_plus1.tsv", sep="/")
  RNA_expression = read.table(RNA_input ,sep = "\t",header = T,check.names = F)
  RNA_expression = RNA_expression[!is.na(RNA_expression$Sample),]
  RNA_expression <- mutate_at(RNA_expression, colnames(RNA_expression), ~replace(., is.na(.), 0))

  ## SNV indel mutation
  non_silent_snvindel = read.table(mut_input,sep = "\t",header = T,check.names = F)

  ## select samples
  non_silent_snvindel = non_silent_snvindel[non_silent_snvindel$Sample %in% RNA_expression$Sample,]
  RNA_expression = RNA_expression[RNA_expression$Sample %in% non_silent_snvindel$Sample,]
  if ("nonhyper" %in% args) {
    print("select nonhyper samples")
    non_silent_snvindel = non_silent_snvindel[non_silent_snvindel$Sample %in% nonhyper_sample,]
    RNA_expression = RNA_expression[RNA_expression$Sample %in% nonhyper_sample,]
  }

  gene_list <- c(input_gene)

  ## test_target_gene_list
  test_target_gene <- paste(input_dir, "FPKM.matrix.trans.log2_plus1.gene_list.tsv", sep="/")
  gene_lines <- readLines(test_target_gene)
  test_target_gene_list <- strsplit(gene_lines, "\n")

  formula_estimate_purity = test_target_Y ~ test_target_X + ESTIMATE_Purity

  ## linear.regression
  omics <- paste(output_prefix, "_RNA_estimate_purity", sep="")
  output_file <- paste(output_prefix, "_RNA_lm_estimate_purity.result.tsv", sep="")
  data = linear.regression(RNA_expression,non_silent_snvindel,gene_list,test_target_gene_list,formula_estimate_purity,omics)

  write.table(data, file = output_file, quote=FALSE, sep = "\t",row.names = F, col.names = TRUE)
}

#######
if ("protein" %in% args) {
  print("calculating correlation versus protein data")
  #### Protein vs SNV_indel 
  protein_input <- paste(input_dir, "protein.trans.log2.tsv", sep="/")
  protein = read.table(protein_input, sep = "\t", header = T, check.names = F)
  protein = protein[!is.na(protein$Sample),]

  ## SNV indel mutation
  non_silent_snvindel = read.table(mut_input, sep = "\t",header = T,check.names = F)

  ## select samples
  non_silent_snvindel = non_silent_snvindel[non_silent_snvindel$Sample %in% protein$Sample,]
  protein = protein[protein$Sample %in% non_silent_snvindel$Sample,]
  if ("nonhyper" %in% args) {
    print("select nonhyper samples")
    non_silent_snvindel = non_silent_snvindel[non_silent_snvindel$Sample %in% nonhyper_sample,]
    protein = protein[protein$Sample %in% nonhyper_sample,]
  }
  gene_list <- c(input_gene)

  ## test target gene list
  test_target_gene <- paste(input_dir, "protein.trans.log2.gene.list", sep="/")
  gene_lines <- readLines(test_target_gene)
  test_target_gene_list <- strsplit(gene_lines, "\n")

  formula_estimate_purity = test_target_Y ~ test_target_X + ESTIMATE_Purity

  omics <- paste(output_prefix, "_protein_estimate_purity", sep="")
  output_file <- paste(output_prefix, "_protein_trans_lm_estimate_purity.result.tsv", sep="")
  data <- linear.regression(protein,non_silent_snvindel,gene_list,test_target_gene_list,formula_estimate_purity,omics)

  write.table(data, file =output_file, quote=FALSE, sep = "\t",row.names = F, col.names = TRUE)

}

# ########
if ("phosphate" %in% args) {
  print("calculating correlation versus phosphorylation data")
  ## Phospho sites vs SNV_indel
  phos_protein_input <- paste(input_dir, "phosSite.trans.log2.tsv", sep="/")
  phosphate = read.table(phos_protein_input,sep= "\t", header = T, check.names = F)
  # ## SNV indel mutation
  non_silent_snvindel = read.table(mut_input,sep = "\t",header = T,check.names = F)

  # ## select samples
  non_silent_snvindel = non_silent_snvindel[non_silent_snvindel$Sample %in% phosphate$Sample,]
  phosphate = phosphate[phosphate$Sample %in% non_silent_snvindel$Sample,]
  if ("nonhyper" %in% args) {
    print("select nonhyper samples")
    non_silent_snvindel = non_silent_snvindel[non_silent_snvindel$Sample %in% nonhyper_sample,]
    phosphate = phosphate[phosphate$Sample %in% nonhyper_sample,]
  }
  ## gene list
 
  gene_list <- c(input_gene)

  ## phosphate sites list
  test_target_gene <- paste(input_dir, "phosSite.trans.log2.gene_phosphate.list", sep="/")
  gene_lines <- readLines(test_target_gene) 
  phosphate_sites_list <- strsplit(gene_lines, "\n")

  formula_estimate_purity = test_target_Y ~ test_target_X + ESTIMATE_Purity

  omics <- paste(output_prefix, "_phosphate_protein_estimate_purity", sep="")
  output_file <- paste(output_prefix, "_phosphate_protein_lm_estimate_purity.result.tsv", sep="")
  data <- linear.regression_phos(phosphate,non_silent_snvindel,gene_list,phosphate_sites_list,formula_estimate_purity, omics)

  # head(data)
  write.table(data, file = output_file, quote=FALSE, sep = "\t",row.names = F, col.names = TRUE)
}

# # ########
if ("methyl" %in% args) {
  print("calculating correlation versus methyl data")
  ## methyl
  methyl_input <- paste(input_dir, "methyaltion_level.tsv", sep="/")
  methyl = read.table(methyl_input,sep= "\t", header = T, check.names = F)
  # ## SNV indel mutation
  non_silent_snvindel = read.table(mut_input,sep = "\t",header = T,check.names = F)

  ## select samples
  non_silent_snvindel = non_silent_snvindel[non_silent_snvindel$Sample %in% methyl$Sample,]
  methyl = methyl[methyl$Sample %in% non_silent_snvindel$Sample,]
  if ("nonhyper" %in% args) {
    print("select nonhyper samples")
    non_silent_snvindel = non_silent_snvindel[non_silent_snvindel$Sample %in% nonhyper_sample,]
    methyl = methyl[methyl$Sample %in% nonhyper_sample,]
  }
  ## gene list
  gene_list <- c(input_gene)

  ## test target gene list
  test_target_gene <- paste(input_dir, "methyaltion_level.list", sep="/")
  gene_lines <- readLines(test_target_gene) 
  test_target_gene_list <- strsplit(gene_lines, "\n")

  ## regression
  formula_estimate_purity = test_target_Y ~ test_target_X + ESTIMATE_Purity
  # formula_facets_purity = test_target_Y ~ test_target_X + Purity + noise
  registerDoParallel(makeCluster(threads))
  omics <- paste(output_prefix, "_methyl_estimate_purity", sep="")
  output_file <- paste(output_prefix, "_methyl_lm_estimate_purity.result.tsv", sep="")

  data = linear.regression(methyl,non_silent_snvindel,gene_list,test_target_gene_list,formula_estimate_purity,omics)
  # head(data)
  write.table(data, file = output_file, quote=FALSE, sep = "\t",row.names = F, col.names = TRUE)

}

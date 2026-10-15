library(ComplexHeatmap)
library(ggplot2)
library(tidyverse)
library(colorRamp2)

## color set
col_fun <- circlize::colorRamp2(
  breaks = c(-2, 0, 2), 
  colors = c("blue", "white", "red")
)

col_fun2 <- circlize::colorRamp2(
  breaks = c(-1.5, 0, 1.5), 
  colors = c("blue", "white", "red")
)

## input data
mut <- read.table("mutation.txt",sep="\t",header=T, check.names = F, stringsAsFactors = F)
mut <- mut[,c("Sample","Driver:KRAS","Driver:TP53","Driver:ACVR2A","Driver:RNF43")]
colnames(mut) = c("Sample","KRAS","TP53","ACVR2A","RNF43")

sign <- read.csv("data_for_plot.csv", stringsAsFactors = F, check.names = F)
sign <- sign %>% filter(Draw == "Y")

iobr <- read.table("immune.txt", sep="\t", header=T, stringsAsFactors = F, check.names = F, row.names = 1)
estimate_cols <- grep("ESTIMATE", rownames(iobr))
iobr <- iobr[(rowSums(as.matrix(iobr) > 0.00001) > (ncol(iobr) * 0.7)) | 
               (rownames(iobr) %in% rownames(iobr)[estimate_cols]),]
dim(iobr)

dat_cluster <- as.matrix(t(iobr) )
dat_cluster <- dat_cluster[-grep("ESTIMATE", colnames(dat_cluster)), ]
dat_cluster <- scale(dat_cluster)
set.seed(123)
kmeans_result <- kmeans(dat_cluster, centers = 4)
print(kmeans_result$centers)

df = data.frame(Cluster = as.factor(kmeans_result$cluster))
subtype <- read.table("cmos_information.txt",header = T,stringsAsFactors = F,sep = '\t')
df$Sample = rownames(df)
subtype <- subtype %>%
  left_join(df, by = "Sample") %>%
  left_join(mut, by = "Sample")

subtype <- subtype[order(subtype$CMOS), ]

iobr <- iobr[rownames(iobr) %in% sign$Subtype, subtype$Sample]
dim(iobr)

tme <- read.table("tme.out.txt", sep="\t", header=T, stringsAsFactors = F, check.names = F, row.names = 1)
tme <- tme[rownames(tme) %in% sign$Subtype, subtype$Sample]
dim(tme)

gsva_rna <- read.table("pathway_gsva.tsv.gz", head=T, stringsAsFactors=F, check.names=F, row.names = 1)
gsva_rna <- gsva_rna[rownames(gsva_rna) %in% sign$Subtype, subtype$Sample]
dim(gsva_rna)

paty <- read.table("progeny_out.txt", sep="\t", header=T, stringsAsFactors = F, check.names = F, row.names = 1)
paty <- paty[rownames(paty) %in% sign$Subtype, subtype$Sample]
dim(paty)

immu_est <- iobr[grep("ESTIMATE", rownames(iobr)), ]

subtype$logTMB = log10(subtype$TMB)
subtype$EMB = subtype$Epimutation_burden
subtype$CMS = subtype$Final_CMS

lables = c("CMOS","CMS","WGD_status","logTMB","SVcount","EMB","PGA","TP53","KRAS","ACVR2A","RNF43")

cols = list(CMS = c(CMS1="#E7A307",CMS2="#0072B1",CMS3="#CD78A7",CMS4="#009F74",Mixed="#999999"),
           CMOS = c(CMOS1="#4372A1",CMOS2="#58A553",CMOS3="#DE9C49",CMOS4="#B43A2F"),
           WGD_status = c(WGD="#67a9cf", noWGD="#f7f7f7"),
           logTMB = colorRamp2(c(0, 3), c("white", "#8e0152")),
           PGA = colorRamp2(c(0, 1), c("white", "#08306b")),
           SVcount = colorRamp2(c(0, 300), c("white", "#8c510a")),
           EMB = colorRamp2(c(0, 0.4), c("white", "#004529")),
           ACVR2A = c("0"="#f0f0f0","1"="black"),
           RNF43 = c("0"="#f0f0f0","1"="black"),
           TP53 = c("0"="#f0f0f0","1"="black"),
           KRAS = c("0"="#f0f0f0","1"="black"))


gsv <- ComplexHeatmap::Heatmap(
  t(scale(t(gsva_rna))),
  name = "z-score",
  col = col_fun,
  cluster_rows = TRUE,
  clustering_distance_rows = "manhattan",
  clustering_method_rows = "ward.D2",
  cluster_columns = FALSE,
  show_column_names = FALSE,
  row_title = "Hallmark",
  column_split = subtype$CMOS,
  row_title_rot = 90, 
  row_title_gp = grid::gpar(fontsize = 10, fontface = "bold"),
  row_names_gp = grid::gpar(fontsize = 7))

sig <- ComplexHeatmap::Heatmap(
  t(scale(t(tme))),
  name = "z-score",
  col = col_fun,
  cluster_rows = TRUE,
  clustering_distance_rows = "manhattan",
  clustering_method_rows = "ward.D2",
  cluster_columns = FALSE,
  show_column_names = FALSE,
  column_split = subtype$CMOS,
  row_title = "TME signature",
  row_title_rot = 90,
  row_title_gp = grid::gpar(fontsize = 10, fontface = "bold"),
  row_names_gp = grid::gpar(fontsize = 7)) 

pro <- ComplexHeatmap::Heatmap(
  t(scale(t(paty))),
  name = "z-score",
  col = col_fun,
  cluster_rows = TRUE,
  clustering_distance_rows = "manhattan",
  clustering_method_rows = "ward.D2",
  cluster_columns = FALSE,
  show_column_names = FALSE,
  column_split = subtype$CMOS,
  row_title = "PROGENy", 
  row_title_rot = 90,
  row_title_gp = grid::gpar(fontsize = 10, fontface = "bold"),
  row_names_gp = grid::gpar(fontsize = 7)) 

immu_rm_est <- iobr[-grep("ESTIMATE", rownames(iobr)), ]
immu_est = iobr[grep("ESTIMATE", rownames(iobr)), ]

immue <- ComplexHeatmap::Heatmap(
  t(scale(t(immu_rm_est))),
  name = "z-score",
  col = col_fun,
  cluster_rows = TRUE,
  clustering_distance_rows = "manhattan",
  clustering_method_rows = "ward.D2",
  cluster_columns = FALSE,
  column_split = subtype$CMOS,
  show_column_names = FALSE,
  top_annotation = HeatmapAnnotation(df = subtype[,lables],
                                     col = cols),
  row_title = "TME",
  row_title_rot = 90, 
  row_title_gp = grid::gpar(fontsize = 10, fontface = "bold"),
  row_names_gp = grid::gpar(fontsize = 7)) 

draw(immue)

est <- ComplexHeatmap::Heatmap(
  t(scale(t(immu_est))),
  name = "z-score",
  col = col_fun,
  cluster_rows = TRUE,
  clustering_distance_rows = "manhattan",
  clustering_method_rows = "ward.D2",
  cluster_columns = FALSE,
  show_column_names = FALSE,
  column_split = subtype$CMOS,
  row_title = "PROGENy",
  row_title_rot = 90,
  row_title_gp = grid::gpar(fontsize = 10, fontface = "bold"),
  row_names_gp = grid::gpar(fontsize = 7)) 

draw(est)

pdf("cmos_feature_heatmap.pdf", h=9, w=9)
draw(immue %v% sig %v% gsv %v% pro %v% est)
dev.off()


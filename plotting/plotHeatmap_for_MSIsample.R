library(ComplexHeatmap)
library(ggplot2)
library(tidyverse)
library(gridExtra)
library(ggpubr)
library(colorRamp2)
library(cowplot)
library(readxl)

dat <- read_excel("MSIsample_information.xlsx", sheet = 1)

msi_df <- dat %>% filter(`MSIstatus_update` == "MSI")

## input mutations
mat <- read.table("mutation.txt",header = T,stringsAsFactors = F,sep = '\t',check.names = F,row.names = 1)
mat[is.na(mat)] <- ""

mat_tran <- as.data.frame(t(mat))
mat_tran[mat_tran == ""] <- "WT"
mat_tran$Sample <- rownames(mat_tran)
mat_tran <- mat_tran %>% rename_with(~ paste0(., "_MUT"), -Sample)

msi <- merge(msi_df, mat_tran, by = "Sample")


## sort
msi <- msi[order(msi$`Germline dMMR gene`, 
                 msi$`MLH1_methylation`,
                 msi$MLH1,
                 msi$MSH2,
                 msi$MLH1_MUT,
                 msi$PMS2_MUT,
                 msi$MSH2_MUT,
                 msi$MSH6_MUT),]

kmat <- mat[, msi$Sample]

col = c("INFRAME" = "#9970ab", 
        "Germline" = "#018571",
        "Missense" = "#018571", 
        "TRUNC" = "#d95f0e", 
        "SV" = "#ae017e",
        "LOH" = "#a6bddb", 
        "HOMDEL" = "#2166ac", 
        "bg" = "#CCCCCC")

col_scale_exp <- colorRamp2(c(-2, 0, 2), c("blue", "white" , "red"))

alter_fun = list(
  background = alter_graphic("rect", fill = "#f0f0f0", col = NA, lwd = 0.01),   
  INFRAME = alter_graphic("rect", fill = col["INFRAME"]),
  TRUNC = alter_graphic("rect", fill = col["TRUNC"]),
  SPLICE = alter_graphic("rect", fill = col["SPLICE"]),
  HOMDEL = alter_graphic("rect", fill = col["HOMDEL"]),
  SV = alter_graphic("rect", fill = "NA", lwd = 2, col = col["SV"]),
  LOH = alter_graphic("rect", height = 0.5, fill = col["LOH"]),
  Germline = alter_graphic("point", pch = "*", col = "#b2182b", size = unit(2, "cm")),
  Missense = alter_graphic("rect", fill = col["Missense"]))

## main heatmap
msi <- msi %>%
  mutate(
    `Germline dMMR gene` = if_else(is.na(`Germline dMMR gene`), "others", `Germline dMMR gene`),
    `MLH1_methylation` = if_else(is.na(`MLH1_methylation`), "others", `MLH1_methylation`)
  )

kmut_main <- oncoPrint(kmat[2:5,], get_type = function(x) strsplit(x, ";")[[1]],
          right_annotation = NULL,
          top_annotation = HeatmapAnnotation(Lynch = msi$`Germline dMMR gene`,
                                             MLH1_met = msi$`MLH1_methylation`,
                                             BRAF = msi$BRAF,
                                             KRAS = msi$KRAS,
                                             NRAS = msi$NRAS,
                                             annotation_name_gp = gpar(fontsize = 10, fontface = "bold"),
                                             na_col = "#d9d9d9",
                                             col = list(MLH1_met = c("HYPER"="#ef8a62","others"="#5ab4ac"),
                                                        Lynch = c("charger > 4"="#d8b365","others"="#a1d76a"),
                                                        KRAS = c("MUT"="black","WT"="gray95"),
                                                        NRAS = c("MUT"="black","WT"="gray95"),
                                                        BRAF = c("V600E"="red"))),
          show_row_names = T,
          row_order = c("MLH1","PMS2","MSH2","MSH6"),
          column_order = msi$Sample,
          row_names_gp = gpar(fontsize = 10),
          alter_fun = alter_fun, 
          border = TRUE,
          gap = unit(0.01, "mm"),
          col = col)

draw(kmut_main)

# Protein
msi_pro <- data.frame(MLH1_Protein = msi$MLH1_Protein,
                      PMS2_Protein = "NA",
                      MSH2_Protein = msi$MSH2_Protein,
                      MSH6_Protein = msi$MSH6_Protein)

colnames(msi_pro) = c("MLH1_Pro","PMS2_Pro","MSH2_Pro","MSH6_Pro")
rownames(msi_pro) = msi$Sample
msi_pro$MLH1_Pro <- as.numeric(msi_pro$MLH1_Pro)
msi_pro$MSH2_Pro <- as.numeric(msi_pro$MSH2_Pro)
msi_pro$MSH6_Pro <- as.numeric(msi_pro$MSH6_Pro)
msi_pro$PMS2_Pro <- as.numeric(msi_pro$PMS2_Pro)
msi_pro_log <- log2(msi_pro)
msi_pro_log <- t(msi_pro_log)

msi_ihc <- msi[,c("MLH1","PMS2","MSH2","MSH6")]
rownames(msi_ihc) = msi$Sample
msi_ihc <- t(msi_ihc)

msi_ihc <- msi_ihc[, msi$Sample]

msi_pro_log_scaled <- scale(msi_pro_log, center = TRUE, scale = FALSE)
msi_pro_log_scaled <- msi_pro_log_scaled[, msi$Sample]

col_scale_pro <- colorRamp2(c(-2, 0, 2), c("#6A2B6F", "white" , "#FF895D"))

kpro <- ComplexHeatmap::Heatmap(
  msi_pro_log_scaled,
  name = "z-score",
  cluster_rows = FALSE,
  cluster_columns = FALSE,
  show_column_names = FALSE,
  col = col_scale_exp,
  na_col = "#d9d9d9",
  border = TRUE, 
  row_title = "Protein",
  row_names_gp = grid::gpar(fontsize = 8),
  cell_fun = function(j, i, x, y, width, height, fill) {
    grid.text(
      msi_ihc[i, j], 
      x = x, y = y,
      gp = gpar(fontsize = 6, col = "white")
    )
  }
)

# expression
msi_exp <- msi[,c("MLH1_Exp","PMS2_Exp","MSH2_Exp","MSH6_Exp")]
msi_exp$MLH1_Exp <- as.numeric(msi_exp$MLH1_Exp)
msi_exp$MSH2_Exp <- as.numeric(msi_exp$MSH2_Exp)
msi_exp$MSH6_Exp <- as.numeric(msi_exp$MSH6_Exp)
msi_exp$PMS2_Exp <- as.numeric(msi_exp$PMS2_Exp)
rownames(msi_exp) = msi$Sample

msi_exp_log <- log2(msi_exp)
rownames(msi_exp_log) = msi$Sample
msi_exp_log <- t(msi_exp_log)
msi_exp_log_scaled <- scale(msi_exp_log, center = TRUE, scale = FALSE)

msi_exp_log_scaled = msi_exp_log_scaled[, msi$Sample]

kexp <- ComplexHeatmap::Heatmap(
  msi_exp_log_scaled,
  name = "z-score",
  cluster_rows = FALSE,
  cluster_columns = FALSE,
  col = col_scale_exp,
  na_col = "#d9d9d9",
  row_title = "Expression",
  row_names_gp = grid::gpar(fontsize = 8),
  border = TRUE, 
  show_column_names = FALSE)

## methylation
msi_met = msi[,c("MLH1_tumor","PMS2_tumor",
                 "MSH2_tumor","MSH6_tumor")]
colnames(msi_met) = c("MLH1_Met","PMS2_Met","MSH2_Met","MSH6_Met")
rownames(msi_met) = msi$Sample
msi_met = t(msi_met)

msi_met = msi_met[, msi$Sample]

col_scale_met <- colorRamp2(
  breaks = c(0, 0.01, 0.2, 0.5),
  colors = c("#0074FE", "#96EBF9", "#FEE900", "#F00003")
)

col_scale_met <- colorRamp2(c(0, 0.05, 0.5), c("white", "#dadaeb", "#3f007d"))
kmet <- ComplexHeatmap::Heatmap(
  msi_met,
  name = "MET",
  cluster_rows = FALSE,
  cluster_columns = FALSE,
  show_column_names = FALSE,
  row_names_gp = grid::gpar(fontsize = 8),
  row_title = "Methylation",
  border = TRUE, 
  col = col_scale_met)

draw(kmet)

ht_list_k = kmut_main %v% kpro %v% kexp %v% kmet
draw(ht_list_k)

p2 <- grid.grabExpr(draw(ht_list_k), wrap = TRUE, warn = FALSE)


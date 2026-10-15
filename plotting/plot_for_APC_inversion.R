library(ComplexHeatmap)
library(ggplot2)
library(tidyverse)
library(gridExtra)
library(ggpubr)
library(colorRamp2)

############################################ color set  ##############################################
col_fun_rna <- circlize::colorRamp2(
  breaks = c(-2, 0, 2), 
  colors = c("blue", "white", "red")
)
## color for mutation type
col = c("Missense" = "#3287bc", 
        "TRUNC" = "#acd7a3",
        "chr5:112815487_A>G" = "#d43f50", 
        "other_splice" = "#f9ac61",
        "LOH" = "#fcde89", 
        "HOMDEL" = "#2166ac", 
        "HETLOSS" = "#92c5de",
        "INV_hotspot" = "#ba2025",
        "BND_hotspot" = "#374ea1")

## color for clinical information
col_age <- c("<60" = "#efedf5",">=60"="#9970ab","UNKNOW"="#e0e0e0")
col_loc <- c("Left_colon" = "#dddb1f","Right_colon"="#7FC9CF","UNKNOW"="#e0e0e0")

exp <- read.table("expression.txt.gz", sep="\t", stringsAsFactors = F, header=T,  check.names = F)
exp <- exp %>%
  mutate(across(where(is.numeric), ~ log2(. + 1))) %>%  
  relocate(names(exp)[1], .before = everything())
exp <- exp %>% rename_with(~ paste0(., "_EXP"), -Sample)

pro <- read.table("protein.txt.gz", sep="\t", stringsAsFactors = F, header=T,  check.names = F)
colnames(pro) <- make.unique(colnames(pro), sep = "_")
pro <- pro %>%
  mutate(across(where(is.numeric), ~ log2(. + 1))) %>%  
  relocate(names(pro)[1], .before = everything())
pro <- pro %>% rename_with(~ paste0(., "_PRO"), -Sample)


status <- read.table("apc_all_mut_clinic.txt", sep="\t", stringsAsFactors = F, header=T,  check.names = F)
status <- status %>% 
  mutate(Agetype=case_when(Age == "UNKNOW" ~ "UNKNOW",
                           Age < 60 ~ "<60",
                           Age >=60 ~ ">=60"))

merge_data_all <- status %>%
  left_join(exp[,c("Sample","APC_EXP")], by = "Sample") %>%
  left_join(pro[,c("Sample","APC_PRO")], by = "Sample")


############################################ heatmap ##############################################
df <- merge_data_all[,c("APC_EXP","APC_PRO")]
df <- as.data.frame(scale(df))
df[df > 2] = 2
df[df < -2] = -2

df$Sample = merge_data_all$Sample
df$Type = merge_data_all$GenomicSubtype2025
df$Location = merge_data_all$Location_2group
df$Age = merge_data_all$Agetype
rownames(df) = df$Sample


## input mutations
mat <- read.table("apc_for_oncoprint.txt",header = T,stringsAsFactors = F,sep = '\t',check.names = F,row.names = 1)
mat[is.na(mat)] <- ""

alter_fun = list(
  background = alter_graphic("rect", fill = "#f0f0f0"),   
  TRUNC = alter_graphic("rect", fill = col["TRUNC"]),
  Missense = alter_graphic("rect", height = 0.5, fill = col["Missense"]),
  other_splice = alter_graphic("rect", fill = col["other_splice"]),
  `chr5:112815487_A>G` = alter_graphic("rect", fill = col["chr5:112815487_A>G"]),
  HOMDEL = alter_graphic("rect", fill = col["HOMDEL"]),
  HETLOSS = alter_graphic("rect", fill = col["HETLOSS"]),
  BND_hotspot = alter_graphic("rect", fill = col["BND_hotspot"]),
  INV_hotspot = alter_graphic("rect", fill = col["INV_hotspot"]),
  LOH = alter_graphic("rect", height = 0.5, fill = col["LOH"]))


## sort
mat_tran <- as.data.frame(t(mat))
mat_tran[mat_tran == ""] <- "WT"
mat_tran$Sample <- rownames(mat_tran)
merge_df <- merge(df, mat_tran, by = "Sample")

merge_df <- merge_df[order(merge_df$APC_SV,
                           merge_df$APC_noncoding,
                           merge_df$APC_coding,
                           merge_df$APC_CNV),]


mat_sort <- mat[,merge_df$Sample]
onco_obj <- oncoPrint(mat_sort, get_type = function(x) strsplit(x, ";")[[1]],
          right_annotation = NULL,
          top_annotation = HeatmapAnnotation(df = merge_df[,c("Location","Age","APC_EXP")],
                                             show_annotation_name = FALSE,
                                             simple_anno_size = unit(0.7, "cm"),
                                             show_legend = FALSE,
                                             col = list(Age=col_age,
                                                        Location = col_loc,
                                                        APC_EXP = col_fun_rna)),
          show_row_names = TRUE,
          column_order = colnames(mat_sort),
          row_order = c("APC_SV","APC_noncoding","APC_coding","APC_CNV"),
          row_names_gp = gpar(fontsize = 10),
          alter_fun = alter_fun, 
          show_heatmap_legend = FALSE,
          col = col)


mutation_legend <- Legend(
  title = "Mutation Types",
  at = names(col),
  labels = names(col),
  legend_gp = gpar(fill = col),
  nrow = 2
)

age_legend <- Legend(
  title = "Age",
  at = c("<60", ">=60", "UNKNOW"),
  legend_gp = gpar(
    fill = col_age
  ),
  nrow = 2
)

loc_legend <- Legend(
  title = "Location",
  at = c("Left_colon", "Right_colon", "UNKNOW"),
  legend_gp = gpar(
    fill = col_loc
  ),
  nrow = 2
)

exp_legend <- Legend(
  title = "APC_EXP",
  col_fun = col_fun_rna,
  direction = "horizontal"
)


packed_legends <- packLegend(
  mutation_legend, exp_legend,
  direction = "horizontal",
  row_gap = unit(5, "mm") 
)

ht <- draw(
  onco_obj,
  heatmap_legend_side = "bottom",
  annotation_legend_list = packed_legends,
  annotation_legend_side = "bottom"
)

p2 <- grid.grabExpr(draw(ht), wrap = TRUE, warn = FALSE)
grid.arrange(p1, p2, ncol = 1, heights = c(1.5, 1))

############################################ boxplot ##############################################

mut_exp <- read.table("apc_mutation_merge_exp.txt", sep="\t", stringsAsFactors = F, header=T,  check.names = F)
mut_exp$APC_EXP = log2(mut_exp$APC_EXP + 1)
mut_exp$MutType <- mut_exp$MutTyp2
colnames(mut_exp)[3] = "MutType"

group_counts <- table(mut_exp$MutType)
groups <- names(group_counts)[group_counts > 5]
filtered_dat <- mut_exp[ mut_exp$MutType %in% groups & mut_exp$MutType != "", ]

filtered_dat <- filtered_dat %>%
  mutate(MutType = factor(
    MutType,
    levels = names(sort(tapply(APC_EXP, MutType, mean, na.rm = TRUE))),
    ordered = TRUE
  ))

non_wt_levels <- setdiff(levels(filtered_dat$MutType), "WT")

box_cols = c("#d43f50","#f9ac61","#4575b4","#af8dc3", "#92c5de","#acd7a3","gray")
ggplot(data = filtered_dat, aes(x = MutType, y = APC_EXP)) + 
  geom_boxplot(width = 0.4,
               outlier.shape = NA,
               color="#636363") + 
  geom_point(aes(color=MutType),
             position = position_jitter(width=0.15, height=0.15),
             alpha=0.5,size=1) +
  labs(y = "APC RNA expression (log2)", x="") + 
  scale_color_manual(values = box_cols) +
  scale_x_discrete(
    name = "",
    labels = function(x) {
      counts <- table(filtered_dat$MutType)[x]
      paste0(x, "\n(n=", counts, ")")
    }
  ) +
  stat_compare_means(
    method = "wilcox.test",
    comparisons = c(
      list(c("HETLOSS", "HETLOSS/INV")),
      list(c("INV/TRUNC", "TRUNC"))),
    size = 3,
    step.increase = 0.1,
    tip.length = 0.02,
    vjust = 0.1
  ) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1, size = 8, colour = "black"),
        axis.line = element_line(colour = "black"),
        panel.grid.major = element_line(colour=NA),
        panel.grid.minor = element_line(colour=NA),
        panel.spacing = unit(0.75, "lines"),
        panel.background = element_rect(fill=NA, color=NA),
        legend.position = "none")


############################################ barplot ##############################################
genomic_data <- read.table("apc_hotspot_inv_for_plot.txt", header = TRUE, sep = "\t", check.names = FALSE)
genomic_data$Sample <- factor(genomic_data$Sample, levels = unique(genomic_data$Sample))

focus_start <- 112682736
focus_end <- 112794558
scale_factor_outside <- 100 

trans_custom <- trans_new(
  name = "focus_zoom",
  transform = function(x) {
    ifelse(x >= focus_start & x <= focus_end,
           x, 
           ifelse(x < focus_start,
                  focus_start - (focus_start - x) / scale_factor_outside, 
                  focus_end + (x - focus_end) / scale_factor_outside)
    )
  },
  inverse = function(x) {
    ifelse(x >= focus_start & x <= focus_end,
           x,
           ifelse(x < focus_start,
                  focus_start - (focus_start - x) * scale_factor_outside, 
                  focus_end + (x - focus_end) * scale_factor_outside) 
    )
  }
)

format_m <- function(x) {
  paste0(round(x / 1e6, 2), "M")
}

ggplot(genomic_data) +
  geom_segment(aes(x = start, 
                   xend = end, 
                   y = Sample, 
                   yend = Sample), 
               size = 2,
               color = "#4393c3") +
  geom_vline(xintercept = c(focus_start, focus_end), 
             linetype = "dashed", 
             color = "red",
             size = 0.8) +
  scale_x_continuous(trans = trans_custom, labels = format_m) +
  labs(x= "Chromosome 5", y = "COACH-CRC (n=26)") +
  theme_bw() +
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        axis.title = element_text(size = 10, face = "bold"))

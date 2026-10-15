library(MOVICS)
library(dplyr)
library(survival)
library(survminer)
library(networkD3)


rate <- 0.1          # top 10% for RNA, protein and phosphosites
met_rate <- 0.05     # top 5% for methylation

read_and_filter <- function(file_path, samples) {
  df <- read.table(file_path, sep = "\t", header = TRUE, row.names = 1, check.names = FALSE)
  df <- df[, intersect(samples, colnames(df)), drop = FALSE]
}

samples <- read.table("input/sample.lst", sep = "\t", header = TRUE, check.names = FALSE)$Sample
exp_raw <- read_and_filter("input/fpkm_mRNA.txt", samples)
cnv_raw <- read_and_filter("input/focal_value.txt", samples)
met_raw <- read_and_filter("input/DMR.txt", samples)
pro_raw <- read_and_filter("input/protein_FOTknn.txt", samples)
pho_raw <- read_and_filter("input/phoSite_knn.txt", samples)


# RNA expression (top 10%)
exp_mad <- apply(exp_raw, 1, mad)
exp_filtered <- exp_raw[rev(order(exp_mad))[1:round(length(exp_mad) * rate)], ]

# Proteins (top 10%)
pro_mad <- apply(pro_raw, 1, mad)
pro_filtered <- pro_raw[rev(order(pro_mad))[1:round(length(pro_mad) * rate)], ]

# Phosphosites (top 10%)
pho_mad <- apply(pho_raw, 1, mad)
pho_filtered <- pho_raw[rev(order(pho_mad))[1:round(length(pho_mad) * rate)], ]

# DMRs (top 5%)
met_mad <- apply(met_raw, 1, mad)
met_filtered <- met_raw[rev(order(met_mad))[1:round(length(met_mad) * met_rate)], ]

# Construct mo5.data list
mo5.data <- list(
  exp = log2(exp_filtered + 1),
  cnv = cnv_raw,
  met = met_filtered,
  pro = pro_filtered,
  pho = pho_filtered
)


common_samples <- Reduce(intersect, lapply(mo5.data, colnames))
mo5.data <- lapply(mo5.data, function(x) x[, common_samples, drop = FALSE])

optk <- getClustNum(
  data        = mo5.data,
  is.binary   = c(FALSE, FALSE, FALSE, FALSE, FALSE),
  try.N.clust = 2:10,
  fig.name    = "omics_cluster"
)

nemo_raw <- getMOIC(
  data        = mo5.data,
  methodslist = list("NEMO"),
  N.clust     = 3,
  type        = rep("gaussian", 5)
)

nemo_k3 <- data.frame(
  samID = rownames(nemo_raw),
  CMOS  = paste0("CMOS", nemo_raw[, "NEMO"] + 1),
  stringsAsFactors = FALSE
)

plotdata <- getStdiz(
  data       = mo5.data,
  halfwidth  = c(2, 1, 2, 2, 2),
  centerFlag = c(TRUE, FALSE, FALSE, TRUE, TRUE),
  scaleFlag  = c(TRUE, FALSE, FALSE, TRUE, TRUE)
)


cnv.col  <- c("#2872A7", "white", "#E21E27")
exp.col  <- c("#00FF00", "#000000", "#FF0000")
met.col  <- c("#0074FE", "#96EBF9", "#FEE900", "#F00003")
pro.col  <- c("#01665e", "white", "#8c510a")
phos.col <- c("#6A2B6F", "white", "#FF895D")

getMoHeatmap(
  data          = plotdata,
  row.title     = c("Expression", "focalCNV", "Methylation", "Proteomics", "Phosphosite"),
  is.binary     = c(FALSE, FALSE, FALSE, FALSE, FALSE),
  legend.name   = c("EXP", "CNV", "MET", "PRO", "PHOS"),
  clust.res     = nemo_k3,
  clust.dend    = NULL,
  show.row.dend = rep(FALSE, 5),
  show.rownames = rep(FALSE, 5),
  show.col.dend = FALSE,
  show.colnames = FALSE,
  color         = list(exp.col, cnv.col, met.col, pro.col, phos.col),
  width         = 10,
  height        = 3,
  fig.name      = "NEMO_heatmap"
)


clin <- read.table("input/clinical_genetic_profile.txt",
                   header = TRUE, stringsAsFactors = FALSE, sep = "\t")

dat <- clin %>%
  left_join(nemo_k3, by = c("Sample" = "samID")) %>%
  mutate(CMOS = if_else(is.na(CMOS), "CMOS1", CMOS))

dat <- dat %>% mutate(CMS = Final_CMS)
sankey_df <- dat %>%
  count(CMOS, CMS) %>%
  mutate(across(c(CMOS, CMS), as.factor))

nodes <- data.frame(node = unique(c(as.character(sankey_df$CMOS),
                                    as.character(sankey_df$CMS))))
sankey_df <- sankey_df %>%
  mutate(
    IDsource = match(CMOS, nodes$node) - 1,
    IDtarget = match(CMS, nodes$node) - 1
  )

sankeyNetwork(
  Links = sankey_df, Nodes = nodes,
  Source = "IDsource", Target = "IDtarget", Value = "n", NodeID = "node",
  iterations = 0, fontSize = 30, nodeWidth = 50
)

dat <- dat %>%
  mutate(
    OS_5year      = as.numeric(OS_5year),
    Death_5year   = as.numeric(Death_5year),
    DFS_5year     = as.numeric(DFS_5year),
    Relapse_Death = as.numeric(Relapse_Death_5year)
  )

# Overall survival
surv_os <- survfit(Surv(OS_5year, Death_5year) ~ CMOS, data = dat)
s1 <- ggsurvplot(
  surv_os, data = dat,
  pval = TRUE, pval.method = FALSE, pval.coord = c(45, 0.5),
  risk.table = TRUE, legend.title = "", legend = c(0.3, 0.5),
  ggtheme = theme_bw() + theme(panel.grid = element_blank()),
  risk.table.pos = 'in',
  xlab = "Months", ylab = "Overall survival",
  font.x = 12, font.y = 12, font.tickslab = c(10, "black")
)

# Disease-free survival
surv_dfs <- survfit(Surv(DFS_5year, Relapse_Death) ~ CMOS, data = dat)
s2 <- ggsurvplot(
  surv_dfs, data = dat,
  pval = TRUE, pval.method = FALSE, pval.coord = c(45, 0.5),
  risk.table = TRUE, legend.title = "", legend = c(0.3, 0.5),
  ggtheme = theme_bw() + theme(panel.grid = element_blank()),
  risk.table.pos = 'in',
  xlab = "Months", ylab = "Disease-free survival",
  font.x = 12, font.y = 12, font.tickslab = c(10, "black")
)

arrange_ggsurvplots(list(s1, s2), print = TRUE, ncol = 2, nrow = 1)

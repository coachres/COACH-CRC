library(caret)
library(randomForest)
library(dplyr)

## differentially expressed genes (DEGs) from pairwise comparisons among CMOS subtypes
feat <- read.table("cmos_significant_rna_gene.txt", header = T, check.names = F, sep = '\t')

## bulk RNA (COACH-CRC)
coach <- read.csv("coach_exp_common_genes.csv.gz", header = T, check.names = F, row.names = 1)
coach = as.data.frame(t(coach))

## single-cell pseudobulk
cell <- read.csv("singleCell_exp_common_genes.csv.gz", header = T, check.names = F, row.names = 1)

## top 10% most variable genes (by MAD) expression data for multi-omics-based molecular classification (CMOS subtype)
mo5.data <- readRDS("omics5_mRNA_mad10.rds")
mad_top_gene = row.names(mo5.data$exp)
union_genes <- union(feat$Gene, mad_top_gene)

## candidate feature set : union(DEGs, top-variable genes), kept only if present in both datasets
union_genes  <- union(feat$Gene, mad_top_gene)
common_genes <- intersect(union_genes, intersect(colnames(coach), colnames(cell)))
message(sprintf("Union genes: %d | Common genes used: %d",
                length(union_genes), length(common_genes)))

preprocess <- function(mat) {
  m <- as.matrix(mat)
  log2(m + 1)
}

cell_exp   <- preprocess(cell[, common_genes, drop = FALSE])
coach_exp  <- preprocess(coach[, common_genes, drop = FALSE])
coach_exp$Sample <- rownames(coach_exp)

subtype <- read.table("cmos_information.txt",header = T, stringsAsFactors = F, sep = '\t')
subtype <- subtype %>% filter(CMOS != "CMOS1")
coach_exp <- subtype[,c("Sample","CMOS")] %>% left_join(coach_exp, by = "Sample")

set.seed(1234)
train_index <- createDataPartition(coach_exp$CMOS, p = 0.7, list = FALSE)
train_data <- coach_exp[train_index, ]
row.names(train_data) = train_data$Sample
train_labels = factor(train_data$CMOS)
names(train_labels) <- train_data$Sample

test_data  <- coach_exp[-train_index, ]
row.names(test_data) = test_data$Sample
test_labels = factor(test_data$CMOS)
names(test_labels) <- test_data$Sample


## randomForest model
mtry_value <- floor(sqrt(ncol(train_data[,-c(1:2)])))
if (mtry_value < 1) mtry_value <- 1

class_freq <- table(train_labels)
class_weights <- 1 / class_freq
class_weights <- class_weights / sum(class_weights)
min_size <- min(class_freq)
samp_sizes <- rep(min_size, length(class_freq))

rf_model <- randomForest(
  x = scale(train_data[,-c(1:2)]),
  y = train_labels,
  ntree = 1000,
  strata = train_labels,
  sampsize = samp_sizes,
  classwt = class_weights,
  mtry = mtry_value,
  nodesize = 3,
  importance = TRUE,
  proximity = TRUE
)

cmos_test_pred <- predict(rf_model, scale(test_data[,-c(1:2)]), type = "response")
cmos_test_pred <- test_data[,-c(1:2)] %>%
  transmute(id = rownames(test_data)) %>%
  cbind(data.frame(CMOS_pred = cmos_test_pred)) %>%
  inner_join(data.frame(id = names(test_labels), CMOS = test_labels), by = "id")

cm = confusionMatrix(cmos_test_pred$CMOS_pred, cmos_test_pred$CMOS)
metrics <- data.frame(
  Feature_Count = dim(test_data[,-c(1:2)])[2],
  Class = gsub("Class: ", "", rownames(cm$byClass)),
  Sensitivity = cm$byClass[, "Sensitivity"],
  Specificity = cm$byClass[, "Specificity"],
  Balanced_Accuracy = cm$byClass[, "Balanced Accuracy"]
)

print(metrics)

test_prob_matrix <- predict(rf_model, scale(test_data[,-c(1:2)]), type = "prob")
test_predictions <- predict(rf_model, scale(test_data[,-c(1:2)]), type = "response")

## probability-based confidence filtering
top2_diff <- function(x) {
  sorted <- sort(x, decreasing = TRUE)
  sorted[1] - sorted[2]
}

test_prob_diff <- apply(test_prob_matrix, 1, top2_diff)
test_max_prob <- apply(test_prob_matrix, 1, max)

cmos_test_pro_pred <- ifelse(test_max_prob > 0.4 & test_prob_diff > 0.1, 
                             as.character(test_predictions), NA)

cmos_test_pro_pred <- test_data[,-c(1:2)] %>%
  transmute(id = rownames(test_data)) %>%
  cbind(data.frame(CMOS_pred = cmos_test_pro_pred)) %>%
  inner_join(data.frame(id = names(test_labels), CMOS = test_labels), by = "id")

cmos_test_pro_pred <- cmos_test_pro_pred[!is.na(cmos_test_pro_pred$CMOS_pred), ]

cmos_test_pro_pred$CMOS_pred <- factor(cmos_test_pro_pred$CMOS_pred)
cmos_test_pro_pred$CMOS <- factor(cmos_test_pro_pred$CMOS)

dim(cmos_test_pro_pred)

cm_pro = confusionMatrix(cmos_test_pro_pred$CMOS_pred, cmos_test_pro_pred$CMOS)

metrics_pro <- data.frame(
  Feature_Count = dim(test_data[,-c(1:2)])[2],
  Class = gsub("Class: ", "", rownames(cm_pro$byClass)),
  Sensitivity = cm_pro$byClass[, "Sensitivity"],
  Specificity = cm_pro$byClass[, "Specificity"],
  Balanced_Accuracy = cm_pro$byClass[, "Balanced Accuracy"]
)

print(metrics_pro)


## Predict single-cell pseudobulk (MSS samples only)
info <- read.csv("singlecell_sample_clinical_info.csv",header = T, check.names = F)
cell_mss_sample = info[!is.na(info$microsatellite_status) & 
                     info$microsatellite_status!="MSI" & 
                     info$microsatellite_status!="MSI-H" & 
                     info$microsatellite_status!="", ]$sample_id

cell_mss_sample = intersect(cell_mss_sample, row.names(cell_exp))
cell_data = cell_exp[cell_mss_sample,common_genes]
dim(cell_data)

cell_prob_matrix <- predict(rf_model, scale(cell_data), type = "prob")
cell_predictions <- predict(rf_model, scale(cell_data), type = "response")

cell_prob_diff <- apply(cell_prob_matrix, 1, top2_diff)
cell_max_prob <- apply(cell_prob_matrix, 1, max)

cell_filtered_pred <- ifelse(cell_max_prob > 0.4 & cell_prob_diff > 0.1, 
                        as.character(cell_predictions), NA)

cell_prob_result <- cbind(cell_prob_matrix, 
                          data.frame(cell_prob_diff = cell_prob_diff,
                                     cell_pred_classes = cell_predictions,
                                     cmos_cell_pro_pred = cell_filtered_pred))

table(cell_prob_result$cmos_cell_pro_pred)
dim(cell_prob_result)
cell_prob_result$sample_id = rownames(cell_prob_result)

cell_prob_result <- cell_prob_result %>%
  left_join(info, by = "sample_id") 

write.csv(cell_prob_result, file = "cell_predict_result.csv", quote = F, row.names = F)

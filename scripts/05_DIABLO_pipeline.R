
# ================================
#         Load & Install Packages
# ================================

# Install BiocManager if not already installed
if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")

# Install mixOmics package
BiocManager::install("mixOmics", update = FALSE, ask = FALSE)

# Set CRAN repository and install additional packages
options(repos = c(CRAN = "https://cran.r-project.org"))
install.packages(c("ragg", "ggplot2", "pheatmap", "igraph", "readxl"), dependencies = TRUE, update = FALSE)

# Load Libraries
library(circlize)

library(ragg)
library(mixOmics)
library(ggplot2)
library(igraph)
library(readxl)
library(Rtsne)
library(readr)
library(pheatmap)
library(cluster)
library(caret)
library(tools)
library(pls)
library(ComplexHeatmap)



# ================================
#              Load Data
# ================================


# Load Data
setwd("C:\\Users\\Abdulaziz\\Desktop\\DIABLO_Data")
metadata <- read_excel("Metad.xlsx")
methylation_data <- read_csv("mmVals_sig.csv", col_names = TRUE)
expression_data <- read_csv("Normalized_count_matrix.csv", col_names = TRUE)

# Extract Feature Names
feature_names_meth <- methylation_data[[1]]  # CpG sites
feature_names_expr <- expression_data[[1]]   # Gene names

# Extract Sample Names from Column Headers
sample_names <- colnames(methylation_data)[-1]       
sample_names_expr <- colnames(expression_data)[-1]    

# Remove Feature Names Column & Convert to Matrix
methylation_matrix <- as.matrix(methylation_data[, -1])
expression_matrix <- as.matrix(expression_data[, -1])

# Assign Correct Row and Column Names BEFORE Transposing
rownames(methylation_matrix) <- feature_names_meth
rownames(expression_matrix) <- feature_names_expr
colnames(methylation_matrix) <- sample_names
colnames(expression_matrix) <- sample_names_expr

# Transpose the Data
T_Methylation_data <- t(methylation_matrix)
T_expression_data <- t(expression_matrix)

# Assign Sample Names as Row Names
rownames(T_Methylation_data) <- sample_names
rownames(T_expression_data) <- sample_names_expr

# Check Row Names
print(rownames(T_Methylation_data))
print(rownames(T_expression_data))

# Convert to numeric WITHOUT losing row names
mode(T_Methylation_data) <- "numeric"
mode(T_expression_data) <- "numeric"



# ================================
#     Sample Matching & Cleaning
# ================================


# Get common samples
common_samples <- intersect(rownames(T_Methylation_data), rownames(T_expression_data))

# Subset data (ONLY if common_samples is not empty)
if (length(common_samples) > 0) {
  T_Methylation_data <- T_Methylation_data[common_samples, ]
  T_expression_data <- T_expression_data[common_samples, ]
} else {
  stop("No matching sample names found between datasets! Check row names before subsetting.")
}

# Verify that both matrices have the same samples
print(dim(T_Methylation_data))
print(dim(T_expression_data))

# Check for NA values and impute if necessary
imputeNA <- function(mat) {
  for (j in 1:ncol(mat)) {
    if (any(is.na(mat[, j]))) {
      mat[is.na(mat[, j]), j] <- median(mat[, j], na.rm = TRUE)
    }
  }
  return(mat)
}
T_Methylation_data <- imputeNA(T_Methylation_data)
T_expression_data <- imputeNA(T_expression_data)

DMS <- T_Methylation_data
RNA_seq <- T_expression_data

# Adjust sample names
rownames(DMS) <- gsub("Patient ", "", rownames(DMS))
rownames(RNA_seq) <- gsub("Patient ", "", rownames(RNA_seq))
metadata$Sample_id_mrna <- gsub("-", ".", metadata$Sample_id_mrna)
metadata$Sample_id_cpg <- gsub("-", ".", metadata$Sample_id_cpg)

# Ensure matching samples
common_samples <- intersect(rownames(DMS), rownames(RNA_seq))
DMS <- DMS[common_samples, , drop = FALSE]
RNA_seq <- RNA_seq[common_samples, , drop = FALSE]
metadata <- metadata[metadata$sample %in% common_samples, ]
metadata <- metadata[match(common_samples, metadata$sample), ]


colnames(DMS) <- paste0("meth_", colnames(DMS))
colnames(RNA_seq) <- paste0("expr_", colnames(RNA_seq))



# ================================
#       Train/Test Split
# ================================


# Prepare omics data
omics_data <- list(Methylation = as.matrix(DMS), Transcriptomics = as.matrix(RNA_seq))
Condition <- factor(metadata$Sample_source_name_ch1_cpg)

# Split Data into Training and Testing Sets
set.seed(123)
train_samples <- sample(rownames(omics_data$Methylation), size = round(0.8 * nrow(omics_data$Methylation)), replace = FALSE)

train_data <- list(
  Methylation = omics_data$Methylation[train_samples, ],
  Transcriptomics = omics_data$Transcriptomics[train_samples, ]
)

test_data <- list(
  Methylation = omics_data$Methylation[-which(rownames(omics_data$Methylation) %in% train_samples), ],
  Transcriptomics = omics_data$Transcriptomics[-which(rownames(omics_data$Transcriptomics) %in% train_samples), ]
)

train_condition <- Condition[rownames(omics_data$Methylation) %in% train_samples]


# ================================
#       Filter Low Variance Genes
# ================================







filtered_RNA <- RNA_seq

filtered_DMS <- DMS

# توحيد العينات
common_samples <- intersect(rownames(filtered_RNA), rownames(filtered_DMS))
filtered_RNA <- filtered_RNA[common_samples, ]
filtered_DMS <- filtered_DMS[common_samples, ]
train_condition <- metadata$Sample_source_name_ch1_cpg
train_condition <- train_condition[match(common_samples, metadata$sample)]
train_condition <- as.factor(train_condition)



# ================================
#          Build DIABLO Model
# ================================

# نحدد عدد الميزات المختارة من كل بلوك
list.keepX <- list(
  Transcriptomics = 50,
  Methylation = 50
)


# مصفوفة التصميم
design <- matrix(0.1, ncol = length(train_data), nrow = length(train_data), 
                 dimnames = list(names(train_data), names(train_data)))
diag(design) <- 0


# بناء الموديل
set.seed(123)
diablo_model <- block.splsda(
  X = train_data,
  Y = train_condition,
  ncomp = 5,
  keepX = list.keepX,
  design = design
)
plotDiablo(diablo_model, ncomp = 1)
plotDiablo(diablo_model, ncomp = 2)
plotDiablo(diablo_model, ncomp = 3)
plotDiablo(diablo_model, ncomp = 4)
plotDiablo(diablo_model, ncomp = 5)


# استخراج الميزات المختارة
selected_features <- selectVar(diablo_model, comp =2)
selected_genes <- selected_features$Transcriptomics$name
selected_cpgs <- selected_features$Methylation$name

# نرجع البيانات الأصلية للميزات المختارة
selected_DMS <- DMS[, colnames(DMS) %in% selected_cpgs]
selected_RNA <- RNA_seq[, colnames(RNA_seq) %in% selected_genes]

# توحيد العينات مرة تانية
common_samples_final <- intersect(rownames(selected_DMS), rownames(selected_RNA))
selected_DMS <- selected_DMS[common_samples_final, ]
selected_RNA <- selected_RNA[common_samples_final, ]

# حساب الارتباط بين الجينات والـ CpGs
cor_matrix <- cor(selected_DMS, selected_RNA, method = "pearson")

# إعداد التدرج اللوني
col_fun <- colorRamp2(c(-1, 0, 1), c("blue", "white", "red"))  # تدرج ألوان من -1 إلى 1

# حفظ الخريطة في ملف PDF
pdf("complex_heatmap.pdf", width = 10, height = 10)

# رسم الخريطة الحرارية
Heatmap(cor_matrix,
        name = "Correlation",                # عنوان شريط الألوان
        col = col_fun,                       # التدرج اللوني
        cluster_rows = TRUE,                 # تجميع الصفوف
        cluster_columns = TRUE,              # تجميع الأعمدة
        show_row_names = TRUE,               # عرض أسماء الصفوف
        show_column_names = TRUE             # عرض أسماء الأعمدة
)

# إغلاق ملف PDF
dev.off()


write.csv(cor_matrix, file = "correlation_matrix.csv")

pdf("circosPlot.pdf", width = 10, height = 10)  # فتح ملف PDF بالحجم المطلوب
# Circos plot
circosPlot(diablo_model, cutoff = 0.8, line = TRUE, 
           comp = 1, title = "Circos Plot: CpGs vs Genes", 
           color.blocks = c("orchid", "lightgreen"))

dev.off()


pdf("plotIndiv.pdf", width = 10, height = 10)  # فتح ملف PDF بالحجم المطلوب
# Sample plot (PCA-like plot)
plotIndiv(diablo_model, comp = c(1,2), group = train_condition,
          ind.names = FALSE, legend = TRUE, 
          title = "DIABLO Sample Plot")
dev.off()


# plot Arrow
pdf("plotArrow.pdf", width = 10, height = 10)  # Open a PDF file
plotArrow(diablo_model, ind.names = FALSE, legend = TRUE, 
          title = 'TCGA, DIABLO comp 1 - 2')
dev.off()

# Loading plot for component 1
pdf("plotLoadings_mRNA.pdf", width = 10, height = 10)  # Open a PDF file
plotLoadings(diablo_model, comp = 1, block = 'Transcriptomics', method = 'mean', contrib = 'max', name.var = FALSE)
dev.off()


# ROC Curve
pdf("roc_curves_blocks.pdf", width = 10, height = 10)
par(mfrow = c(1, 2))  # لو عندك 2 omics blocks

# ROC for Transcriptomics
auroc(diablo_model, roc.block = "Transcriptomics", roc.comp = 1:2, print = TRUE)

# ROC for Methylation
auroc(diablo_model, roc.block = "Methylation", roc.comp = 1:2, print = TRUE)

dev.off()
roc_results <- auroc(diablo_model, roc.comp = 1:2, print = FALSE)
saveRDS(roc_results, "roc_auc_values.rds")



# DMS
pdf("plotLoadings_DMS.pdf", width = 8, height = 8)  # Open a PDF file
plotLoadings(diablo_model, comp = 1, block = 'Methylation', method = 'mean', contrib = 'max', name.var = TRUE, size.name = 0.5)
dev.off()


# Biplot
pdf("biplot.pdf", width = 8, height = 8)  # Open a PDF file
biplot(diablo_model, comp = 1:2, group = train_condition, 
       col = c("lightgreen", "orchid"), 
       title = "Biplot: DIABLO Components vs Samples")
dev.off()


































plotIndiv(diablo_model, comp = c(1,2), 
          group = train_condition, 
          legend = TRUE, 
          title = "Individual Plot (comp 1 & 2)")


plotArrow(diablo_model, comp = 1:2, 
          group = train_condition, 
          legend = TRUE, 
          title = "Arrow Plot")

# مشكلة 
plotVar(diablo_model, comp = 2, 
        var.names = TRUE, 
        style = "graphics", 
        legend = TRUE)


circosPlot(diablo_model, 
           comp = 1, 
           cutoff = 0.85,  # You can adjust based on correlation strength
           line = TRUE, 
           color.blocks = c("blue", "red"), 
           title = "Circos Plot: Genes vs CpGs")



# RNA
plotLoadings(diablo_model, comp = 1, block = 'Transcriptomics', method = 'mean', contrib = 'max', name.var = FALSE)
# DMS
plotLoadings(diablo_model, comp = 1, block = 'Methylation', method = 'mean', contrib = 'max', name.var = TRUE, size.name = 0.5)


perf_diablo <- perf(diablo_model, validation = "Mfold", folds = 5, nrepeat = 10, progressBar = TRUE)

auc_results <- perf_diablo$AUC


auc_results <- auc(diablo_model, roc.comp = 1)
print(auc_results)



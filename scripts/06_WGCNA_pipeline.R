#Step 1: Install and Load Required Packages

# Install Bioconductor manager (if not installed)
if (!requireNamespace("BiocManager", quietly = TRUE)) install.packages("BiocManager")

# Install necessary packages
BiocManager::install(c("WGCNA", "flashClust", "dynamicTreeCut", "GO.db"))

# Load the packages
library(WGCNA)
library(flashClust)
library(dynamicTreeCut)

# Allow large dataset handling
options(stringsAsFactors = FALSE)

# Set working directory (change this to your dataset location)
setwd("your/directory/path")



#Step 2: Load Your transcriptomics and methylation Data

# Load methylation and transcriptomics data
meth_data <- read.csv("methylation_data.csv", row.names = 1, header = TRUE)
expr_data <- read.csv("expression_data.csv", row.names = 1, header = TRUE)

# Check dimensions
dim(meth_data)   # Should return (genes, samples)
dim(expr_data)   # Should return (genes, samples)

# Find common sample names
common_samples <- intersect(colnames(meth_data), colnames(expr_data))
print(common_samples)  # Make sure the names match exactly

# Subset both datasets to have the same samples
expr_data <- expr_data[, common_samples]
meth_data <- meth_data[, common_samples]

# Convert to numeric (ensure data is numerical)
expr_data <- as.data.frame(lapply(expr_data, as.numeric))
meth_data <- as.data.frame(lapply(meth_data, as.numeric))

# Check new dimensions
dim(expr_data)
dim(meth_data)


# Convert to matrix format
expr_matrix <- as.matrix(expr_data)
meth_matrix <- as.matrix(meth_data)

# Check structure
dim(expr_matrix)
dim(meth_matrix)

# Remove genes/CpGs with near-zero variance
expr_matrix <- expr_matrix[apply(expr_matrix, 1, var) > 0, ]
meth_matrix <- meth_matrix[apply(meth_matrix, 1, var) > 0, ]

# Check new dimensions
dim(expr_matrix)
dim(meth_matrix)

combined_data <- rbind(expr_matrix, meth_matrix)

# Check final dataset
dim(combined_data)

library(WGCNA)

#Step 3: 

# Select soft-thresholding power
powers <- c(1:20)
sft <- pickSoftThreshold(expr_data, powerVector = powers, verbose = 5)

# Plot results
par(mfrow = c(1, 2))
plot(sft$fitIndices[, 1], -sign(sft$fitIndices[, 3]) * sft$fitIndices[, 2], 
     xlab = "Soft Threshold (power)", ylab = "Scale-Free Topology Model Fit", type = "b", col = "blue")
plot(sft$fitIndices[, 1], sft$fitIndices[, 5], 
     xlab = "Soft Threshold (power)", ylab = "Mean Connectivity", type = "b", col = "red")

# Choose power based on where R² > 0.9
softPower <- sft$powerEstimate
print(softPower)


#Step 4: Choose Soft-Thresholding Power

powers <- c(1:20)  # Range of power values
sft <- pickSoftThreshold(meth_data, powerVector = powers, verbose = 5)

# Plot scale-free topology fit
plot(sft$fitIndices[,1], -sign(sft$fitIndices[,3]) * sft$fitIndices[,2], type = "b",
     xlab = "Soft Threshold (power)", ylab = "Scale Free Topology Model Fit",
     main = "Choosing the Soft-Thresholding Power")

# Select power where model fit reaches 0.9
softPower <- sft$powerEstimate

#Step 5: Construct WGCNA Network

# Construct adjacency matrix
adjacency <- adjacency(expr_data, power = softPower)

expr_data_sub <- expr_data[,1:1000 ]  # Subset to the first 1000 genes, for example
adjacency <- adjacency(expr_data_sub, power = softPower)


# Convert adjacency matrix into a topological overlap matrix (TOM)
TOM <- TOMsimilarity(adjacency)
dissTOM <- 1 - TOM  # Dissimilarity measure

# Hierarchical clustering
geneTree <- hclust(as.dist(dissTOM), method = "average")

# Define modules using dynamic tree cut
dynamicMods <- cutreeDynamic(dendro = geneTree, distM = dissTOM,
                             deepSplit = 2, pamRespectsDendro = FALSE)

# Convert module labels into colors
moduleColors <- labels2colors(dynamicMods)

# Plot dendrogram with module colors
plotDendroAndColors(geneTree, moduleColors, "Dynamic Tree Cut",
                    main = "WGCNA Module Detection", dendroLabels = FALSE)

library(WGCNA)

# Step 1: Calculate the dissimilarity matrix and perform hierarchical clustering
dissTOM <- 1 - cor(t(expr_data))  # Calculate dissimilarity matrix (correlation-based)
geneTree <- hclust(as.dist(dissTOM), method = "average")  # Hierarchical clustering

# Step 2: Apply dynamic tree cutting with the hybrid method (providing the dissimilarity matrix)
dynamicTreeCut <- cutreeDynamic(dendro = geneTree, distM = dissTOM, method = "hybrid", deepSplit = 2)

# Check if modules were detected (should show the number of modules)
table(dynamicTreeCut)

# Step 3: Merge small modules based on minModuleSize (after dynamic tree cutting)
# Use the module labels (dynamicTreeCut$colors)
dynamicTreeCutMerged <- mergeCloseModules(expr_data, dynamicTreeCut, cutHeight = 0.25, minModuleSize = 30)

# Step 4: Convert labels to colors
moduleColors <- labels2colors(dynamicTreeCutMerged$colors)

# Step 5: Calculate module eigengenes
MEList <- moduleEigengenes(expr_data, colors = moduleColors)

#Step 6: Identify Disease-Associated Modules

# Compute module eigengenes (representative of each module)
MEList <- moduleEigengenes(expr_data, colors = moduleColors)
MEs <- MEList$eigengenes

# Correlation between module eigengenes and disease status
traitData <- as.numeric(factor(your_MS_status_vector))  # Replace with your disease labels (0 = control, 1 = MS)
moduleTraitCor <- cor(MEs, traitData, use = "p")

# Visualize correlation heatmap
pheatmap(moduleTraitCor, cluster_rows = FALSE, cluster_cols = FALSE,
         display_numbers = TRUE, main = "Module-Disease Correlation")

#Step 7: Integrate Methylation Data with WGCNA

# Compute module eigengenes for methylation data
MEList_meth <- moduleEigengenes(meth_matrix, colors = moduleColors)
MEs_meth <- MEList_meth$eigengenes

# Compute correlation between gene expression and methylation modules
corExprMeth <- cor(MEs, MEs_meth, use = "p")

# Plot correlation heatmap
pheatmap(corExprMeth, cluster_rows = TRUE, cluster_cols = TRUE, display_numbers = TRUE,
         main = "Correlation Between Expression and Methylation Modules")

#Step 8: Extract Key Hub Genes

# Identify hub genes in the most disease-associated module
selectedModule <- "turquoise"  # Change based on most correlated module from previous steps
moduleGenes <- colnames(expr_matrix)[moduleColors == selectedModule]

# Find top hub genes by connectivity
geneModuleMembership <- cor(expr_matrix[moduleGenes, ], MEs[selectedModule], use = "p")
hubGenes <- names(sort(geneModuleMembership, decreasing = TRUE))[1:10]  # Select top 10

print("Top Hub Genes:")
print(hubGenes)

#Step 9: Export Network for Visualization in Cytoscape

# Export network for Cytoscape
TOM_selected <- TOM[moduleColors == selectedModule, moduleColors == selectedModule]

# Convert to edge list format
edges <- which(TOM_selected > 0.1, arr.ind = TRUE)  # Adjust threshold as needed
edgeList <- data.frame(from = rownames(TOM_selected)[edges[,1]],
                       to = colnames(TOM_selected)[edges[,2]],
                       weight = TOM_selected[edges])

# Save to file
write.csv(edgeList, "WGCNA_network_for_Cytoscape.csv", row.names = FALSE)


#Step 10: Pathway Enrichment Analysis for Disease Modules

if (!require("clusterProfiler")) install.packages("clusterProfiler")
library(clusterProfiler)

# Perform GO enrichment analysis on selected module genes
ego <- enrichGO(gene = moduleGenes, OrgDb = org.Hs.eg.db, keyType = "SYMBOL",
                ont = "BP", pAdjustMethod = "BH", qvalueCutoff = 0.05)

# Plot results
dotplot(ego, showCategory = 10, title = "GO Enrichment of MS-Associated Module")

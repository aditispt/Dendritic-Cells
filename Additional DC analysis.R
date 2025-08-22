# Set the file paths for sample ##Ctrl-1
# Load the matrix file (expression matrix) using readMM, not load
matrix_file <- "~/Desktop/Mor Lab/DC analysis from other database/GSE183837_RAW (1)/GSM5572238_Ctrl-1.matrix.mtx"
matrix_data <- readMM(matrix_file)

# Load the gene file using fread
gene_file <- "~/Desktop/Mor Lab/DC analysis from other database/GSE183837_RAW (1)/GSM5572238_Ctrl-1.genes.tsv"
gene_data <- fread(gene_file, header = FALSE)  # header = FALSE assumes no column names in the file

rownames(matrix_data) <- make.unique(gene_data$V2)  # This will add a suffix (.1, .2, etc.) to duplicated gene names

## Create the Seurat object for Sample 1 without barcodes
seurat_object1 <- CreateSeuratObject(counts = matrix_data, project = "Ctrl-1")
#Downstream Seurat Processing
seurat_object1<-NormalizeData(seurat_object1)
seurat_object1<-FindVariableFeatures(seurat_object1)
seurat_object1 <- ScaleData(seurat_object1)
seurat_object1 <- RunPCA(seurat_object1)  

ElbowPlot(seurat_object1)

seurat_object1 <- FindNeighbors(seurat_object1, dims=1:15)
seurat_object1 <- FindClusters(seurat_object1, resolution = 0.5)
seurat_object1 <- RunUMAP(seurat_object1, dims=1:15)

DimPlot(seurat_object1, reduction = "umap", pt.size = 0.5, label = TRUE)
FeaturePlot(seurat_object1, features = c("LYZ", "IRF8"), pt.size = 0.5, label = TRUE)

##Extract the DC ##Cluster 7
DC_ctrl1 <- subset(seurat_object1, idents = "7")

# Set the file paths for sample ##Ctrl-2
# Load the matrix file (expression matrix) using readMM, not load
matrix_file <- "~/Desktop/Mor Lab/DC analysis from other database/GSE183837_RAW (1)/GSM5572239_Ctrl-2.matrix.mtx"
matrix_data <- readMM(matrix_file)

# Load the gene file using fread
gene_file <- "~/Desktop/Mor Lab/DC analysis from other database/GSE183837_RAW (1)/GSM5572239_Ctrl-2.genes.tsv"
gene_data <- fread(gene_file, header = FALSE)  # header = FALSE assumes no column names in the file

rownames(matrix_data) <- make.unique(gene_data$V2) # Assuming the gene names are in the first column (V1)

## Create the Seurat object for Sample 2 without barcodes
seurat_object2 <- CreateSeuratObject(counts = matrix_data, project = "Ctrl-2")

#Downstream Seurat Processing
seurat_object2<-NormalizeData(seurat_object2)
seurat_object2<-FindVariableFeatures(seurat_object2)
seurat_object2 <- ScaleData(seurat_object2)
seurat_object2 <- RunPCA(seurat_object2)  

ElbowPlot(seurat_object2)

seurat_object2 <- FindNeighbors(seurat_object2, dims=1:15)
seurat_object2 <- FindClusters(seurat_object2, resolution = 0.5)
seurat_object2 <- RunUMAP(seurat_object2, dims=1:15)

DimPlot(seurat_object2, reduction = "umap", pt.size = 0.5, label = TRUE)
FeaturePlot(seurat_object2, features = c("LYZ", "IRF8"), pt.size = 0.5, label = TRUE)

##Extract the DC ##Cluster 9
DC_ctrl2 <- subset(seurat_object2, idents = "9")


# Set the file paths for sample ##Ctrl-3
# Load the matrix file (expression matrix) using readMM, not load
matrix_file <- "~/Desktop/Mor Lab/DC analysis from other database/GSE183837_RAW (1)/GSM5572240_Ctrl-3.matrix.mtx"
matrix_data <- readMM(matrix_file)

# Load the gene file using fread
gene_file <- "~/Desktop/Mor Lab/DC analysis from other database/GSE183837_RAW (1)/GSM5572240_Ctrl-3.genes.tsv"
gene_data <- fread(gene_file, header = FALSE)  # header = FALSE assumes no column names in the file

rownames(matrix_data) <- make.unique(gene_data$V2)  # Assuming the gene names are in the first column (V1)

## Create the Seurat object for Sample 3 without barcodes
seurat_object3 <- CreateSeuratObject(counts = matrix_data, project = "Ctrl-3")

#Downstream Seurat Processing
seurat_object3<-NormalizeData(seurat_object3)
seurat_object3<-FindVariableFeatures(seurat_object3)
seurat_object3 <- ScaleData(seurat_object3)
seurat_object3 <- RunPCA(seurat_object3)  

ElbowPlot(seurat_object3)

seurat_object3 <- FindNeighbors(seurat_object3, dims=1:15)
seurat_object3 <- FindClusters(seurat_object3, resolution = 0.5)
seurat_object3 <- RunUMAP(seurat_object3, dims=1:15)

DimPlot(seurat_object3, reduction = "umap", pt.size = 0.5, label = TRUE)
FeaturePlot(seurat_object3, features = c("LYZ", "IRF8"), pt.size = 0.5, label = TRUE)

##Extract the DC ##Cluster 7
DC_ctrl3 <- subset(seurat_object3, idents = "7")

# Set the file paths for sample ##RIF-1
# Load the matrix file (expression matrix) using readMM, not load
matrix_file <- "~/Desktop/Mor Lab/DC analysis from other database/GSE183837_RAW (1)/GSM5572241_RIF-1.matrix.mtx"
matrix_data <- readMM(matrix_file)

# Load the gene file using fread
gene_file <- "~/Desktop/Mor Lab/DC analysis from other database/GSE183837_RAW (1)/GSM5572241_RIF-1.genes.tsv"
gene_data <- fread(gene_file, header = FALSE)  # header = FALSE assumes no column names in the file

rownames(matrix_data) <- make.unique(gene_data$V2)  # Assuming the gene names are in the first column (V1)

## Create the Seurat object for RIF-1 without barcodes
seurat_object4 <- CreateSeuratObject(counts = matrix_data, project = "RIF-1")

##Downstream Seurat Processing
seurat_object4<-NormalizeData(seurat_object4)
seurat_object4<-FindVariableFeatures(seurat_object4)
seurat_object4 <- ScaleData(seurat_object4)
seurat_object4 <- RunPCA(seurat_object4)  

ElbowPlot(seurat_object4)

seurat_object4 <- FindNeighbors(seurat_object4, dims=1:15)
seurat_object4 <- FindClusters(seurat_object4, resolution = 0.5)
seurat_object4 <- RunUMAP(seurat_object4, dims=1:15)

DimPlot(seurat_object4, reduction = "umap", pt.size = 0.5, label = TRUE)
FeaturePlot(seurat_object4, features = c("LYZ", "IRF8"), pt.size = 0.5, label = TRUE)

##Extract the DC ##Cluster 11
DC_rif1 <- subset(seurat_object4, idents = "11")

# Set the file paths for sample ##RIF-2
# Load the matrix file (expression matrix) using readMM, not load
matrix_file <- "~/Desktop/Mor Lab/DC analysis from other database/GSE183837_RAW (1)/GSM5572242_RIF-2.matrix.mtx"
matrix_data <- readMM(matrix_file)

# Load the gene file using fread
gene_file <- "~/Desktop/Mor Lab/DC analysis from other database/GSE183837_RAW (1)/GSM5572242_RIF-2.genes.tsv"
gene_data <- fread(gene_file, header = FALSE)  # header = FALSE assumes no column names in the file

rownames(matrix_data) <- make.unique(gene_data$V2)  # Assuming the gene names are in the first column (V1)

## Create the Seurat object for RIF-2 without barcodes
seurat_object5 <- CreateSeuratObject(counts = matrix_data, project = "RIF-2")

#Downstream Seurat Processing
seurat_object5<-NormalizeData(seurat_object5)
seurat_object5<-FindVariableFeatures(seurat_object5)
seurat_object5 <- ScaleData(seurat_object5)
seurat_object5 <- RunPCA(seurat_object5)  

ElbowPlot(seurat_object5)

seurat_object5 <- FindNeighbors(seurat_object5, dims=1:15)
seurat_object5 <- FindClusters(seurat_object5, resolution = 0.5)
seurat_object5 <- RunUMAP(seurat_object5, dims=1:15)

DimPlot(seurat_object5, reduction = "umap", pt.size = 0.5, label = TRUE)
FeaturePlot(seurat_object5, features = c("LYZ", "IRF8"), pt.size = 0.5, label = TRUE)

##Extract the DC ##Cluster 10
DC_rif2 <- subset(seurat_object5, idents = "10")

# Set the file paths for sample ##RIF-3
# Load the matrix file (expression matrix) using readMM, not load
matrix_file <- "~/Desktop/Mor Lab/DC analysis from other database/GSE183837_RAW (1)/GSM5572243_RIF-3.matrix.mtx"
matrix_data <- readMM(matrix_file)

# Load the gene file using fread
gene_file <- "~/Desktop/Mor Lab/DC analysis from other database/GSE183837_RAW (1)/GSM5572243_RIF-3.genes.tsv"
gene_data <- fread(gene_file, header = FALSE)  # header = FALSE assumes no column names in the file

rownames(matrix_data) <- make.unique(gene_data$V2)  # Assuming the gene names are in the first column (V1)

## Create the Seurat object for RIF-3 without barcodes
seurat_object6 <- CreateSeuratObject(counts = matrix_data, project = "RIF-3")

#Downstream Seurat Processing
seurat_object6<-NormalizeData(seurat_object6)
seurat_object6<-FindVariableFeatures(seurat_object6)
seurat_object6 <- ScaleData(seurat_object6)
seurat_object6 <- RunPCA(seurat_object6)  

ElbowPlot(seurat_object6)

seurat_object6 <- FindNeighbors(seurat_object6, dims=1:15)
seurat_object6 <- FindClusters(seurat_object6, resolution = 0.5)
seurat_object6 <- RunUMAP(seurat_object6, dims=1:15)

DimPlot(seurat_object6, reduction = "umap", pt.size = 0.5, label = TRUE)
FeaturePlot(seurat_object6, features = c("LYZ", "IRF8"), pt.size = 0.5, label = TRUE)

##Extract the DC ##Cluster 11
DC_rif3 <- subset(seurat_object6, idents = "11")


# Set the file paths for sample ##RIF-4
# Load the matrix file (expression matrix) using readMM, not load
matrix_file <- "~/Desktop/Mor Lab/DC analysis from other database/GSE183837_RAW (1)/GSM5572244_RIF-4.matrix.mtx"
matrix_data <- readMM(matrix_file)

# Load the gene file using fread
gene_file <- "~/Desktop/Mor Lab/DC analysis from other database/GSE183837_RAW (1)/GSM5572244_RIF-4.genes.tsv"
gene_data <- fread(gene_file, header = FALSE)  # header = FALSE assumes no column names in the file

rownames(matrix_data) <- make.unique(gene_data$V2)  # Assuming the gene names are in the first column (V1)

## Create the Seurat object for RIF-4 without barcodes
seurat_object7 <- CreateSeuratObject(counts = matrix_data, project = "RIF-4")

#Downstream Seurat Processing
seurat_object7<-NormalizeData(seurat_object7)
seurat_object7<-FindVariableFeatures(seurat_object7)
seurat_object7 <- ScaleData(seurat_object7)
seurat_object7 <- RunPCA(seurat_object7)  

ElbowPlot(seurat_object7)

seurat_object7 <- FindNeighbors(seurat_object7, dims=1:15)
seurat_object7 <- FindClusters(seurat_object7, resolution = 0.5)
seurat_object7 <- RunUMAP(seurat_object7, dims=1:15)

DimPlot(seurat_object7, reduction = "umap", pt.size = 0.5, label = TRUE)
FeaturePlot(seurat_object7, features = c("LYZ", "IRF8"), pt.size = 0.5, label = TRUE)

##Extract the DC ##Cluster 11
DC_rif4 <- subset(seurat_object7, idents = "11")

# Set the file paths for sample ##RIF-5
# Load the matrix file (expression matrix) using readMM, not load
matrix_file <- "~/Desktop/Mor Lab/DC analysis from other database/GSE183837_RAW (1)/GSM5572245_RIF-5.matrix.mtx"
matrix_data <- readMM(matrix_file)

# Load the gene file using fread
gene_file <- "~/Desktop/Mor Lab/DC analysis from other database/GSE183837_RAW (1)/GSM5572245_RIF-5.genes.tsv"
gene_data <- fread(gene_file, header = FALSE)  # header = FALSE assumes no column names in the file

rownames(matrix_data) <- make.unique(gene_data$V2)  # Assuming the gene names are in the first column (V1)

## Create the Seurat object for RIF-5 without barcodes
seurat_object8 <- CreateSeuratObject(counts = matrix_data, project = "RIF-5")

##Downstream Seurat Processing
seurat_object8<-NormalizeData(seurat_object8)
seurat_object8<-FindVariableFeatures(seurat_object8)
seurat_object8 <- ScaleData(seurat_object8)
seurat_object8 <- RunPCA(seurat_object8)  

ElbowPlot(seurat_object8)

seurat_object8 <- FindNeighbors(seurat_object8, dims=1:15)
seurat_object8 <- FindClusters(seurat_object8, resolution = 0.5)
seurat_object8 <- RunUMAP(seurat_object8, dims=1:15)

DimPlot(seurat_object8, reduction = "umap", pt.size = 0.5, label = TRUE)
FeaturePlot(seurat_object8, features = c("LYZ", "IRF8"), pt.size = 0.5, label = TRUE)

##Extract the DC ##Cluster 8
DC_rif5 <- subset(seurat_object8, idents = "8")

# Set the file paths for sample ##RIF-6
# Load the matrix file (expression matrix) using readMM, not load
matrix_file <- "~/Desktop/Mor Lab/DC analysis from other database/GSE183837_RAW (1)/GSM5572246_RIF-6.matrix.mtx"
matrix_data <- readMM(matrix_file)

# Load the gene file using fread
gene_file <- "~/Desktop/Mor Lab/DC analysis from other database/GSE183837_RAW (1)/GSM5572246_RIF-6.genes.tsv"
gene_data <- fread(gene_file, header = FALSE)  # header = FALSE assumes no column names in the file

rownames(matrix_data) <- make.unique(gene_data$V2)  # Assuming the gene names are in the first column (V1)

## Create the Seurat object for RIF-6 without barcodes
seurat_object9 <- CreateSeuratObject(counts = matrix_data, project = "RIF-6")

#Downstream Seurat Processing
seurat_object9<-NormalizeData(seurat_object9)
seurat_object9<-FindVariableFeatures(seurat_object9)
seurat_object9 <- ScaleData(seurat_object9)
seurat_object9 <- RunPCA(seurat_object9)  

ElbowPlot(seurat_object9)

seurat_object9 <- FindNeighbors(seurat_object9, dims=1:15)
seurat_object9 <- FindClusters(seurat_object9, resolution = 0.5)
seurat_object9 <- RunUMAP(seurat_object9, dims=1:15)

DimPlot(seurat_object9, reduction = "umap", pt.size = 0.5, label = TRUE)
FeaturePlot(seurat_object9, features = c("LYZ", "IRF8"), pt.size = 0.5, label = TRUE)

##Extract the DC ##Cluster 9
DC_rif6 <- subset(seurat_object9, idents = "9")

##End of creating separate seurat objects

##Create a big seurat with all DC clusters
##Merging again correctly
combined_DC <- merge(DC_ctrl1, 
                     y = c(DC_ctrl2, DC_ctrl3, DC_rif1, DC_rif2, DC_rif3, DC_rif4, DC_rif5, DC_rif6),
                     add.cell.ids = c("Ctrl1", "Ctrl2", "Ctrl3", "RIF1", "RIF2", "RIF3", "RIF4", "RIF5", "RIF6"),
                     project = "DC_Combined")

# Assign orig.ident based on the cell names (barcodes)
combined_DC$orig.ident <- sapply(strsplit(colnames(combined_DC), "_"), `[`, 1)

# Check the result
head(combined_DC@meta.data)

# Save the Seurat object as an RDS file
saveRDS(combined_DC, file = "combined_DC_all Dendritic Cells extratced after seurat processing of all 9 samples.rds")

# This will combine all your DC Seurat objects into one.
##QC
#Mitochondrial contamination

#Find the Mitochondrial PCT%information and plot a graph
combined_DC[["percent.mt"]]<- PercentageFeatureSet(combined_DC, pattern = "^MT-")

#Add the mitochondrial PCT% information to the Metadata
combined_DC$percent.mt <- PercentageFeatureSet(combined_DC, pattern = "^MT-")

#Make a Violin plot of the QC columns
VlnPlot(combined_DC, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)

#MT genes found in the given dataset
combined_DC <- subset(combined_DC, subset = nFeature_RNA > 200 & nFeature_RNA < 3000 & percent.mt < 40)

VlnPlot(combined_DC, features = c("nFeature_RNA", "percent.mt"), ncol = 2)
FeatureScatter(combined_DC, feature1 = "nCount_RNA", feature2 = "percent.mt")

# Normalization and scaling for the combined object
combined_DC <- NormalizeData(combined_DC)
combined_DC <- FindVariableFeatures(combined_DC)
combined_DC <- ScaleData(combined_DC)

# Run PCA and UMAP
combined_DC <- RunPCA(combined_DC)

ElbowPlot(combined_DC)

combined_DC <- FindNeighbors(combined_DC, dims=1:12)
combined_DC <- FindClusters(combined_DC, resolution = 0.5)
combined_DC <- RunUMAP(combined_DC, dims=1:12)

DimPlot(combined_DC, reduction = "umap", pt.size = 0.5, label = TRUE)

DimPlot(combined_DC, reduction = "umap", pt.size = 0.5, group.by = "orig.ident", label = TRUE)


FeaturePlot(combined_DC, features = c ("XCR1","CLEC9A"), pt.size = 0.5, label = TRUE)

##These were done on not filtered data (not filtered for 40% mito contamination)

##Cluster 11 showed marker genes for the XCR1 and CLEC9A.
cluster_11_cells <- subset(combined_DC, idents = "11")

# Extract metadata
cluster_11_metadata <- cluster_11_cells@meta.data

# Count the number of cells from each sample (orig.ident)
cell_counts <- table(cluster_11_metadata$orig.ident)

# Convert to a data frame for easy plotting
cell_counts_df <- as.data.frame(cell_counts)
colnames(cell_counts_df) <- c("Sample", "Cell_Count")

library(ggplot2)

# Plotting the histogram
ggplot(cell_counts_df, aes(x = Sample, y = Cell_Count)) +
  geom_bar(stat = "identity", fill = "skyblue") +
  theme_minimal() +
  labs(title = "Cell Counts in Cluster 11 by Sample",
       x = "Sample",
       y = "Number of Cells") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

##NUmber of cells from each samples for each cluster
# Extract metadata
metadata <- combined_DC@meta.data

# Make sure that the cluster identities are stored in a column named "seurat_clusters"
metadata$cluster <- Idents(combined_DC)

library(dplyr)

# Create a table with counts of cells for each cluster and sample
cell_counts_all <- metadata %>%
  group_by(cluster, orig.ident) %>%
  summarise(Cell_Count = n()) %>%
  as.data.frame()

library(ggplot2)

# Get unique cluster IDs (from 0 to 11 in this case)
unique_clusters <- unique(metadata$cluster)

# Loop over each cluster and plot the histogram
for (cluster_id in unique_clusters) {
  # Subset the data for the current cluster
  cluster_data <- cell_counts_all[cell_counts_all$cluster == cluster_id, ]
  
  # Plot the histogram
  p <- ggplot(cluster_data, aes(x = orig.ident, y = Cell_Count)) +
    geom_bar(stat = "identity", fill = "skyblue") +
    theme_minimal() +
    labs(title = paste("Cell Counts in Cluster", cluster_id, "by Sample"),
         x = "Sample",
         y = "Number of Cells") +
    theme(axis.text.x = element_text(angle = 45, hjust = 1))
  
  # Print the plot for each cluster
  print(p)
}

##DotPlot ##Too many markers, so reduced n=3
# Create an empty list to store marker dataframes for each cluster
cluster_markers_list <- list()

# Loop through clusters 0 to 11
for (cluster_id in 0:11) {
  cat("Finding markers for Cluster", cluster_id, "\n")
  
  # Find markers for the current cluster
  cluster_markers <- FindMarkers(combined_DC, ident.1 = cluster_id, only.pos = TRUE,
                                 logfc.threshold = 0.5, min.pct = 0.25)
  
  # Keep the top 5 markers based on p-value
  top_markers <- cluster_markers %>%
    top_n(5, -log10(p_val))
  
  # Store the top markers in the list
  cluster_markers_list[[paste0("Cluster", cluster_id, "_markers")]] <- top_markers
}

##Make Bubble Plots
bubble_features = rownames(bind_rows(cluster_markers_list, .id = "column_label"))

DotPlot(object = combined_DC, features = bubble_features) +
  theme(axis.text.x = element_text(angle = 90, hjust = 1))

#The layers are separated so need to run join layers
combined_DC <- JoinLayers(combined_DC)
DefaultAssay(combined_DC) <- "RNA"

# Create an empty list to store marker dataframes for each cluster
cluster_markers_list <- list()

# Loop through clusters 0 to 11
for (cluster_id in 0:11) {
  cat("Finding markers for Cluster", cluster_id, "\n")
  
  # Find markers for the current cluster
  cluster_markers <- FindMarkers(combined_DC, ident.1 = cluster_id, only.pos = TRUE,
                                 logfc.threshold = 0.5, min.pct = 0.25)
  
  # Keep the top 5 markers based on p-value
  top_markers <- cluster_markers %>%
    top_n(5, -log10(p_val))
  
  # Store the top markers in the list
  cluster_markers_list[[paste0("Cluster", cluster_id, "_markers")]] <- top_markers
}

##Make Bubble Plots
bubble_features = rownames(bind_rows(cluster_markers_list, .id = "column_label"))

DotPlot(object = combined_DC, features = bubble_features) +
  theme(axis.text.x = element_text(angle = 90, hjust = 1))



##A big seurat without the DC clusters extratced individually to create the combined object
##This Big seurat contains all cell types.
##Create a Big seurat with all seurats combined
# Merge all Seurat objects into a single object
DC_seurat <- merge(seurat_object1, 
                   y = list(seurat_object2, seurat_object3, seurat_object4, seurat_object5, 
                            seurat_object6, seurat_object7, seurat_object8, seurat_object9), 
                   add.cell.ids = c("Ctrl-1", "Ctrl-2", "Ctrl-3", "RIF-1", "RIF-2", "RIF-3", "RIF-4", "RIF-5", "RIF-6"), 
                   project = "Combined_Project")

# Check the combined Seurat object
print(DC_seurat)

#Find the Mitochondrial PCT%information and plot a graph
DC_seurat[["percent.mt"]]<- PercentageFeatureSet(DC_seurat, pattern = "^MT-")

#Add the mitochondrial PCT% information to the Metadata
DC_seurat$percent.mt <- PercentageFeatureSet(DC_seurat, pattern = "^MT-")

#Make a Violin plot of the QC columns
VlnPlot(DC_seurat, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)

#MT genes found in the given dataset
DC_seurat <- subset(DC_seurat, subset = nFeature_RNA > 200 & nFeature_RNA < 2500)

#Log-transform the counts
DC_seurat <- NormalizeData(DC_seurat)

#Find variable Features
DC_seurat <- FindVariableFeatures(DC_seurat)

#Scale the Data
DC_seurat <- ScaleData(DC_seurat)

#Run PCA
DC_seurat <- RunPCA(DC_seurat)

#Choose the number of principle components to keep
ElbowPlot(DC_seurat)

#Find nearest neighbors and construct the graph
DC_seurat <- FindNeighbors(DC_seurat, dims=1:12)

#Find the clusters
DC_seurat <- FindClusters(DC_seurat, resolution = 0.5)

#Get the UMAP embedding (Dims is the number of PCA plots you want the R to consider (PCA #number can be found from Elbow plot))
DC_seurat <- RunUMAP(DC_seurat, dims=1:12)

#Plot the UMAP with all the clusters
DimPlot(DC_seurat, reduction = "umap", label = TRUE)

FeaturePlot(DC_seurat, features = c("LYZ","IRF8","CD14"), label = TRUE)
FeaturePlot(DC_seurat, features = c("CD11C","ITGAX","CD11B"), label = TRUE)

combined_seurat <- DC_seurat

DimPlot(combined_seurat, reduction = "umap", label = TRUE)
DimPlot(combined_seurat, reduction = "umap", split.by = "active.ident", label = TRUE, repel = TRUE)


# Save the Seurat object as an RDS file
saveRDS(combined_seurat, file = "combined_seurat_all_endometrial_cells_processed.rds")

# Subsetting cluster 9
DC_seurat <- subset(DC_seurat, idents = 9)

#Log-transform the counts
DC_seurat <- NormalizeData(DC_seurat)

#Find variable Features
DC_seurat <- FindVariableFeatures(DC_seurat)

#Scale the Data
DC_seurat <- ScaleData(DC_seurat)

#Run PCA
DC_seurat <- RunPCA(DC_seurat)

#Choose the number of principle components to keep
ElbowPlot(DC_seurat)

#Find nearest neighbors and construct the graph
DC_seurat <- FindNeighbors(DC_seurat, dims=1:12)

#Find the clusters
DC_seurat <- FindClusters(DC_seurat, resolution = 0.5)

#Get the UMAP embedding (Dims is the number of PCA plots you want the R to consider (PCA #number can be found from Elbow plot))
DC_seurat <- RunUMAP(DC_seurat, dims=1:12)

#Plot the UMAP with all the clusters
DimPlot(DC_seurat, reduction = "umap", label = TRUE)

FeaturePlot(DC_seurat, features = c("LYZ","IRF8","CD14"), label = TRUE)

FeaturePlot(DC_seurat, features = c("CD11C","ITGAX","CD11B"), label = TRUE)

# Split UMAP by Ctrl and RIF using the "split.by" argument in DimPlot
DimPlot(DC_seurat, reduction = "umap", split.by = "active.ident", label = TRUE, repel = TRUE)

# Add active.ident to the metadata as a new column
DC_seurat$condition <- Idents(DC_seurat)

# Now split the UMAP by the new 'condition' column
DimPlot(DC_seurat, reduction = "umap", split.by = "condition", label = TRUE, repel = TRUE)


##Global normalization
# Normalize the combined data
DC_seurat <- NormalizeData(DC_seurat, normalization.method = "LogNormalize", scale.factor = 10000)

# Find variable features
DC_seurat <- FindVariableFeatures(DC_seurat, selection.method = "vst", nfeatures = 2000)

# Scale the data
DC_seurat <- ScaleData(DC_seurat)

# Perform PCA
DC_seurat <- RunPCA(DC_seurat, features = VariableFeatures(object = DC_seurat))

# Find neighbors and clusters
DC_seurat <- FindNeighbors(DC_seurat, dims = 1:20)
DC_seurat <- FindClusters(DC_seurat, resolution = 0.5)

# Run UMAP
DC_seurat <- RunUMAP(DC_seurat, dims = 1:20)

# Plot UMAP
DimPlot(DC_seurat, reduction = "umap", label = TRUE, repel = TRUE)


# --------------------------------------------------------------------------------------------------------
# --------------------------------------------------------------------------------------------------------
############ Single cell analysis ################
# --------------------------------------------------------------------------------------------------------
# --------------------------------------------------------------------------------------------------------

# USER-DEFINED OUTPUT FOLDER
output_dir <- "Path to file"

#load libraries
library(Seurat)
library(dplyr)
library(tidyr)
library(ggplot2)
library(stringr)
library(svglite)

## load single cell datasets and filter to retain only primary tumor samples ##

# Ebert
obj_Ebert <- readRDS("Path to file/obj_Ebert.rds")
Idents(obj_Ebert) <- "combined_clusters_annot"

# Abdelfattah
obj_Abdelfattah <- readRDS("Path to file/obj_Abdelfattah.rds")
table(obj_Abdelfattah$Patient)
Idents(obj_Abdelfattah) <- "Patient"
DimPlot(obj_Abdelfattah)
#remove recurrent samples
obj_Abdelfattah <- subset(obj_Abdelfattah, idents = c("ndGBM-01", "ndGBM-02", "ndGBM-04", "ndGBM-05", "ndGBM-06", "ndGBM-07", "ndGBM-08", "ndGBM-09", "ndGBM-10", "ndGBM-11"))
Idents(obj_Abdelfattah) <- "Assignment"

# Couturier 
obj_Couturier <- readRDS("Path to file/obj_Couturier.rds")
Idents(obj_Couturier) <- "cellkb"

# Nomura
obj_Nomura <- readRDS("Path to file/obj_Nomura.rds")
#remove recurrent samples
Idents(obj_Nomura) <- "PvsR"
obj_Nomura <- subset(obj_Nomura, idents = grep("Primary$", levels(Idents(obj_Nomura)), value = TRUE))
Idents(obj_Nomura) <- "CellType"

# --------------------------------------------------------------------------------------------------------
# Align author cell annotations by adding new metadata column "Pooled_Cell_Type" to harmonize nomenclature#
# --------------------------------------------------------------------------------------------------------

# List of Seurat objects and their corresponding cell type column names
seurat_list <- list(obj_Ebert, obj_Abdelfattah, obj_Couturier, obj_Nomura)
cell_type_columns <- c("combined_clusters_annot", "Assignment", "cellkb", "CellType")

# Define cell type mappings for each Seurat object
cell_type_mappings <- list(
  # Mapping for obj_Ebert
  c(
    "Endothelial cells_2" = "Endothelial",
    "Macrophages_5" = "Myeloid",
    "Pericyte cells_4" = "Pericyte",
    "T cells_1" = "Lymphocyte",
    "Tumor cells_0" = "Malignant",
    "Tumor/prolif. cells_3" = "Myeloid"
  ),
  # Mapping for obj_Abdelfattah
  c(
    "BCells" = "Lymphocyte",
    "Endo" = "Endothelial",
    "Glioma" = "Malignant",
    "Oligo" = "Oligodendrocyte",
    "Pericytes" = "Pericyte",
    "TCells" = "Lymphocyte",
    "Other" = "Other"
  ),
  # Mapping for obj_Couturier
  c(
    "Tumor" = "Malignant",
    "Myeloid" = "Myeloid",
    "Endothelial" = "Endothelial",
    "T cell" = "Lymphocyte"
  ),
  # Mapping for obj_Nomura
  c(
    "Endothel" = "Endothelial",
    "Excitatory neuron" = "Neuron",
    "Inhibitory neuron" = "Neuron",
    "Lymphocyte" = "Lymphocyte",
    "Oligodendrocyte" = "Oligodendrocyte",
    "OPC" = "OPC_non_malignant",
    "Pericyte" = "Pericyte",
    "Astrocyte" = "Astrocyte",
    "TAM" = "Myeloid",
    "Malignant" = "Malignant"
  )
)

# Loop through each Seurat object to add the new annotations
for (i in 1:length(seurat_list)) {
  # Get the current Seurat object, cell type column name, and mapping
  seurat_obj <- seurat_list[[i]]
  cell_type_column <- cell_type_columns[i]
  cell_type_mapping <- cell_type_mappings[[i]]
  
  # Recode the cell type annotations to the new pooled annotations
  seurat_obj@meta.data$Pooled_Cell_Type <- recode(
    seurat_obj@meta.data[[cell_type_column]], !!!cell_type_mapping
  )
  
  # Update the original Seurat object explicitly
  if (i == 1) {
    obj_Ebert <- seurat_obj
  } else if (i == 2) {
    obj_Abdelfattah <- seurat_obj
  } else if (i == 3) {
    obj_Couturier <- seurat_obj
  } else if (i == 4) {
    obj_Nomura <- seurat_obj
  }
}

# Example: View the updated metadata for obj_Ebert
table(obj_Ebert$Pooled_Cell_Type, useNA = "ifany")


# --------------------------------------------------------------------------------------------------------
# Heatmap of genes of interest across cell types #
# --------------------------------------------------------------------------------------------------------

# Genes of interest
genes_of_interest <- c("CD2", "CD3E", "CD3G", "LAG3", "TNFRSF9", "CD8B", "CD48",
  "CD247", "CD6", "LCK", "CD5", "TIGIT", "CD3D", "CD27", "CD70","ICOS","CD28", "CTLA4", "ZAP70",
  "GZMB", "PRF1", "BTLA", "IFNG", "CD40LG", "IL2", "CD160", "FASLG", "PDCD1"
)

# Function to calculate average expression and Z-score per dataset
process_dataset <- function(seurat_obj, dataset_name) {
  if (!"Pooled_celltype" %in% colnames(seurat_obj@meta.data)) {
    stop(paste("Metadata column 'Pooled_celltype' not found in dataset:", dataset_name))
  }
  genes_present <- intersect(genes_of_interest, rownames(seurat_obj))
  if (length(genes_present) == 0) {
    stop(paste("No genes of interest found in the dataset:", dataset_name))
  }
  assay_data <- GetAssayData(seurat_obj, layer = "data")
  assay_data <- assay_data[genes_present, , drop = FALSE]
  assay_data <- as.matrix(assay_data)
  expression_df <- as.data.frame(t(assay_data))
  expression_df$cell_type <- seurat_obj@meta.data$Pooled_celltype
  avg_expression <- expression_df %>%
    group_by(cell_type) %>%
    summarise_at(
      vars(all_of(genes_present)),
      mean,
      na.rm = TRUE
    ) %>%
    pivot_longer(
      -cell_type,
      names_to = "gene",
      values_to = "avg_expression"
    )
  plot_data <- avg_expression %>%
    group_by(gene) %>%
    mutate(z_score = as.numeric(scale(avg_expression))) %>%
    ungroup() %>%
    mutate(
      cell_type_dataset = paste(cell_type, dataset_name, sep = "_")
    )
  return(plot_data)
}

# Process each dataset
plot_data_Ebert <- process_dataset(obj_Ebert, "Ebert")
plot_data_Adelfattah <- process_dataset(obj_Abdelfattah, "Adelfattah")
plot_data_Couturier <- process_dataset(obj_Couturier, "Couturier")
plot_data_Nomura <- process_dataset(obj_Nomura, "Nomura")

# Combine all datasets
combined_plot_data <- bind_rows(
  plot_data_Couturier,
  plot_data_Adelfattah,
  plot_data_Ebert,
  plot_data_Nomura
)

# Specify gene order
custom_gene_order <- c(
  "CD2", "CD3E", "CD3G", "LAG3", "TNFRSF9", "CD8B", "CD48",
  "CD247", "CD6", "LCK", "CD5", "TIGIT", "CD3D", "CD27", "CD70","ICOS","CD28", "CTLA4", "ZAP70",
  "GZMB", "PRF1", "BTLA", "IFNG", "CD40LG", "IL2", "CD160", "FASLG", "PDCD1"
)

# Specify cell type order
custom_celltype_order <- c(
  "Lymphocyte",
  "Malignant",
  "Endothelial",
  "Pericyte",
  "Myeloid",
  "Oligodendrocyte",
  "Other"
)

# Set gene order
combined_plot_data$gene <- factor(
  combined_plot_data$gene,
  levels = custom_gene_order
)

# Set cell type/dataset order
combined_plot_data$cell_type_dataset <- factor(
  combined_plot_data$cell_type_dataset,
  levels = unlist(
    lapply(custom_celltype_order, function(cell_type) {
      grep(
        paste0("^", cell_type, "_"),
        unique(combined_plot_data$cell_type_dataset),
        value = TRUE
      )
    })
  )
)

# Exclude "Other"
library(stringr)

combined_plot_data_filtered <- combined_plot_data %>%
  filter(!str_detect(cell_type_dataset, "^Other_"))

# Plot heatmap
ggplot(
  combined_plot_data_filtered,
  aes(x = gene, y = cell_type_dataset, fill = z_score)
) +
  geom_tile() +
  scale_fill_gradient2(
    low = "blue",
    mid = "white",
    high = "red",
    midpoint = 0,
    limits = c(-2, 2),
    oob = scales::squish
  ) +
  labs(
    x = "Gene",
    y = "Cell Type (Dataset)",
    fill = "Z-score"
  ) +
  theme_minimal() +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1),
    panel.grid = element_blank()
  )

# Create heatmap-style matrix/dataframe:
# rows = celltype_dataset
# columns = genes
heatmap_csv <- combined_plot_data_filtered %>%
  select(cell_type_dataset, gene, z_score) %>%
  mutate(
    cell_type_dataset = factor(
      cell_type_dataset,
      levels = levels(combined_plot_data$cell_type_dataset)
    ),
    gene = factor(
      gene,
      levels = custom_gene_order
    )
  ) %>%
  arrange(cell_type_dataset, gene) %>%
  pivot_wider(
    names_from = gene,
    values_from = z_score
  ) %>%
  arrange(cell_type_dataset)

# Write CSV
write.csv(heatmap_csv, "Path to file/cell_type_gene_heatmap.csv", row.names=FALSE)

# --------------------------------------------------------------------------------------------------------
# Dotplot of genes of interest across Nomura metaprograms #
# --------------------------------------------------------------------------------------------------------

# Add Neftel annotations to Nomura metaprograms metadata
# To add malignant Neftel cell state annotations metadata column ("Malignant_subpop") see Github: https://github.com/chloe-shard/Houdova_et_al_2026 - 
# please use the following R script: Houdova_et_al_R_script_12_scRNAseq_dataset_Malignant_subpopulation_annotations

# Copy AssignedState (Nomura annotation of meta programs) into a new metadata column
obj_Nomura$AssignedState_Neftel <- as.character(obj_Nomura$AssignedState)

# Identify malignant cells
malignant_cells <- obj_Nomura$AssignedState == "Malignant"

# Replace Malignant with the corresponding Neftel cell state
obj_Nomura$AssignedState_Neftel[malignant_cells] <- paste0(
  "Malignant_",
  obj_Nomura$Malignant_subpop[malignant_cells]
)

# Convert back to factor
obj_Nomura$AssignedState_Neftel <- factor(obj_Nomura$AssignedState_Neftel)

# Check the result
table(obj_Nomura$AssignedState_Neftel)

# Generate dotplot
p4<- DotPlot(
  obj_Nomura,
  group.by = "AssignedState_Neftel",
  features = c(
    "CD2", "CD3E", "CD3G", "LAG3", "TNFRSF9", "CD8B", "CD48",
    "CD247", "CD6", "LCK", "CD5", "TIGIT", "CD3D", "CD27", "CD70","ICOS","CD28", "CTLA4", "ZAP70",
    "GZMB", "PRF1", "BTLA", "IFNG", "CD40LG", "IL2", "CD160", "FASLG", "PDCD1"
  )
) +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1)
  )

# Save as SVG
ggsave(
  filename = "Path to file/Nomura_subtype_dotplot.svg",
  plot = p4,
  device = svglite,
  width = 14,
  height = 10,
  units = "in"
)



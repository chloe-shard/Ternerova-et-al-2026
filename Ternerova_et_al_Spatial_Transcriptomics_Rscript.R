# --------------------------------------------------------------------------------------------------------
# --------------------------------------------------------------------------------------------------------
############# Spatial Analysis ####################
# --------------------------------------------------------------------------------------------------------
# --------------------------------------------------------------------------------------------------------

# USER-DEFINED OUTPUT FOLDER
output_dir <- "Path to file"

# Load libraries
library(Seurat)
library(dplyr)
library(tidyr)
library(ggplot2)
library(svglite)

# Load spatial tissues
path_to_files <- "Path to files"

# List of Seurat object names (Vessel enriched)
seurat_names <- c( "T243", "T248", "T255",
                   "T259", "T260", "T262",
                   "GBM_ZH881T1", "GBM_ZH916T1", "GBM_ZH1007nec",
                   "GBM_ZH8811Abulk", "GBM_ZH8812bulk",
                   "GBM_ZH1007inf", "GBM_ZH1019inf")

# Load each Seurat object if the file exists
for (seurat_name in seurat_names) {
  file_path <- file.path(path_to_files, paste0(seurat_name, ".rds"))
  if (file.exists(file_path)) {
    assign(seurat_name, readRDS(file_path))
    message(paste("Loaded:", seurat_name))
  } else {
    warning(paste("File not found:", file_path))
  }
}

# Note FAP+PC and broad cell types (e.g. Lymphocytes) were annotated using RCTD analysis and vessel spots ('Vascular_broad') 
# were annotated by summing pericyte and endothelial scores and thresholding (mean + 1SD) - for more information on this and 
# how to select vessel enriched tissues please see the following R scripts in Github: https://github.com/chloe-shard/Houdova_et_al_2026
# 1) Houdova_et_al_R_script_13_spatial_transcriptomics_RCTD_cell_type_annotation
# 2) Houdova_et_al_R_script_14_spatial_transcriptomics_spot_enrichment_annotation

# --------------------------------------------------------------------------------------------------------
# Annotate FAP+PC positive and FAP+PC negative perivascular regions (only in vessel enriched tissues)#
# --------------------------------------------------------------------------------------------------------

# create new metadata column called "Perivascular_annotation"
seurat_names <- c(
  "T243", "T248", "T255",
  "T259", "T260", "T262",
  "GBM_ZH881T1", "GBM_ZH916T1", "GBM_ZH1007nec",
  "GBM_ZH8811Abulk", "GBM_ZH8812bulk",
  "GBM_ZH1007inf", "GBM_ZH1019inf"
)

required_columns <- c(
  "FAP_PCs_thresholded",
  "Vascular_broad_thresholded")

# Function to annotate perivascular regions
for (i in seq_along(seurat_names)) {
  object_name <- seurat_names[i]
  message(
    "[", i, "/", length(seurat_names), "] Processing: ",
    object_name
  )
  if (!exists(object_name, envir = .GlobalEnv)) {
    warning("Object not found: ", object_name)
    next
  }
  seurat_obj <- get(object_name, envir = .GlobalEnv)
  missing_columns <- setdiff(
    required_columns,
    colnames(seurat_obj@meta.data)
  )
  
  if (length(missing_columns) > 0) {
    warning(
      object_name,
      " is missing: ",
      paste(missing_columns, collapse = ", ")
    )
    next
  }
  # Create the new metadata annotation
  seurat_obj$Perivascular_annotation <- case_when(
    
    seurat_obj$FAP_PCs_thresholded == "positive" &
      seurat_obj$Vascular_broad_thresholded == "positive" ~
      "FAP_PC_pos_ves",
    
    seurat_obj$FAP_PCs_thresholded == "negative" &
      seurat_obj$Vascular_broad_thresholded == "positive" ~
      "FAP_PC_neg_ves",
    
    seurat_obj$FAP_PCs_thresholded == "positive" &
      seurat_obj$Vascular_broad_thresholded == "negative" ~
      "FAP_PC_pos_only",
    
    TRUE ~ "negative"
  )
  seurat_obj$Perivascular_annotation <- factor(
    seurat_obj$Perivascular_annotation,
    levels = c(
      "negative",
      "FAP_PC_pos_only",
      "FAP_PC_neg_ves",
      "FAP_PC_pos_ves"
    )
  )
  assign(
    object_name,
    seurat_obj,
    envir = .GlobalEnv
  )
  message("Completed: ", object_name)
}

#check
p2<- SpatialDimPlot(GBM_ZH8811Abulk, group.by = ("Perivascular_annotation"), pt.size.factor = 3)

# Generate spatial plot of Lymphocytes
p1<- SpatialFeaturePlot(
  GBM_ZH8811Abulk,
  features = "Lymphocyte",
  pt.size.factor = 3,
  max.cutoff = 0.04
)

p3 <- p1+p2

# Save as SVG
ggsave(
  filename = "Path to file/Lymphocyte_spatial_plot.svg",
  plot = p3,
  device = svglite,
  width = 12,
  height = 6,
  units = "in"
)

# ----------------------------------------------------------------------------------------
# Quantify expression of genes in local neighbourhood of FAP_PC_pos_ves vs FAP_PC_neg_ves #
# ----------------------------------------------------------------------------------------

# Genes of interest
genes_of_interest <- c( "CD2", "CD3E", "CD3G", "LAG3", "TNFRSF9", "CD8B", "CD48",
                        "CD247", "CD6", "LCK", "CD5", "TIGIT", "CD3D", "CD27", "CD70","ICOS","CD28", "CTLA4", "ZAP70",
                        "GZMB", "PRF1", "BTLA", "IFNG", "CD40LG", "IL2", "CD160", "FASLG", "PDCD1")

# Metadata column and annotations to compare
metadata_column <- "Perivascular_annotation"

group1 <- "FAP_PC_pos_ves"
group2 <- "FAP_PC_neg_ves"

# List of vessel-enriched Seurat objects
seurat_names <- c(
  "T243", "T248", "T255",
  "T259", "T260", "T262",
  "GBM_ZH881T1", "GBM_ZH916T1", "GBM_ZH1007nec",
  "GBM_ZH8811Abulk", "GBM_ZH8812bulk",
  "GBM_ZH1007inf", "GBM_ZH1019inf"
)

# list to store results
all_results <- list()

# Calculate average neighbouring-spot expression
calculate_fold_change_expression <- function(
    seurat_obj,
    obj_name,
    metadata_column,
    group1,
    group2,
    gene,
    number_of_spots = 7
) {
  if (!metadata_column %in% colnames(seurat_obj@meta.data)) {
    stop(
      paste0(
        "Metadata column '",
        metadata_column,
        "' was not found."
      )
    )
  }
  if (!gene %in% rownames(seurat_obj)) {
    stop(
      paste0(
        "Gene '",
        gene,
        "' was not found in the Seurat object."
      )
    )
  }
  available_annotations <- unique(
    as.character(seurat_obj@meta.data[[metadata_column]])
  )
  missing_groups <- setdiff(
    c(group1, group2),
    available_annotations
  )
  if (length(missing_groups) > 0) {
    stop(
      paste0(
        "The following annotations were not found in ",
        metadata_column,
        ": ",
        paste(missing_groups, collapse = ", ")
      )
    )
  }
  gene_data <- FetchData(
    object = seurat_obj,
    vars = gene,
    layer = "data"
  )
  coords <- GetTissueCoordinates(seurat_obj)
  common_spots <- Reduce(
    intersect,
    list(
      rownames(gene_data),
      rownames(coords),
      rownames(seurat_obj@meta.data)
    )
  )
  if (length(common_spots) == 0) {
    stop("No matching spots were found between expression data and coordinates.")
  }
  gene_data <- gene_data[common_spots, , drop = FALSE]
  coords <- coords[common_spots, , drop = FALSE]
  metadata <- seurat_obj@meta.data[common_spots, , drop = FALSE]
  # Calculate mean expression across all spots, including zeros
  gene_all_spot_mean <- mean(
    gene_data[[gene]],
    na.rm = TRUE
  )
  if (is.na(gene_all_spot_mean) || gene_all_spot_mean == 0) {
    stop(
      paste0(
        "Gene '",
        gene,
        "' has a mean expression of zero across all spots in ",
        obj_name,
        "."
      )
    )
  }
  
  # Calculate neighbourhood expression for one annotation
  get_avg_expression <- function(annotation) {
    annotation_spots <- rownames(metadata)[
      !is.na(metadata[[metadata_column]]) &
        metadata[[metadata_column]] == annotation
    ]
    
    if (length(annotation_spots) == 0) {
      return(
        list(
          mean_expression = NA_real_,
          spot_expressions = numeric(0),
          number_of_spots = 0
        )
      )
    }
    avg_expressions <- numeric(length(annotation_spots))
    for (i in seq_along(annotation_spots)) {
      spot <- annotation_spots[i]
      spot_coords <- as.numeric(
        coords[spot, , drop = TRUE]
      )
      distances <- sqrt(
        rowSums(
          sweep(
            as.matrix(coords),
            MARGIN = 2,
            STATS = spot_coords,
            FUN = "-"
          )^2
        )
      )
      number_to_select <- min(
        number_of_spots,
        length(distances)
      )
      closest_neighbors <- names(
        sort(distances, decreasing = FALSE)
      )[seq_len(number_to_select)]
      neighbor_gene_expression <- gene_data[
        closest_neighbors,
        gene,
        drop = TRUE
      ]
      avg_expressions[i] <- mean(
        neighbor_gene_expression,
        na.rm = TRUE
      )
    }
    list(
      mean_expression = mean(
        avg_expressions,
        na.rm = TRUE
      ),
      spot_expressions = avg_expressions,
      number_of_spots = length(annotation_spots)
    )
  }
  
  
  # Calculate expression around each annotation
  group1_results <- get_avg_expression(group1)
  group2_results <- get_avg_expression(group2)
  
  avg_expression_group1 <- group1_results$mean_expression
  avg_expression_group2 <- group2_results$mean_expression
  
  # Normalise to the mean expression among gene-positive spots
  normalized_expression_group1 <-
    avg_expression_group1 / gene_all_spot_mean
  
  normalized_expression_group2 <-
    avg_expression_group2 / gene_all_spot_mean
  
  
  # Return one row per tissue and gene
  result_df <- data.frame(
    Object = obj_name,
    Gene = gene,
    Group1 = group1,
    Group2 = group2,
    Number_FAP_PC_pos_ves_Spots =
      group1_results$number_of_spots,
    Number_FAP_PC_neg_ves_Spots =
      group2_results$number_of_spots,
    Avg_Expression_FAP_PC_pos_ves =
      avg_expression_group1,
    Avg_Expression_FAP_PC_neg_ves =
      avg_expression_group2,
    Normalized_Expression_FAP_PC_pos_ves =
      normalized_expression_group1,
    Normalized_Expression_FAP_PC_neg_ves =
      normalized_expression_group2,
    stringsAsFactors = FALSE
  )
  
  return(result_df)
}


# Process each Seurat object and gene
for (obj_name in seurat_names) {
  
  cat("Processing", obj_name, "\n")
  
  if (!exists(obj_name, envir = .GlobalEnv)) {
    message("Object not found: ", obj_name)
    next
  }
  
  seurat_obj <- get(
    obj_name,
    envir = .GlobalEnv
  )
  
  for (gene in genes_of_interest) {
    
    tryCatch({
      
      result <- calculate_fold_change_expression(
        seurat_obj = seurat_obj,
        obj_name = obj_name,
        metadata_column = metadata_column,
        group1 = group1,
        group2 = group2,
        gene = gene,
        number_of_spots = 7
      )
      
      all_results[[paste(obj_name, gene, sep = "_")]] <- result
      
    }, error = function(e) {
      
      message(
        "Error processing object ",
        obj_name,
        " and gene ",
        gene,
        ": ",
        e$message
      )
      
    })
  }
}


# Combine all tissue-level results
if (length(all_results) == 0) {
  stop("No results were generated.")
}

merged_results <- bind_rows(all_results)

# Combine all tissue-level results
if (length(all_results) == 0) {
  stop("No results were generated.")
}

merged_results <- bind_rows(all_results)


# Perform one paired two-tailed t-test per gene across matched tissues
statistical_results <- lapply(
  genes_of_interest,
  function(gene) {
    
    gene_data <- merged_results %>%
      filter(Gene == gene) %>%
      filter(
        is.finite(Avg_Expression_FAP_PC_pos_ves),
        is.finite(Avg_Expression_FAP_PC_neg_ves)
      )
    
    # Require at least 2 matched tissues
    if (nrow(gene_data) < 2) {
      return(
        data.frame(
          Gene = gene,
          Number_of_Tissues = nrow(gene_data),
          Mean_Expression_FAP_PC_pos_ves = mean(
            gene_data$Avg_Expression_FAP_PC_pos_ves,
            na.rm = TRUE
          ),
          Mean_Expression_FAP_PC_neg_ves = mean(
            gene_data$Avg_Expression_FAP_PC_neg_ves,
            na.rm = TRUE
          ),
          P_Value = NA_real_
        )
      )
    }
    
    # Paired two-tailed t-test
    test_result <- tryCatch(
      {
        t.test(
          gene_data$Avg_Expression_FAP_PC_pos_ves,
          gene_data$Avg_Expression_FAP_PC_neg_ves,
          paired = TRUE,
          alternative = "two.sided"
        )
      },
      error = function(e) NULL
    )
    
    data.frame(
      Gene = gene,
      Number_of_Tissues = nrow(gene_data),
      
      Mean_Expression_FAP_PC_pos_ves = mean(
        gene_data$Avg_Expression_FAP_PC_pos_ves,
        na.rm = TRUE
      ),
      
      Mean_Expression_FAP_PC_neg_ves = mean(
        gene_data$Avg_Expression_FAP_PC_neg_ves,
        na.rm = TRUE
      ),
      
      P_Value = if (is.null(test_result)) {
        NA_real_
      } else {
        test_result$p.value
      }
    )
  }
) %>%
  bind_rows() %>%
  mutate(
    Adjusted_P_Value = p.adjust(
      P_Value,
      method = "BH"
    )
  )


# Keep genes significant by p-value < 0.05
significant_genes <- statistical_results %>%
  filter(
    !is.na(P_Value),
    P_Value < 0.05
  ) %>%
  arrange(P_Value)


# Print significant genes
print(significant_genes)

# Export tissue-level expression results
write.csv(
  merged_results,
  file = file.path(
    output_dir,
    "final_results_gene_expression_FAP_PC_pos_vs_neg_vessels.csv"
  ),
  row.names = FALSE
)

# Export gene-level statistical results
write.csv(
  statistical_results,
  file = file.path(
    output_dir,
    "statistical_results_FAP_PC_pos_vs_neg_vessels.csv"
  ),
  row.names = FALSE
)


# ------------------------------------------------------------
# Create boxplots of normalised expression
# ------------------------------------------------------------

# Select genes to plot
genes_to_plot <- c(
  "CD2",
  "CD3E",
  "CD3G",
  "LAG3",
  "TNFRSF9",
  "CD8B",
  "CD48"
)


# Convert normalised-expression columns to long format
normalized_boxplot_data <- merged_results %>%
  filter(Gene %in% genes_to_plot) %>%
  select(
    Object,
    Gene,
    Normalized_Expression_FAP_PC_pos_ves,
    Normalized_Expression_FAP_PC_neg_ves
  ) %>%
  pivot_longer(
    cols = c(
      Normalized_Expression_FAP_PC_pos_ves,
      Normalized_Expression_FAP_PC_neg_ves
    ),
    names_to = "Annotation",
    values_to = "Normalized_Expression"
  ) %>%
  mutate(
    Gene = factor(
      Gene,
      levels = genes_to_plot
    ),
    Annotation = recode(
      Annotation,
      Normalized_Expression_FAP_PC_pos_ves = "FAP_PC_pos_ves",
      Normalized_Expression_FAP_PC_neg_ves = "FAP_PC_neg_ves"
    )
  )


# Plot boxplots
normalized_boxplot <- ggplot(
  normalized_boxplot_data,
  aes(
    x = Gene,
    y = Normalized_Expression,
    fill = Annotation
  )
) +
  geom_hline(
    yintercept = 1,
    linetype = "dashed"
  ) +
  geom_boxplot(
    position = position_dodge(width = 0.8),
    outlier.shape = NA
  ) +
  geom_point(
    aes(group = Annotation),
    position = position_jitterdodge(
      jitter.width = 0.15,
      dodge.width = 0.8
    ),
    size = 1.5,
    alpha = 0.6
  ) +
  labs(
    title = "Normalised gene expression across tissues",
    subtitle = "Each point represents one tissue",
    x = "Gene",
    y = "Normalised expression",
    fill = "Annotation"
  ) +
  theme_minimal() +
  theme(
    axis.text.x = element_text(
      angle = 45,
      hjust = 1
    )
  )

normalized_boxplot

# Save as an SVG file
ggsave(
  filename = file.path(
    output_dir,
    "normalized_expression_boxplot_FAP_PC_pos_ves_vs_neg.svg"
  ),
  plot = normalized_boxplot,
  device = svglite::svglite,
  width = 10,
  height = 6,
  units = "in"
)


# -----------------------------------------------------------------------------------------
# Lymphocyte neighbourhood composition score analysis
# -----------------------------------------------------------------------------------------

# Seurat object names
seurat_names <- c(
  "T243", "T248", "T255",
  "T259", "T260", "T262",
  "GBM_ZH881T1", "GBM_ZH916T1", "GBM_ZH1007nec",
  "GBM_ZH8811Abulk", "GBM_ZH8812bulk",
  "GBM_ZH1007inf", "GBM_ZH1019inf"
)

# Metadata column containing the perivascular annotations
perivascular_metadata <- "Perivascular_annotation"

# Perivascular annotations to compare
perivascular_annotations <- c(
  "FAP_PC_pos_ves",
  "FAP_PC_neg_ves"
)

# Neighbor cell types to analyse
neighbor_thresholded_cell_types <- c(
  "Lymphocyte_thresholded"
)

# Initialise a dataframe to store the results
results <- data.frame(
  SeuratObject = character(),
  PerivascularAnnotation = character(),
  NeighborCellType = character(),
  AverageFractionPositive = numeric(),
  stringsAsFactors = FALSE
)

# Function to process each Seurat object
process_spatial_object <- function(
    spatial_obj,
    perivascular_metadata,
    perivascular_annotation,
    neighbor_cell_type
) {
  
  if (!perivascular_metadata %in% colnames(spatial_obj@meta.data)) {
    warning(
      paste(
        "Annotation column",
        perivascular_metadata,
        "not found in object:",
        deparse(substitute(spatial_obj))
      )
    )
    return(NA)
  }
  
  if (!neighbor_cell_type %in% colnames(spatial_obj@meta.data)) {
    warning(
      paste(
        "Neighbour annotation",
        neighbor_cell_type,
        "not found in object:",
        deparse(substitute(spatial_obj))
      )
    )
    return(NA)
  }
  
  coords <- GetTissueCoordinates(spatial_obj)
  
  # Select spots matching the chosen Perivascular_annotation category
  positive_spots <- WhichCells(
    spatial_obj,
    cells = rownames(
      spatial_obj@meta.data[
        spatial_obj@meta.data[[perivascular_metadata]] ==
          perivascular_annotation,
        ,
        drop = FALSE
      ]
    )
  )
  
  if (length(positive_spots) == 0) {
    warning(
      paste(
        "No spots annotated as",
        perivascular_annotation,
        "found in object:",
        deparse(substitute(spatial_obj))
      )
    )
    return(NA)
  }
  
  distance_matrices <- list()
  
  calculate_distances <- function(spot_coords, all_coords) {
    
    spot_vector <- as.numeric(spot_coords)
    
    distances <- sqrt(
      rowSums(
        (t(t(all_coords) - spot_vector))^2
      )
    )
    
    return(distances)
  }
  
  for (spot in positive_spots) {
    
    spot_coords <- coords[
      spot,
      ,
      drop = FALSE
    ]
    
    distances <- calculate_distances(
      spot_coords,
      coords
    )
    
    distance_df <- data.frame(
      SpotID = rownames(coords),
      Distance = distances
    )
    
    distance_df <- distance_df[
      order(distance_df$Distance),
    ]
    
    distance_matrices[[spot]] <- distance_df
  }
  
  neighbor_matrix <- matrix(
    NA,
    nrow = length(positive_spots),
    ncol = 7
  )
  
  rownames(neighbor_matrix) <- positive_spots
  
  colnames(neighbor_matrix) <- paste0(
    "Neighbor_",
    1:7
  )
  
  for (i in seq_along(positive_spots)) {
    
    spot <- positive_spots[i]
    distance_df <- distance_matrices[[spot]]
    
    neighbor_matrix[i, 1:7] <-
      distance_df$SpotID[1:7]
  }
  
  neighbor_matrix_annotated <- neighbor_matrix
  
  for (i in seq_len(nrow(neighbor_matrix))) {
    
    for (j in seq_len(ncol(neighbor_matrix))) {
      
      spot_id <- neighbor_matrix[i, j]
      
      if (!is.na(spot_id)) {
        
        neighbor_matrix_annotated[i, j] <-
          spatial_obj@meta.data[
            spot_id,
            neighbor_cell_type
          ]
      }
    }
  }
  
  positive_fractions <- apply(
    neighbor_matrix_annotated,
    1,
    function(row) {
      
      positive_count <- sum(
        row == "positive",
        na.rm = TRUE
      )
      
      total_count <- sum(
        !is.na(row)
      )
      
      fraction_positive <-
        positive_count / total_count
      
      return(fraction_positive)
    }
  )
  
  average_fraction_positive <- mean(
    positive_fractions,
    na.rm = TRUE
  )
  
  return(average_fraction_positive)
}

# Process each Seurat object, perivascular annotation and neighbour cell type
for (obj_name in seurat_names) {
  
  spatial_obj <- get(
    obj_name,
    envir = .GlobalEnv
  )
  
  for (perivascular_annotation in perivascular_annotations) {
    
    for (neighbor_cell_type in neighbor_thresholded_cell_types) {
      
      if (
        perivascular_metadata %in%
        colnames(spatial_obj@meta.data)
      ) {
        
        if (
          neighbor_cell_type %in%
          colnames(spatial_obj@meta.data)
        ) {
          
          avg_fraction_positive <-
            process_spatial_object(
              spatial_obj = spatial_obj,
              perivascular_metadata =
                perivascular_metadata,
              perivascular_annotation =
                perivascular_annotation,
              neighbor_cell_type =
                neighbor_cell_type
            )
          
          if (!is.na(avg_fraction_positive)) {
            
            results <- rbind(
              results,
              data.frame(
                SeuratObject = obj_name,
                PerivascularAnnotation =
                  perivascular_annotation,
                NeighborCellType =
                  neighbor_cell_type,
                AverageFractionPositive =
                  avg_fraction_positive,
                stringsAsFactors = FALSE
              )
            )
          }
          
        } else {
          
          warning(
            paste(
              "Neighbour annotation",
              neighbor_cell_type,
              "not found in object:",
              obj_name
            )
          )
        }
        
      } else {
        
        warning(
          paste(
            "Annotation column",
            perivascular_metadata,
            "not found in object:",
            obj_name
          )
        )
      }
    }
  }
}

# Order results
results <- results[
  order(
    results$PerivascularAnnotation,
    -results$AverageFractionPositive
  ),
]

# Print results
print(results)

# Clean up names for plotting
results_plot <- results %>%
  mutate(
    NeighbourPopulation = gsub(
      "_thresholded",
      "",
      NeighborCellType
    ),
    PerivascularAnnotation = factor(
      PerivascularAnnotation,
      levels = perivascular_annotations
    )
  )

# Calculate mean neighbourhood composition score for each group
mean_results <- results_plot %>%
  group_by(
    PerivascularAnnotation,
    NeighbourPopulation
  ) %>%
  summarise(
    Mean_NCS = mean(
      AverageFractionPositive,
      na.rm = TRUE
    ),
    .groups = "drop"
  )

# Plot both perivascular annotations as boxplots
neighbourhood_plot <- ggplot(
  results_plot,
  aes(
    x = PerivascularAnnotation,
    y = AverageFractionPositive
  )
) +
  geom_boxplot(
    outlier.shape = NA
  ) +
  geom_jitter(
    width = 0.15,
    size = 3,
    alpha = 0.75,
    aes(colour = SeuratObject)
  ) +
  facet_wrap(
    ~NeighbourPopulation
  ) +
  labs(
    x = "Perivascular annotation",
    y = "Lymphocyte neighbourhood composition score",
    colour = "Tissue",
    title = paste0(
      "Lymphocyte neighbourhood composition around ",
      "FAP_PC_pos_ves and FAP_PC_neg_ves spots"
    )
  ) +
  theme_classic(base_size = 14) +
  theme(
    axis.text.x = element_text(
      angle = 45,
      hjust = 1
    ),
    legend.position = "right"
  )

print(neighbourhood_plot)

# Reshape so each tissue has matched values for both annotations
paired_results <- results_plot %>%
  select(
    SeuratObject,
    NeighbourPopulation,
    PerivascularAnnotation,
    AverageFractionPositive
  ) %>%
  pivot_wider(
    names_from = PerivascularAnnotation,
    values_from = AverageFractionPositive
  ) %>%
  filter(
    !is.na(FAP_PC_pos_ves),
    !is.na(FAP_PC_neg_ves)
  )

# Two-tailed paired t-test
paired_t_test <- t.test(
  paired_results$FAP_PC_pos_ves,
  paired_results$FAP_PC_neg_ves,
  paired = TRUE,
  alternative = "two.sided"
)

print(paired_t_test)

# Save as an SVG file
ggsave(
  filename = file.path(
    output_dir,
    "NCS_lymphocyte_FAP_PC_pos_ves_vs_neg.svg"
  ),
  plot = neighbourhood_plot,
  device = svglite::svglite,
  width = 6,
  height = 10,
  units = "in"
)

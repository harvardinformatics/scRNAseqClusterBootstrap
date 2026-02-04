args <- commandArgs(trailingOnly = TRUE)
seurat_rds <- args[1]
output <- args[2]
library("tidyverse")
library("Seurat")
library("glmGamPoi")
options(future.globals.maxSize = 2 * 1024^3) 

## define functions

JaccardSimilarity <- function(set1, set2) {
  intersect_length <- length(intersect(set1, set2))
  union_length <- length(set1) + length(set2) - intersect_length
  intersect_length / union_length
}


RandomSubsetData<- function(object, rate, random.subset.seed = NULL, ...){
        ncells<- nrow(object@meta.data)
        ncells.subsample<- round(ncells * rate)

        set.seed(random.subset.seed)

        selected.cells<- sample(colnames(object), ncells.subsample)
        object<- subset(object, cells =  selected.cells,
                            ...)
        return(object)
}

SubSampleReSCTSeuratObject <- function(seurat_obj,subrate){
  subsampled_seurat <- RandomSubsetData(seurat_obj, rate = subrate)
  subsampled_seurat$presub_clusters <- seurat_obj$seurat_clusters
  subsampled_seurat <- SCTransform(subsampled_seurat, vars.to.regress = "percent.mt", verbose = FALSE)
  subsampled_seurat <-RunPCA(subsampled_seurat, verbose = FALSE)
  subsampled_seurat <- RunUMAP(subsampled_seurat, dims = 1:30)
  subsampled_seurat <- FindNeighbors(subsampled_seurat, dims = 1:30)
  subsampled_seurat <- FindClusters(subsampled_seurat)
  return(subsampled_seurat)
}


GetJaccardMaxByCluster <- function(seurat_obj,bootstrap) {
  jaccard_max_stats <- tibble::tibble(
    clusterid   = factor(),
    max_jaccard = numeric(),
    bootstrap_number = integer()
  )

  for (original_cluster in unique(seurat_obj$presub_clusters)) {

    barcodes <- rownames(
      subset(seurat_obj@meta.data,presub_clusters == original_cluster))

    dat <- tibble::tibble(
      cell_id = names(seurat_obj@active.ident),
      cluster = seurat_obj$seurat_clusters
    ) %>%
      tidyr::nest(data = -cluster) %>%
      dplyr::arrange(cluster)

    maxstat <- dat %>%
      dplyr::mutate(
        jaccard = purrr::map(data, ~ JaccardSimilarity(barcodes, .x$cell_id))
      ) %>%
      dplyr::pull(jaccard) %>%
      unlist() %>%
      max()

    jaccard_max_stats <- jaccard_max_stats %>%
      tibble::add_row(clusterid = original_cluster, max_jaccard = maxstat,bootstrap_number=bootstrap)
  }

  return(jaccard_max_stats)
}


## load data
seurat_obj <-readRDS(seurat_rds)

## do bootstrapping analysis
jaccard_max_stats <-NULL
jaccard_max_stats <- tibble::tibble(
    clusterid   = factor(),
    max_jaccard = numeric(),
    bootstrap_number = integer())



for (replicate in seq(1,100)){
  subsampled_obj <- SubSampleReSCTSeuratObject(seurat_obj,0.8)
  replicate_stats <- GetJaccardMaxByCluster(subsampled_obj,replicate)
  jaccard_max_stats <-bind_rows(jaccard_max_stats,replicate_stats)
}

write_tsv(jaccard_max_stats, output)

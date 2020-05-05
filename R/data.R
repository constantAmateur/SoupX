#' PBMC 4K meta data
#'
#' Collection of bits of meta data relating to the 10X PBMC 4K data.
#' 
#' This data set pertains to the 10X demonstration PBMC 4K data and includes metadata about it in the \code{data.frame} named \code{PBMC_DR}.
#'
#' \code{PBMC_DR} was created using Seurat (v2) to calculate a tSNE representation of the data and cluster cells with these commands.
#' \itemize{
#'   \item \code{set.seed(1)}
#'   \item \code{srat = CreateSeuratObject(sc$toc)}
#'   \item \code{srat = NormalizeData(srat)}
#'   \item \code{srat = ScaleData(srat)}
#'   \item \code{srat = FindVariableGenes(srat)}
#'   \item \code{srat = RunPCA(srat,pcs.compute=30)}
#'   \item \code{srat = RunTSNE(srat,dims.use=seq(30))}
#'   \item \code{srat = FindClusters(srat,dims.use=seq(30),resolution=1)}
#'   \item \code{PBMC_DR = as.data.frame(srat@dr$tsne@cell.embeddings)}
#'   \item \code{colnames(PBMC_DR) = c('RD1','RD2')}
#'   \item \code{PBMC_DR$Cluster = factor(srat@meta.data[rownames(PBMC_DR),'res.1'])}
#'   \item \code{PBMC_DR$Annotation = factor(c('7'='B','4'='B','1'='T_CD4','2'='T_CD4','3'='T_CD8','5'='T_CD8','6'='NK','8'='NK','0'='MNP','9'='MNP','10'='MNP','11'='?')[as.character(PBMC_DR$Cluster)])}
#' }
#' 
#' @format \code{PBMC_DR} is a data.frame with 4 columns: RD1, RD2, Cluster, and Annotation.
#' @usage data(PBMC_DR)
#' @name PBMC_DR
#' @docType data
#' @source \url{https://support.10xgenomics.com/single-cell-gene-expression/datasets/2.1.0/pbmc4k}
"PBMC_DR"

#' PBMC 4K cellular barcodes
#'
#' A vector of cellular barcodes
#' 
#' This vector was created by running \code{PBMC_cellBarcodes = colnames(Seurat::Read10X('filtered_gene_bc_matrices/GRCh38/'))}.
#'
#' @format \code{PBMC_cellBarcodes} is a vector containing the barcodes for droplets containing cells.
#' @usage data(PBMC_cellBarcodes)
#' @name PBMC_cellBarcodes
#' @docType data
#' @source \url{https://support.10xgenomics.com/single-cell-gene-expression/datasets/2.1.0/pbmc4k}
"PBMC_cellBarcodes"

#' tSNE coordinates for PBMC data
#' 
#' tSNE coordinates for PBMC 4k data used in the vignette.  tSNE have been calculated using standard Seurat settings.  Also includes a column indicating which Seurat cluster each cell belongs to.
#' 
#' The full list of commands used to create this object (assuming \code{scl} is a \code{SoupChannelList} object containing the PBMC4k data):
#' \code{set.seed(1)}
#' \code{srat = CreateSeuratObject(scl$toc)}
#' \code{srat = NormalizeData(srat)}
#' \code{srat = ScaleData(srat)}
#' \code{srat = FindVariableGenes(srat)}
#' \code{srat = RunPCA(srat,pcs.compute=30)}
#' \code{srat = RunTSNE(srat,dims.use=seq(30))}
#' \code{srat = FindClusters(srat,dims.use=seq(30),resolution=1)}
#' \code{PBMC_DR = as.data.frame(srat@dr$tsne@cell.embeddings)}
#' \code{colnames(PBMC_DR) = c('RD1','RD2')}
#' \code{PBMC_DR$Cluster = factor(srat@meta.data[rownames(PBMC_DR),'res.1'])}
#'
#' @name PBMC_DR
#' @usage data(PBMC_DR)
#' @docType data
#' @format Note that cells have been named assuming the PBMC4k channel is named "Channel1".  A data frame with 4,340 rows (one per cell) and 3 variables:
#' \describe{
#'   \item{RD1}{The first tSNE coordinate}
#'   \item{RD2}{The second tSNE coordinate}
#'   \item{Cluster}{A factor giving cluster assignment}
#' }
#' @source \url{https://support.10xgenomics.com/single-cell-gene-expression/datasets/2.1.0/pbmc4k}
NULL

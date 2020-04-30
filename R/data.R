#' PBMC data
#' 
#' Collection of bits of data relating to the 10X PBMC 4K data.
#' 
#' This collection includes \code{PBMC_DR}, \code{PBMC_cellBarcodes}, and \code{PBMC_tod}.  See individual help for more details.
#'
#' @name PBMC
#' @usage data(PBMC)
#' @docType data
#' @source \url{https://support.10xgenomics.com/single-cell-gene-expression/datasets/2.1.0/pbmc4k}
NULL



#' tSNE coordinates for PBMC data
#' 
#' tSNE coordinates for PBMC 4k data used in the vignette.  tSNE have been calculated using standard Seurat(v2) settings.  Also includes a column indicating which Seurat cluster each cell belongs to and a basic cell type annotation.
#' 
#' The full list of commands used to create this object (assuming \code{sc} is a \code{SoupChannel} object containing the PBMC4k data):
#' \code{set.seed(1)}
#' \code{srat = CreateSeuratObject(sc$toc)}
#' \code{srat = NormalizeData(srat)}
#' \code{srat = ScaleData(srat)}
#' \code{srat = FindVariableGenes(srat)}
#' \code{srat = RunPCA(srat,pcs.compute=30)}
#' \code{srat = RunTSNE(srat,dims.use=seq(30))}
#' \code{srat = FindClusters(srat,dims.use=seq(30),resolution=1)}
#' \code{PBMC_DR = as.data.frame(srat@dr$tsne@cell.embeddings)}
#' \code{colnames(PBMC_DR) = c('RD1','RD2')}
#' \code{PBMC_DR$Cluster = factor(srat@meta.data[rownames(PBMC_DR),'res.1'])}
#' \code{PBMC_DR$Annotation = factor(c('7'='B','4'='B','1'='T_CD4','2'='T_CD4','3'='T_CD8','5'='T_CD8','6'='NK','8'='NK','0'='MNP','9'='MNP','10'='MNP','11'='?')[as.character(PBMC_DR$Cluster)])}
#'
#' @name PBMC_DR
#' @usage data(PBMC)
#' @docType data
#' @format Note that cells have been named assuming the PBMC4k channel is named "Channel1".  A data frame with 4,340 rows (one per cell) and 3 variables:
#' \describe{
#'   \item{RD1}{The first tSNE coordinate}
#'   \item{RD2}{The second tSNE coordinate}
#'   \item{Cluster}{A factor giving cluster assignment}
#'   \item{Annotation}{A factor giving a cell type allocation for each cell}
#' }
#' @source \url{https://support.10xgenomics.com/single-cell-gene-expression/datasets/2.1.0/pbmc4k}
NULL


#' Raw PBMC data
#' 
#' Raw PBMC data used in the vignette.  This is the full matrix including all droplets in sparse format.  Created from the PBMC 4k data using:
#' \code{PBMC_tod = Read10X('raw_gene_bc_matrices/GRCh38/')}
#' 
#' @name PBMC_tod
#' @usage data(PBMC)
#' @docType data
#' @format Sparse matrix with rows representing genes, columns cells and entries number of unique molecular identifiers detected.
#' @source \url{https://support.10xgenomics.com/single-cell-gene-expression/datasets/2.1.0/pbmc4k}
NULL


#' Cell barcodes for PBMC data
#' 
#' Which barcodes contain cells (as determined by cellranger) in the PBMC data.  Created by running:
#' \code{PBMC_cellBarcodes = colnames(Read10X('filtered_gene_bc_matrices/GRCh38/'))}
#' 
#' @name PBMC_cellBarcodes
#' @usage data(PBMC)
#' @docType data
#' @format Vector of cellular barcodes.
#' @source \url{https://support.10xgenomics.com/single-cell-gene-expression/datasets/2.1.0/pbmc4k}
NULL

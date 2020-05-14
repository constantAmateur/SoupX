#' PBMC 4K meta data
#'
#' Collection of bits of meta data relating to the 10X PBMC 4K data.
#' 
#' This data set pertains to the 10X demonstration PBMC 4K data and includes metadata about it in the \code{data.frame} named \code{PBMC_metaData}.
#'
#' \code{PBMC_metaData} was created using Seurat (v2) to calculate a tSNE representation of the data and cluster cells with these commands.
#' \itemize{
#'   \item \code{set.seed(1)}
#'   \item \code{srat = CreateSeuratObject(sc$toc)}
#'   \item \code{srat = NormalizeData(srat)}
#'   \item \code{srat = ScaleData(srat)}
#'   \item \code{srat = FindVariableGenes(srat)}
#'   \item \code{srat = RunPCA(srat,pcs.compute=30)}
#'   \item \code{srat = RunTSNE(srat,dims.use=seq(30))}
#'   \item \code{srat = FindClusters(srat,dims.use=seq(30),resolution=1)}
#'   \item \code{PBMC_metaData = as.data.frame(srat@dr$tsne@cell.embeddings)}
#'   \item \code{colnames(PBMC_metaData) = c('RD1','RD2')}
#'   \item \code{PBMC_metaData$Cluster = factor(srat@meta.data[rownames(PBMC_metaData),'res.1'])}
#'   \item \code{PBMC_metaData$Annotation = factor(c('7'='B','4'='B','1'='T_CD4','2'='T_CD4','3'='T_CD8','5'='T_CD8','6'='NK','8'='NK','0'='MNP','9'='MNP','10'='MNP','11'='?')[as.character(PBMC_metaData$Cluster)])}
#' }
#' 
#' @format \code{PBMC_metaData} is a data.frame with 4 columns: RD1, RD2, Cluster, and Annotation.
#' @usage data(PBMC_metaData)
#' @name PBMC_metaData
#' @docType data
#' @source \url{https://support.10xgenomics.com/single-cell-gene-expression/datasets/2.1.0/pbmc4k}
"PBMC_metaData"

#' SoupChannel from PBMC data
#'
#' \code{\link{SoupChannel}} created from 10X demonstration PBMC 4k data.
#' 
#' \code{PBMC_sc} was created by running the following commands.
#' \itemize{
#'   \item \code{set.seed(1137)}
#'   \item \code{tmpDir = tempdir(check=TRUE)}
#'   \item \code{download.file('http://cf.10xgenomics.com/samples/cell-exp/2.1.0/pbmc4k/pbmc4k_raw_gene_bc_matrices.tar.gz',destfile=file.path(tmpDir,'tod.tar.gz'))}
#'   \item \code{download.file('http://cf.10xgenomics.com/samples/cell-exp/2.1.0/pbmc4k/pbmc4k_filtered_gene_bc_matrices.tar.gz',destfile=file.path(tmpDir,'toc.tar.gz'))}
#'   \item \code{untar(file.path(tmpDir,'tod.tar.gz'),exdir=tmpDir)}
#'   \item \code{untar(file.path(tmpDir,'toc.tar.gz'),exdir=tmpDir)}
#'   \item \code{library(SoupX)}
#'   \item \code{PBMC_sc = load10X(tmpDir)}
#' }
#' 
#' @format \code{PBMC_sc} is a \code{SoupChannel} object with 33,694 genes and 4,340 cells.
#' @usage data(PBMC_sc)
#' @name PBMC_sc
#' @docType data
#' @source \url{https://support.10xgenomics.com/single-cell-gene-expression/datasets/2.1.0/pbmc4k}
"PBMC_sc"


#' Toy SoupChanel object
#'
#' A \code{\link{SoupChannel}} object created from the toy data used in examples.
#' 
#' The toy data is created from a modified version of the extremely reduced \code{Seurat} \code{pbmc_small} dataset.  It includes clusters, tSNE coordinates and a flat estimate of 0.1 contamination.  It includes data for only 226 genes and 62 cells and should not be used for anything other than testing functions as it is not representative of real data in any way.
#'
#' @format \code{scToy} is a \code{SoupChannel} object.
#' @usage data(scToy)
#' @name scToy
#' @docType data
"scToy"


#' Calculate which cells genuinely do not express a particular gene or set of genes
#' 
#' Given a list of correlated genes (e.g. Haemoglobin genes, Immunoglobulin genes, etc.), make an attempt to estimate which cells genuinely do not express each of these gene sets in turn.  The central idea is that in cells that are not genuinely producing a class of mRNAs (such as haemoglobin genes), any observed expression of these genes must be due to ambient RNA contamination.  As such, if we can identify these cells, we can use the observed level of expression of these genes to estimate the level of contamination.
#' 
#' The ideal way to do this would be to have a prior annotation of your data indicating which cells are (for instance) red blood cells and genuinely expression haemoglobin genes, and which do not and so only express haemoglobin genes due to contamination.  If this is your circumstance, there is no need to run this function, you can instead pass a matrix encoding which cells are haemoglobin expressing and which are not to \code{\link{calculateContaminationFraction}} via the \code{useToEst} parameter.
#' 
#' This function will use a conservative approach to excluding cells that it thinks may express one of your gene sets.  This is because falsely including a cell in the set of non-expressing cells may erroneously inflate your estimated contamination, whereas failing to include a geunine non-expressing cell in this set has no significant effect.
#'
#' To this end, this function will exclude any cluster of cells in which any cell is deemed to have geunine expression of a gene set.  Clustering of data is beyond the scope of this package, but can be performed by the user.  In the case of 10X data mapped using cellranger and loaded using \code{\link{load10X}}, the cellranger graph based clustering is automatically loaded and used.
#'
#' To decide if a cell is genuinely expressing a set of genes, a Poisson test is used.  This tests whether the observed expression is greater than \code{maximumContamination} times the expected number of counts for a set of genes, if the cell were assumed to be derived wholy from the background.  This process can be made less conservative (i.e., excluding fewer cells/clusters) by either decreasing the value of the maximum contamination the user believes is plausible (\code{maximumContamination}) or making the significance threshold for the test more strict (by reducing \code{FDR}).
#'
#' @export
#' @param sc A SoupChannel object. 
#' @param nonExpressedGeneList A list containing sets of genes which will be used to estimate the contamination fraction.
#' @param clusters A named vector indicating how to cluster cells.  Names should be cell IDs, values cluster IDs.  If NULL, we will attempt to load it from sc$metaData$clusters.  If set to FALSE, each cell will be considered individually.
#' @param maximumContamination The maximum contamination fraction that you would reasonably expect.  The lower this value is set, the more aggressively cells are excluded from use in estimation.
#' @param FDR A Poisson test is used to identify cells to exclude, this is the false discevore rate to use.  Higher FDR = more aggressive exclusion.
#' @seealso calculateContaminationFraction plotMarkerMap
#' @return A matrix indicating which cells to be used to estimate contamination for each set of genes.  Typically passed to the \code{useToEst} parameter of \code{\link{calculateContaminationFraction}} or \code{\link{plotMarkerMap}}.
estimateNonExpressingCells = function(sc,nonExpressedGeneList,clusters=NULL,maximumContamination=1.0,FDR=0.05){
  if(!is(sc,'SoupChannel'))
    stop("sc is not a valid SoupChannel object")
  #Get clusters if they exist, if they don't, set to individual cells
  if(is.null(clusters)){
    if('clusters' %in% colnames(sc$metaData)){
      clusters = setNames(as.character(sc$metaData$clusters),rownames(sc$metaData))
    }
  }
  #Using each cell as its own cluster
  if(is.null(clusters) || (length(clusters)==1 && clusters==FALSE)){
      message("No clusters found or supplied, using every cell as its own cluster.")
      clusters = setNames(rownames(sc$metaData),rownames(sc$metaData))
  }
  #Check we have coverage of everything
  if(!all(colnames(sc$toc) %in% names(clusters)))
    stop("Invalid cluster specification.  clusters must be a named vector with all column names in the table of counts appearing.")
  #Convert gene list to genuine list if vector
  if(!is.list(nonExpressedGeneList))
    stop("nonExpressedGeneList must be a list of sets of genes.  e.g. list(HB = c('HBB','HBA2'))")
  #Now work out which clusters to use which genes on 
  tgtGns = unique(unlist(nonExpressedGeneList))
  dat = sc$toc[tgtGns,,drop=FALSE]
  cnts = do.call(rbind,lapply(nonExpressedGeneList,function(e) colSums(dat[e,,drop=FALSE])))
  #Work out how many counts we'd expect if the cell were maximally contaminated and all expression came from the contamination
  exp = outer(sc$soupProfile[tgtGns,'est'],sc$metaData$nUMIs*maximumContamination)
  colnames(exp) = colnames(cnts)
  rownames(exp) = tgtGns
  exp = do.call(rbind,lapply(nonExpressedGeneList,function(e) colSums(exp[e,,drop=FALSE])))
  #Identify those cells where a gene is definitely expressed
  s = split(names(clusters),clusters)
  clustExp = ppois(cnts-1,exp,lower.tail=FALSE)
  clustExp = t(apply(clustExp,1,p.adjust,method='BH'))
  clustExp = do.call(rbind,lapply(s,function(e) apply(clustExp[,e,drop=FALSE],1,min)))
  clustExp = clustExp>=FDR
  #Expand it out into a full cell matrix
  clustExp = clustExp[match(clusters,rownames(clustExp)),,drop=FALSE]
  rownames(clustExp) = names(clusters)
  #Check that we actually got som
  if(sum(clustExp)==0){
    warning("No non-expressing cells identified.  Consider setting clusters=FALSE, increasing maximumContamination and/or FDR")
  }
  #A small number found
  if(sum(clustExp)>0 && sum(clustExp)<100){
    warning("Fewer than 100 non-expressing cells identified.  The estimation of the contamination fraction may be inaccurate.  Consider setting clusters=FALSE, increasing maximumContamination and/or FDR")
  }
  return(clustExp)
}

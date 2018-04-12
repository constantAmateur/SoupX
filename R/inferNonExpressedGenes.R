#' Plots the genes most likely to be non-expressed in some fraction of cells
#'
#' Ideally, prior biological knowledge can be used to select a set of genes that are known to be non-expressed in at least some cells, that can be used to estimate the contamination fraction.  For solid tissues, blood is a ubiquitous and unavoidable contaminant, making Haemoglobin genes ideal for this estimation (as they are genuinely expressed only in red blood cells, which can be trivially identified).  
#'
#' The purpose of this function is to attempt to identify genes that are highly cell specific and so can be used to estimate background contamination in those cells that do not express them.  Although it would be possible to select these algorithmically, this is deliberately not done as selection of inappropriate genes can lead to over estimation of the contamination fraction and over correction of the expression profiles.  Instead, this function produces a diagnostic plot and a series of tables that the user should use to guide together with their biological knowledge to select sensible candidates.  Usually 1 or 2 highly expressed genes are sufficient to get a good channel level estimation of the contamination fraction, which is preferable to using genes of uncertain validity.
#'
#' If \code{soupProfile} is not given, \code{tod} must be given as it will be used to estimate the soup profile from the empty droplets.
#'
#' The names of markers can be passed to \code{\link{plotMarkerDistribution}} using the parameters \code{nonExpressedGeneList} to visualise candidate genes.
#'
#' @export
#' @param toc Table of UMIs for all cells.
#' @param soupProfile The soup profile for this channel (output of \code{\link{estimateSoup}}).  One of \code{soupProfile} or \code{tod} must be specified.
#' @return A table summarising the suitability of each gene as a candidate for estimating contamination.  The columns in this table are respectively: the number of cells expressing this gene, the number of cells with expression less than the soup, the fraction that are low, an "extremity score" (mean squared log-ratio), a "centrality score" (mean of 1/(1+x^2)) and an indicator if this gene has passed the criteria for being potentially useful.
inferNonExpressedGenes = function(toc,soupProfile){
  #The usual thing
  nUMIs = colSums(toc)
  #Construct expression fractions from table of counts
  rat = t(t(toc)/nUMIs)
  #Convert to useful format
  rat = as(rat,'dgTMatrix')
  #Now convert to ratio
  rat@x = rat@x/soupProfile[rat@i+1,'est']
  #Exclude the useless rows
  toKeep = unique(rat@i+1)
  toKeep = toKeep[soupProfile[toKeep,'est']>0]
  rat = rat[toKeep,]
  #Get summary stats
  nCells = table(factor(rownames(rat)[rat@i+1],levels=rownames(rat)))
  lowCount = table(factor(rownames(rat)[rat@i[rat@x<1]+1],levels=rownames(rat)))
  tmp = split(rat@x,rownames(rat)[rat@i+1])[names(nCells)]
  extremity = sapply(tmp,function(e) sum(log10(e)**2)/length(e))[names(lowCount)]
  centrality = sapply(tmp,function(e) sum(1/(1+log10(e)**2))/length(e))[names(lowCount)]
  minFrac = log10(sapply(tmp,min))
  #Construct targets data frame
  targets = data.frame(nCells = as.integer(nCells),
                       lowCount = as.integer(lowCount),
                       lowFrac = as.integer(lowCount)/as.integer(nCells),
                       extremity = extremity,
                       centrality = centrality,
                       minFrac = minFrac,
                       isUseful = as.integer(lowCount)>0.1*ncol(rat))
  #Order by something
  targets = targets[order(targets$isUseful,targets$extremity,decreasing=TRUE),]
  return(targets)
}

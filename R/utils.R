#' Expands soup counts calculated at the cluster level to the cell level
#'
#' Given a clustering of cells and soup counts calculated for each of those clusters, determines a most likely allocation of soup counts at the cell level.
#'
#' @param clustSoupCnts Matrix of genes (rows) by clusters (columns) where counts are number of soup counts for that gene/cluster combination.
#' @param cellObsCnts Matrix of genes (rows) by cells (columns) giving the observed counts
#' @param clusters Mapping from cells to clusters.
#' @param cellWeights Weighting to give to each cell when distributing counts.  This would usually be set to the number of expected soup counts for each cell.
#' @param verbose Integer giving level of verbosity.  0 = silence, 1 = Basic information, 2 = Very chatty, 3 = Debug.
#' @return A matrix of genes (rows) by cells (columns) giving the number of soup counts estimated for each cell.  Non-integer values possible.
expandClusters = function(clustSoupCnts,cellObsCnts,clusters,cellWeights,verbose=1){
  ws = cellWeights
  #Do one cluster at a time
  if(verbose>0)
    message(sprintf("Expanding counts from %d clusters to %d cells.",ncol(clustSoupCnts),ncol(cellObsCnts)))
  #lapply instead of loop is a hold-over from when mclapply was an option
  out = lapply(seq(ncol(clustSoupCnts)),
                 function(j) {
    if(verbose>1)
      message(sprintf("Expanding cluster %s",colnames(clustSoupCnts)[j]))
    #Which cells
    wCells = which(clusters==colnames(clustSoupCnts)[j])
    #How should they be weighted
    ww = ws[wCells]/sum(ws[wCells])
    #What is the limits
    lims = cellObsCnts[,wCells,drop=FALSE]
    #And how many soup
    nSoup = clustSoupCnts[,j]
    #Create the output object
    expCnts = as(lims,'dgTMatrix')
    #Most cases are easily dealt with.  In rough order of frequency.
    #1. No soup for gene - set to zero
    #2. All counts for gene are soup - set to lims
    #3. Not all counts are soup, but every entry is a 0 or 1 so no iteration needed.
    #4. Some iteration needed.
    #Deal with case 1
    expCnts@x[(expCnts@i+1) %in% which(nSoup==0)]=0
    #Case 2 is dealt with by construction
    #Get set of genes for cases 3 and 4
    wGenes = which(nSoup>0 & nSoup<rowSums(lims))
    #And deal with them as appropriate.  Save time by looking only at non-zero entries
    w = which((expCnts@i+1) %in% wGenes)
    w = split(w,expCnts@i[w]+1)
    tmp = lapply(w,function(e) alloc(nSoup[expCnts@i[e[1]]+1],expCnts@x[e],ww[expCnts@j[e]+1]))
    expCnts@x[unlist(w,use.names=FALSE)] = unlist(tmp,use.names=FALSE)
    return(expCnts)
  })
  out = do.call(cbind,out)
  out = out[,colnames(cellObsCnts)]
  return(out)
}

#' Create Seurat style progress bar
#'
#' Creates progress bar that won't ruin log files and shows progress towards 100%.
#'
#' @param min Minimum value of parameter.
#' @param max Maximum value of parameter.
#' @param ... Passed to \code{\link{txtProgressBar}}
#' @return A txtProgressBar object to use updating progress.
initProgBar = function(min,max,...){
  message('0%   10   20   30   40   50   60   70   80   90   100%')
  message('|----|----|----|----|----|----|----|----|----|----|')
  pb=txtProgressBar(min=min,max=max,style=1,width=51,char='*',...)
  return(pb)
}

#' Allocate values to "buckets" subject to weights and constraints
#'
#' Allocates \code{tgt} of something to \code{length(bucketLims)} different "buckets" subject to the constraint that each bucket has a maximum value of \code{bucketLims} that cannot be exceeded.  By default counts are distributed equally between buckets, but weights can be provided using \code{ws} to have the redistribution prefer certain buckets over others.
#'
#' @param tgt Value to distribute between buckets.
#' @param bucketLims The maximum value that each bucket can take.  Must be a vector of positive values.
#' @param ws Weights to be used for each bucket.  Default value makes all buckets equally likely.
#' @return A vector of the same length as \code{bucketLims} containing values distributed into buckets.
#' @examples
#' set.seed(1137)
#' ws = abs(rnorm(10))
#' tops = round(runif(10)*3)
#' #Simple case where the bucket limits change nothing
#' SoupX:::alloc(1,tops,ws)
#' #Case where some buckets get full
#' SoupX:::alloc(8,tops,ws)
alloc = function(tgt,bucketLims,ws=rep(1/length(bucketLims),length(bucketLims))){
  #Normalise weights
  ws = ws/sum(ws)
  #Save time in line
  if(all(tgt*ws<=bucketLims))
    return(tgt*ws)
  #Need to order things in the order they'll be removed as the tgt increases
  o = order(bucketLims/ws)
  w = ws[o]
  y = bucketLims[o]
  #The formula for number removed at entry i is
  #k_i = \frac{y_i}{w_i} (1- \sum_j=0^{i-1} w_j) + \sum_j=0^{i-1} y_j
  cw = cumsum(c(0,w[-length(w)]))
  cy = cumsum(c(0,y[-length(y)]))
  k = y/w* (1 - cw) + cy
  #Handle zero-weights appropriately
  k[w==0] = Inf
  #Everything that has k<=tgt will be set to y
  b = (k<=tgt)
  #We then need to work out how many counts to distribute we have left over and distribute them according to re-normalised weights
  resid = tgt-sum(y[b])
  w = w/(1-sum(w[b]))
  out = ifelse(b,y,resid*w)
  #Need to reverse sort
  return(out[order(o)])
}

#' Expands soup counts calculated at the cluster level to the cell level
#'
#' Given a clustering of cells and soup counts calculated for each of those clusters, determines a most likely allocation of soup counts at the cell level.
#'
#' @param clustSoupCnts Matrix of genes (rows) by clusters (columns) where counts are number of soup counts for that gene/cluster combination.
#' @param cellObsCnts Matrix of genes (rows) by cells (columns) giving the observed counts
#' @param clusters Mapping from cells to clusters.
#' @param cellWeights Weighting to give to each cell when distributing counts.  This would usually be set to the number of expected soup counts for each cell.
#' @param verbose Integer giving level of verbosity.  0 = silence, 1 = Basic information, 2 = Very chatty, 3 = Debug.
#' @param nCores Number of cores to use.  Defaults to all cores.
#' @param ... Passed to mclapply
#' @return A matrix of genes (rows) by cells (columns) giving the number of soup counts estimated for each cell.  Non-integer values possible.
#' @importFrom parallel mclapply
expandClusters = function(clustSoupCnts,cellObsCnts,clusters,cellWeights,verbose=1,nCores=getOption('mc.cores',1),...){
  ws = cellWeights
  #Do one cluster at a time
  if(verbose>0)
    message(sprintf("Expanding counts from %d clusters to %d cells.",ncol(clustSoupCnts),ncol(cellObsCnts)))
  #Do the expansion if needed
  out = mclapply(seq(ncol(clustSoupCnts)),
                 mc.cores = nCores,
                 FUN = function(j) {
    if(nCores>1 & verbose>1)
      verbose=1
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
    #And one gene at a time
    if(length(wGenes)>0 & verbose>1)
      pb = initProgBar(1,length(wGenes))
    ii=1
    for(i in wGenes){
      ii=ii+1
      #Initialise
      toDo = nSoup[i]
      wLocal = which(expCnts@i+1 ==i)
      outLocal = rep(0,length(wLocal))
      tgtLocal = seq_along(wLocal)
      limLocal = expCnts@x[wLocal]
      weightsLocal = ww[expCnts@j[wLocal]+1]
      while(TRUE){
        #Adjust guess
        outLocal[tgtLocal] = outLocal[tgtLocal] + toDo*ww[tgtLocal]/sum(ww[tgtLocal])
        #Update counter of remaining allocation budget
        w = outLocal[tgtLocal]-limLocal[tgtLocal]
        toDo = sum(w[w>0])
        if(verbose>2)
          message(sprintf("%d at or above limit, with %g unallocated",sum(w>0),toDo))
        if(toDo==0)
          break
        #If we have work to do, truncate.
        outLocal[tgtLocal[w>0]] = limLocal[tgtLocal[w>0]]
        #Update definition of local
        tgtLocal = tgtLocal[w<0]
      }
      expCnts@x[wLocal] = outLocal
      if(verbose>1)
        setTxtProgressBar(pb,ii)
    }
    if(length(wGenes)>0 & verbose>1)
      close(pb)
    return(expCnts)
  },...)
  ##  out[[j]] = expCnts
  ##}
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

#' Remove background contamination from count matrix
#'
#' After the level of background contamination has been estimated or specified for a channel, calculate the resulting corrected count matrix with background contamination removed.  
#'
#' This essentially subtracts off the mean expected background counts for each gene, then redistributes any "unused" counts.  A count is unused if its subtraction has no effect.  For example, subtracting a count from a gene that has zero counts to begin with.
#'
#' As expression data is highly sparse at the single cell level, it is highly recommended that clustering information be provided to allow the subtraction method to share information between cells.  Without grouping cells into clusters, it is difficult (and usually impossible) to tell the difference between a count of 1 due to background contamination and a count of 1 due to endogenous expression.  This ambiguity is removed at the cluster level where counts can be aggregated across cells.  This information can then be propogated back to the individual cell level to provide a more accurate removal of contaminating counts.
#' 
#' To provide clustering information, either set clustering on the SoupChannel object with \code{\link{setClusters}} or explicitly passing the \code{clusters} paramater.  
#'
#' If \code{roundToInt=TRUE}, this function will round the result to integers.  That is, it will take the floor of the connected value and then round back up with probability equal to the fractional part of the number.
#' 
#' The \code{method} parameter controls how the removal of counts in performed.  This should almost always be left at the default ('subtraction'), which iteratively subtracts counts from all genes as described above.  The 'soupOnly' method will use a p-value based estimation procedure to identify those genes that can be confidently identified as having endogenous expression and removes everything else (described in greater detail below).  Because this method either removes all or none of the expression for a gene in a cell, the correction procedure is much faster.  Finally, the 'multinomial' method explicitly maximises the multinomial likelihood for each cell.  This method gives essentially identical results as 'subtraction' and is considerably slower.
#'
#' In greater detail, the 'soupOnly' method is done by sorting genes within each cell by their p-value under the null of the expected soup fraction using a Poisson model.  So that genes that definitely do have a endogenous contribution are at the end of the list with p=0.  Those genes for which there is poor evidence of endogenous cell expression are removed, until we have removed approximately nUMIs*rho molecules.  The cut-off to prevent removal of genes above nUMIs*rho in each cell is achieved by calculating a separate p-value for the total number of counts removed to exceed nUMIs*rho, again using a Poisson model.  The two p-values are combined using Fisher's method and the cut-off is applied to the resulting combined p-value calculated using a chi-squared distribution with 4 degrees of freedom.
#'
#' @export
#' @param sc A SoupChannel object.
#' @param clusters A vector of cluster IDs, named by cellIDs.  If NULL clusters auto-loaded from \code{sc}.  If FALSE, no clusters are used.  See details.
#' @param method Method to use for correction.  See details.  One of 'multinomial', 'soupOnly', or 'subtraction'
#' @param roundToInt Should the resulting matrix be rounded to integers?
#' @param verbose Integer giving level of verbosity.  0 = silence, 1 = Basic information, 2 = Very chatty, 3 = Debug.
#' @param tol Allowed deviation from expected number of soup counts.  Don't change this.
#' @param pCut The p-value cut-off used when \code{method='soupOnly'}.
#' @param ... Passed to expandClusters.  Of particular interest is the nCores parameter which can be used to parallelise the calculation.
#' @return A modified version of the table of counts, with background contamination removed.
#' @importFrom Matrix sparseMatrix Matrix
#' @importFrom stats rbinom pchisq
adjustCounts = function(sc,clusters=NULL,method=c('subtraction','soupOnly','multinomial'),roundToInt=FALSE,verbose=1,tol=1e-3,pCut=0.01,...){
  #####################
  # Parameter checking
  method = match.arg(method)
  if(!is(sc,'SoupChannel'))
    stop("sc must be an object of type SoupChannel")
  if(!'rho' %in% colnames(sc$metaData))
    stop("Contamination fractions must have already been calculated/set.")
  #Check the clusters parameter. If it's null, try and auto-fetch
  if(is.null(clusters)){
    if('clusters' %in% colnames(sc$metaData)){
      clusters = setNames(as.character(sc$metaData$clusters),rownames(sc$metaData))
    }else{
      warning("Clustering data not found.  Adjusting counts at cell level.  You will almost certainly get better results if you cluster data first.")
      clusters=FALSE
    }
  }
  ################################################
  # Recursive application for when using clusters
  if(clusters[1]!=FALSE){
    #Check we have coverage of everything
    if(!all(colnames(sc$toc) %in% names(clusters)))
      stop("Invalid cluster specification.  clusters must be a named vector with all column names in the table of counts appearing.")
    #OK proceed
    s = split(colnames(sc$toc),clusters[colnames(sc$toc)])
    tmp = sc
    #Create by cluster table of counts.  Must be sparse matrix
    tmp$toc = do.call(cbind,lapply(s,function(e) rowSums(sc$toc[,e,drop=FALSE])))
    tmp$toc = Matrix(tmp$toc,sparse=TRUE)
    tmp$metaData = data.frame(nUMIs = sapply(s,function(e) sum(sc$metaData[e,'nUMIs'])),
                              rho = sapply(s,function(e) sum(sc$metaData[e,'rho']*sc$metaData[e,'nUMIs'])/sum(sc$metaData[e,'nUMIs'])))
    #Recursively apply.  Could be done more neatly I'm sure.  Don't round to integer until we've re-expanded.
    out = adjustCounts(tmp,clusters=FALSE,method=method,roundToInt=FALSE,verbose=verbose,tol=tol,pCut=pCut)
    #This gives the corrected table, we want the soup counts
    out = tmp$toc - out
    #Finally re-expand the results back to single cell level
    out = expandClusters(out,sc$toc,clusters,sc$metaData$nUMIs*sc$metaData$rho,verbose=verbose,...)
    #And convert back to a corrected table of counts
    out = sc$toc - out
  }else{
    ##############################
    # Actual adjustment of counts
    if(method=='multinomial'){
      #Quick initial guess at best fit
      if(verbose>1)
        message("Initialising with subtraction method.")
      fitInit = sc$toc - adjustCounts(sc,clusters=FALSE,method='subtraction',roundToInt=TRUE)
      ps = sc$soupProfile$est
      out = list()
      #Loop over cells
      if(verbose>0){
        message(sprintf("Fitting multinomial distribution to %d cells/clusters.",ncol(sc$toc)))
        pb=initProgBar(1,ncol(sc$toc))
      }
      for(i in seq(ncol(sc$toc))){
        if(verbose>0)
          setTxtProgressBar(pb,i)
        #How many soup molecules do we expect for this cell?
        tgtN = round(sc$metaData$rho[i]*sc$metaData$nUMIs[i])
        #And what are the observational limits on which genes they can be assigned to
        lims = sc$toc[,i]
        #Initialise 
        fit = fitInit[,i]
        while(TRUE){
          #Work out which we can increase
          increasable = fit<lims
          decreasable = fit>0
          #And what the likelihood gain/cost for changing them is
          delInc = log(ps[increasable]) - log(fit[increasable]+1)
          delDec = -log(ps[decreasable]) +log(fit[decreasable])
          #Decide which swap(s) would lead to best likelihood gain
          wInc = wIncAll = which(increasable)[which(delInc==max(delInc))]
          wDec = wDecAll = which(decreasable)[which(delDec==max(delDec))]
          nInc = length(wIncAll)
          nDec = length(wDecAll)
          if(nInc>1)
            wInc = sample(wIncAll,1)
          if(nDec>1)
            wDec = sample(wDecAll,1)
          #How many do we have
          curN = sum(fit)
          if(curN<tgtN){
            if(verbose>2)
              message(sprintf("# choices: nInc=%d nDec=%d, Under-allocated (%d of %d), increasing...",nInc,nDec,curN,tgtN))
            fit[wInc] = fit[wInc]+1
          }else if(curN>tgtN){
            if(verbose>2)
              message(sprintf("# choices: nInc=%d nDec=%d, Over-allocated (%d of %d), decreasing...",nInc,nDec,curN,tgtN))
            fit[wDec] = fit[wDec]-1
          }else{
            #Check if the swap will increase likelihood
            delTot = max(delInc)+max(delDec)
            #Three possibilites
            #1. del>0, keep going, we're improving the fit
            #2. del<0, stop, we won't get any better
            #3. del==0, accumulate all the ones that can go up/down.  As del==0 this is reversable, so want to share out "excess" counts between this group.
            if(verbose>2)
              message(sprintf("# choices: nInc=%d nDec=%d, Total log likelihood difference %s",nInc,nDec,delTot))
            if(delTot==0){
              #As the difference is zero, all movements between wInc and wDec are reversible.  So want to distribute evenly the "available" counts between those in the ambiguous set.
              #Take them away from those that presently have them
              fit[wDecAll] = fit[wDecAll]-1
              #And share them equally amongst those in bucket
              zeroBucket = unique(c(wIncAll,wDecAll))
              fit[zeroBucket] = fit[zeroBucket] + length(wDecAll)/length(zeroBucket)
              if(verbose>2)
                message(sprintf("Ambiguous final configuration. Shared %d reads between %d equally likely options",length(wDecAll),length(zeroBucket)))
              break
            }else if(delTot<0){
              #Minimum reached.
              if(verbose>2)
                message("Unique final configuration.")
              break
            }else{
              #Improvements to be made.
              fit[wInc] = fit[wInc]+1
              fit[wDec] = fit[wDec]-1
            }
          }
        }
        out[[i]]=fit
      }
      if(verbose>0)
        close(pb)
      out = do.call(cbind,out)
      out = as(out,'dgTMatrix')
      rownames(out) = rownames(sc$toc)
      colnames(out) = colnames(sc$toc)
      out = sc$toc - out
    }else if(method=='soupOnly'){
      if(verbose>0)
        message("Identifying and removing genes likely to be pure contamination in each cell.")
      #Start by calculating the p-value against the null of soup.
      out = as(sc$toc,'dgTMatrix')
      if(verbose>1)
        message("Calculating probability of each gene being soup")
      p = ppois(out@x-1,sc$metaData$nUMIs[out@j+1]*sc$soupProfile$est[out@i+1]*sc$metaData$rho[out@j+1],lower.tail=FALSE)
      #Order them by cell, then by p-value
      o = order(-(out@j+1),p,decreasing=TRUE)
      #Get the running total for removal.  Could probably make this faster with some tricks.
      if(verbose>1)
        message("Calculating probability of the next count being soup")
      s = split(o,out@j[o]+1)
      rTot = unlist(lapply(s,function(e) cumsum(out@x[e])),use.names=FALSE)
      #Now we need to get the soup probability vector as well.
      pSoup = ppois(rTot-out@x[o]-1,sc$metaData$nUMIs[out@j[o]+1]*sc$metaData$rho[out@j[o]+1],lower.tail=FALSE)
      if(verbose>1)
        message("Filtering table of counts")
      #Now construct the final probability vector
      pp = p[o]*pSoup
      q = pchisq(-2*log(pp),4,lower.tail=FALSE)
      #And we then want to drop anything meeting our termination criteria
      w = which(q<pCut)
      #Keep a list of dropped genes
      dropped = data.frame(cell=colnames(out)[out@j[o[-w]]+1],
                           gene=rownames(out)[out@i[o[-w]]+1],
                           cnt = out@x[o[-w]])
      if(verbose>2){
        message(sprintf("Most removed genes are:"))
        x = sort(table(dropped$gene)/ncol(out),decreasing=TRUE)
        print(x[seq_len(min(length(x),100))])
      }
      #Construct the corrected count matrix
      out = sparseMatrix(i=out@i[o[w]]+1,
                         j=out@j[o[w]]+1,
                         x=out@x[o[w]],
                         dims=dim(out),
                         dimnames=dimnames(out),
                         giveCsparse=FALSE)
    }else if(method=='subtraction'){
      #Create the final thing efficiently without making a big matrix
      out = as(sc$toc,'dgTMatrix')
      expSoupCnts = sc$metaData$nUMIs * sc$metaData$rho
      soupFrac = sc$soupProfile$est
      #Distribute counts according to the soup profile.  Could be made faster by not considering zeros, but eh.
      out = out - do.call(cbind,lapply(seq(ncol(out)),function(e) alloc(expSoupCnts[e],out[,e],soupFrac)))
      out = as(out,'dgTMatrix')
      #Iteratively remove until we've removed the expected number of counts
      #How many counts should we have left when we're done?
      #tgts = sc$metaData$nUMIs - expSoupCnts
      #Which genes do we still need to bother trying to remove counts from in which cells
      #toAdjust = seq_along(out@i)
      #if(verbose>0)
      #  message("Subtracting contaminating counts")
      #while(TRUE){
      #  #How many left to do?
      #  toDo = colSums(out)-tgts
      #  #Get the soup frac correction factors.
      #  tmp = rep(1,length(toDo))
      #  soupFracSum = sapply(split(soupFrac[out@i[toAdjust]+1],out@j[toAdjust]+1),sum)
      #  tmp[as.numeric(names(soupFracSum))]=soupFracSum
      #  toDo = toDo/tmp
      #  #Do the adjustment
      #  out@x[toAdjust] = out@x[toAdjust]-soupFrac[out@i[toAdjust]+1]*toDo[out@j[toAdjust]+1]
      #  #Only keep updating those that need it
      #  toAdjust = toAdjust[out@x[toAdjust]>0]
      #  out@x[out@x<0]=0
      #  if(verbose>1)
      #    print(quantile(colSums(out)-tgts))
      #  if(max(colSums(out)-tgts)<tol)
      #    break
      #}
      ##This is the clearer but slower version of above
      #out@x = pmax(0,out@x - soupFrac[out@i+1]*expSoupCnts[out@j+1])
      #while(max(colSums(out)-tgts)>tol){
      #  toDo = colSums(out)-tgts
      #  out@x = pmax(0,out@x - soupFrac[out@i+1]*toDo[out@j+1])
      #  print(quantile(colSums(out)-tgts))
      #}
      #Fix the object internals
      w = which(out@x>0)
      out = sparseMatrix(i=out@i[w]+1,
                         j=out@j[w]+1,
                         x=out@x[w],
                         dims=dim(out),
                         dimnames=dimnames(out),
                         giveCsparse=FALSE)
    }else{
      stop("Impossible!")
    }
  }
  #Do stochastic rounding to integers if needed
  if(roundToInt){
    if(verbose>1)
      message("Rounding to integers.")
    #Round to integer below, then probabilistically bump back up
    out@x = floor(out@x)+rbinom(length(out@x),1,out@x-floor(out@x))
  }
  return(out)
}

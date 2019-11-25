#' Remove background contamination from cell expression profile
#'
#' After the level of background contamination has been estimated or specified for a channel, calculate the resulting corrected count matrix with background contamination removed.  
#'
#' This essentially subtracts off the mean expected background counts for each gene, then redistributes any "unused" counts.  A count is unused if its subtraction has no effect.  For example, subtracting a count from a gene that has zero counts to begin with.
#'
#' If \code{dropWholeCountsOnly=TRUE} then an aggressive alternative strategy is used.  **This is primarily kept for historical reasons and is not recommended.**  For each gene in each cell, this strategy either does nothing, or sets the counts to 0.  This is very effective at removing the smattering of low-level contamination, but leaves any gene that has a contribution from both background and the cell itself unchanged.  This estimation is done by sorting genes within each cell by their p-value under the null of the expected soup fraction.  So that genes that definitely do have a endogenous contribution are at the end of the list with p=0.  Those genes for which there is poor evidence of endogenous cell expression are removed, until we have removed approximately nUMIs*rho molecules.
#' 
#' If \code{roundToInt=TRUE}, this function will round the result to integers.  That is, it will take the floor of the connected value and then round back up with probability equal to the fractional part of the number.
#'
#' @export
#' @param sc A SoupChannel object.
#' @param roundToInt Should the resulting matrix be rounded to integers?
#' @param verbose Should function be verbose?
#' @param tol Allowed deviation from expected number of soup counts.  Don't change this.
#' @param dropWholeCountsOnly Don't set this to TRUE unless you know what you're doing.  See details.
#' @param pCut The p-value cut-off used when \code{dropWholeCountsOnly=TRUE}.
#' @return A modified version of the table of counts, with background contamination removed.
#' @importFrom Matrix sparseMatrix
#' @importFrom stats rbinom
adjustCounts = function(sc,roundToInt=FALSE,verbose=FALSE,tol=1e-3,dropWholeCountsOnly=FALSE,pCut=0.01){
  if(!is(sc,'SoupChannel'))
    stop("sc must be an object of type SoupChannel")
  if(!'rho' %in% colnames(sc$metaData))
    stop("Contamination fractions must have already been calculated/set.")
  if(!dropWholeCountsOnly){
    #Create the final thing efficiently without making a big matrix
    out = as(sc$toc,'dgTMatrix')
    expSoupCnts = sc$metaData$nUMIs * sc$metaData$rho
    soupFrac = sc$soupProfile$est
    #Iteratively remove until we've removed the expected number of counts
    #How many counts should we have left when we're done?
    tgts = sc$metaData$nUMIs - expSoupCnts
    #Which genes do we still need to bother trying to remove counts from in which cells
    toAdjust = seq_along(out@i)
    if(verbose)
      message("Adjusting count matrix")
    while(TRUE){
      #How many left to do?
      toDo = colSums(out)-tgts
      #Get the soup frac correction factors.
      tmp = rep(1,length(toDo))
      soupFracSum = sapply(split(soupFrac[out@i[toAdjust]+1],out@j[toAdjust]+1),sum)
      tmp[as.numeric(names(soupFracSum))]=soupFracSum
      toDo = toDo/tmp
      #Do the adjustment
      out@x[toAdjust] = out@x[toAdjust]-soupFrac[out@i[toAdjust]+1]*toDo[out@j[toAdjust]+1]
      #Only keep updating those that need it
      toAdjust = toAdjust[out@x[toAdjust]>0]
      out@x[out@x<0]=0
      if(verbose)
        print(quantile(colSums(out)-tgts))
      if(max(colSums(out)-tgts)<tol)
        break
    }
    ##This is the clearer but slower version of above
    #out@x = pmax(0,out@x - soupFrac[out@i+1]*expSoupCnts[out@j+1])
    #while(max(colSums(out)-tgts)>tol){
    #  toDo = colSums(out)-tgts
    #  out@x = pmax(0,out@x - soupFrac[out@i+1]*toDo[out@j+1])
    #  print(quantile(colSums(out)-tgts))
    #}
  }else{
    #Retain this approach for backwards compatability
    #Start by calculating the p-value against the null of soup.
    out = as(sc$toc,'dgTMatrix')
    if(verbose)
      message("Calculating probability of each gene being soup")
    p = pbinom(out@x-1,sc$metaData$nUMIs[out@j+1],sc$soupProfile$est[out@i+1]*sc$metaData$rho[out@j+1],lower.tail=FALSE)
    #Order them by cell, then by p-value
    o = order(-(out@j+1),p,decreasing=TRUE)
    #Get the running total for removal.  Could probably make this faster with some tricks.
    if(verbose)
      message("Calculating probability of the next count being soup")
    s = split(o,out@j[o]+1)
    rTot = unlist(lapply(s,function(e) cumsum(out@x[e])),use.names=FALSE)
    #Now we need to get the soup probability vector as well.
    pSoup = pbinom(rTot-out@x[o]-1,sc$metaData$nUMIs[out@j[o]+1],sc$metaData$rho[out@j[o]+1],lower.tail=FALSE)
    if(verbose)
      message("Filtering table of counts")
    #Now construct the final probability vector
    pp = p[o]*pSoup
    #And we then want to drop anything meeting our termination criteria
    w = which(pp<pCut)
    #Keep a list of dropped genes
    dropped = data.frame(cell=colnames(out)[out@j[o[-w]]+1],
                         gene=rownames(out)[out@i[o[-w]]+1],
                         cnt = out@x[o[-w]])
    if(verbose){
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
  }
  #Do stochastic rounding to integers if needed
  if(roundToInt){
    if(verbose)
      message("Rounding to integers.")
    #Round to integer below, then probabilistically bump back up
    out@x = floor(out@x)+rbinom(length(out@x),1,out@x-floor(out@x))
  }
  #Fix the object internals
  w = which(out@x>0)
  out = sparseMatrix(i=out@i[w]+1,
                     j=out@j[w]+1,
                     x=out@x[w],
                     dims=dim(out),
                     dimnames=dimnames(out))
  return(out)
}

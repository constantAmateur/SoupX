#' Adjust counts to remove soup
#'
#' Explicitly drop counts from the count matrix in order of most likely to be soup to least likely until we have removed as many counts as we expect there to be soup molecules in each cell.  The actual removal criteria is to calculate the probability of each count being from the soup in each cell, pMol, and the probability of removing another soup molecule from this cell, pSoup, and keep sequentially removing counts until pSoup*pMol<pCut.  That is, pSoup makes sure we don't throw out many more than we expect to be present in each cell and pMol makes sure we don't throw out any molecule that is unlikely to be soup.
#'
#' @export
#' @param scl A SoupChannel or SoupChannelList object for which rho has been calculated.
#' @param pCut The threshold for excluding counts.
#' @param verbose By default, this method prints out the 100 most removed genes by channel.  If this is set to FALSE, these are not printed.
#' @seealso \code{\link{strainCells}}
#' @return A modified version of \code{scl} with the corrected count matrix stored in atoc.
#' @importFrom Matrix sparseMatrix
adjustCounts = function(scl,pCut=0.01,verbose=TRUE){
  #If we're just doing the one channel...
  if(is(scl,'SoupChannel')){
    scl = adjustCounts(SoupChannelList(list(scl)),pCut=pCut,verbose=verbose)
    sc = scl$channels[[1]]
    sc$atoc = scl$atoc
    return(sc)
  }
  scl$toc = as(scl$toc,'dgTMatrix')
  #Construct a global rho vector
  rhos = unlist(lapply(scl$channels,function(e) e$rhos),use.names=FALSE)
  names(rhos) = unlist(lapply(scl$channels,function(e) names(e$rhos)),use.names=FALSE)
  #Should be redundant, but do it anyway
  rhos = rhos[colnames(scl$toc)]
  #The fast version.  The logic here is that when we remove a count, the adjustment to the probability is given by:
  #p_new = 1-sum(dbinom(k,n-1,rho*f_gs)) 
  # = 1-(1-rho*f_gs)**-1 * sum((1-(k/n))*dbinom(k,n,rho*f_gs))
  #Which in the limit n>>k and rho*f_gs<<1, which are both usually good assumptions
  # = 1-sum(dbinom(k,n,rho*f_gs)) = p_old
  #As such, we should basically be fine to just calculate the probability vectors once and then traverse them until the termination criteria is met.
  #Which channel is each cell from
  if(verbose)
    message("Calculating probability of each gene being soup")
  cMap = match(gsub('___.*','',colnames(scl$toc)),colnames(scl$soupMatrix))
  p  = pbinom(scl$toc@x-1,scl$nUMIs[scl$toc@j+1],rhos[scl$toc@j+1]*scl$soupMatrix[cbind(scl$toc@i+1,cMap[scl$toc@j+1])],lower.tail=FALSE)
  #Order them by cell, then by p
  o = order(-(scl$toc@j+1),p,decreasing=TRUE)
  #Get the running total for removal.  Could probably make this faster with some tricks.
  if(verbose)
    message("Calculating probability of the next count being soup")
  s = split(o,scl$toc@j[o]+1)
  rTot = unlist(lapply(s,function(e) cumsum(scl$toc@x[e])),use.names=FALSE)
  #Now we need to get the soup probability vector as well.
  pSoup = pbinom(rTot-scl$toc@x[o]-1,scl$nUMIs[scl$toc@j[o]+1],rhos[scl$toc@j[o]+1],lower.tail=FALSE)
  if(verbose)
    message("Filtering table of counts")
  #Now construct the final probability vector
  pp = p[o]*pSoup
  #And we then want to drop anything meeting our termination criteria
  w = which(pp<pCut)
  #Keep a list of dropped genes
  dropped = data.frame(cell=colnames(scl$toc)[scl$toc@j[o[-w]]+1],
                       gene=rownames(scl$toc)[scl$toc@i[o[-w]]+1],
                       channel = colnames(scl$soupMatrix)[cMap[scl$toc@j[o[-w]]+1]])
  if(verbose){
    for(channel in colnames(scl$soupMatrix)){
      message(sprintf("Most removed genes for channel %s are:",channel))
      x = sort(table(dropped$gene[dropped$channel==channel])/sum(colnames(scl$soupMatrix)[cMap]==channel),decreasing=TRUE)
      print(x[seq_len(min(length(x),100))])
    }
  }
  #Construct the corrected count matrix
  scl$atoc = sparseMatrix(i=scl$toc@i[o[w]]+1,j=scl$toc@j[o[w]]+1,x=scl$toc@x[o[w]],dims=dim(scl$toc),dimnames = dimnames(scl$toc))
  return(scl)
  #repeat{
  #  message('Calculating probabilities')
  #  #Calculate everyone's probability at once
  #  p = pbinom(toc@x-1,nLeft[w[toc@j+1]],rhos[w[toc@j+1]]*sc$soupProfile$est[toc@i+1],lower.tail=FALSE)
  #  message('Getting best')
  #  o = order(toc@j+1,p,decreasing=TRUE)
  #  m = !duplicated((toc@j+1)[o])
  #  #Rather than do it a million times
  #  m = o[m]
  #  #Get the max probability by cell
  #  pMax = p[m]
  #  #And which gene it is
  #  gMax = (toc@i+1)[m]
  #  #And how many counts it is. Cheat a little here.  Because if we remove one count, it's almost always favourable to remove all the counts (as p increases hugely with each removal), just remove them all in one go.
  #  xMax = toc@x[m]
  #  #How likely are we to have another soup molecule for the remaining cells
  #  pSoup = pbinom(nCnts[w]-1,sc$nUMIs[w],rhos[w],lower.tail=FALSE)
  #  message('Updating counters')
  #  #And set the counts we just removed to zero
  #  toc@i = toc@i[-m]
  #  toc@j = toc@j[-m]
  #  toc@x = toc@x[-m]
  #  #Have we finished any of them?
  #  toDrop = which(pMax*pSoup<pCut)
  #  #Check for termination criteria
  #  if(length(toDrop)==length(w))
  #    break
  #  if(length(toDrop)>0){
  #    toc = toc[,-toDrop]
  #    w = w[-toDrop]
  #    pMax = pMax[-toDrop]
  #    gMax = gMax[-toDrop]
  #    xMax = xMax[-toDrop]
  #    pSoup = pSoup[-toDrop]
  #  }
  #  #Store the newly dropped ones
  #  cnts$i = c(cnts$i,gMax)
  #  cnts$j = c(cnts$j,w)
  #  cnts$x = c(cnts$x,xMax)
  #  #And update counters
  #  nCnts[w] = nCnts[w] + xMax
  #  nLeft[w] = nLeft[w] - xMax
  #  print(toDrop)
  #  print(quantile(nLeft))
  #  print(quantile(nCnts))
  #  print(quantile(rhos*sc$nUMIs-nCnts))
  #  print(quantile((pMax*pSoup)))
  #}
}


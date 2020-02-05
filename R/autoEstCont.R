#' Automatically calculate the contamination fraction
#'
#' The idea of this method is that genes that are highly expressed in the soup and are marker genes for some population can be used to estimate the background contamination.  Marker genes are identified using the tfidf method (see \code{\link{quickMarkers}}).  The contamination fraction is then calculated at the cluster level for each of these genes and clusters are then aggressively pruned to remove those that give implausible estimates.  
#'
#' The pruning of implausible clusters is based on several parameters.  Firstly, clusters are excluded for estimation on a gene-by-gene basis using \code{\link{estimateNonExpressingCells}} as would be done had the genes been specified by the user.  Any cluster whose estimation for rho is not significantly less than \code{rhoMax} (as measured by a Poisson test) is also thrown out (with false discovery rate \code{FDR}).
#'
#' Finally, clusters are excluded if they number of observed counts is less than \code{minCnts}.  The default (where \code{minCnts=1}) essentially prevents an excess of estimates of zero contamination skewing the estimation.
#' 
#' An estimate of the global contamination rate is then provided based on the value of rho with the highest density in the estimates using the pruned cluster/gene pairs.  The logic behind this is that  clusters that truly estimate the contamination fraction will cluster around the true value, while erroneous estimates will be spread out across the range (0,1) without a 'preferred value'.
#' 
#' Note that the interpretation of the tf-idf cut-off in this context is that a cut-off t implies that a marker gene has the property that geneFreqGlobal < exp(-t/geneFreqInClust).
#'
#' @export
#' @param sc The SoupChannel object.
#' @param nMarks How many marker genes to use. These must also pass the tfidf cut-off and be highly expressed in the soup, so there may end up being fewer than this number.
#' @param tfidfMin Minimum value of tfidf to accept for a marker gene.
#' @param soupMin Minimum value of expression in the background to accept for marker.
#' @param rhoMax What contamination fraction is the maximum that is plausible.  Must be a value between 0 and 1.
#' @param FDR False discovery rate for statistical tests.
#' @param minCnts Don't use estimates with fewer than this many counts.
#' @param doPlot Create a plot summarising the estimation?
#' @param ... Passed to density.
autoEstCont = function(sc,nMarks = 50,tfidfMin=0.8,soupMin=1e-4,rhoMax=1.0,FDR=0.1,minCnts=1,doPlot=TRUE,...){
  #First collapse by cluster
  s = split(rownames(sc$metaData),sc$metaData$clusters)
  tmp = do.call(cbind,lapply(s,function(e) rowSums(sc$toc[,e,drop=FALSE])))
  ssc = sc 
  ssc$toc = tmp
  ssc$metaData = data.frame(nUMIs = colSums(tmp),row.names=colnames(tmp))
  ###################
  # Get best markers
  #Get the top N soup Genes
  tgts = ssc$soupProfile[order(ssc$soupProfile$est,decreasing=TRUE),]
  tgts = rownames(tgts)[tgts$est>soupMin]
  #Refine this to the best markers we can manage
  mrks = quickMarkers(sc$toc,sc$metaData$clusters,N=Inf)
  #Keep only the ones that are high in soup
  mrks = mrks[mrks$gene %in% tgts,]
  #And only the most specific entry for each gene
  mrks = mrks[order(mrks$gene,-mrks$tfidf),]
  mrks = mrks[!duplicated(mrks$gene),]
  #Order by tfidif maxness
  mrks = mrks[order(-mrks$tfidf),]
  #Apply tf-idf cut-off
  mrks = mrks[mrks$tfidf > tfidfMin,]
  tgts = head(mrks$gene,nMarks)
  if(length(tgts)==0)
    stop("No plausible marker genes found.  Try reducing tfidfMin.")
  if(length(tgts)<10)
    warning("Fewer than 10 marker genes found.  Consider reducing tfidfMin")
  ############################
  # Get estimates in clusters
  #Get which ones we'd use and where with canonical method
  tmp = as.list(tgts)
  names(tmp) = tgts
  ute = estimateNonExpressingCells(sc,tmp,maximumContamination=rhoMax,FDR=FDR)
  m = rownames(sc$metaData)[match(rownames(ssc$metaData),sc$metaData$clusters)]
  ute = t(ute[m,])
  colnames(ute) = rownames(ssc$metaData)
  #Now calculate the observed and expected counts for each cluster for 
  expCnts = outer(ssc$soupProfile$est,ssc$metaData$nUMIs)
  rownames(expCnts) = rownames(ssc$soupProfile)
  colnames(expCnts) = rownames(ssc$metaData)
  expCnts = expCnts[tgts,]
  #And the observed ones
  obsCnts = ssc$toc[tgts,]
  #Filter out the shite
  #Get the p-value for this being less than 1
  pp = ppois(obsCnts,expCnts*rhoMax,lower.tail=TRUE)
  qq = p.adjust(pp,method='BH')
  qq = matrix(qq,nrow=nrow(pp),ncol=ncol(pp),dimnames=dimnames(pp))
  #Get the cluster level ratio
  rhos = obsCnts/expCnts
  #Index in range
  rhoIdx = t(apply(rhos,1,function(e) order(order(e))))
  #Make a data.frame with everything
  dd = data.frame(gene = rep(rownames(ute),ncol(ute)),
                  passNonExp = as.vector(ute),
                  rho = as.vector(rhos),
                  rhoIdx = as.vector(rhoIdx),
                  obsCnt = as.vector(obsCnts),
                  expCnt = as.vector(expCnts),
                  qVal = as.vector(qq)
                  )
  dd$pass = dd$obsCnt >= minCnts & 
    dd$qVal < FDR & 
    #dd$rhoIdx <= min(clustPerGene,floor(ncol(rhoIdx)*maxClustFrac)) & 
    dd$passNonExp
  #I think the best way to do this is based on the density.
  tmp = density(dd$rho[dd$pass],...)
  dd$rhoEst = tmp$x[which.max(tmp$y)]
  if(doPlot){
    plot(0,0,
         xlim=c(0,1),
         ylim=c(0,max(tmp$y)),
         type='n',
         frame.plot=FALSE,
         xlab='Contamination Fraction',
         ylab='Density'
         )
    lines(tmp$x,tmp$y)
    abline(v=dd$rhoEst[1])
  }
  return(dd)
}

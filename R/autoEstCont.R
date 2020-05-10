#' Automatically calculate the contamination fraction
#'
#' The idea of this method is that genes that are highly expressed in the soup and are marker genes for some population can be used to estimate the background contamination.  Marker genes are identified using the tfidf method (see \code{\link{quickMarkers}}).  The contamination fraction is then calculated at the cluster level for each of these genes and clusters are then aggressively pruned to remove those that give implausible estimates.  
#'
#' This set of marker genes is filtered to include only those with tf-idf value greater than \code{tfidfMin}.  A higher tf-idf value implies a more specific marker.  Specifically a cut-off t implies that a marker gene has the property that geneFreqGlobal < exp(-t/geneFreqInClust).  See \code{\link{quickMarkers}}.  It may be necessary to decrease this value for data sets with few good markers.
#'
#' This set of marker genes is filtered down to include only the genes that are highly expressed in the soup, controlled by the \code{soupQuantile} parameter.  Genes highly expressed in the soup provide a more precise estimate of the contamination fraction.
#'
#' The pruning of implausible clusters is based on a call to \code{\link{estimateNonExpressingCells}}.  The parameters \code{maximumContamination} and \code{rhoMaxFDR} are passed to this function.  The defaults set here are calibrated to aggressively prune anything that has even the weakest of evidence that it is genuinely expressed. 
#'
#' For each cluster/gene pair the posterior distribution of the contamination fraction is calculated (based on gamma prior, controlled by \code{priorRho} and \code{priorRhoStdDev}).  These posterior distributions are aggregated to produce a final estimate of the contamination fraction. The logic behind this is that estimates from clusters that truly estimate the contamination fraction will cluster around the true value, while erroneous estimates will be spread out across the range (0,1) without a 'preferred value'.  The most probable value of the contamination fraction is then taken as the final global contamination fraction.
#' 
#' @export
#' @param sc The SoupChannel object.
#' @param topMarkers A data.frame giving marker genes.  Must be sorted by decreasing specificity of marker and include a column 'gene' that contains the gene name.  If set to NULL, markers are estimated using \code{\link{quickMarkers}}.
#' @param tfidfMin Minimum value of tfidf to accept for a marker gene.
#' @param soupQuantile Only use genes that are at or above this expression quantile in the soup.  This prevents inaccurate estimates due to using genes with poorly constrained contribution to the background.
#' @param maxMarkers If we have heaps of good markers, keep only the best \code{maxMarkers} of them.
#' @param maximumContamination What contamination fraction is the maximum that is plausible.  Must be a value between 0 and 1.  Passed to \code{\link{estimateNonExpressingCells}}.
#' @param rhoMaxFDR False discovery rate passed to \code{\link{estimateNonExpressingCells}}, to test if rho is less than \code{maximumContamination}.
#' @param priorRho Mode of gamma distribution prior on contamination fraction.
#' @param priorRhoStdDev Standard deviation of gamma distribution prior on contamination fraction.
#' @param doPlot Create a plot showing the density of estimates?
#' @param forceAccept Passed to \code{\link{setContaminationFraction}}.
#' @param verbose Be verbose?
#' @seealso quickMarkers
#' @return A modified SoupChannel object where the global contamination rate has been set.  Information about the estimation is also stored in the slot \code{fit}
#' @examples
#' #Use less specific markers
#' scToy = autoEstCont(scToy,tfidfMin=0.8)
#' #Allow large contamination fractions to be allocated
#' scToy = autoEstCont(scToy,forceAccept=TRUE)
#' #Be quiet
#' scToy = autoEstCont(scToy,verbose=FALSE,doPlot=FALSE)
#' @importFrom stats dgamma qgamma 
#' @importFrom graphics abline lines legend plot
autoEstCont = function(sc,topMarkers=NULL,tfidfMin=1.0,soupQuantile=0.90,maxMarkers=100,maximumContamination=0.8,rhoMaxFDR=0.2,priorRho=0.05,priorRhoStdDev=0.10,doPlot=TRUE,forceAccept=FALSE,verbose=TRUE){
  if(!'clusters' %in% colnames(sc$metaData))
    stop("Clustering information must be supplied, run setClusters first.")
  #First collapse by cluster
  s = split(rownames(sc$metaData),sc$metaData$clusters)
  tmp = do.call(cbind,lapply(s,function(e) rowSums(sc$toc[,e,drop=FALSE])))
  ssc = sc 
  ssc$toc = tmp
  ssc$metaData = data.frame(nUMIs = colSums(tmp),row.names=colnames(tmp))
  ###################
  # Get best markers
  #Get the top N soup Genes
  soupProf = ssc$soupProfile[order(ssc$soupProfile$est,decreasing=TRUE),]
  soupMin = quantile(soupProf$est,soupQuantile)
  #Find or load markers.
  if(is.null(topMarkers)){
    #Refine this to the best markers we can manage
    mrks = quickMarkers(sc$toc,sc$metaData$clusters,N=Inf)
    #And only the most specific entry for each gene
    mrks = mrks[order(mrks$gene,-mrks$tfidf),]
    mrks = mrks[!duplicated(mrks$gene),]
    #Order by tfidif maxness
    mrks = mrks[order(-mrks$tfidf),]
    #Apply tf-idf cut-off
    mrks = mrks[mrks$tfidf > tfidfMin,]
  }else{
    mrks = topMarkers
  }
  #Filter to include only those that exist in soup 
  tgts = rownames(soupProf)[soupProf$est>soupMin]
  #And get the ones that pass our tfidf cut-off
  filtPass = mrks[mrks$gene %in% tgts,]
  tgts = head(filtPass$gene,n=maxMarkers)
  if(verbose)
    message(sprintf("%d genes passed tf-idf cut-off and %d soup quantile filter.  Taking the top %d.",nrow(mrks),nrow(filtPass),length(tgts)))
  #mrks = mrks[mrks$gene %in% tgts,]
  #tgts = head(mrks$gene,nMarks)
  if(length(tgts)==0){
    stop("No plausible marker genes found.  Reduce tfidfMin or soupQuantile")
  }
  if(length(tgts)<10){
    warning("Fewer than 10 marker genes found.  Consider reducing tfidfMin or soupQuantile")
  }
  ############################
  # Get estimates in clusters
  #Get which ones we'd use and where with canonical method
  tmp = as.list(tgts)
  names(tmp) = tgts
  ute = estimateNonExpressingCells(sc,tmp,maximumContamination=maximumContamination,FDR=rhoMaxFDR)
  m = rownames(sc$metaData)[match(rownames(ssc$metaData),sc$metaData$clusters)]
  ute = t(ute[m,,drop=FALSE])
  colnames(ute) = rownames(ssc$metaData)
  #Now calculate the observed and expected counts for each cluster for 
  expCnts = outer(ssc$soupProfile$est,ssc$metaData$nUMIs)
  rownames(expCnts) = rownames(ssc$soupProfile)
  colnames(expCnts) = rownames(ssc$metaData)
  expCnts = expCnts[tgts,,drop=FALSE]
  #And the observed ones
  obsCnts = ssc$toc[tgts,,drop=FALSE]
  #We're done, but record some extra data for fun and profit
  #Filter out the shite
  #Get the p-value for this being less than 1
  pp = ppois(obsCnts,expCnts*maximumContamination,lower.tail=TRUE)
  qq = p.adjust(pp,method='BH')
  qq = matrix(qq,nrow=nrow(pp),ncol=ncol(pp),dimnames=dimnames(pp))
  #Get the cluster level ratio
  rhos = obsCnts/expCnts
  #Index in range
  rhoIdx = t(apply(rhos,1,function(e) order(order(e))))
  #Make a data.frame with everything
  dd = data.frame(gene = rep(rownames(ute),ncol(ute)),
                  passNonExp = as.vector(ute),
                  rhoEst = as.vector(rhos),
                  rhoIdx = as.vector(rhoIdx),
                  obsCnt = as.vector(obsCnts),
                  expCnt = as.vector(expCnts),
                  isExpressedFDR = as.vector(qq)
                  )
  dd$geneIdx = match(dd$gene,mrks$gene)
  dd$tfidf = mrks$tfidf[dd$geneIdx]
  dd$soupIdx = match(dd$gene,rownames(soupProf))
  dd$soupExp = soupProf$est[dd$soupIdx]
  dd$useEst = #dd$obsCnt >= minCnts & 
    #dd$isExpressedFDR < rhoMaxFDR & 
    #dd$rhoIdx <= min(clustPerGene,floor(ncol(rhoIdx)*maxClustFrac)) & 
    dd$passNonExp
  #The logic of piling up desity around the true value gets wonky if the number of estimates is low
  if(sum(dd$useEst)<10)
    warning("Fewer than 10 independent estimates, rho estimation is likely to be unstable.  Consider reducing tfidfMin or increasing SoupMin.")
  if(verbose)
    message(sprintf("Using %d independent estimates of rho.",sum(dd$useEst)))
  #Now aggregate the posterior probabilities for the ones we're including
  p.L = function(x,alpha){if(x==0){0}else{qgamma(alpha,x)}}
  p.U = function(x,alpha){qgamma(1-alpha,x+1)}
  alpha=0.95
  alpha=(1-alpha)/2
  dd$rhoHigh=sapply(seq(nrow(dd)),function(e) p.U(dd$obsCnt[e],alpha)/dd$expCnt[e])
  dd$rhoLow=sapply(seq(nrow(dd)),function(e) p.L(dd$obsCnt[e],alpha)/dd$expCnt[e])
  rhoProbes=seq(0,1,.001)
  #Using 95% confidence intervals
  #tmp = sapply(rhoProbes,function(e) {w=which(dd$useEst & dd$tfidf<1.5);sum(e>=dd$rhoLow[w] & e<=dd$rhoHigh[w])/length(w)})
  #Do a posterior estimation instead.  Use gamma prior defined by mode (priorRho) and standard deviation (priorRhoStdDev), which yields a posterior distribution for gamma of the form dgamma(rho,obsCnt+k,scale=theta/(1+theta*expCnts)). Where k and theta are the parameters for prior distribution derived using the above constraints.
  v2 = (priorRhoStdDev/priorRho)**2
  k = 1 +v2**-2/2*(1+sqrt(1+4*v2))
  theta = priorRho/(k-1)
  tmp = sapply(rhoProbes,function(e) {
                 tmp = dd[dd$useEst,]
                 mean(dgamma(e,k+tmp$obsCnt,scale=theta/(1+theta*tmp$expCnt)))
                  })
  #Calculate prior curve
  xx=dgamma(rhoProbes,k,scale=theta)
  #Get estimates
  rhoEst = rhoProbes[which.max(tmp)]
  rhoFWHM = range(rhoProbes[which(tmp>=(max(tmp)/2))])
  contEst = rhoEst
  if(verbose)
    message(sprintf("Estimated global rho of %.2f",rhoEst))
  ##I think the best way to do this is based on the density.
  #tmp = density(dd$rhoEst[dd$useEst],...)
  #contEst = tmp$x[which.max(tmp$y)]
  if(doPlot){
    plot(rhoProbes,tmp,'l',
         xlim=c(0,1),
         ylim=c(0,max(c(xx,tmp))),
         frame.plot=FALSE,
         xlab='Contamination Fraction',
         ylab='Probability Density')
    #Add prior
    lines(rhoProbes,xx,lty=2)
    abline(v=rhoProbes[which.max(tmp)],col='red')
    legend(x='topright',
           legend=c(sprintf('prior rho %g(+/-%g)',priorRho,priorRhoStdDev),
                    sprintf('post rho %g(%g,%g)',rhoEst,rhoFWHM[1],rhoFWHM[2]),
                    'rho max'),
           lty=c(2,1,1),
           col=c('black','black','red'),
           bty='n')
    #plot(0,
    #     xlim=c(0,1),
    #     ylim=c(0,max(tmp$y)),
    #     type='n',
    #     frame.plot=FALSE,
    #     xlab='Contamination Fraction',
    #     ylab='Density'
    #     )
    #lines(tmp$x,tmp$y)
    #abline(v=contEst,col='red')
  }
  sc$fit = list(dd=dd,
                priorRho=priorRho,
                priorRhoStdDev=priorRhoStdDev,
                posterior = tmp,
                rhoEst = rhoEst,
                rhoFWHM = rhoFWHM
                )
  #Set the contamination fraction
  sc = setContaminationFraction(sc,contEst,forceAccept=forceAccept)
  return(sc)
}

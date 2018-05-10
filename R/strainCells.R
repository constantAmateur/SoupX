#' Remove background contamination from cell expression profile
#'
#' Given a table of counts composed of cells from multiple channels, estimates of the background contamination in each cell and the soup expression profile for each channel, calculates the true cell specific expression profile for each cell or a set of cells.  If background correction is being performed on cells all from the same channel, \code{soupProfiles} can be the output of \code{estimateSoup} for a single channel.
#' 
#' @export
#' @param scl SoupChannelList object.
#' @param groupVec Vector named by unique droplet IDs of the format channelName___dropletBarcode that are unique across channels, with the entries indicating the name of the group each cell belongs to.  Cells in the same group are given an aggregated expression profile after background correction.  If NULL, all cells are corrected individually.  Any unmentioned cells are also corrected individually.
#' @param retainSparsity The MLE optimisation is much faster for non-sparse matrices when \code{groupVec} is not NULL.  By default the code will coerce a sparse matrix to be non-sparse in this case unless this flag is set to TRUE.
#' @seealso \code{\link{adjustCnts}}
#' @return A modified version of \code{scl} that has an extra matrix, scl$strainedExp, that contains the corrected expression matrix.
strainCells = function(scl,groupVec=NULL,retainSparsity=FALSE){
  if(!is(scl,'SoupChannelList'))
    stop("scl must be an object of type SoupChannelList")
  if(any(sapply(scl$channels,function(e) is.null(e$rhos))))
    stop("Cell specific contamination fractions must be calculated for all channels.")
  #Get the big global list of rhos
  rhos = unlist(lapply(scl$channels,function(e) e$rhos),use.names=FALSE)
  names(rhos) = unlist(lapply(scl$channels,function(e) names(e$rhos)),use.names=FALSE)
  rhos = rhos[colnames(scl$toc)]
  #This is a useful thing to use
  cMap = gsub('___.*','',colnames(scl$toc))
  #If we're just doing single cells, the estimation is easy
  if(is.null(groupVec)){
    fgc = t(t(scl$toc)/scl$nUMIs)
    fgc = as(fgc,'dgTMatrix')
    fgc@x = (fgc@x - scl$soupMatrix[cbind(fgc@i+1,match(cMap[fgc@j+1],colnames(scl$soupMatrix)))]*rhos[fgc@j+1])/(1-rhos[fgc@j+1])
    #fgc = t((t(scl$toc)/scl$nUMIs-t(scl$soupMatrix[,cMap,drop=FALSE])*rhos)/(1-rhos))
    #fgc = t((t(toc)/nUMIs-t(soupProfiles[,channelMap])*rhos)/(1-rhos))
  }else{
    if(!all(names(groupVec) %in% colnames(scl$toc)))
      stop("groupVec must be named by unique droplet identifier.")
    #Create complete groupVec
    groupVec = factor(groupVec)
    groupNoms = levels(groupVec)
    groupVec = setNames(as.numeric(groupVec),groupNoms)
    #Add in non-mentioned ones
    w = colnames(scl$toc)[!(colnames(scl$toc)%in%names(groupVec))]
    groupVec[w] = seq_len(w)+max(groupVec)+1
    groupNoms = c(groupNoms,w)
    groupVec = setNames(groupNoms[groupVec],names(groupvec))
    #Re-order as toc
    groupVec = groupVec[colnames(scl$toc)]
    #Create the output matrix
    fgc = matrix(0,nrow=nrow(scl$toc),ncol=length(unique(groupVec)),dimnames=list(rownames(scl$toc),unique(groupVec)))
    #Work out which ones are singletons and do them quickly
    w = (!(groupVec %in% unique(groupVec[duplicated(groupVec)])))
    fgc[,w,drop=FALSE] = t((t(scl$toc[,w,drop=FALSE])/scl$nUMIs[w]-t(scl$soupMatrix[,cMap[w],drop=FALSE])*rhos[w])/(1-rhos[w]))
    #No closed form for the MLE in this case, have to explicitly optimise
    #To make this time efficient, we have to sacrifice sparsity.  If this isn't an option be prepared to wait for a long time.
    toc = scl$toc
    if(is(toc,'sparseMatrix')){
      if(retainSparsity){
        message("Retaining sparse matrix to preserve memory.  The resulting calculation will be very slow.")
      }else{
        toc = as.matrix(toc)
      }
    }
    grpNoms = unique(groupVec[!w])
    for(grpNom in grpNoms){
      message(sprintf('Calculating cell expression profile via MLE for group of cells %s',grpNom))
      w = which(grpNom==groupVec)
      est = rowSums(toc[,w])/sum(scl$nUMIs[w])
      ww = rowSums(toc[,w])==0
      pb = txtProgressBar(min=1,max=sum(!ww),initial=1,style=3)
      cnt=1
      for(i in which(!ww)){
        #Pre-calculate for speed
        dat = toc[i,w]
        rs = rhos[w]
        tots = scl$nUMIs[w]
        soups = soupMatrix[i,cMap[w]]
        A = tots*rhos*soups
        B = tots*(1-rhos)
        optFun = function(par) { -1*sum(dpois(dat,A+B*par,log=TRUE))
        }
        #grr = function(par){
        #  sum((1-rhos)*tots*(1-dat/(1e-8 + tots*(rhos*soups+(1-rhos)*par))))
        #}
        fit = optimise(optFun,c(0,1),tol=1e-10)
        fgc[i,grpNom]=fit$minimum
        cnt=cnt+1
        setTxtProgressBar(pb,cnt)
        #print(cnt)
      }
      close(pb)
    }
  }
  #Zero anything that's too low
  fgc[fgc<0]=0
  #Re-normalise 
  fgc = t(t(fgc)/colSums(fgc))
  scl$strainedExp = fgc
  return(scl)
}

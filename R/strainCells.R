#' Remove background contamination from cell expression profile
#'
#' Given a table of counts composed of cells from multiple channels, estimates of the background contamination in each cell and the soup expression profile for each channel, calculates the true cell specific expression profile for each cell or a set of cells.  If background correction is being performed on cells all from the same channel, \code{soupProfiles} can be the output of \code{estimateSoup} for a single channel.
#' 
#' @export
#' @param toc Table of UMIs for all droplets containing cells.
#' @param rhos Vector giving contamination fraction for each cell.
#' @param soupProfiles Matrix giving the soup expression profile for each channel.  rows are genes and columns are channels, named as in \code{channelMap}.  If only one channel is used (i.e., if \code{channelMap} is NULL) can be the output of \code{\link{estimateSoup}}.
#' @param channelMap Vector giving the name of the channel from which this cell is derived.  Must correspond to a column in \code{soupProfiles}.  If NULL, cells are assumed to all be from the one channel.
#' @param groupVec Vector indicating how to group cells when calculating expression profile.  If NULL, a profile for each cell is calculated.
#' @param zeroCut Set any estimated value less than this to zero.
#' @param retainSparsity The MLE optimisation is much faster for non-sparse matrices when \code{groupVec} is not NULL.  By default the code will coerce a sparse matrix to be non-sparse in this case unless this flag is set to TRUE.
#' @return A matrix containing background corrected expression profiles for each group specified.
strainCells = function(toc,rhos,soupProfiles,channelMap=NULL,groupVec=NULL,zeroCut=1/sum(toc),retainSparsity=FALSE){
  nUMIs = colSums(toc)
  #Handle the special case of just one channel.
  if(is.null(channelMap)){
    channelMap = rep('Channel',ncol(toc))
    if(is.null(ncol(soupProfiles)))
      soupProfiles = matrix(soupProfiles,ncol=1,dimnames=list(names(soupProfiles),'Channel'))
    if(ncol(soupProfiles)>1)
      soupProfiles = matrix(soupProfiles$est,ncol=1,dimnames=list(rownames(soupProfiles),'Channel'))
  }
  #Cell level estimation?
  if(is.null(groupVec)){
    #MLE is simple in this case...
    fgc = t((t(toc)/nUMIs-t(soupProfiles[,channelMap])*rhos)/(1-rhos))
  }else{
    #No closed form for the MLE in this case, have to explicitly optimise
    if(!is.factor(groupVec))
      groupVec = factor(groupVec)
    #Initialise estimation matrix
    fgc = matrix(0,nrow=nrow(toc),ncol=length(levels(clMap)),dimnames=list(rownames(toc),levels(groupVec)))
    #To make this time efficient, we have to sacrifice sparsity.  If this isn't an option be prepared to wait for a long time.
    if(is(toc,'sparseMatrix')){
      if(retainSparsity){
        message("Retaining sparse matrix to preserve memory.  The resulting calculation will be very slow.")
      }else{
        toc = as.matrix(toc)
      }
    }
    for(grpNom in levels(groupVec)){
      message(sprintf('Calculating cell expression profile via MLE for group of cells %s',grpNom))
      w = which(grpNom==groupVec)
      est = rowSums(toc[,w])/sum(nUMIs[w])
      ww = rowSums(toc[,w])==0
      pb = txtProgressBar(min=1,max=sum(!ww),initial=1,style=3)
      cnt=1
      for(i in which(!ww)){
        #Pre-calculate for speed
        dat = toc[i,w]
        rs = rhos[w]
        tots = nUMIs[w]
        soups = soupProfiles[i,channelMap[w]]
        A = tots*rhos*soups
        B = tots*(1-rhos)
        optFun = function(par) { -1*sum(dpois(dat,A+B*par,log=TRUE))
        }
        #grr = function(par){
        #  sum((1-rhos)*tots*(1-dat/(1e-8 + tots*(rhos*soups+(1-rhos)*par))))
        #}
        fit = optimise(optFun,c(0,1),tol=1e-10)
        fgc[i,clNom]=fit$minimum
        cnt=cnt+1
        setTxtProgressBar(pb,cnt)
        #print(cnt)
      }
      close(pb)
    }
  }
  #Zero anything that's too low
  fgc[fgc<zeroCut]=0
  #Re-normalise 
  fgc = t(t(fgc)/colSums(fgc))
  return(fgc)
}

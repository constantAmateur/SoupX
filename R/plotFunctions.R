#' Plot contamination fraction for channel
#'
#' Plots the contamination fraction as a function of nUMIs in a droplet.  Global estimate is shown in red and lowess curve (usually used for interpolation) in green.
#'
#' @export
#' @param sc A \code{SoupChannel} object on which \code{\link{calculateContaminationFraction}} has been run.  Alternatively, a \code{SoupChannelList} object containing such a channel.
#' @param channelName The name of the channel to use if \code{sc} is a \code{SoupChannelList}
#' @param showErrorBars Should error bars showing the 95 percent confidence interval on the estimate of rho be included?
#' @param ... Extra parameters passed to lowess smoother.
#' @return A ggplot2 object containing the plot.
plotChannelContamination = function(sc,channelName,showErrorBars=TRUE,...){
  if(is(sc,'SoupChannelList'))
    sc = sc$channels[[channelName]]
  if(!is(sc,'SoupChannel'))
    stop("sc not a valid SoupChannel or SoupChannelList object.")
  if(is.null(sc$rhoGrouped))
    stop("calculateContaminationFraction not run on channel.")
  rhos = sc$rhoGrouped
  #Separate out global estimate
  globRho = rhos['Global',]
  rhos = rhos[rownames(rhos)!='Global',]
  gg = ggplot(rhos,aes(log10(nUMIs),est)) +
    geom_point() +
    geom_line()
  if(showErrorBars){
    gg = gg +
      geom_errorbar(aes(ymin=lower,ymax=upper))
  }
  gg = gg +
    geom_hline(aes(yintercept=ifelse(isLogged,log10(globRho$est),globRho$est)),linetype=1,colour='red') +
    geom_hline(aes(yintercept=ifelse(isLogged,log10(globRho$lower),globRho$lower)),linetype=2,colour='red') +
    geom_hline(aes(yintercept=ifelse(isLogged,log10(globRho$upper),globRho$upper)),linetype=2,colour='red') +
    ylab('Contamination Fraction') +
    xlab('log10(nUMIs/cell)')
  x = rhos$nUMIs
  y = rhos$est
  w = is.finite(y) & !is.na(y) & x>0
  x = x[w]
  y = y[w]
  fit = lowess(y~x,...)
  gg = gg + geom_line(data=data.frame(x=fit$x,y=fit$y),
                      aes(log10(x),y),colour='green')
  gg$df = rhos
  return(gg)
}

#' Plot correlation of expression profiles of soup and aggregated cells
#' 
#' Calculates an expression profile by aggregating counts across all cells and plots this (on a log10 scale) against the expression profile of the soup.
#'
#' @export
#' @param sc A SoupChannel or SoupChannelList object.
#' @param channelName The name of the channel to use if \code{sc} is a \code{SoupChannelList}
#' @return A ggplot2 object containing the plot.
plotSoupCorrelation = function(sc,channelName){
  if(is(sc,'SoupChannelList'))
    sc = sc$channels[[channelName]]
  if(!is(sc,'SoupChannel'))
    stop("sc not a valid SoupChannel or SoupChannelList object.")
  #Calculate the cell profile
  cellProfile = rowSums(sc$toc)
  cellProfile = (cellProfile/sum(cellProfile))
  df = as.data.frame(cbind(cellProfile,soupProfile=sc$soupProfile$est))
  gg = ggplot(df,aes(log10(cellProfile),log10(sc$soupProfile))) +
    geom_point(alpha=1/3) +
    geom_abline(intercept=0,slope=1) +
    ylab('log10(Soup Expression)')+
    xlab('log10(Aggregate cell Expression)')
  gg$df = df
  return(gg)
}

#' Plot non-expressed genes for channel
#' 
#' A helper function to plot the top N non-expressed genes for a channel.
#' 
#' @export
#' @param sc A SoupChannel or SoupChannelList object.
#' @param channelName The name of the channel to plot.  Ignored if \code{sc} is of class \code{SoupChannel}.
#' @param n Number of candidate genes to plot.
#' @param ... Parameters passed to \code{\link{plotMarkerDistribution}}.
#' @return A ggplot2 object containing the plot.
plotCandidateMarkerGenes = function(sc,channel,n=20,...){
  if(is(sc,'SoupChannelList')){
    plotCandidateMarkerGenes(sc$channels[[channel]],n=n,...)
  }else if(is(sc,'SoupChannel')){
    if(is.null(sc$nonExpressedGenes))
      stop("Run inferNonExpressedGenes first!")
    plotMarkerDistribution(sc,nonExpressedGeneList = rownames(sc$nonExpressedGenes)[seq(n)],...)
  }else{
    stop("sc must be a SoupChannel or SoupChannelList object")
  }
}


#' Plots the distribution of the observed to expected expression for marker genes
#'
#' If each cell were made up purely of background reads, the expression fraction would equal that of the soup.  This plot compares this expectation of pure background to the observed expression fraction in each cell, for each of the groups of genes in \code{nonExpressedGeneList}.  For each group of genes, the distribution of this ratio is plotted across all cells.  A value significantly greater than 1 (0 on log scale) can only be obtained if some of the genes in each group are genuinely expressed by the cell.  
#'
#' This plot is a useful diagnostic for the assumption that a list of genes is non-expressed in most cell types.  For non-expressed cells, the ratio should cluster around the contamination fraction, while for expressed cells it should be elevated.  The most useful non-expressed gene sets are those for which the genes are either strongly expressed, or not expressed at all.  Such groups of genes will show up in this plot as a bimodal distribution, with one mode containing the cells that do not express these genes around the contamination fraction for this channel and another around a value at some value equal to or greater than 0 (1 on non-log scale) for the expressed cells.
#'
#' The red line shows the global estimate of rho for each group of markers.  This is usually lower than the low mode of the distribution as there will typically be a non-negligible number of cells with 0 observed counts (and hence -infinity log ratio).
#'
#' @export
#' @param sc A SoupChannel or SoupChannelList object.
#' @param channelName The name of the channel to use if \code{sc} is a \code{SoupChannelList}
#' @param nonExpressedGeneList Which genes to use to estimate soup (see \code{\link{calculateContaminationFraction}}).  Can also be a vector of gene names, in which case each gene is treated as a separate gene set.
#' @param maxCells Randomly plot only this many cells to prevent over-crowding.
#' @param minRho A statistical test is performed to identify all those cells with more expression than would be possible from the soup alone.  This is done by testing against a null hypothesis of all expression coming from the soup, with the soup contamination fraction being less than or equal to \code{minRho}.
#' @param qCut The FDR cut-off for the hypothesis test described above.
#' @importFrom stats setNames
#' @return A ggplot2 object containing the plot.
plotMarkerDistribution = function(sc,channelName,nonExpressedGeneList,maxCells=150,minRho=1.0,qCut=0.05){
  if(is(sc,'SoupChannelList'))
    sc = sc$channels[[channelName]]
  if(!is(sc,'SoupChannel'))
    stop("sc not a valid SoupChannel or SoupChannelList object.")
  #Get nonExpressedGeneList from sc if missing
  if(missing(nonExpressedGeneList))
    nonExpressedGeneList = rownames(sc$nonExpressedGenes)[seq(min(20,nrow(sc$nonExpressedGenes)))]
  #Make non-lists into lists
  if(!is.list(nonExpressedGeneList))
    nonExpressedGeneList = as.list(setNames(nonExpressedGeneList,nonExpressedGeneList))
  #Calculate the ratio to the soup for each marker group in each cell
  obsProfile = t(t(sc$toc)/sc$nUMIs)
  #Get the ratio
  tst = lapply(nonExpressedGeneList,function(e) colSums(obsProfile[e,,drop=FALSE])/sum(sc$soupProfile[e,'est']))
  #Unlist the thing
  df = data.frame(MarkerGroup = rep(names(tst),lengths(tst)),
                  Barcode=unlist(lapply(tst,names),use.names=FALSE),
                  Values=unlist(tst,use.names=FALSE))
  #Work out which cells to over-plot
  keep = sample(colnames(sc$toc),min(ncol(sc$toc),maxCells))
  #Calculate p-value for each being over some cut-off 
  qVals = do.call(rbind,lapply(nonExpressedGeneList,function(e) p.adjust(pbinom(colSums(sc$toc[e,,drop=FALSE])-1,sc$nUMIs,minRho*sum(sc$soupProfile[e,'est']),lower.tail=FALSE),method='BH')))
  df$qVals = qVals[cbind(match(df[,1],rownames(qVals)),match(df[,2],colnames(qVals)))]
  df$nUMIs = colSums(sc$toc)[df$Barcode]
  #Get the expected number of counts
  expCnts = do.call(rbind,lapply(nonExpressedGeneList,function(e) sc$nUMIs*sum(sc$soupProfile[e,'est'])))
  df$expCnts = expCnts[cbind(match(df[,1],rownames(expCnts)),match(df[,2],colnames(expCnts)))]
  #Set order of marker group as in input
  df$MarkerGroup = factor(df$MarkerGroup,levels=names(nonExpressedGeneList))
  #Add a line estimating the global rho from each group
  tmp = df[df$qVals>0.05,]
  globRhos = sapply(split(tmp,tmp$MarkerGroup),function(e) sum(sc$toc[nonExpressedGeneList[[unique(e$MarkerGroup)]],e$Barcode])/sum(e$nUMIs)/sum(sc$soupProfile[nonExpressedGeneList[[unique(e$MarkerGroup)]],'est']))
  globRhos = data.frame(MarkerGroup=factor(names(globRhos),levels=names(nonExpressedGeneList)),
                        rho = log10(globRhos))
  #Now turn it into a bunch of violin plots
  gg = ggplot(df,aes(MarkerGroup,log10(Values))) +
    geom_violin() +
    geom_jitter(data=df[df$Barcode %in% keep,],aes(colour=qVals<qCut,size=log10(expCnts)),height=0,width=0.3,alpha=1/2) +
    geom_line(data=globRhos,aes(MarkerGroup,rho,group=1),colour='red') +
    scale_colour_manual(values=c('TRUE'='red','FALSE'='black')) + 
    theme(axis.text.x = element_text(angle = 90, hjust = 1)) +
    labs(colour='expressed\nby cell')+
    ylab('log10(observed/expected)') +
    xlab('Marker group')
  gg$df = df
  return(gg)
}

#' Plot ratio of observed to expected counts on reduced dimension map
#'
#' Given some reduced dimensional representation of the data (such as UMAP or tSNE) that has been calculated however you would like, this provides a way to visualise how likely a set of genes are to be soup derived on that map.  That is, given a set of genes, this function calculates how many counts would be expected if that droplet were nothing but soup and compares that to the observed count.  This is done via a log2 ratio of the two values.  A Poisson test is performed and points that have a statistically significant enrichment over the background (at 1% FDR) are marked.
#'
#' @export
#' @param scl SoupChannelList object.
#' @param geneSet A vector with the names of the genes to aggregate and plot evidence for.
#' @param DR A data.frame, with rows named by unique cell IDs (i.e., <ChannelName>___<Barcode>) the first two columns of which give the coordinates of each cell in some reduced dimension representation of the data.
#' @return A ggplot2 containing the plot.
plotMarkerMap = function(scl,geneSet,DR,ratLims=c(-2,2)){
  #Make sure DR is sensible
  DR = as.data.frame(DR)
  if(ncol(DR)<2)
    stop("Need at least two reduced dimensions.")
  if(!(all(rownames(DR) %in% colnames(scl$toc))))
    stop("rownames of DR need to match column names of scl$toc")
  #Get the ratio of observed to expected
  obs = colSums(scl$toc[geneSet,,drop=FALSE])
  exp = colSums(scl$soupMatrix[geneSet,,drop=FALSE])[gsub('___.*','',colnames(scl$toc))]*scl$nUMIs
  expRatio = obs/exp
  #Add it to the dimensions
  DR$geneRatio = expRatio[rownames(DR)]
  colnames(DR)[1:2] = c('RD1','RD2')
  #Sanitise the values a little
  tgtScale = c(ratLims[1],0,ratLims[2])
  #Rescale to be between zero and 1
  rescaled = (tgtScale-tgtScale[1])/(max(tgtScale)-tgtScale[1])
  DR$logRatio = log10(DR$geneRatio)
  #Keep -Inf as NA as we're not really interested in those that have zero expression
  DR$logRatio[DR$logRatio < ratLims[1]] = ratLims[1]
  DR$logRatio[DR$logRatio > ratLims[2]] = ratLims[2]
  DR$logRatio[DR$geneRatio==0]=NA
  #Calculate the corrected p-value of 
  DR$qVals = p.adjust(ppois(obs-1,exp,lower.tail=FALSE),method='BH')[rownames(DR)]
  #Create the plot
  gg = ggplot(DR,aes(RD1,RD2)) +
    #Stick NAs underneath
    geom_point(data=DR[is.na(DR$logRatio),],aes(colour=qVals<.01),size=0.25) +
    geom_point(data=DR[!is.na(DR$logRatio),],aes(fill=logRatio,colour=qVals<0.01),size=2.0,shape=21,stroke=0.5) +
    scale_colour_manual(values=c(`FALSE`='black',`TRUE`='#009933'))+
    xlab('ReducedDim1') +
    ylab('ReducedDim2') +
    scale_fill_gradientn(colours = c('blue','white','red'),
                        values = rescaled,
                        guide='colorbar',
                        limits=ratLims
                        )
    gg$df = DR
    gg
}

#' Plot maps comparing corrected/raw expression
#'
#' Given some reduced dimensional representation of the data (such as UMAP or tSNE) that has been calculated however you would like, this provides a way to visualise how the expression of a geneSet changes after soup correction.
#'
#' @export
#' @param scl SoupChannelList object.
#' @param geneSet A vector with the names of the genes to aggregate and plot evidence for.
#' @param DR A data.frame, with rows named by unique cell IDs (i.e., <ChannelName>___<Barcode>) the first two columns of which give the coordinates of each cell in some reduced dimension representation of the data.
#' @param dataType How should data be represented.  Binary sets each cell to expressed or not, counts converts everything to counts, expression converts everything to expression.  Ratio plots the ratio of everything to uncorrected (and uncorrected is shown as z-scaled expression).
#' @param log Should we log data?
#' @param includePanels Which of 'Uncorrected','CorrectedExpression' and 'CorrectedCounts' to include.  Also sets plot ordering.
#' @return A ggplot2 containing the plot.
plotChangeMap = function(scl,geneSet,DR,dataType=c('binary','counts','expression','ratio'),logData=TRUE,includePanels=c('Uncorrected','CorrectedExpression','CorrectedCounts')){
  dataType = match.arg(dataType)
  if(dataType=='binary')
    logData=FALSE
  if(dataType=='ratio')
    logData=TRUE
  if(is.null(scl$strainedExp) & is.null(scl$atoc))
    stop("No adjusted expression matrix found.")
  #Make sure DR is sensible
  DR = as.data.frame(DR)
  if(ncol(DR)<2)
    stop("Need at least two reduced dimensions.")
  if(!(all(rownames(DR) %in% colnames(scl$toc))))
    stop("rownames of DR need to match column names of scl$toc")
  dfs=list()
  #Get the raw data
  df = DR
  df$correction = 'Uncorrected'
  if(dataType=='binary'){
    df$data = colSums(scl$toc[geneSet,rownames(df),drop=FALSE])>0
  }else if(dataType=='counts'){
    df$data = colSums(scl$toc[geneSet,rownames(df),drop=FALSE])
  }else{
    df$data = colSums(scl$toc[geneSet,rownames(df),drop=FALSE])/scl$nUMIs[rownames(df)]
  }
  #Force z-scale
  if(dataType=='ratio'){
    df$data = scale(df$data)[,1]
  }else{
    if(logData)
      df$data = log10(df$data)
  }
  dfs[['raw']]=df
  #Get the strained expression if present
  if(!is.null(scl$strainedExp)){
    df = DR
    df$correction = 'CorrectedExpression'
    if(dataType=='binary'){
      df$data = colSums(scl$strainedExp[geneSet,rownames(df),drop=FALSE])>0
    }else if(dataType=='counts'){
      df$data = colSums(scl$strainedExp[geneSet,rownames(df),drop=FALSE])*1e4
    }else if(dataType=='ratio'){
      df$data = colSums(scl$strainedExp[geneSet,rownames(df),drop=FALSE])/(colSums(scl$toc[geneSet,rownames(df),drop=FALSE])/scl$nUMIs[rownames(df)])
    }else{
      df$data = colSums(scl$strainedExp[geneSet,rownames(df),drop=FALSE])
    }
    if(logData)
      df$data = log10(df$data)
    dfs[['correctedExpression']]=df
  }
  if(!is.null(scl$strainedExp)){
    df = DR
    df$correction = 'CorrectedCounts'
    if(dataType=='binary'){
      df$data = colSums(scl$atoc[geneSet,rownames(df),drop=FALSE])>0
    }else if(dataType=='counts'){
      df$data = colSums(scl$atoc[geneSet,rownames(df),drop=FALSE])
    }else if(dataType=='ratio'){
      df$data = (colSums(scl$atoc[geneSet,rownames(df),drop=FALSE])/colSums(scl$atoc[,rownames(df),drop=FALSE]))/(colSums(scl$toc[geneSet,rownames(df),drop=FALSE])/scl$nUMIs[rownames(df)])
    }else{
      df$data = colSums(scl$atoc[geneSet,rownames(df),drop=FALSE])/colSums(scl$atoc[,rownames(df),drop=FALSE])
    }
    if(logData)
      df$data = log10(df$data)
    dfs[['correctedCounts']]=df
  }
  dfs = do.call(rbind,dfs)
  #Define order to plot
  dfs = dfs[dfs$correction %in% includePanels,]
  lvls = includePanels
  dfs$correction = factor(dfs$correction,levels=lvls[lvls%in%dfs$correction])
  #Truncate things on ratio plot to make limits clear
  if(dataType=='ratio'){
    dfs$data[dfs$data< -1]=-1
    dfs$data[dfs$data> 1]=1
  }
  #Now make the plot
  gg = ggplot(dfs,aes(RD1,RD2)) +
    geom_point(aes(colour=data),size=0.5) +
    xlab('ReducedDim1') +
    ylab('ReducedDim2') +
    labs(colour='geneSet') + 
    ggtitle('Comparison of before and after correction') +
    facet_wrap(~correction)
  if(dataType=='ratio')
    gg = gg + scale_colour_gradientn(colours=rainbow(50),limits=c(-1,1))
  #Store the input in case the user wants to make their own plot
  gg$df = dfs
  gg
}


####' Plots a summary of the contamination for a channel
####' 
####' Takes the output of \code{\link{strainChannel}} (or \code{\link{strain10X}}) and plots a summary.  This just pieces together the plots produced by \code{\link{plotChannelSummary}}, \code{\link{plotSoupCorrelation}}, \code{\link{plotMarkerDistribution}} and \code{\link{inferNonExpressedGenes}}.
####'
####' @export
####' @param cleanedChannel Output of \code{\link{strainChannel}} or \code{\link{strain10X}}.
####' @return A ggplot2 object containing the plot.
####' @importFrom cowplot ggdraw draw_plot draw_plot_label
###plotChannelSummary = function(cleanedChannel){
###  #Make the plots
###  gg_cFrac = plotChannelContamination(cleanedChannel$contaminationFractions)
###  gg_soupCorrelation = plotSoupCorrelation(cleanedChannel$rawCounts,cleanedChannel$soupProfile)
###  gg_markers = plotMarkerDistribution(cleanedChannel$rawCounts,cleanedChannel$soupProfile,cleanedChannel$nonExpressedGeneList)
###  tmp = inferNonExpressedGenes(cleanedChannel$rawCounts,cleanedChannel$soupProfile)
###  gg_bimod = plotMarkerDistribution(cleanedChannel$rawCounts,cleanedChannel$soupProfile,rownames(tmp)[seq(min(20,nrow(tmp),sum(tmp$isUseful)))])
###  #Plaster them together
###  ggdraw() + 
###    draw_plot(gg_soupCorrelation,x=0,y=0.5,width=0.5,height=0.5) +
###    draw_plot(gg_cFrac,x=0.5,y=0.5,width=0.5,height=0.5) +
###    draw_plot(gg_markers,x=0,y=0,width=0.5,height=0.5) +
###    draw_plot(gg_bimod,x=0.5,y=0,width=0.5,height=0.5) +
###    draw_plot_label(label=c('A','B','C','D'),size=15,x=c(0,0.5,0,0.5),y=c(1,1,0.5,0.5))
###}
 

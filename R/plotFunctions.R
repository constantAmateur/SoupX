#' Plot contamination fraction for channel
#'
#' Plots the contamination fraction as a function of nUMIs in a droplet.
#'
#' @export
#' @param rho Output from \code{\link{calculateContaminationFraction}}.
#' @param showErrorBars Should error bars showing the 95 percent confidence interval on the estimate of rho be included?
#' @return A ggplot2 object containing the plot.
plotChannelContamination = function(rhos,showErrorBars=TRUE){
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
  return(gg)
}

#' Plot correlation of expression profiles of soup and aggregated cells
#' 
#' Calculates an expression profile by aggregating counts across all cells and plots this (on a log10 scale) against the expression profile of the soup.
#'
#' @export
#' @param toc Table of UMIs for all cells.
#' @param soupProfile The expression profile of the soup (output from \code{\link{estimateSoup}})
#' @return A ggplot2 object containing the plot.
plotSoupCorrelation = function(toc,soupProfile){
  #Calculate the cell profile
  cellProfile = rowSums(toc)
  cellProfile = (cellProfile/sum(cellProfile))
  df = as.data.frame(cbind(cellProfile,soupProfile=soupProfile$est))
  gg = ggplot(df,aes(log10(cellProfile),log10(soupProfile))) +
    geom_point(alpha=1/3) +
    geom_abline(intercept=0,slope=1) +
    ylab('log10(Soup Expression)')+
    xlab('log10(Aggregate cell Expression)')
  return(gg)
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
#' @param toc Table of UMIs for all cells.
#' @param soupProfile The expression profile of the soup (output from \code{\link{estimateSoup}})
#' @param nonExpressedGeneList Which genes to use to estimate soup (see \code{\link{calculateContaminationFraction}}).  Can also be a vector of gene names, in which case each gene is treated as a separate gene set.
#' @param maxCells Randomly plot only this many cells to prevent over-crowding.
#' @param minRho A statistical test is performed to identify all those cells with more expression than would be possible from the soup alone.  This is done by testing against a null hypothesis of all expression coming from the soup, with the soup contamination fraction being less than or equal to \code{minRho}.
#' @param qCut The FDR cut-off for the hypothesis test described above.
#' @importFrom stats setNames
#' @return A ggplot2 object containing the plot.
plotMarkerDistribution = function(toc,soupProfile,nonExpressedGeneList,maxCells=150,minRho=1.0,qCut=0.05){
  #Make non-lists into lists
  if(!is.list(nonExpressedGeneList))
    nonExpressedGeneList = as.list(setNames(nonExpressedGeneList,nonExpressedGeneList))
  #Calculate the ratio to the soup for each marker group in each cell
  obsProfile = t(t(toc)/colSums(toc))
  #Get the ratio
  tst = lapply(nonExpressedGeneList,function(e) colSums(obsProfile[e,,drop=FALSE])/sum(soupProfile[e,'est']))
  #Unlist the thing
  df = data.frame(MarkerGroup = rep(names(tst),lengths(tst)),
                  Barcode=unlist(lapply(tst,names),use.names=FALSE),
                  Values=unlist(tst,use.names=FALSE))
  #Work out which cells to over-plot
  keep = sample(colnames(toc),min(ncol(toc),maxCells))
  #Calculate p-value for each being over some cut-off 
  qVals = do.call(rbind,lapply(nonExpressedGeneList,function(e) p.adjust(pbinom(colSums(toc[e,,drop=FALSE])-1,nUMIs,minRho*sum(soupProfile[e,'est']),lower.tail=FALSE),method='BH')))
  df$qVals = qVals[cbind(match(df[,1],rownames(qVals)),match(df[,2],colnames(qVals)))]
  df$nUMIs = colSums(toc)[df$Barcode]
  #Get the expected number of counts
  expCnts = do.call(rbind,lapply(nonExpressedGeneList,function(e) nUMIs*sum(soupProfile[e,'est'])))
  df$expCnts = expCnts[cbind(match(df[,1],rownames(expCnts)),match(df[,2],colnames(expCnts)))]
  #Set order of marker group as in input
  df$MarkerGroup = factor(df$MarkerGroup,levels=names(nonExpressedGeneList))
  #Add a line estimating the global rho from each group
  tmp = df[df$qVals>0.05,]
  globRhos = sapply(split(tmp,tmp$MarkerGroup),function(e) sum(toc[nonExpressedGeneList[[unique(e$MarkerGroup)]],e$Barcode])/sum(e$nUMIs)/sum(soupProfile[nonExpressedGeneList[[unique(e$MarkerGroup)]],'est']))
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
  return(gg)
}

#' Plots a summary of the contamination for a channel
#' 
#' Takes the output of \code{\link{strainChannel}} (or \code{\link{strain10X}}) and plots a summary.  This just pieces together the plots produced by \code{\link{plotChannelSummary}}, \code{\link{plotSoupCorrelation}}, \code{\link{plotMarkerDistribution}} and \code{\link{inferNonExpressedGenes}}.
#'
#' @export
#' @param cleanedChannel Output of \code{\link{strainChannel}} or \code{\link{strain10X}}.
#' @return A ggplot2 object containing the plot.
#' @importFrom cowplot ggdraw draw_plot draw_plot_label
plotChannelSummary = function(cleanedChannel){
  #Make the plots
  gg_cFrac = plotChannelContamination(cleanedChannel$contaminationFractions)
  gg_soupCorrelation = plotSoupCorrelation(cleanedChannel$rawCounts,cleanedChannel$soupProfile)
  gg_markers = plotMarkerDistribution(cleanedChannel$rawCounts,cleanedChannel$soupProfile,cleanedChannel$nonExpressedGeneList)
  tmp = inferNonExpressedGenes(cleanedChannel$rawCounts,cleanedChannel$soupProfile)
  gg_bimod = plotMarkerDistribution(cleanedChannel$rawCounts,cleanedChannel$soupProfile,rownames(tmp)[seq(min(20,nrow(tmp),sum(tmp$isUseful)))])
  #Plaster them together
  ggdraw() + 
    draw_plot(gg_soupCorrelation,x=0,y=0.5,width=0.5,height=0.5) +
    draw_plot(gg_cFrac,x=0.5,y=0.5,width=0.5,height=0.5) +
    draw_plot(gg_markers,x=0,y=0,width=0.5,height=0.5) +
    draw_plot(gg_bimod,x=0.5,y=0,width=0.5,height=0.5) +
    draw_plot_label(label=c('A','B','C','D'),size=15,x=c(0,0.5,0,0.5),y=c(1,1,0.5,0.5))
}
 

#' Plot correlation of expression profiles of soup and aggregated cells
#' 
#' Calculates an expression profile by aggregating counts across all cells and plots this (on a log10 scale) against the expression profile of the soup.
#'
#' @param sc A SoupChannel object.
#' @return A ggplot2 object containing the plot.
plotSoupCorrelation = function(sc){
  if(!is(sc,'SoupChannel'))
    stop("sc not a valid SoupChannel object.")
  #Calculate the cell profile
  cellProfile = rowSums(sc$toc)
  cellProfile = (cellProfile/sum(cellProfile))
  df = data.frame(cellProfile,soupProfile=sc$soupProfile$est)
  gg = ggplot(df,aes(log10(cellProfile),log10(soupProfile))) +
    geom_point(alpha=1/3) +
    geom_abline(intercept=0,slope=1) +
    ylab('log10(Soup Expression)')+
    xlab('log10(Aggregate cell Expression)')
  return(gg)
}

#' Plots the distribution of the observed to expected expression for marker genes
#'
#' If each cell were made up purely of background reads, the expression fraction would equal that of the soup.  This plot compares this expectation of pure background to the observed expression fraction in each cell, for each of the groups of genes in \code{nonExpressedGeneList}.  For each group of genes, the distribution of this ratio is plotted across all cells.  A value significantly greater than 1 (0 on log scale) can only be obtained if some of the genes in each group are genuinely expressed by the cell.  That is, the assumption that the cell is pure background does not hold for that gene.
#'
#' This plot is a useful diagnostic for the assumption that a list of genes is non-expressed in most cell types.  For non-expressed cells, the ratio should cluster around the contamination fraction, while for expressed cells it should be elevated.  The most useful non-expressed gene sets are those for which the genes are either strongly expressed, or not expressed at all.  Such groups of genes will show up in this plot as a bimodal distribution, with one mode containing the cells that do not express these genes around the contamination fraction for this channel and another around a value at some value equal to or greater than 0 (1 on non-log scale) for the expressed cells.
#'
#' The red line shows the global estimate of the contamination for each group of markers.  This is usually lower than the low mode of the distribution as there will typically be a non-negligible number of cells with 0 observed counts (and hence -infinity log ratio).
#'
#' If \code{nonExpressedGeneList} is missing, this function will try and find genes that are very specific to different clusters, as these are often the most useful in estimating the contamination fraction.   This is meant only as a heuristic, which can hopefully provide some inspiration as to a class of genes to use to estimation the contamination for your experiment.  Please do **NOT** blindly use the top N genes found in this way to estimate the contamination.  That is, do not feed this list of genes into \code{\link{calculateContaminationFraction}} without any manual consideration or filtering as this *will over-estimate your contamination* (often by a large amount).  For this reason, these gene names are not returned by the function. 
#'
#' @export
#' @param sc A SoupChannel object.
#' @param nonExpressedGeneList Which sets of genes to use to estimate soup (see \code{\link{calculateContaminationFraction}}).
#' @param maxCells Randomly plot only this many cells to prevent over-crowding.
#' @param tfidfMin Minimum specificity cut-off used if finding marker genes (see \code{\link{quickMarkers}}).
#' @param ... Passed to \code{\link{estimateNonExpressingCells}}
#' @importFrom stats setNames
#' @return A ggplot2 object containing the plot.
#' @examples
#' gg = plotMarkerDistribution(scToy,list(CD7='CD7',LTB='LTB'))
plotMarkerDistribution = function(sc,nonExpressedGeneList,maxCells=150,tfidfMin=1,...){
  if(!is(sc,'SoupChannel'))
    stop("sc not a valid SoupChannel object.")
  #Get nonExpressedGeneList algorithmically if missing...
  if(missing(nonExpressedGeneList)){
    message("No gene lists provided, attempting to find and plot cluster marker genes.")
    #Get marker genes instead.  Obviously this requires clustering
    if(!'clusters' %in% colnames(sc$metaData))
      stop("Failed as no clusters found!  Clustering must be set via 'setClusters' to find marker genes.")
    #Get top markers
    mrks = quickMarkers(sc$toc,sc$metaData$clusters,N=Inf)
    #And only the most specific entry for each gene
    mrks = mrks[order(mrks$gene,-mrks$tfidf),]
    mrks = mrks[!duplicated(mrks$gene),]
    #Order by tfidif maxness
    mrks = mrks[order(-mrks$tfidf),]
    #Apply tf-idf cut-off
    mrks = mrks[mrks$tfidf > tfidfMin,]
    message(sprintf("Found %d marker genes",nrow(mrks)))
    #Of the top markers, order by soup 
    mrks = mrks[order(sc$soupProfile[mrks$gene,'est'],decreasing=TRUE),]
    #And keep the top 20
    nonExpressedGeneList = mrks$gene[seq(min(nrow(mrks),20))]
    nonExpressedGeneList = setNames(as.list(nonExpressedGeneList),nonExpressedGeneList)
  }
  #Make non-lists into lists
  if(!is.list(nonExpressedGeneList))
    stop("nonExpressedGeneList must be a list of sets of genes.  e.g. list(HB = c('HBB','HBA2'))")
  #Get the non-expressing matrix
  nullMat = estimateNonExpressingCells(sc,nonExpressedGeneList,...)
  #Calculate the ratio to the soup for each marker group in each cell
  obsProfile = t(t(sc$toc)/sc$metaData$nUMIs)
  #Get the ratio
  tst = lapply(nonExpressedGeneList,function(e) colSums(obsProfile[e,,drop=FALSE])/sum(sc$soupProfile[e,'est']))
  #Unlist the thing
  df = data.frame(MarkerGroup = rep(names(tst),lengths(tst)),
                  Barcode=unlist(lapply(tst,names),use.names=FALSE),
                  Values=unlist(tst,use.names=FALSE))
  #Work out which cells to over-plot
  keep = sample(colnames(sc$toc),min(ncol(sc$toc),maxCells))
  #Calculate p-value for each being over some cut-off 
  #qVals = do.call(rbind,lapply(nonExpressedGeneList,function(e) p.adjust(pbinom(colSums(sc$toc[e,,drop=FALSE])-1,sc$metaData$nUMIs,minRho*sum(sc$soupProfile[e,'est']),lower.tail=FALSE),method='BH')))
  #df$qVals = qVals[cbind(match(df[,1],rownames(qVals)),match(df[,2],colnames(qVals)))]
  df$nUMIs = sc$metaData[df$Barcode,'nUMIs']
  #Get the expected number of counts
  expCnts = do.call(rbind,lapply(nonExpressedGeneList,function(e) sc$metaData$nUMIs*sum(sc$soupProfile[e,'est'])))
  colnames(expCnts) = rownames(sc$metaData)
  df$expCnts = expCnts[cbind(match(df$MarkerGroup,rownames(expCnts)),match(df$Barcode,colnames(expCnts)))]
  #Set order of marker group as in input
  df$MarkerGroup = factor(df$MarkerGroup,levels=names(nonExpressedGeneList))
  #Add a line estimating the global rho from each group
  #First calculate global rho using nullMat
  globalRhos=c()
  for(i in seq_along(nonExpressedGeneList)){
    if(sum(nullMat[,i])>0){
      tmp = suppressMessages(calculateContaminationFraction(sc,nonExpressedGeneList[i],nullMat[,i,drop=FALSE],forceAccept=TRUE))
      globalRhos = c(globalRhos,tmp$metaData$rho[1])
    }else{
      globalRhos = c(globalRhos,NA)
    }
  }
  names(globalRhos) = names(nonExpressedGeneList)
  globalRhos = data.frame(MarkerGroup = factor(names(globalRhos),levels=names(nonExpressedGeneList)),
                          rho = log10(globalRhos))
  #tmp = df[df$qVals>0.05,]
  #globRhos = sapply(split(tmp,tmp$MarkerGroup),function(e) sum(sc$toc[nonExpressedGeneList[[unique(e$MarkerGroup)]],e$Barcode])/sum(e$nUMIs)/sum(sc$soupProfile[nonExpressedGeneList[[unique(e$MarkerGroup)]],'est']))
  #globRhos = data.frame(MarkerGroup=factor(names(globRhos),levels=names(nonExpressedGeneList)),
  #                      rho = log10(globRhos))
  #Now turn it into a bunch of violin plots
  gg = ggplot(df,aes(MarkerGroup,log10(Values))) +
    geom_violin() +
    geom_jitter(data=df[df$Barcode %in% keep,],aes(size=log10(expCnts)),height=0,width=0.3,alpha=1/2) +
    geom_line(data=globalRhos,aes(MarkerGroup,rho,group=1),colour='red') +
    geom_point(data=globalRhos,aes(MarkerGroup,rho,group=1),colour='red',shape=2) +
    scale_colour_manual(values=c('TRUE'='red','FALSE'='black')) + 
    theme(axis.text.x = element_text(angle = 90, hjust = 1)) +
    labs(colour='expressed\nby cell')+
    ylab('log10(observed/expected)') +
    xlab('Marker group')
  return(gg)
}

#' Plot ratio of observed to expected counts on reduced dimension map
#'
#' Given some reduced dimensional representation of the data (such as UMAP or tSNE) that has been calculated however you would like, this provides a way to visualise how likely a set of genes are to be soup derived on that map.  That is, given a set of genes, this function calculates how many counts would be expected if that droplet were nothing but soup and compares that to the observed count.  This is done via a log2 ratio of the two values.  A Poisson test is performed and points that have a statistically significant enrichment over the background (at 5% FDR) are marked.
#'
#' @export
#' @param sc SoupChannel object.
#' @param geneSet A vector with the names of the genes to aggregate and plot evidence for.
#' @param DR A data.frame, with rows named by unique cell IDs (i.e., <ChannelName>_<Barcode>) the first two columns of which give the coordinates of each cell in some reduced dimension representation of the data.  Try and fetch automatically if missing.
#' @param ratLims Truncate log ratios at these values.
#' @param FDR False Discovery Rate for statistical test of enrichment over background.
#' @param useToEst A vector (usually obtained from \code{\link{estimateNonExpressingCells}}), that will be used to mark cells instead of the usual Poisson test.
#' @param pointSize Size of points
#' @param pointShape Shape of points
#' @param pointStroke Stroke size for points
#' @param naPointSize Point size for NAs.
#' @return A ggplot2 containing the plot.
#' @examples
#' gg = plotMarkerMap(scToy,'CD7')
plotMarkerMap = function(sc,geneSet,DR,ratLims=c(-2,2),FDR=0.05,useToEst=NULL,pointSize=2.0,pointShape=21,pointStroke=0.5,naPointSize=0.25){
  if(!is(sc,'SoupChannel'))
    stop("sc not a valid SoupChannel object.")
  #Try and get DR if missing
  if(missing(DR))
    DR = sc$metaData[,sc$DR]
  #Make sure DR is sensible
  DR = as.data.frame(DR)
  if(ncol(DR)<2)
    stop("Need at least two reduced dimensions.")
  if(!(all(rownames(DR) %in% colnames(sc$toc))))
    stop("rownames of DR need to match column names of sc$toc")
  #Get the ratio of observed to expected
  obs = colSums(sc$toc[geneSet,,drop=FALSE])
  exp = sc$metaData$nUMIs*sum(sc$soupProfile[geneSet,'est'])
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
  colVal = 'qVals<FDR'
  if(!is.null(useToEst)){
    DR$useToEst = useToEst
    colVal='useToEst'
  }
  #Create the plot
  gg = ggplot(DR,aes(RD1,RD2)) +
    #Stick NAs underneath
    geom_point(data=DR[is.na(DR$logRatio),],aes_string(colour=colVal),size=naPointSize) +
    geom_point(data=DR[!is.na(DR$logRatio),],aes_string(fill='logRatio',colour=colVal),size=pointSize,shape=pointShape,stroke=pointStroke) +
    scale_colour_manual(values=c(`FALSE`='black',`TRUE`='#009933'))+
    xlab('ReducedDim1') +
    ylab('ReducedDim2') +
    scale_fill_gradientn(colours = c('blue','white','red'),
                        values = rescaled,
                        guide='colorbar',
                        limits=ratLims
                        )
    gg
}

#' Plot maps comparing corrected/raw expression
#'
#' Given some reduced dimensional representation of the data (such as UMAP or tSNE) that has been calculated however you would like, this provides a way to visualise how the expression of a geneSet changes after soup correction.
#'
#' @export
#' @param sc SoupChannel object.
#' @param cleanedMatrix A cleaned matrix to compare against the raw one.  Usually the output of \code{\link{adjustCounts}}.
#' @param geneSet A vector with the names of the genes to aggregate and plot evidence for.
#' @param DR A data.frame, with rows named by unique cell IDs (i.e., <ChannelName>_<Barcode>) the first two columns of which give the coordinates of each cell in some reduced dimension representation of the data.
#' @param dataType How should data be represented.  Binary sets each cell to expressed or not, counts converts everything to counts, soupFrac plots the fraction of the observed counts that are identified as contamination (i.e., (old-new)/old) for each cell and is the default.
#' @param logData Should we log the thing we plot?
#' @param pointSize Size of points
#' @return A ggplot2 containing the plot.
#' @examples
#' out = adjustCounts(scToy)
#' gg = plotChangeMap(scToy,out,'S100A9')
plotChangeMap = function(sc,cleanedMatrix,geneSet,DR,dataType=c('soupFrac','binary','counts'),logData=FALSE,pointSize=0.5){
  dataType = match.arg(dataType)
  if(dataType=='binary')
    logData=FALSE
  #Try and get DR if missing
  if(missing(DR))
    DR = sc$metaData[,sc$DR]
  #Make sure DR is sensible
  DR = as.data.frame(DR)
  if(ncol(DR)<2)
    stop("Need at least two reduced dimensions.")
  if(!(all(rownames(DR) %in% colnames(sc$toc))))
    stop("rownames of DR need to match column names of sc$toc")
  colnames(DR)[1:2] = c('RD1','RD2')
  #Simple one panel version
  if(dataType=='soupFrac'){
    df = DR
    old = colSums(sc$toc[geneSet,rownames(df),drop=FALSE])
    new = colSums(cleanedMatrix[geneSet,rownames(df),drop=FALSE])
    relChange = (old-new)/old
    df$old = old
    df$new = new
    df$relChange=relChange
    nom = 'SoupFrac'
    if(logData){
      df$relChange = log10(df$relChange)
      nom = paste0('log10(',nom,')')
      zLims=c(-2,0)
    }else{
      zLims=c(0,1)
    }
    df = df[order(!is.na(df$relChange)),]
    #Truncate to prevent -Inf -> NA coloured dots
    df$relChange[which(df$relChange<zLims[1])] = zLims[1]
    df$relChange[which(df$relChange>zLims[2])] = zLims[2]
    gg = ggplot(df,aes(RD1,RD2)) +
      geom_point(aes(col=relChange),size=0.5) +
      xlab('ReducedDim1') +
      ylab('ReducedDim2') +
      labs(colour=nom) + 
      ggtitle('Change in expression due to soup correction')
  }else{
    dfs=list()
    #Get the raw data
    df = DR
    df$correction = 'Uncorrected'
    if(dataType=='binary'){
      df$data = colSums(sc$toc[geneSet,rownames(df),drop=FALSE])>0
    }else if(dataType=='counts'){
      df$data = colSums(sc$toc[geneSet,rownames(df),drop=FALSE])
    }
    if(logData)
      df$data = log10(df$data)
    dfs[['raw']]=df
    #And the corrected expression
    df = DR
    df$correction = 'Corrected'
    if(dataType=='binary'){
      df$data = colSums(cleanedMatrix[geneSet,rownames(df),drop=FALSE])>0
    }else if(dataType=='counts'){
      df$data = colSums(cleanedMatrix[geneSet,rownames(df),drop=FALSE])
    }
    if(logData)
      df$data = log10(df$data)
    dfs[['correctedExpression']]=df
    dfs = do.call(rbind,dfs)
    #Define order to plot
    lvls = c('Uncorrected','Corrected')
    dfs$correction = factor(dfs$correction,levels=lvls[lvls%in%dfs$correction])
    #Stick the NAs at the bottom
    dfs = dfs[order(!is.na(dfs$data)),]
    zLims=c(NA,NA)
    #Now make the plot
    gg = ggplot(dfs,aes(RD1,RD2)) +
      geom_point(aes(colour=data),size=pointSize) +
      xlab('ReducedDim1') +
      ylab('ReducedDim2') +
      labs(colour='geneSet') + 
      ggtitle('Comparison of before and after correction') +
      facet_wrap(~correction)
  }
  #Make less ugly colour scheme
  if(dataType!='binary')
    gg = gg + scale_colour_gradientn(colours=c('#fff5eb','#fee6ce','#fdd0a2','#fdae6b','#fd8d3c','#f16913','#d94801','#a63603','#7f2704'),limits=zLims)
  gg
}

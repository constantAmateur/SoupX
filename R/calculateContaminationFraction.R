#' Calculate the contamination fraction
#'
#' This function computes the contamination fraction using two user-provided bits of information.  Firstly, a list of sets of genes that can be biologically assumed to be absent in at least some cells in your data set.  For example, these might be haemoglobin genes or immunoglobulin genes, which should not be expressed outside of erythroyctes and antibody producing cells respectively.
#'
#' Secondly, this function needs to know which cells definitely do not express the gene sets described above.  Continuing with the haemoglobin example, which are the erythrocytes that are producing haemoglobin mRNAs and which are non-erythrocytes that we can safely assume produce no such genes.  The assumption made is any expression from a gene set in cell marked as a "non-expressor" for that gene set, must be derived from the soup.  Therefore, the level of contamination present can be estimated from the amount of expression of these genes seen in these cells.
#'
#' Most often, the genesets are user supplied based on your knowledge of the experiment and the cells in which they are genuinely expressed is estimated using \code{\link{estimateNonExpressingCells}}.  However, they can also be supplied directly if other information is available.
#'
#' Usually, there is very little variation in the contamination fraction within a channel and very little power to detect the contamination accurately at a single cell level.  As such, the default mode of operation simply estimates one value of the contamination fraction that is applied to all cells in a channel.  However, to allow for circumstances where there is good power to detect fine grained contamination differences, a cell-specific estimation procedure can be invoked by setting \code{cellSpecificEstimates=TRUE}.
#'
#' The global model fits a simple Poisson glm to the aggregated count data across all cells.  The cell-specific model uses hierarchical Bayes estimation using STAN to share information between cells.  Note that this takes much longer and should only be used if you have a convincing reason.  In the most cases, a global estimate will do a good job.
#'
#' Finally, note that if you are not able to find a reliable set of genes to use for contamination estimation, or you do not trust the values produced, the contamination fraction can be manually set by the user using \code{\link{setContaminationFraction}}.
#'
#' @export 
#' @param sc A SoupChannel object.
#' @param nonExpressedGeneList A list containing sets of genes which can be assumed to be non-expressed in a subset of cells (see details).
#' @param useToEst A boolean matrix of dimensions ncol(toc) x length(nonExpressedGeneList) indicating which gene-sets should not be assumed to be non-expressed in each cell.  Row names must correspond to the names of \code{nonExpressedGeneList}.  Usually produced by \code{\link{estimateNonExpressingCells}}.
#' @param verbose Print best estimate.
#' @param cellSpecificEstimates Allow the estimated contamination fraction to differ at the individual cell level.  Note that this approach is MUCH slower.
#' @param stanParams Extra parameters passed to \code{rstan::stan}, when \code{cellSpecificEstimates=TRUE}.
#' @param forceAccept Passed to \code{\link{setContaminationFraction}}.
#' @return A modified version of \code{sc} with estimates of the contamination (rho) added to the metaData table.
#' @examples
#' \dontrun{
#' geneList = list(HB=c('HBB','HBA2'))
#' ute = estimateNonExpressingCells(sc,geneList)
#' sc = calculateContaminationFraction(sc,geneList,ute)
#' #Use hierarchical bayes method to estimate cell specific contamination (requires rstan)
#' sc = calculateContaminationFraction(sc,geneList,ute,cellSpecificEstimates=TRUE)
#' }
#' @importFrom stats coef confint glm poisson quantile
calculateContaminationFraction = function(sc,nonExpressedGeneList,useToEst,verbose=TRUE,cellSpecificEstimates=FALSE,stanParams=list(chains=1,warmup=8000,iter=24000,cores=1),forceAccept=FALSE){
  if(!is(sc,'SoupChannel')){
    stop("sc must be a SoupChannel object")
  }
  #Check that you've provided the genes in the right format
  if(!is.list(nonExpressedGeneList))
    stop("nonExpressedGeneList must be a list of sets of genes.  e.g. list(HB = c('HBB','HBA2'))")
  #Check we can estimate
  if(sum(useToEst)==0)
    stop("No cells specified as acceptable for estimation.  useToEst must not be all FALSE")
  #Construct the data.frame to perform inferance on
  df = list()
  for(i in seq_along(nonExpressedGeneList)){
    tgts = nonExpressedGeneList[[i]]
    #Get the soup fraction for this set
    sFrac = sum(sc$soupProfile[tgts,'est'])
    w = rownames(useToEst)[useToEst[,i]]
    if(length(w)>0){
      #Get the counts
      cnts = as.matrix(sc$toc[tgts,w,drop=FALSE])
      df[[i]] = data.frame(row.names=NULL,
                           cells=colnames(cnts),
                           geneSet=i,
                           soupFrac = sFrac,
                           counts=colSums(cnts),
                           stringsAsFactors=FALSE)
    }
  }
  df = do.call(rbind,df)
  df$nUMIs = sc$metaData[df$cells,'nUMIs']
  df$expSoupCnts = df$nUMIs * df$soupFrac
  #Make cells a factor, but preserve ordering in sc.  This ensures the the STAN index is the index in sc$metaData
  df$cells = factor(df$cells,levels=rownames(sc$metaData))
  #Now decide what model to fit
  if(cellSpecificEstimates){
    #It just won't seem to work without this...
    if(!requireNamespace("rstan", quietly = TRUE))
      stop("Package \"rstan\" needed for cell specific estimation. Please install it.",call. = FALSE)
    #This needs to be loaded explicitly or we get rstan errors
    require('rstan',quietly=TRUE)
    mod="
data{
    int<lower=1> N;
    int<lower=1> N_cells;
    int counts[N];
    int cells[N];
    real expSoupCnts[N];
}
parameters{
    real rho;
    real dist_mean;
    real<lower=0> dist_sd;
    vector[N_cells] rho_cells;
}
model{
    vector[N] mu;
    dist_mean ~ normal(0,0.5);
    dist_sd ~ normal(0,1);
    rho_cells ~ normal( dist_mean , dist_sd );
    rho ~ normal( -4 , 1 );
    for ( i in 1:N ) {
        mu[i] = rho + rho_cells[cells[i]] + log(expSoupCnts[i]);
    }
    counts ~ poisson_log( mu );
}
generated quantities{
    vector[N] mu;
    real dev;
    dev = 0;
    for ( i in 1:N ) {
        mu[i] = rho + rho_cells[cells[i]] + log(expSoupCnts[i]);
    }
    dev = dev + (-2)*poisson_log_lpmf( counts | mu );
}"
    dat = list(N=nrow(df),
               N_cells = nlevels(df$cells),
               counts = df$counts,
               cells = as.numeric(df$cells),
               nUMIs = df$nUMIs,
               expSoupCnts = df$expSoupCnts
               )
    stanParams$model_code = mod
    stanParams$data = dat
    sc$fit = do.call(rstan::stan,stanParams)
    #Match the cells
    ee = rstan::extract(sc$fit)
    #Construct cell specific estimates
    ee = as.numeric(ee$rho) + as.matrix(ee$rho_cells)
    sc = setContaminationFraction(sc,setNames(exp(colMeans(ee)),rownames(sc$metaData)),forceAccept=forceAccept)
    sc$metaData$rhoLow = exp(apply(ee,2,quantile,0.025))
    sc$metaData$rhoHigh = exp(apply(ee,2,quantile,0.975))
  }else{
    #Fit Poisson GLM with log-link to get global rho
    sc$fit = glm(counts ~ 1,family=poisson(),offset=log(expSoupCnts),data=df)
    #Finally, add it to the meta-data
    sc = setContaminationFraction(sc,exp(coef(sc$fit)),forceAccept=forceAccept)
    tmp=suppressMessages(confint(sc$fit))
    sc$metaData$rhoLow = exp(tmp[1])
    sc$metaData$rhoHigh = exp(tmp[2])
    if(verbose)
      message(sprintf("Estimated global contamination fraction of %0.2f%%",100*exp(coef(sc$fit)))) 
  }
  return(sc)
}

#' Gets top N markers for each cluster
#'
#' Uses tf-idf ordering to get the top N markers of each cluster.  For each cluster, either the top N or all genes passing the hypergeometric test with the FDR specified, whichever list is smallest.
#' 
#' Term Frequency - Inverse Document Frequency is used in natural language processing to identify terms specific to documents.  This function uses the same idea to order genes within a group by how predictive of that group they are.  The main advantage of this is that it is extremely fast and gives reasonable results.
#'
#' To do this, gene expression is binarised in each cell so each cell is either considered to express or not each gene.  That is, we replace the counts with \code{toc > zeroCut}.  The frequency with which a gene is expressed within the target group is compared to the global frequency to calculate the tf-idf score.  We also calculate a multiple hypothesis corrected p-value based on a hypergeometric test, but this is extremely permissive.
#'
#' @export
#' @param toc Table of counts.  Must be a sparse matrix.
#' @param clusters Vector of length \code{ncol(toc)} giving cluster membership.
#' @param N Number of marker genes to return per cluster.
#' @param FDR False discover rate to use. 
#' @param expressCut Value above which a gene is considered expressed.
#' @return data.frame with top N markers (or all that pass the hypergeometric test) and their statistics for each cluster.
#' @examples
#' #Calculate markers of clusters in toy data
#' mrks = quickMarkers(scToy$toc,scToy$metaData$clusters)
#' \dontrun{
#' #Calculate markers from Seurat (v3) object
#' mrks = quickMarkers(srat@assays$RNA@count,srat@active.ident)
#' }
quickMarkers = function(toc,clusters,N=10,FDR=0.01,expressCut=0.9){
  #Convert to the more manipulable format
  #toc = as(toc,'dgTMatrix')
  toc = as(as(as(toc, "dMatrix"), "generalMatrix"), "TsparseMatrix")
  w = which(toc@x>expressCut)
  #Get the counts in each cluster
  clCnts = table(clusters)
  nObs = split(factor(rownames(toc))[toc@i[w]+1],clusters[toc@j[w]+1])
  nObs = sapply(nObs,table)
  #Calculate the observed and total frequency
  nTot = rowSums(nObs)
  tf = t(t(nObs)/as.integer(clCnts[colnames(nObs)]))
  ntf = t(t(nTot - nObs)/as.integer(ncol(toc)-clCnts[colnames(nObs)]))
  idf = log(ncol(toc)/nTot)
  score = tf*idf
  #Calculate p-values
  qvals = lapply(seq_len(ncol(nObs)),function(e)
                 p.adjust(phyper(nObs[,e]-1,nTot,ncol(toc)-nTot,clCnts[colnames(nObs)[e]],lower.tail=FALSE),method='BH'))
  qvals = do.call(cbind,qvals)
  colnames(qvals) = colnames(nObs)
  #Get gene frequency of second best
  sndBest = lapply(seq_len(ncol(tf)),function(e) apply(tf[,-e,drop=FALSE],1,max))
  sndBest = do.call(cbind,sndBest)
  colnames(sndBest) = colnames(tf)
  #And the name
  sndBestName = lapply(seq_len(ncol(tf)),function(e) (colnames(tf)[-e])[apply(tf[,-e,drop=FALSE],1,which.max)])
  sndBestName = do.call(cbind,sndBestName)
  colnames(sndBestName) = colnames(tf)
  rownames(sndBestName) = rownames(tf)
  #Now get the top N for each group
  w = lapply(seq_len(ncol(nObs)),function(e){
             o = order(score[,e],decreasing=TRUE)
             if(sum(qvals[,e]<FDR)>=N){
               o[seq(N)]
             }else{
               o[qvals[o,e]<FDR]
             }
                 })
  #Now construct the data.frame with everything
  ww = cbind(unlist(w,use.names=FALSE),rep(seq_len(ncol(nObs)),lengths(w)))
  out = data.frame(gene = rownames(nObs)[ww[,1]],
                   cluster = colnames(nObs)[ww[,2]],
                   geneFrequency = tf[ww],
                   geneFrequencyOutsideCluster = ntf[ww],
                   geneFrequencySecondBest = sndBest[ww],
                   geneFrequencyGlobal = nTot[ww[,1]]/ncol(toc),
                   secondBestClusterName = sndBestName[ww],
                   tfidf = score[ww],
                   idf = idf[ww[,1]],
                   qval = qvals[ww],
                   stringsAsFactors=FALSE)
  return(out)
}

#' Sets clustering for SoupChannel
#' 
#' Adds or updates clustering information to meta-data table in SoupChannel object.
#'
#' @export
#' @param sc A SoupChannel object.
#' @param clusters A named vector, where entries are the cluster IDs and names are cellIDs.  If no names are provided, the order is assumed to match the order in \code{sc$metaData}.
#' @return An updated SoupChannel object with clustering information stored.
setClusters = function(sc,clusters){
  if(!all(colnames(sc$toc) %in% names(clusters))){
    if(length(clusters)!=nrow(sc$metaData)){
      stop("Invalid cluster specification.  See help.")
    }else{
      sc$metaData$clusters = clusters
    }
  }else{
    sc$metaData$clusters = clusters[rownames(sc$metaData)]
  }
  return(sc)
}

#' Manually set contamination fraction
#'
#' Manually specify the contamination fraction.
#'
#' @export
#' @param sc A SoupChannel object.
#' @param contFrac The contamination fraction.  Either a constant, in which case the same value is used for all cells, or a named vector, in which case the value is set for each cell.
#' @param forceAccept If this function returns an extremely high contamination estimate, an error is usually raised.  Setting this to TRUE forces this estimate to be used. 
#' @return A modified SoupChannel object for which the contamination (rho) has been set.
setContaminationFraction = function(sc,contFrac,forceAccept=FALSE){
  #Do some basic checks
  if(any(contFrac>1)){
    stop("Contamination fraction greater than 1 detected.  This is impossible and likely represents a failure in the estimation procedure used.")
  }else if(any(contFrac>0.5)){
    if(forceAccept){
      warning(sprintf("Extremely high contamination estimated (%.2g).  This likely represents a failure in estimating the contamination fraction.  Consider if you really want to use this value.",max(contFrac)))
    }else{
      stop(sprintf("Extremely high contamination estimated (%.2g).  This likely represents a failure in estimating the contamination fraction.  Set forceAppcept=TRUE to proceed with this value.",max(contFrac)))
    }
  }else if(any(contFrac>0.3)){
    warning(sprintf("Estimated contamination is very high (%.2g).  Check diagnostic plot for smaller peak closer to 0.",max(contFrac)))
  }
  if(length(contFrac)==1){
    sc$metaData$rho=contFrac
  }else{
    if(!all(names(contFrac) %in% rownames(sc$metaData)))
      stop("contFrac must be either of length 1 or a named vector with names matching the rows of sc$metaData")
    sc$metaData$rho[match(names(contFrac),rownames(sc$metaData))] = contFrac
  }
  return(sc)
}

#' Manually set dimension reduction for a channel
#'
#' Manually specify the dimension reduction 
#'
#' @export
#' @param sc A SoupChannel object.
#' @param DR The dimension reduction coordinates (e.g., tSNE).  This must be a data.frame, with two columns giving the two dimension reduction coordinates.  The data.frame must either have row names matching the row names of sc$metaData, or be ordered in the same order as sc$metaData.
#' @param reductName What to name the reduction (defaults to column names provided).
#' @return A modified SoupChannel object for which the dimension reduction has been set.
setDR = function(sc,DR,reductName=NULL){
  #If more than two columns, keep the first two
  if(ncol(DR)>2){
    warning(sprintf("DR has %d columns where 2 were expected.  Using first two.",ncol(DR)))
    DR = DR[,1:2]
  }
  #Check if the rownames match the metadata
  m = match(rownames(sc$metaData),rownames(DR))
  if(any(is.na(m))){
    #Can't use row-names, so the number better match
    if(nrow(DR)!=nrow(sc$metaData)){
      stop(sprintf("Rownames present in metaData not found in DR and row numbers differ (%d in metaData, %d in DR).  Each cell must have a corresponding entry in DR.",nrow(sc$metaData),nrow(DR)))
    }
    m = seq(nrow(DR))
  }
  #Should we change the names?
  if(!is.null(reductName))
    colnames(DR) = paste0(reductName,'_',1:2)
  #Add the entries in
  sc$metaData = cbind(sc$metaData,DR[m,])
  #And point them in the right place
  sc$DR = colnames(DR)
  return(sc)
}

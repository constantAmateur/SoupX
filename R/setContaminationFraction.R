#' Manually set contamination fraction
#'
#' Manually specify the contamination fraction.
#'
#' @export
#' @param sc A SoupChannel object.
#' @param contFrac The contamination fraction.  Either a constant, in which case the same value is used for all cells, or a named vector, in which case the value is set for each cell.
#' @return A modified SoupChannel object for which the contamination (rho) has been set.
setContaminationFraction = function(sc,contFrac){
  if(length(contFrac)==1){
    sc$metaData$rho=contFrac
  }else{
    if(!all(names(contFrac) %in% rownames(sc$metaData)))
      stop("contFrac must be either of length 1 or a named vector with names matching the rows of sc$metaData")
    sc$metaData$rho[match(names(contFrac),rownames(sc$metaData))] = contFrac
  }
  return(sc)
}

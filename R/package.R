#' SoupX: Profile, quantify and remove ambient RNA expression from droplet based RNA-seq
#'
#' This package implements the method described in REF.  First a few notes about nomenclature:
#' soup - Used a shorthand to refer to the ambient RNA which is contained in the input solution to droplet based RNA-seq experiments and ends up being sequenced along with the cell endogenous RNAs that the experiment is aiming to quantify.
#' channel - This refers to a single run input into a droplet based sequencing platform.  For Chromium 10X 3' sequencing there are currently 8 "channels" per run of the instrument.  Because the profile of the soup depends on the input solution, this is the minimal unit on which the soup should be estimated and subtracted.
#'
#' The essential step in performing background correction is deciding which genes are not expressed in a reasonable fraction of cells.  This is because SoupX estimates the contamination fraction by comparing the expression of these non-expressed genes in droplets containing cells to the soup defined from empty droplets.  For solid tissue, the set of Haemoglobin genes usually works well.  The key properties a gene should have are:
#' - it should be easy to identify when it is truly expressed (i.e., when it's expressed, it should be highly expressed) 
#' - it should be highly specific to a certain cell type or group of cell types so that when the expression level is low, you can be confident that the expression is coming from the soup and not a very low level of expression from the cell
#'
#' Spike-in RNAs are the best case scenario.  In the case where you do not have spike-ins and haemoglobin genes are not viable estimators, the user should begin by using the \link{plotMarkerDistribution} function to plot those genes with bi-modal distributions that have a pattern of expression across cells that is consistent with high cell-type specificity.  The user should then select a set of genes that can be used for estimation from this list.  One or two high quality genes is usually sufficient to obtain a good estimate for the average contamination level of a channel.
#'
#' @docType package
#' @name SoupX
#' @import ggplot2
#' @importFrom Matrix colSums rowSums t
#' @importFrom stats lowess spline approx
#' @importFrom grDevices rainbow
#' @importFrom stats approx dpois optimise p.adjust pbinom phyper ppois qbeta
#' @importFrom utils data read.delim setTxtProgressBar txtProgressBar head
#' @importFrom methods as is
NULL
utils::globalVariables(c('RD1','RD2','nUMIs','est','lower','upper','isLogged','MarkerGroup','Values','rho','qVals','logRatio','expSoupCnts','soupProfile'))

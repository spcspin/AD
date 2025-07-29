#' The GV-average function
#'
#' This function is designed to calculate GV-average.
#'
#' @param KA A numeric matrix of additive kinship.
#' @param KD A numeric matrix of dominant kinship.
#' @param sigmaA A numeric value of additive variance.
#' @param sigmaD A numeric value of dominant variance.
#' 
#' @return This function will return each individuals' GV-average values.
#'
#' @export
#' @examples
#' data(wheat)
#' GVaverage(KA,KD)
#
GVaverage <- function(KA,KD, sigmaA=1, sigmaD=1){
  Nc <- length(KA[,1])
  K <- sigmaA*KA+sigmaD*KD
  GVaverage <- diag(K)
  ID_GVaverage <- as.data.frame(cbind(GVaverage,colnum=seq(1,Nc)))
  ID_GVaverage$colnum <- as.numeric(ID_GVaverage$colnum)
  ID_GVaverage$GVaverage <- as.numeric(ID_GVaverage$GVaverage)
  return(ID_GVaverage)
}
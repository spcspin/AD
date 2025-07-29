#' The CDmean(v2) function
#'
#' This function is designed to calculate CDmean(v2).
#'
#' @param KA A numeric matrix of additive kinship.
#' @param KD A numeric matrix of dominant kinship.
#' @param alphaA A numeric value of additive alpha.
#' @param alphaD A numeric value of dominant alpha.
#' @param train An integer vector of which individuals are in the training set.
#' @param test An integer vector of which individuals are in the test set.
#'
#' @return This function will return the CDmean(v2) value.
#'
#' @export
#' @examples
#' data(wheat)
#' CDmeanv2(KA, KD,1:50,1:100)
#
CDmeanv2 <- function(KA, KD,  train, test, alphaA=1, alphaD=1){
  KAt <- KA[train,train]
  KA0 <- KA[test,test]
  KAt0 <- KA[train,test]
  KA0t <- KA[test,train]
  KDt <- KD[train,train]
  KD0 <- KD[test,test]
  KDt0 <- KD[train,test]
  KD0t <- KD[test,train]
  n0 <- length(test)
  nt <- length(train)
  I <- diag(nt)
  Jbar <- matrix(1/nt,nt,nt)
  M <- I-Jbar
  Gt <- alphaA*KAt+alphaD*KDt
  G0t <- alphaA*KA0t+alphaD*KD0t
  Gt <- alphaA*KAt+alphaD*KDt
  G0 <- alphaA*KA0+alphaD*KD0
  A <- diag(G0t%*%solve(M%*%Gt+I)%*%M%*%t(G0t))
  B <- diag(G0)
  CDv2 <- sum(A/B)
  return(CDv2)
}
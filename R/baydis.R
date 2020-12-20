#' @title  Bayesian discrimination of two populations
#' @description Bayesian discrimination of two populations using R
#' @param X1 a matrix or data frame that contains the explanatory variables and belongs to the first population 
#' @param X2 a matrix or data frame that contains the explanatory variables and belongs to the second population 
#' @param X  a matrix or data frame containing the explanatory variables to be classified
#' @param R1 prior probability ratio p2/p1
#' @param R2 the ratio of the two losses L(1|2)/L(2|1)
#' @param eq.var If eq.var = True, the variance is equal. Otherwise, the variances are not equal.
#' @return the classification result 
#' @examples
#' \dontrun{
#' library(MASS)
#' X1 <- iris[1:50,1:4]
#' X2 <- iris[51:100,1:4]
#' res <- baydis(X1, X2)
#' }
#' @export

baydis<-function(X1, X2, X = NULL, R1 = 1, R2 = 1, eq.var = TRUE){
  if (is.matrix(X1) != TRUE) X1<-as.matrix(X1)
  if (is.matrix(X2) != TRUE) X2<-as.matrix(X2)
  mu1 <- colMeans(X1)
  mu2 <- colMeans(X2) 
  R <- R1*R2
  if (is.null(X) == TRUE) X<-rbind(X1,X2)
  if (is.matrix(X) != TRUE) X<-as.matrix(X)
  n<-nrow(X)  
  grouping<-matrix(nrow=1,ncol = n, byrow=TRUE, 
                dimnames=list("class", 1:n))
  if (eq.var == TRUE  || eq.var == T){
    S <- var(rbind(X1,X2))
    beta <- 2*log(R)
    d1 <- mahalanobis(X, mu1, S)
    d2 <- mahalanobis(X, mu2, S)
  }
  else{
    S1 <- var(X1)
    S2 <- var(X2)
    beta <- 2*log(R)+log(det(S1)/det(S2))
    d1 <- mahalanobis(X, mu1, S1)
    d2 <- mahalanobis(X, mu2, S2)
  }
  W <- d2-d1
  for (i in 1:n){
    if (W[i] > beta)
      grouping[i] <- 1
    else
      grouping[i] <- 2
  }
  return(grouping)
}
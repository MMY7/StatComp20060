#' @title  A random sample from Pareto(a,b)
#' @description Generate random samples from Pareto(a,b) using R
#' @param n the number of random samples
#' @param a a
#' @param b b
#' @return a random sample 
#' @examples
#' \dontrun{
#' x <- rPareto(10,1,2)
#' }
#' @export 

rPareto<-function(n,a,b){
  u<-runif(n)
  x<-b/u^(1/a)
  return(x)
}



#' @title  Random variates from the rescaled Epanechnikov kernel
#' @description Generate random variates from the rescaled Epanechnikov kernel using R
#' @param n number of random numbers generated
#' @return n random variates from the rescaled Epanechnikov kernel 
#' @examples
#' \dontrun{
#' x <- rEk(1)
#' }
#' @export 

rEk<-function(n){
  u1<-runif(n,-1,1)
  u2<-runif(n,-1,1)
  u3<-runif(n,-1,1)
  x<-rep(NA, n)
  for (i in 1:n) {
    if(abs(u3[i])>abs(u1[i])&abs(u3[i])>abs(u2[i])){
      x[i]<-u2[i]
    }else{
      x[i]<-u3[i]
    }
  }
  return(x)
}




#' @title  Count Five test
#' @description Performs the Count Five test on vectors of data.
#' @param x a (non-empty) numeric vector of data values.
#' @param y an optional (non-empty) numeric vector of data values.
#' @return 1 (reject H0) or 0 (do not reject H0). 
#' @examples
#' \dontrun{
#' x <- rnorm(100, 0, 1)
#' y <- rnorm(100, 0, 0.5)
#' I <- CFtest(x,y)
#' }
#' @export 

CFtest <- function(x, y) {
  X <- x - mean(x)
  Y <- y - mean(y)
  outx <- sum(X < min(Y))+sum(X > max(Y))
  outy <- sum(Y < min(X))+sum(Y > max(X))
  I<-as.integer(max(c(outx, outy)) > 5)
  return(I)
}



#' @title  The multivariate skewness statistic
#' @description Calculate the multivariate skewness statistic using R
#' @param x samples
#' @return the multivariate skewness statistic 
#' @export 

smsk <- function(x) {
  xbar <- colMeans(x)
  xcov <- cov(x)
  n <- nrow(x)
  xSigma <- xcov*(n-1)/n
  m<-numeric(n)
  b<- sum(((t(t(x)-xbar))%*%solve(xSigma)%*%(t(x)-xbar))^3)/n^2
}


#' @title  A random walk Metropolis sampler for generating the standard Laplace distribution
#' @description A random walk Metropolis sampler for generating the standard Laplace distribution using R
#' @param x0 initial value
#' @param N length of the chain
#' @param sigma variance
#' @return a list with two components:x and acc.rate give the chain and its acceptance rate
#' @examples
#' \dontrun{
#'  rw<- Metro.rw(25,3000,2)
#' }
#' @export 

Metro.rw <- function(x0, N,sigma){
  x <- numeric(N)
  x[1] <- x0
  k <- 0
  #U(0,1)
  u <- runif(N)
  
  for (i in 2:N) {
    y <- rnorm(1, x[i-1], sigma)
    alpha <- exp(-1*abs(y))/exp(-1*abs(x[i-1]))
    if (u[i]<=alpha)
      x[i] <- y 
    else{
      x[i] <- x[i-1]
      #reject
      k <- k + 1
    }
  }
  return(list(x=x, acc.rate=1-k/N))
}



#' @title   Scalar summary statistic 
#' @description Calculate the scalar summary statistic using R
#' @param X samples
#' @return the scalar summary statistic 
#' @export 

psi.st <- function(X)
{
  Y <- t(apply(X, 1, cumsum))
  for (i in 1:nrow(X))
    Y[i,] <- Y[i,] / (1:ncol(X))
  return(Y)
}



#' @title Gelman-Rubin statistic 
#' @description Calculate Gelman-Rubin statistic using R
#' @param psi Scalar summary statistic
#' @return Gelman-Rubin statistic
#' @export 

Gelman.Rubin <- function(psi){ 
  k <- nrow(psi)
  N <- ncol(psi)
  psi.means <- rowMeans(psi) 
  B <- N * var(psi.means)
  W <- mean(apply(psi, 1, "var")) 
  Var.hat <- (N-1)*W/N + B/N
  R.hat <- Var.hat/W 
  return(R.hat)
}



#' @title  The intersection points A(k) in $(0, sqrt k)$ of the curves $S_{k-1}(a)$ and $S_k(a)$
#' @description Find the intersection points A(k) in $(0, sqrt k)$ of the curves $S_{k-1}(a)$ = P (t(k-1) > sqrt {frac{a^2(k-1)}{k-a^2}})$ and $S_k(a) = P(t(k) > frac{a^2k }{k+1-a^2})$
#' @param k k degrees of freedom
#' @return a sequence with three values,the intersection point, the value of the function evaluated at it and the intersection points
#' @examples
#' \dontrun{
#'  res <- f.intersection_find(10)
#' }
#' @export 

f.intersection_find <- function(k){
  #S_{k}(a)
  S_k = function(a) 
    1-pt(sqrt(k*a^2 /(k + 1 - a^2)), df=k)
  #S_{k-1}(a)
  S_k1 <- function(a) 
    1-pt(sqrt(  (k - 1)*a^2 / (k - a^2)), df=k-1)
  #S_{k}(a)-S_{k-1}(a)
  f = function(a) 
    S_k(a) - S_k1(a)
  #open interval
  eps <- .Machine$double.eps^0.5
  out <- uniroot(f, interval = c(0+eps, sqrt(k)-eps),tol = .Machine$double.eps^0.25)
  return(unlist(out)[1:3])
}
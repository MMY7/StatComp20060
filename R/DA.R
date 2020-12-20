#' @title  Linear Discriminant Analysis
#' @description Doing linear discriminant analysis using R 
#' @param x a matrix or data frame or Matrix containing the explanatory variables 
#' @param y a factor specifying the class for each sample
#' @param prior prior probability for each class
#' @return a list with three components:beta0 and beta1 give parameters of discrimination function and label gives the levels of classification
#' @examples
#' \dontrun{
#' n0 <- sample(1 : 150, 100)
#' iris.train <- iris[n0, ]
#' iris.lda <- lda.method(iris.train[, -5], iris.train[, 5],c(1, 1, 1)/3)
#'} 
#' @export

lda.method <- function(x, y, prior){
  p <- ncol(x)
  #Number of classes
  nclass <- length(levels(y))
  res <- mle.est(x, y)
  beta0 <- vector()
  beta1 <- matrix(rep(0, nclass*p), nrow = nclass, ncol = p)
  omega <- rowSums(prior * res$Sigma, dims = 2)
  #Inversion of omega
  Sigma.temp <- solve(omega)
  beta0 <- -1/2 * diag(res$mu %*% Sigma.temp %*% t(res$mu)) + log(prior)
  beta1 <- t(Sigma.temp %*% t(res$mu))
  list(beta0 = beta0, beta1 = beta1, label = unique(y))
}



#' @title Classify one Multivariate observation by Linear Discrimination 
#' @description Classify one Multivariate observation in conjunction with lda.method and also project data onto the linear discriminants
#' @param beta0 beta0 in discrimination function
#' @param beta1 beta1 in discrimination function
#' @param x a vector to be classified
#' @param label levels of classification
#' @return the classification result
#' @export

one.predict.lda <- function(beta0, beta1, x, label){
  nclass <- length(beta0)
  res.class <- sapply(1 : nclass, function(i){beta0[i] + t(beta1[i, ]) %*% x})
  label[which.max(res.class)]
}




#' @title  Classify Multivariate observations by Linear Discrimination
#' @description Classify multivariate observations in conjunction with lda.method and also project data onto the linear discriminants
#' @param struct the result of lda.method function
#' @param newx a matrix or data frame or Matrix to be classified
#' @return the classification result
#' @examples
#' \dontrun{
#' n0 <- sample(1 : 150, 100)
#' iris.train <- iris[n0, ]
#' iris.test <- iris[-n0, ]
#' iris.lda <- lda.method(iris.train[, -5], iris.train[, 5],c(1, 1, 1)/3)
#' pre.lda <- lda.predict(iris.lda, iris.test[, 1 : 4])
#'} 
#' @export

lda.predict <- function(struct, newx){
  beta0 <- struct$beta0
  beta1 <- struct$beta1
  label <- struct$label
  newx <- as.matrix(newx)
  newn <- nrow(newx)
  newy <- sapply(1 : newn, function(i){one.predict.lda(beta0, beta1, newx[i, ], label)})
  newy <- factor(newy)
  levels(newy) <- levels(label)
  newy
}






#' @title  Compute mu and Sigma
#' @description Compute mle of lda/qda model, return mle of mu and Sigma using R
#' @param x a matrix or data frame or Matrix containing the explanatory variables 
#' @param y a factor specifying the class for each sample
#' @return a list with two components: mu and Sigma 
#' @examples
#' \dontrun{
#' n0 <- sample(1 : 150, 100)
#' iris.train <- iris[n0, ]
#' mle.est(iris.train[, -5], iris.train[, 5])
#'} 
#' @export

mle.est <- function(x, y){
  x <- as.matrix(x)
  nclass <- length(levels(y))
  #Number of samples for each class
  class.temp <- as.data.frame(table(y))
  n <- nrow(x)
  p <- ncol(x)
  mu <- matrix(rep(0, nclass*p), nrow = nclass, ncol = p)
  Sigma <- array(data = rep(0, p * p * nclass), dim = c(p, p, nclass))
  for (i in 1 : n) {
    temp <- which(class.temp[,1] == y[i])
    mu[temp, ] <- mu[temp, ] + x[i, ]
  }
  mu <- mu/matrix(rep(class.temp[,2], p), ncol = p)
  for (i in 1 : n) {
    temp <- which(class.temp[,1] == y[i])
    Sigma[, , temp] <- Sigma[, , temp] + ((x[i, ] - mu[temp, ]) %*% t(x[i, ] - mu[temp, ]))/class.temp[temp, 2] }
  list(mu = mu, Sigma = Sigma)
}







#' @title Quadratic discriminant analysis
#' @description Doing Quadratic discriminant analysis using R
#' @param x a matrix or data frame or Matrix containing the explanatory variables
#' @param y a factor specifying the class for each sample
#' @param prior prior probability for each class
#' @return a list with four components:beta0 and beta1 give parameters of discrimination function, omega and label give the inverse of the covariance matrix*(-1/2) for each level and the levels of classification
#' @examples
#' \dontrun{
#' n0 <- sample(1 : 150, 100)
#' iris.train <- iris[n0, ]
#' iris.qda <- qda.method(iris.train[, -5], iris.train[, 5],c(1, 1, 1)/3)
#'} 
#' @export

qda.method <- function(x, y, prior){
  p <- ncol(x)
  nclass <- length(unique(levels(y)))
  res <- mle.est(x, y)
  beta0 <- vector()
  beta1 <- matrix(rep(0, p * nclass), nrow = nclass, ncol = p)
  omega <- array(data = rep(0, p * p * nclass), dim = c(p, p, nclass))
  for (i in 1 : nclass) {
    mu.temp <- res$mu[i, ]
    Sigma.temp <- solve(res$Sigma[, , i])
    beta0[i] <- -1/2 * t(mu.temp) %*% Sigma.temp %*% mu.temp + log(prior[i])
    beta1[i, ] <- Sigma.temp %*% mu.temp
    omega[, , i] <- -1/2 * Sigma.temp
  }
  list(beta0 = beta0, beta1 = beta1, omega = omega, label = unique(y))
}




#' @title Classify one multivariate observation from Quadratic Discriminant Analysis  
#' @description Classify one multivariate observation in conjunction with qda.method
#' @param beta0 beta0 in discrimination function
#' @param beta1 beta1 in discrimination function
#' @param omega the inverse of the covariance matrix*(-1/2) for each level
#' @param x a vector to be classified
#' @param label levels of classification
#' @return the classification result
#' @export

one.predict.qda <- function(beta0, beta1, omega, x, label){
  nclass <- length(beta0)
  res.class <- vector()
  res.class <- sapply(1 : nclass, function(i){beta0[i] + t(beta1[i, ]) %*% x + t(x) %*%omega[, , i] %*% x})
  label[which.max(res.class)]
}


  

#' @title Classify from Quadratic Discriminant Analysis 
#' @description Classify multivariate observations in conjunction with qda.method
#' @param struct the result of qda.method function
#' @param newx a matrix or data frame or Matrix to be classified
#' @return the classification result
#' @export

qda.predict <- function(struct, newx){
  beta0 <- struct$beta0
  beta1 <- struct$beta1
  omega <- struct$omega
  label <- struct$label
  newx <- as.matrix(newx)
  newn <- nrow(newx)
  newy <- sapply(1 : newn, function(i){one.predict.qda(beta0, beta1, omega, newx[i, ],label)})
  newy <- factor(newy)
  levels(newy) <- levels(label)
  newy
}

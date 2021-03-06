% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/DA.R
\name{lda.predict}
\alias{lda.predict}
\title{Classify Multivariate observations by Linear Discrimination}
\usage{
lda.predict(struct, newx)
}
\arguments{
\item{struct}{the result of lda.method function}

\item{newx}{a matrix or data frame or Matrix to be classified}
}
\value{
the classification result
}
\description{
Classify multivariate observations in conjunction with lda.method and also project data onto the linear discriminants
}
\examples{
\dontrun{
n0 <- sample(1 : 150, 100)
iris.train <- iris[n0, ]
iris.test <- iris[-n0, ]
iris.lda <- lda.method(iris.train[, -5], iris.train[, 5],c(1, 1, 1)/3)
pre.lda <- lda.predict(iris.lda, iris.test[, 1 : 4])
} 
}

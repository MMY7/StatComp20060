## ----setup, include = FALSE---------------------------------------------------
knitr::opts_chunk$set(
  collapse = TRUE,
  comment = "#>"
)

## -----------------------------------------------------------------------------
library(StatComp20060)
#LDA
n0 <- sample(1 : 150, 100)
iris.train <- iris[n0, ]
iris.test <- iris[-n0, ]
iris.lda <- lda.method(iris.train[, -5], iris.train[, 5],c(1, 1, 1)/3)
pre.lda <- lda.predict(iris.lda, iris.test[, 1 : 4])
#QDA
n0 <- sample(1 : 150, 100)
iris.train <- iris[n0, ]
iris.qda <- qda.method(iris.train[, -5], iris.train[, 5],c(1, 1, 1)/3)
pre.qda <- qda.predict(iris.qda, iris.test[, 1 : 4]) 

## -----------------------------------------------------------------------------
library(StatComp20060)
pi <- piSugar(1000)


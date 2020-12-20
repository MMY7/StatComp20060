## ----setup, include = FALSE---------------------------------------------------
knitr::opts_chunk$set(
  collapse = TRUE,
  comment = "#>"
)

## -----------------------------------------------------------------------------
library(StatComp20060)
library(MASS)
library(caret)
library(ggplot2)
set.seed(1234)
n0 <- sample(1 : 150, 50)
iris.train <- iris[-n0, ]
iris.test <- iris[n0, ]
ptr <- ggplot(data = iris.train, aes(Petal.Length, Petal.Width, color = Species))
ptr <- ptr + geom_point() + theme_bw() + theme(panel.grid = element_blank())
ptr
##LDA
iris.lda <- lda.method(iris.train[, -5], iris.train[, 5],c(1, 1, 1)/3)
pre.lda <- lda.predict(iris.lda, iris.test[, 1 : 4])
print(iris.lda)
#misclassification error
sum(pre.lda != iris.test[, 5])/length(iris.test[, 5])
confusionMatrix(pre.lda, iris.test[, 5])$table
iris.LDA <- lda(Species~., data = iris.train, prior = c(1,1,1)/3)
pre.LDA <- predict(iris.LDA, newdata = iris.test)
#compare with lda
sum(pre.lda != pre.LDA$class)
##QDA
iris.qda <- qda.method(iris.train[, -5], iris.train[, 5],c(35/100, 35/100, 30/100))
pre.qda <- qda.predict(iris.qda, iris.test[, 1 : 4])
#misclassification error
sum(pre.qda != iris.test[, 5])/length(iris.test[, 5])
confusionMatrix(pre.qda, iris.test[, 5])$table
#compare with qda
iris.QDA <- qda(Species~., data = iris.train)
pre.QDA <- predict(iris.QDA, newdata = iris.test)
sum(pre.qda != pre.QDA$class)

## -----------------------------------------------------------------------------
library(StatComp20060)
pi <- piSugar(1000)
print(pi)


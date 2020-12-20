## ----setup, include = FALSE---------------------------------------------------
knitr::opts_chunk$set(
  collapse = TRUE,
  comment = "#>"
)

## -----------------------------------------------------------------------------
library("datasets")
data("trees")
hist(trees$Height, main = "Histogram of Height",xlab = "Height")
plot(density(trees$Height),main ="Density of Height")
plot(Volume~Girth,data=trees,pch=20,col='blue')
model=lm(Volume~Girth,data=trees)
abline(model,lty=3)

## -----------------------------------------------------------------------------
library(DAAG)
attach(leafshape17)
f <- function(a, v) {
  #Andrews curve f(a) for a data vector v in R^3
  v[1]/sqrt(2) + v[2]*sin(a) + v[3]*cos(a)
}
#scale data to range [-1, 1]
x <- cbind(bladelen, petiole, bladewid)
n <- nrow(x)
mins <- apply(x, 2, min) #column minimums
maxs <- apply(x, 2, max) #column maximums
r <- maxs - mins #column ranges
y <- sweep(x, 2, mins) #subtract column mins
y <- sweep(y, 2, r, "/") #divide by range
x <- 2 * y - 1 #now has range [-1, 1]
#set up plot window, but plot nothing yet
plot(0, 0, xlim = c(-pi, pi), ylim = c(-3,3),
     xlab = "t", ylab = "Andrews Curves",
     main = "", type = "n")
#now add the Andrews curves for each observation
#line type corresponds to leaf architecture
#0=orthotropic, 1=plagiotropic
a <- seq(-pi, pi, len=101)
dim(a) <- length(a)
for (i in 1:n) {
  g <- arch[i] + 1
  y <- apply(a, MARGIN = 1, FUN = f, v = x[i,])
  lines(a, y, lty = g)
}
legend(3, c("Orthotropic", "Plagiotropic"), lty = 1:2)

## -----------------------------------------------------------------------------
data(swiss)
model=lm(Fertility ~.,data=swiss)
library(knitr)
kable(summary(model)$coef,digits=3)

## -----------------------------------------------------------------------------
P <- c(29.6, 24.3, 28.5, 32)
T <- c(27.3, 32.6, 30.8, 34.8)
S <- c(5.8, 6.2, 11, 8.3)
E <- c(21.6, 17.4, 18.3, 19)
C <- c(29.2, 32.8, 25, 24.2)
#glue the columns together in a data frame
x <- data.frame(P, T, S, E, C)
#now stack the data for ANOVA
y <- stack(x)
names(y) <- c("Binding", "Antibiotic")

## -----------------------------------------------------------------------------
#check the default formula
print(formula(y)) #default formula is right one

## -----------------------------------------------------------------------------
lm(y)
anova(lm(y))
kable(anova(lm(y)),digits=3)

## -----------------------------------------------------------------------------
#the function to generate a random sample from Pareto(a,b)
rPareto<-function(n,a,b){
  u<-runif(n)
  x<-b/u^(1/a)
  return(x)
}
#
set.seed(123)
n<-100
a<-2
b<-2
x<-rPareto(100,2,2)
#the density histogram of the sample from Pareto(2,2) with the theoretical density superimposed
hist(x,prob=TRUE,xlab="x",breaks = 100,xlim = c(0, 100),main="Pareto(2,2)")
y <- seq(2,100,0.01)
lines(y,b^a*a*y^(-(a+1)))


## -----------------------------------------------------------------------------
#the function to generate random variates from f_e
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
#generate 10000 random observations
set.seed(246)
n<-100
x<-rEk(100)
#the histogram density estimate
hist(x,prob = TRUE,xlab ="x",breaks = 40,xlim = c(-1.2, 1.2),main = expression(f(x)==3*(1-x^2)/4))
lines(density(x))
y <- seq(-1,1,0.01)
lines(y,3*(1-y^2)/4, col="red")


## -----------------------------------------------------------------------------
#generate 1000 random observations from the mixture 
set.seed(137)
n<-100
r<-4
beta<-2
lambda<-rgamma(n, r, beta)
x<-rexp(n, lambda)
#the density histogram of the sample with the Pareto density curve superimposed
hist(x,prob = TRUE,xlab ="x",breaks = 40,xlim = c(0, 20),main = expression(f(x)==64/(2+y)^5))
y<-seq(0,20,0.01)
lines(y,64/(2+y)^5)


## -----------------------------------------------------------------------------
set.seed(1234)
x<-runif(100,0,pi/3)
#Monte Carlo estimate
sinMC<-mean(pi*sin(x)/3) 
# Exact value
sinEV<-1-cos(pi/3) 
print(c(sinMC,sinEV))

## -----------------------------------------------------------------------------
f<-function(x){
  y<-exp(x)
  return(y)
}
theta_MC<-function(n,duiou = TRUE)
{
  x1<-runif(n/2)
  if(duiou == TRUE)
    x2<-1-x1
  else
    x2<-runif(n/2)
  uniform<-c(x1,x2)
  MC<-mean(f(uniform))
}
m<-10000
set.seed(134)
theta_MC0=theta_MC(m,duiou = FALSE)
set.seed(134)
theta_MC1=theta_MC(m,duiou = TRUE)
#estimators and exact value
print(c(theta_MC0,theta_MC1,exp(1)-1))
#variance
set.seed(1220)
num<-100
eMC0<-numeric(num)
eMC1<-numeric(num)
for(i in 1:num)
{
  eMC0[i]<-theta_MC(m,duiou = FALSE)
  eMC1[i]<-theta_MC(m,duiou = TRUE)
}
print(c(var(eMC0),var(eMC1)))
#percent reduction
print(100*(1-var(eMC1)/var(eMC0)))
#compare
print(c(100*(1-var(eMC1)/var(eMC0)),100-100*(-3*exp(2)+10*exp(1)-5)/(2*exp(2)-2-4*(exp(1)-1)^2)*2))


## -----------------------------------------------------------------------------
#plot the graphs
x=seq(1,4,0.01)
g=x^2*exp(-x^2/2)/(2*pi)
plot(x,g,type="l", lwd=3,col="black",ylim=c(0,1),ylab = "function")
f1=2*exp(-(x-1)^2/2)/sqrt(2*pi)
lines(x,f1,lty=1,col="blue")
f2=exp(-x+1)
lines(x,f2,lty=1,col="green")
f3=1/(1 + x^2)/(pi/4)
lines(x,f3,lty=1,col="red")
legend(3.6, 0.95, c("g","f1","f2","f3"), cex=0.9, col=c("black","blue","green","red"), pch=21:24, lty=c(1,1,1,1), bg="grey", title.col="black")


x=seq(1,10,0.01)
g=x^2*exp(-x^2/2)/(2*pi)
f1=2*exp(-(x-1)^2/2)/sqrt(2*pi)
plot(x,g/f1,type="l", lwd=2,col="blue",ylim=c(0,0.8),ylab = "function")
f2=exp(-x+1)
lines(x,g/f2,lty=1,col="green")
f3=1/(1 + x^2)/(pi/4)
lines(x,f3,lty=1,col="red")
legend(8.2, 0.8, c("g/f1","g/f2","g/f3"), cex=0.9, col=c("blue","green","red"), pch=21:23, lty=c(1,1,1), bg="grey", title.col="black")


#The estimates and the corresponding standard errors
set.seed(520)
m <- 100
theta.hat <- se <- numeric(3)
g <- function(x)
{ x^2*exp(-x^2/2)/(2*pi)*(x>1)}  
f1 <- function(x)
{2*exp(-(x-1)^2/2)/sqrt(2*pi)}
f2 <- function(x)
{exp(-x+1)}
f3 <- function(x)
{1/(1 + x^2)/(pi/4)}
#f1 normal transform
u1<-rnorm(m)
x1<-abs(u1)+1
gf1<- g(x1)/f1(x1)
theta.hat[1] <- mean(gf1)
se[1] <- sd(gf1)

#f2 inverse 
u <-runif(m)
x2<-1-log(1-u)
gf2 <- g(x2)/f2(x2)
theta.hat[2] <- mean(gf2)
se[2] <- sd(gf2)

#f3 inverse 
u <- runif(m) 
x3 <- tan(pi*(u+1)/4)
gf3 <- g(x3)/f3(x3)
theta.hat[3] <- mean(gf3)
se[3] <- sd(gf3)
print(rbind(theta.hat,se))


## -----------------------------------------------------------------------------
##Stratified Importance Sampling
set.seed(124)
n=100
k=5
X<-matrix(data=NA, nrow = n/k, ncol = k)
thetai<-matrix(data=NA, nrow = n/k, ncol = k)
thetasi<-Varsi<-numeric(k)
  for(j in 1:k)
  {
    #g/f on each subinterval
    gfj<-function(x)
    {(exp((-j+1)/5) - exp((-j)/5))/(1+x^2)}
    u <- runif(n/k)
    #randoms on each subinterval
    X[,j] <--log(exp((-j+1)/5) - u*(exp((-j+1)/5) - exp((-j)/5)))
    # Integration on each subinterval
    thetai[,j]<- gfj(X[,j])
    Varsi[j]<-var(thetai[,j])
    thetasi[j]<-mean(thetai[,j])
  }
 theta<-sum(thetasi)
 Var<-(k/n)*sum(Varsi)
 se<-sqrt(Var)
 print(c(theta,se))

## -----------------------------------------------------------------------------
##Importance Sampling
set.seed(137)
n=100#sample size
Moni=50#number of times to repeat the estimation
theta.hat1<-theta.hat2<-numeric(Moni)
gf<-function(x)
{(1-exp(-1))/(1+x^2)}
for(mo in 1:Moni)
{U<-runif(n)
x<--log(1-U*(1-exp(-1)))
theta.hat1[mo]<- mean(gf(x))
}
theta1<-mean(theta.hat1)
se1=sd(theta.hat1)
Var1=var(theta.hat1)
##Stratified Importance Sampling
k=5
X<-matrix(data=NA, nrow = n/k, ncol = k)
theta.si<-matrix(data=NA, nrow = Moni, ncol = k)
for(mo in 1:Moni)
{
  for(j in 1:k)
  {
    #g/f on each subinterval
    gfj<-function(x)
    {(exp((-j+1)/5) - exp((-j)/5))/(1+x^2)}
    u <- runif(n/k)
    #randoms on each subinterval
    X[,j] <--log(exp((-j+1)/5) - u*(exp((-j+1)/5) - exp((-j)/5)))
    # Integration on each subinterval
    theta.si[mo,j]<- mean(gfj(X[,j]))
    theta.hat2[mo]<-sum(theta.si[mo,])
  }
}
theta2<-mean(theta.hat2)
se2=sd(theta.hat2)
Var2=var(theta.hat2) #variance of integrated items
rv=(Var1-Var2)/Var1
theta<-c(theta1,theta2)
Var<-c(Var1,Var2)
se<-c(se1,se2)
print(rbind(theta,Var,se))
print(rv)

## -----------------------------------------------------------------------------
set.seed(676)
## confidence interval for mu
#  sample size
n<-50
# Times to repeat MC experiment
Tm<-100
#100(1-a)% confidence interval
a<-0.05
u<-5
I.ci<-numeric(Tm)
for (i in 1:Tm) 
{
  #random sample
  x<-rlnorm(n,u,1)
  #the confidence interval
  CI<-c(qt(a/2,n-1)*sd(log(x))/sqrt(n),-qt(a/2,n-1)*sd(log(x))/sqrt(n))+mean(log(x))
  I.ci[i]<- u>CI[1]&&u<CI[2]
}
#confidence level
p.ci<-mean(I.ci)
print(p.ci)


## -----------------------------------------------------------------------------
set.seed(77)
## t-interval for mean
# Required sample size
n<-20 
# Times to repeat MC experiment
Tm<-100
#100(1-a)% confidence interval
a<-0.05
I1<-numeric(Tm)
for (i in 1:Tm) 
{
  #random sample
  x<-rchisq(n,2)
  #the confidence interval
  CI<-c(qt(a/2,n-1)*sd(x)/sqrt(n),-qt(a/2,n-1)*sd(x)/sqrt(n))+mean(x)
  I1[i]<- 2>CI[1]&&2<CI[2]
}
#confidence level
p1.ci<-mean(I1)

## Confidence interval for  variance in example 6.4
I2<-numeric(Tm)
for (i in 1:Tm){
  x<-rnorm(n,0,2)
  UCL<-(n-1)*var(x)/qchisq(a,n-1)
  I2[i]<- 4<UCL
}
p2.ci<-mean(I2)
cbind(p1.ci,p2.ci)

## -----------------------------------------------------------------------------
## Generate sample from Beta and t
# Sample size
set.seed(18)
n<-100 
# replication
m<-100  
# Sample skewness function
ssk<-function(x){
  q3<-mean((x-mean(x))^3)
  q2<-mean((x-mean(x))^2)
  y<-q3/q2^1.5
  return(y)
}
#parameters
theta<-c(seq(1,10,0.5)) 
N <- length(theta)
power1<-power2<-numeric(N)
#the quantile
z<-qnorm(0.975,0,sqrt(6/n))
for (i in 1:N)
{#theta
  I1<-I2<-numeric(m)
  for (j in 1:m)
  {#replicate experiment
    x1<-rbeta(n,theta[i],theta[i])
    x2<-rt(n,theta[i])
    ske1<-ssk(x1)
    ske2<-ssk(x2)
    I1[j]<- as.integer(abs(ske1)>=z)
    I2[j]<-as.integer(abs(ske2)>=z)
  }
  power1[i]<-mean(I1)
  power2[i]<-mean(I2)
}
se1 <- sqrt(power1*(1-power1)/m)
se2 <- sqrt(power2*(1-power2)/m)
# power curve
par(mfrow=c(1,2))
plot(theta,power1,xlab="α",ylab="power1",lwd=2,col="blue",main="Beta(α,α)",type="l")
lines(theta, power1+se1, lty = 4)
lines(theta, power1-se1, lty = 4)
plot(theta,power2,xlab="ν",ylab="power2",lwd=2,col="blue",main="t(ν)",type="l")
lines(theta, power2+se2, lty = 3)
lines(theta, power2-se2, lty = 3)

## -----------------------------------------------------------------------------
#small,medium and large sample
set.seed(77)
size<-c(10,100,1000)
N<-length(size)
m<-100
#parameters
sigma1 <- 1
sigma2 <- 1.5
#the quantile
a<-0.055
CFtest <- function(x, y) {
  X <- x - mean(x)
  Y <- y - mean(y)
  outx <- sum(X < min(Y))+sum(X > max(Y))
  outy <- sum(Y < min(X))+sum(Y > max(X))
  I<-as.integer(max(c(outx, outy)) > 5)
  return(I)
}
power1<-power2<-numeric(N)
for(i in 1:N)
{ 
  I1<-I2<-numeric(m)
  for (j in 1:m)
  {
   x <- rnorm(size[i], 0, sigma1)
   y <- rnorm(size[i], 0, sigma2)
   F.test<-var.test(x,y)$p.value
   I1[j]<-CFtest(x, y)
   I2[j]<-as.integer(F.test<=a)
  }
  power1[i]<-mean(I1)
  power2[i]<-mean(I2)  
} 
print(rbind(power1,power2))

## -----------------------------------------------------------------------------
set.seed(777)
library(MASS)
size <- c(10, 20, 30, 50, 100,500)
N <- length(size)
m <- 200
# Sample multivariate  skewness function
smsk <- function(x) {
  xbar <- colMeans(x)
  xcov <- cov(x)
  n <- nrow(x)
  xSigma <- xcov*(n-1)/n
  m<-numeric(n)
  b<- sum(((t(t(x)-xbar))%*%solve(xSigma)%*%(t(x)-xbar))^3)/n^2
}
#estimates of Type I error rate
TIer <- numeric(N)
for (i in 1:N) 
{
  I <- numeric(m)
  for (j in 1:m) 
  {
    x <- mvrnorm(size[i], rep(0, 2), diag(2))
    #test reject or not reject
    I[j] <- as.integer(size[i]*abs(smsk(x))/6 >= qchisq(0.95,4))
  }
  TIer[i] <- mean(I) 
}
print(TIer)
plot(size,TIer,xlab="size",ylab="power",lwd=2,col="black",type="l")
points(size,TIer,type="b")

## -----------------------------------------------------------------------------
set.seed(21)
library(MASS)
n <- 20
m <- 100
e <- c(seq(0, .2, .02), seq(.2, 1, .1))
N <- length(e)
smsk <- function(x) {
  xbar <- colMeans(x)
  xcov <- cov(x)
  n <- nrow(x)
  xSigma <- xcov*(n-1)/n
  m<-numeric(n)
  b<- sum(((t(t(x)-xbar))%*%solve(xSigma)%*%(t(x)-xbar))^3)/n^2
}
power <- numeric(N)
for (i in 1:N) 
{#epsilon
  I1 <- numeric(m)
  for (j in 1:m) 
  {#replicate experiment
   c<- sample(c(1,100), replace = TRUE,size = n, prob = c(1-e[i], e[i]))
   x<-matrix(nr=n,nc=2) 
     for(s in 1:n)
     {x[s,] <- mvrnorm(1, rep(0, 2), c[s]*diag(2))}
   I1[j] <- as.integer(n*abs(smsk(x))/6 >= qchisq(0.9,4))
  }
  power[i] <- mean(I1)
}
#plot power against ε-contaminated multivariate normal
plot(e,power,xlab="ε",ylim = c(0,1),ylab="power",lwd=2,col="blue",type="l")
se <- sqrt(power * (1-power) / m) #add standard errors
lines(e, power+se, lty = 4)
lines(e, power-se, lty = 4)

## -----------------------------------------------------------------------------
library(bootstrap)
data("law")
law<-as.matrix(law)
law82<-as.matrix(law82)
law82<-law82[,2:3]
R.hat<-cor(law[,1],law[,2])
#bias Jackknife
n<-nrow(law)
R.jack<-numeric(n)
for (i in 1:n)
{R.jack[i]<-cor(law[-i,1],law[-i,2])}
#estimate of bias
ebias.jack<-(n-1)*(mean(R.jack)-R.hat)
#standard error jackknife
#estimate of standard error
ese.jack<-sqrt((n-1)*mean((R.jack-mean(R.jack))^2))
out<-data.frame(ebias.jack,ese.jack)
knitr::kable(out)

## -----------------------------------------------------------------------------
#Bootstrap estimate of bias of a ratio estimate
set.seed(187)
library(bootstrap)
library(boot)
data("aircondit")
aircondit<-as.matrix(aircondit)
#the statistic function
theta.boot <- function(Data, index) 
{
  x<-Data[index]
  mean(x) 
}
boot.result <- boot(aircondit, statistic = theta.boot, R = 100)
print(boot.result)
#Compute 95% bootstrap confidence intervals
print(boot.ci(boot.result,type = c("norm", "basic", "perc","bca")))

## -----------------------------------------------------------------------------
set.seed(1234)
# Sample size
n1 <- 20
n2 <- 30
N <- n1+n2
m <- 100
pv <-pv1<- numeric(m)
#Confidence level
a <- 0.05
R<-999
ex.piont <- function(x,y) {
  X.c <- x-mean(x)
  Y.c <- y-mean(y)
  outx <- sum(X.c>max(Y.c))+sum(X.c<min(Y.c))
  outy <- sum(Y.c>max(X.c))+sum(Y.c<min(X.c))
  return(max(c(outx,outy)))
}
pv <- replicate(m,expr = {#Calculate pvalue
   x <- rnorm(n1,0,1)
   y <- rnorm(n2,0,1)
   z <- c(x,y)
   Cs0 <- ex.piont(x,y)
   Cs <- numeric(R)
   for (i in 1:R)
   {#Permutation
   index <- sample(1:N,n1)
   x.p <- z[index]
   y.p <- z[-index]
   Cs[i] <- ex.piont(x.p,y.p)
   }
   p <- mean(c(Cs0,Cs)>=Cs0)
   return(p)
})
#Estimate tpye 1 error
t1er <- mean(pv<= a)
print(t1er)

## ----eval=FALSE---------------------------------------------------------------
#  library(RANN)
#  library(boot)
#  library(energy)
#  library(Ball)
#  library(MASS)
#  n1 <- n2 <- 30
#  n <- n1 + n2
#  d <- 2
#  m <- 100
#  N <- c(n1, n2)
#  R <- 999
#  J <- 3
#  alpha <- 0.05
#  NNJ <- function(z, index, sizes,J)
#  {
#    n1 <- sizes[1]
#    n2 <- sizes[2]
#    n<- n1 + n2
#    if(is.vector(z))
#      z <- data.frame(z, 0)
#    z <- z[index, ]
#    NN <- nn2(z, k = J+1)
#    block1 <- NN$nn.idx[1:n1, -1]
#    block2 <- NN$nn.idx[(n1+1):n,-1 ]
#    I1 <- sum(block1 < n1 + .5)
#    I2 <- sum(block2 > n1 + .5)
#    return((I1 + I2) / (J * n))
#  }
#  eqdist.nn <- function(z,sizes,J)
#  {
#    boot.obj <- boot(data = z, statistic = NNJ,
#                         sim = "permutation", R = 999, sizes = sizes,J=J)
#    TS <- c(boot.obj$t0,boot.obj$t)
#    p.value <- mean(TS >= boot.obj$t0)
#    list(statistic = boot.obj$t0,p.value = p.value)
#  }
#  #Unequal variances and equal expectations
#  p.values <- matrix(NA,m,3)
#  for(i in 1:m)
#  {
#    x <- mvrnorm(n = n1, rep(0, 2), diag(1,d))
#    y <- mvrnorm(n = n2, rep(0, 2), diag(2.5,d))
#    z<-rbind(x,y)
#    p.values[i,1] <- eqdist.nn(z,N,J)$p.value
#    p.values[i,2] <- eqdist.etest(z,sizes = N,R = R)$p.value
#    p.values[i,3] <- bd.test(x = x,y = y,num.permutations = 999,seed = i*125)$p.value
#  }
#  pwr1 <- colMeans(p.values<alpha)
#  print(pwr1)
#  #Unequal variances and unequal expectations
#  p.values <- matrix(NA,m,3)
#  for(i in 1:m)
#  {
#    x <- mvrnorm(n=n1, c(1, 0.5), diag(1,d))
#    y <- mvrnorm(n=n2, c(1, 0), diag(0.5,d))
#    z <- rbind(x,y)
#    p.values[i,1] <- eqdist.nn(z,N,J)$p.value
#    p.values[i,2] <- eqdist.etest(z,sizes = N,R = R)$p.value
#    p.values[i,3] <- bd.test(x = x,y = y,num.permutations=999,seed = i*145)$p.value
#  }
#  pwr2 <- colMeans(p.values<alpha)
#  print(pwr2)
#  #t distribution with 1 df, bimodel distribution
#  p.values <- matrix(NA,m,3)
#  for(i in 1:m)
#  {
#    x <- matrix(rt(n1*d,1),ncol=d);
#    mix <- sample(c(0, 1), size = n2*d, prob = c(0.3, 0.7), replace = T)
#    y <- rnorm(n2*d, mean = ifelse(mix == 0, 0, 1.5), sd = ifelse(mix == 0, 2, 4))
#    y <- matrix(y,ncol=d)
#    z <- rbind(x,y)
#    p.values[i,1] <- eqdist.nn(z,N,J)$p.value
#    p.values[i,2] <- eqdist.etest(z,sizes = N,R = R)$p.value
#    p.values[i,3] <- bd.test(x = x,y = y,num.permutations=999,seed = i*145)$p.value
#  }
#  pwr3 <- colMeans(p.values<alpha)
#  print(pwr3)
#  #Unbalanced samples(1:10)
#  n1 <- 10
#  n2 <- 100
#  n <- n1 + n2
#  N <- c(n1, n2)
#  p.values <- matrix(NA,m,3)
#  for(i in 1:m)
#  {
#    x <- mvrnorm(n=n1, c(1, 0.5), diag(1,d))
#    y <- mvrnorm(n=n2, c(1, 0), diag(0.5,d))
#    z <- rbind(x,y)
#    p.values[i,1] <- eqdist.nn(z,N,J)$p.value
#    p.values[i,2] <- eqdist.etest(z,sizes = N,R = R)$p.value
#    p.values[i,3] <- bd.test(x = x,y = y,num.permutations=999,seed = i*125)$p.value
#  }
#  pwr4 <- colMeans(p.values<alpha)
#  print(pwr4)
#  out<-rbind(pwr1,pwr2,pwr3,pwr4)
#  colnames(out)<-c("NN","Energy","Ball")
#  knitr::kable(out)

## -----------------------------------------------------------------------------
set.seed(127)
#x0 initial value 
Metro.rw <- function(x0, N,sigma)
{
  x <- numeric(N)
  x[1] <- x0
  k <- 0
  #U(0,1)
  u <- runif(N)
 
  for (i in 2:N) 
  {
    y <- rnorm(1, x[i-1], sigma)
    alpha <- exp(-1*abs(y))/exp(-1*abs(x[i-1]))
    if (u[i]<=alpha)
      x[i] <- y 
    else
    {
        x[i] <- x[i-1]
        #reject
        k <- k + 1
    }
  }
  return(list(x=x, acc.rate=1-k/N))
}
x0 <- 20
#length N
N <- 300
#n chains
n <- 4 
sigma <- c(0.05, 0.5, 2, 16)
Metro.rw1 <- Metro.rw(x0, N, sigma[1])
Metro.rw2 <- Metro.rw(x0, N, sigma[2])
Metro.rw3 <- Metro.rw(x0, N, sigma[3])
Metro.rw4 <- Metro.rw(x0, N, sigma[4])
#accept
acc.rate<-c(Metro.rw1$acc.rate, Metro.rw2$acc.rate, Metro.rw3$acc.rate, Metro.rw4$acc.rate)
acc.rate<-round(acc.rate,3)
ACC<-data.frame(sigma=sigma,Acceptance_rate=acc.rate)
knitr::kable(ACC)

## -----------------------------------------------------------------------------
#par(mfrow=c(2,2))  
Metro.rw = cbind(Metro.rw1$x, Metro.rw2$x, Metro.rw3$x,  Metro.rw4$x)
for (i in 1:4)
  plot(Metro.rw[,i], type="l",sub=bquote(sigma == .(sigma[i])),xlab = "n",ylab="X", ylim=range(Metro.rw[,i]))

## -----------------------------------------------------------------------------
#Q1
set.seed(123)
#k chains
k <- 4 
#length N
N <- 150
#burn-in length
N0 <- 10
#initial values
x0 <- c(-10, -5, 5, 10)
sigma <- c(0.05, 0.5, 2, 16)
#generate chains
Metro.rw <- function(x0, N,sigma)
{
  x <- numeric(N)
  x[1] <- x0
  k <- 0
  #U(0,1)
  u <- runif(N)
 
  for (i in 2:N) 
  {
    y <- rnorm(1, x[i-1], sigma)
    alpha <- exp(-1*abs(y))/exp(-1*abs(x[i-1]))
    if (u[i]<=alpha)
      x[i] <- y 
    else
    {
        x[i] <- x[i-1]
        #reject
        k <- k + 1
    }
  }
  return(list(x=x, acc.rate=1-k/N))
}
#psi statistics
psi.st <- function(X)
{
  Y <- t(apply(X, 1, cumsum))
  for (i in 1:nrow(X))
    Y[i,] <- Y[i,] / (1:ncol(X))
  return(Y)
}
#Gelman Rubin statistic
Gelman.Rubin <- function(psi) 
{ 
  k <- nrow(psi)
  N <- ncol(psi)
  psi.means <- rowMeans(psi) 
  B <- N * var(psi.means)
  W <- mean(apply(psi, 1, "var")) 
  Var.hat <- (N-1)*W/N + B/N
  R.hat <- Var.hat/W 
  return(R.hat)
}
#calculate Rhat
rhat <- function(sigma,x0,N,N0)
{
  X <- matrix(NA, nrow=k, ncol=N)
  for (i in 1:k)
  X[i, ] <- Metro.rw(x0[i], N, sigma)$x
  psi <- psi.st(X)
  #par(mfrow=c(2,2))
  for (i in 1:k)
  plot(psi[i, (N0+1):N], type="l",sub=bquote(sigma == .(sigma)),xlab=i, ylab=bquote(psi))
  R.hat <- numeric(N)
  for (j in (N0+1):N)
  R.hat[j] <- Gelman.Rubin(psi[,1:j])
  R.hat <- R.hat[(N0+1):N]
  return(R.hat)
}
#number of sigma is 4
R <- matrix(NA,nrow=4, ncol=N-N0)
for(i in 1:4)
  R[i,] <- rhat(sigma[i],x0,N,N0)
#par(mfrow=c(2,2))
for(i in 1:4)
{
  plot(R[i,], type="l", sub=bquote(sigma == .(sigma[i])),xlab="", ylab="R")
  abline(h=1.2, lty=2,col="red")
}

## -----------------------------------------------------------------------------
k<-c(4:25, 100, 500, 1000)
f.intersection_find <- function(k) 
{
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
res <- sapply(k, function (k) {f.intersection_find(k)})
intersection_find<-data.frame(k=k,intersection=res[1,],f.root=res[2,],iter=res[3,])
knitr::kable(intersection_find)


## -----------------------------------------------------------------------------
library(nloptr)
#Maximum likelihood estimation
#objective function
f.ob = function(x_mstep,x_estep,n_A,n_B,nOO,nAB)
{
  
  #E-step
  r_estep = 1-sum(x_estep)
  nAA = n_A*x_estep[1]^2/(x_estep[1]^2+2*x_estep[1]*r_estep)
  nBB = n_B*x_estep[2]^2/(x_estep[2]^2+2*x_estep[2]*r_estep)
  #M-step
  r_mstep = 1-sum(x_mstep)
  l1 <- 2*nAA*log(x_mstep[1])+2*nBB*log(x_mstep[2])+2*nOO*log(r_mstep)
  l2 <- (n_A-nAA)*log(2*x_mstep[1]*r_mstep)+(n_B-nBB)*log(2*x_mstep[2]*r_mstep)+nAB*log(2*x_mstep[1]*x_mstep[2])
  #log-likelihood function
  ll <- l1+l2
  return(-ll)
}
#Constraint condition
f.cc <- function(x_mstep,x_estep,n_A,n_B,nOO,nAB) 
{
  return(sum(x_mstep)-0.9999999)
}
para <- matrix(0,1,2)
#p0 and q0
para <- rbind(para,c(0.25,0.15))
EM <- NULL
i <- 2
while (mean(abs(para[i,]-para[i-1,]))>1e-8) 
{
  res <- nloptr( x0=c(0.32,0.12),eval_f = f.ob,
                 lb = c(0,0), ub = c(1,1), 
                 eval_g_ineq = f.cc, 
                 ##a fractional tolerance xtol_rel
                 opts = list("algorithm"="NLOPT_LN_COBYLA","xtol_rel"=1.0e-8),
                 x_estep=para[i,],n_A=444,n_B=132,nOO=361,nAB=63)
  para <- rbind(para,res$solution)
  EM <- c(EM,res$objective)
  i <- i+1
}
#EM algorithm estimate
print(para[-1,])
#Objective function value
plot(EM,type = 'l')

## -----------------------------------------------------------------------------
library(base)
formulas <- list(#four models
  mpg ~ disp,
  mpg ~ I(1 / disp),
  mpg ~ disp + wt,
  mpg ~ I(1 / disp) + wt
)
data("mtcars")
attach(mtcars)
#using loop
Models1 <-  vector("list", length(formulas))
for (i in 1:length(formulas)) 
  Models1[[i]] <- lm(formulas[[i]], data = mtcars)
print(Models1)
# using lapply
Models2 <- lapply(X = formulas, FUN = function(model){ y <- lm(model, data = mtcars)})
# print the models in list
print(Models2)

## -----------------------------------------------------------------------------
###12.3
set.seed(123)
trials <- replicate(
  #generate 100 p values
  100,
  #use 10+7 Poisson random values to generate one p value
  t.test(rpois(10, 10), rpois(7, 10)),
  simplify = FALSE
)
p.value1 <- sapply(1:length(trials), function(index) trials[[index]]$p.value)
p.value2 <- sapply(trials, function(X) X[['p.value']])
# [[ is a function return trials[[]]
p.value3 <- sapply(trials, '[[', 'p.value')
#Returns the first parts (10)
print(p.value1)
print(p.value2)
print(p.value3)

## -----------------------------------------------------------------------------
#a variant of lapply() 
vari.lapply<-function (f,n,data_type, ...)
{
#to match as function
f<-match.fun(FUN = f, descend = TRUE)
f1 <- Map(f, ...)
#logical < integer < double < complex < character 
if(data_type == "character")     res <- vapply(f1,cbind,FUN.VALUE = character(n))
else if(data_type == "complex")  res <- vapply(f1,cbind,FUN.VALUE = complex(n))
else if(data_type == "numeric")  res <- vapply(f1,cbind,FUN.VALUE = numeric(n))
else if(data_type == "logical")  res <- vapply(f1,cbind,FUN.VALUE = logical(n))
return(res)
}
#Example(self-defining)
#matrix
set.seed(123)
f.xy <- function(x,y) {x>y}
x1<-list(runif(5,-1, 1), runif(5,-2,2))
x2<-list(rnorm(5), rnorm(5,0,1.5))
matrix1 <-vari.lapply(f.xy,n=5,data_type="logical",x1,x2)
print(matrix1)
#vector
set.seed(123)
y <- list(rnorm(5), rnorm(10))
vector1 <- vari.lapply(mean,n=1,data_type="numeric",y)
print(vector1)

## ----eval=FALSE---------------------------------------------------------------
#  #the standard Laplace distribution*2
#  pdf.Laplace_R <- function(x){
#    return(exp(-1*abs(x)))
#  }
#  #
#  rw.Metro_R <-function(sigma,N,x_initial){
#   u <- runif(N)
#   x <- numeric(N)
#   x[1] <- x_initial
#   k <- 0
#   for (i in 2:N){
#      y <- rnorm(1,x[i-1],sigma)
#      if(u[i] <= (pdf.Laplace_R(y)/pdf.Laplace_R(x[i-1]))){
#       x[i] <- y
#      }
#      else{
#       x[i] <- x[i-1]
#       k <- k+1
#   }
#   return (list(x=x,k))
#  }

## ----eval=FALSE---------------------------------------------------------------
#  #include <Rcpp.h>
#  using namespace Rcpp;
#  // [[Rcpp::export]]
#  double pdfLaplaceC(double x)
#  {
#    return exp(-1*abs(x));
#   }
#  // [[Rcpp::export]]
#  List rwMetroC(double sigma, int N, double xinitial){
#    NumericVector u = runif(N);
#    NumericVector x(N);
#    x[0] = xinitial;
#    int k = 0;
#    NumericVector y;
#    for(int i = 1;i <= (N-1);i++){
#      y = rnorm(1, x[i-1], sigma);
#      if (u[i] <= (pdfLaplaceC(y[0])/pdfLaplaceC(x[i-1]))) {
#        x[i] = y[0];
#      }
#      else {
#        x[i] = x[i-1];
#        k++;
#      }
#    }
#    List res = List::create(x,k);
#    return (res);
#  }
#  

## -----------------------------------------------------------------------------
#the standard Laplace distribution*2 
pdf.Laplace_R <- function(x){
  return(exp(-1*abs(x)))
}
#
rw.Metro_R <-function(sigma,N,x_initial){
 u <- runif(N)
 x <- numeric(N) 
 x[1] <- x_initial
 k <- 0
 for (i in 2:N){
    y <- rnorm(1,x[i-1],sigma) 
    if(u[i] <= (pdf.Laplace_R(y)/pdf.Laplace_R(x[i-1]))){
     x[i] <- y
    }
    else{
      x[i] <- x[i-1]
      k <- k+1
    }
 }
 return (list(x=x,k))
}
library(Rcpp)
library(microbenchmark)
sourceCpp(
  code='
  #include<Rcpp.h>
  using namespace Rcpp;
  // [[Rcpp::export]]
  double pdfLaplaceC(double x) 
  {
    return exp(-1*abs(x));
  }
  // [[Rcpp::export]]
  List rwMetroC(double sigma, int N, double xinitial){
    NumericVector u = runif(N);
    NumericVector x(N); 
    x[0] = xinitial;
    int k = 0;
    NumericVector y;
    for(int i = 1;i <= (N-1);i++){
      y = rnorm(1, x[i-1], sigma);
      if (u[i] <= (pdfLaplaceC(y[0])/pdfLaplaceC(x[i-1]))) {
        x[i] = y[0];
      }
      else {
        x[i] = x[i-1];
  	    k++;
  	}
    }
    List res = List::create(x,k);
    return (res);
  } 
  '
)
# Campare the computation time when sigma = 2
runtime<-microbenchmark(list_R=rw.Metro_R(2,200,25),list_C=rwMetroC(2,200,25))
summary(runtime)

## -----------------------------------------------------------------------------
set.seed(12345)
n <- 4
sigma<-c(0.05,0.5,2,16)
N <- 200
x0 <- 25
X_R <- X_C<- matrix(nrow = N ,ncol = n)
rej.rateR <- rej.rateC <- numeric(n)
#par(mfrow=c(2,2))
for (i in 1:n){
  rw_R <- rw.Metro_R(sigma[i],N,x0)
  rw_C <- rwMetroC(sigma[i],N,x0)
  X_R[,i]<- rw_R[[1]]
  X_C[,i]<- rw_C[[1]]
  qqplot(rw_R[[1]],rw_R[[1]],xlab = "rw_R",ylab = "rw_C",main = paste("σ=",sigma[i]))
  abline(a = 0,b = 1,col = 2, lwd = 2)
  rej.rateR[i]<-rw_R[[2]]/N
  rej.rateC[i]<-rw_C[[2]]/N
}

## -----------------------------------------------------------------------------
library(knitr)
kable(data.frame(sigma = sigma,reject.rate_R = rej.rateR,reject.rate_C = rej.rateC))


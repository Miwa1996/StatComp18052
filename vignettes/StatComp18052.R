## ------------------------------------------------------------------------
a=rnorm(10)
a

## ------------------------------------------------------------------------
matrix(1:6, 2, 3)

## ------------------------------------------------------------------------
x=rnorm(10)
y=rnorm(10)
plot(x,y)

## ----Code chunk1, echo = FALSE-------------------------------------------
x <- 0:4; p <- c(0.1,0.2,0.2,0.2,0.3)
cp <- cumsum(p); m <- 1e3; r <- numeric(m)
r <- x[findInterval(runif(m),cp)+1]
ct <- as.vector(table(r)); ct/sum(ct) 

## ----Code chunk2, echo = FALSE-------------------------------------------
n <- 1e3;j<-k<-0;y <- numeric(n)
while (k < n) {
u <- runif(1)
j <- j + 1
x <- runif(1) #random variate from g
if (27/4*x^2 * (1-x) > u) {
#we accept x
k <- k + 1
y[k] <- x
}
}
hist(y,breaks=20,freq=F)
f<- function(x)  return(12*x^2*(1-x))
curve(f,add=T)

## ----Code chunk3, echo = FALSE-------------------------------------------
n <- 1e3; r <- 4; beta <- 3
lambda <- rgamma(n, shape=r, scale=beta)
y <- rexp(n, lambda)
hist(y,breaks=20,freq=F)

## ----Code chunk4, echo = FALSE-------------------------------------------
Pbeta<-function(y){
if(y<=0) value=0
if(y>1) value=1
if(y>=0&y<1){
m <- 1e4; x <- runif(m, min=0, max=y)
value <- mean(30*x^2*(1-x)^2)*(y-0)
}
value
}
sapply((1:9)/10,Pbeta)
pbeta1<-function(x){
  pbeta(x,3,3)
}
sapply((1:9)/10,pbeta1)

## ----Code chunk5_1,echo = FALSE------------------------------------------
MC.Phi<-function(x,sigma,R= 10000, antithetic = TRUE) {
u <- runif(R/2)
if (!antithetic)
v<- runif(R/2)
else 
v<-1-u
w<-c(u, v)
cdf <- numeric(length(x))
for (i in 1:length(x)) {
g <-x[i]*w/(sigma)^2*exp((-(x[i]*w)^2/(2*sigma*sigma)))*x[i]
cdf[i] <- mean(g)
}
return(cdf)
}
x <- seq(.1,2.5,length=5);
pRayleigh<-function(x,sigma){
  s<-sigma;
  p<-numeric(length(x));
  intergrand<-function(x){
    x/(s^2)*exp((-x^2/(2*(s^2))))};
  for(i in 1:length(x)){
  p[i]<-integrate(intergrand,0,x[i])$value;
  }
  return(p)
}
Phi<-pRayleigh(x,sigma=2)
set.seed(123)
MC1<- MC.Phi(x,sigma=2,anti=FALSE) #for (X1+X2)/2 which X1,X2 is independent
set.seed(123)
MC2<- MC.Phi(x,sigma=2,anti=TRUE)  #for antithetic variables (X+X')/2
print(round(rbind(x, MC1, MC2, Phi),5))

## ----Code chunk5_2,echo = FALSE------------------------------------------
m <- 1000
MC1 <- MC2 <- numeric(m)
x <- 1.95
for (i in 1:m) {
MC1[i] <- MC.Phi(x,2,R = 1000, anti = FALSE)
MC2[i] <- MC.Phi(x,2,R = 1000,anti=TRUE)
}
print(sd(MC1))

## ----Code chunk5_3,echo = FALSE------------------------------------------
print(sd(MC2))

## ----Code chunk5_4,echo = FALSE------------------------------------------
print((var(MC1) - var(MC2))/var(MC1))

## ----Code chunk6, echo = FALSE-------------------------------------------
x <- seq(1,5,length.out = 20)
    g <- exp(-x^2/2)*x^2/sqrt(2*pi)
    f1 <- exp(-x+1)
    f2 <- 1/2*(x-1)^2 *exp(-x+1)
#figure (a)
plot(g~x,type = "l",col=1)
lines(f1~x,col=2)
lines(f2~x,col=3)
legend("topright", legend =c("g", "f1", "f2"),
           lty = 1:3, lwd = 2, inset = 0.02,col=1:3)
#figure (b)
plot(g/f1~x,type = "l",col=1)
lines(g/f2~x,col=2)
legend("topright", legend =c("f1", "f2"),
           lty = 1:2, lwd = 2, inset = 0.02,col=1:2)
m <- 10000
  theta.hat <- se <- numeric(2)
  g <- function(x) {
  exp(-x^2/2)*x^2/sqrt(2*pi) * (x > 1)
  }
x <- rexp(m, rate= 1)+1 #using f1
  fg <- g(x)/exp(-x+1)
  theta.hat[1] <- mean(fg)
  se[1] <- sd(fg)
x <- rgamma(m, shape=3, rate = 1)+1 #using f2
  fg <- g(x)/(1/2*(x-1)^2 *exp(-x+1))
  theta.hat[2] <- mean(fg)
  se[2] <- sd(fg)
  res <- rbind(theta=round(theta.hat,3), se=round(se,3))
  colnames(res) <- paste0('f',1:2)
  knitr::kable(res,align='c')

## ----Code chunk7, echo = FALSE-------------------------------------------
m <- 10000
  theta.hat <- se <- numeric(1)
  g <- function(x) {
  exp(-x^2/2)*x^2/sqrt(2*pi) * (x > 1)
  }
x <- rexp(m, rate= 1)+1 #using f
  fg <- g(x)/exp(-x+1)
  theta.hat[1] <- mean(fg)
  se[1] <- sd(fg)
  res <- rbind(theta=round(theta.hat,3), se=round(se,3))
  colnames(res) <- 'f'
  knitr::kable(res,align='c')

## ----Code chunk8, echo = FALSE-------------------------------------------
m <- 1000 # Number of Monte Carlo trials
n <- 100
set.seed(1)
G.hat1 <- numeric(m) # Storage for test statistics from the MC trials
g1 <- numeric(n)
Id <- 1:n
## Start the simulation
for (i in 1:m){
  x1 <- rlnorm(n) # x1 generated from standard lognormal distribution
  mu_hat1 <- mean(x1) # The estimation of mu
  x1_order <- sort(x1) # x1_order is the order statistic of x1
  G.hat1[i] <- (2 * Id -n-1)%*%x1_order/(n^2 *mu_hat1)# Estimate the value of G
}
print(c(mean(G.hat1),median(G.hat1))) # the mean and median of G.hat1
print(quantile(G.hat1,probs=seq(0.1,1,0.1))) # the deciles of G.hat1
hist(G.hat1,prob = TRUE)

m <- 1000 # Number of Monte Carlo trials
n <- 100
set.seed(12)
G.hat2 <- numeric(m) # Storage for test statistics from the MC trials
g2 <- numeric(n)
Id <- 1:n
## Start the simulation
for (i in 1:m){
  x2 <- runif(n) # x2 generated from uniform distribution
  mu_hat2 <- mean(x2) # The estimation of mu
  x2_order <- sort(x2) # x2_order is the order statistic of x2
  G.hat2[i] <- (2 * Id -n-1)%*%x2_order/(n^2 *mu_hat2)# Estimate the value of G
}
print(c(mean(G.hat2),median(G.hat2))) # the mean and median of G.hat2
print(quantile(G.hat2,probs=seq(0.1,1,0.1))) # the deciles of G.hat2
hist(G.hat2,prob = TRUE)

m <- 1000 # Number of Monte Carlo trials
n <- 100
set.seed(123)
G.hat3 <- numeric(m) # Storage for test statistics from the MC trials
g3 <- numeric(n)
Id <- 1:n
## Start the simulation
for (i in 1:m){
  x3 <- rbinom(n,1,0.1) # x3 generated from Bernoulli(0.1) distribution
  mu_hat3 <- mean(x3) # The estimation of mu
  x3_order <- sort(x3) # x3_order is the order statistic of x3
  G.hat3[i] <- (2 * Id -n-1)%*%x3_order/(n^2 *mu_hat3)# Estimate the value of G
}
print(c(mean(G.hat3),median(G.hat3))) # the mean and median of G.hat3
print(quantile(G.hat3,probs=seq(0.1,1,0.1))) # the deciles of G.hat3
hist(G.hat3,prob = TRUE)


## ----Code chunk9_1,echo = FALSE------------------------------------------
# function to calculate the confidence interval
compu.interval <- function(a,b){
m<-1e3
G<-numeric(m)
I<-2*c(1:m)-m-1
set.seed(123)
for(i in 1:m){
  x<-rlnorm(m,a,b) #generate random numbers
  x<-sort(x) #sorting x
  mu=mean(x)
  G[i]<-1/m^2/mu*(t(I)%*%x) #compute G
}
CI<-c(mean(G)-1.96*sd(G)/sqrt(m),mean(G)+1.96*sd(G)/sqrt(m))#compute confidence interval
return(CI)
}

## ----Code chunk9_2,echo = FALSE------------------------------------------
#approximate Coverage probability(ECP) of confidence interval
N<-100
bar<-numeric(N)
k<-0
a <- 0
b <- 1
I<-2*c(1:m)-m-1
G.true<-numeric(m)
CI <- compu.interval(a,b)
set.seed(1234)
for(j in 1:N){
  for(i in 1:m){
    x<-rlnorm(m,0,1) 
    x<-sort(x) 
    mu<-mean(x)
    G.true[i]<-1/m^2/mu*(t(I)%*%x)
  }
  bar[j]<-mean(G.true)
  if(bar[j]>CI[1]&bar[j]<CI[2]){
    k<-k+1}
}
k/N

## ----Code chunk11, echo = FALSE------------------------------------------
set.seed(12345)
LSAT=c(576,635,558,578,666,580,555,661,651,605,653,575,545,572,594)
GPA=c(339,330,281,303,344,307,300,343,336,313,312,274,276,288,296)
x=cbind(LSAT,GPA)
n=15
b.cor <- function(x,i) cor(x[i,1],x[i,2])
theta.hat <- b.cor(x,1:n)
theta.jack <- numeric(n)
for(i in 1:n){
theta.jack[i] <- b.cor(x,(1:n)[-i])
}
bias.jack <- (n-1)*(mean(theta.jack)-theta.hat)
se.jack <- sqrt((n-1)*mean((theta.jack-theta.hat)^2))
bias.jack
se.jack

## ----Code chunk12_1,echo = FALSE-----------------------------------------
set.seed(12345)
library(boot)
x=c(3,5,7,18,43,85,91,98,100,130,230,487)
boot.mean=function(x,i) mean(x[i])
de=boot(data=x,statistic=boot.mean,R=1024)
ci=boot.ci(de,type=c("norm","basic","perc","bca"))
ci


## ----Code chunk13, echo = FALSE------------------------------------------
library(bootstrap)
data=scor
u=c(mean(data[,1]),mean(data[,2]),mean(data[,3]),mean(data[,4]),mean(data[,5]))
m=matrix(0,5,5)
for (i in 1:88) m=m+(as.numeric(data[i,])-u)%*%t(as.numeric(data[i,])-u)  
m=m/88   ##MLE of covariance matrix
lambda=eigen(m)$values  ## the eigenvalues
theta.hat=lambda[1]/sum(lambda)
theta.jack=numeric(5)
for (i in 1:5) theta.jack[i]=lambda[1]/sum(lambda[-i])
bias.jack <- (n-1)*(mean(theta.jack)-theta.hat)
se.jack <- sqrt((n-1)*mean((theta.jack-theta.hat)^2))
bias.jack
se.jack

## ----Code chunk14_1,echo = FALSE-----------------------------------------
##leave-one-out (n-fold) cross validation
library(DAAG); attach(ironslag)
n <- length(magnetic) #in DAAG ironslag
e1 <- e2 <- e3 <- e4 <- numeric(n)
# for n-fold cross validation
# fit models on leave-one-out samples
for (k in 1:n) {
y <- magnetic[-k]
x <- chemical[-k]
J1 <- lm(y ~ x)
yhat1 <- J1$coef[1] + J1$coef[2] * chemical[k]
e1[k] <- magnetic[k] - yhat1
J2 <- lm(y ~ x + I(x^2))
yhat2 <- J2$coef[1] + J2$coef[2] * chemical[k] +
J2$coef[3] * chemical[k]^2
e2[k] <- magnetic[k] - yhat2
J3 <- lm(log(y) ~ x)
logyhat3 <- J3$coef[1] + J3$coef[2] * chemical[k]
yhat3 <- exp(logyhat3)
e3[k] <- magnetic[k] - yhat3
J4 <- lm(log(y) ~ log(x))
logyhat4 <- J4$coef[1] + J4$coef[2] * log(chemical[k])
yhat4 <- exp(logyhat4)
e4[k] <- magnetic[k] - yhat4
}
c(mean(e1^2), mean(e2^2), mean(e3^2), mean(e4^2))
##According to the prediction error criterion, Model 2, the quadratic model,
##would be the best fit for the data.
##leave-two-out cross validation
n <- length(magnetic) #in DAAG ironslag
e1 <- e2 <- e3 <- e4 <- numeric(n)
E1 <- E2 <- E3 <- E4 <-rep(0,2)
subscript<-t(combn(n,2))
for (k in 1:choose(n,2)) {
K<-subscript[k,]
y <- magnetic[-K]
x <- chemical[-K]
J1 <- lm(y ~ x)
yhat1 <- J1$coef[1] + J1$coef[2] * chemical[K]
E1<- magnetic[K] - yhat1
e1[k]<-sum(abs(E1))
J2 <- lm(y ~ x + I(x^2))
yhat2 <- J2$coef[1] + J2$coef[2] * chemical[K] +J2$coef[3] * chemical[K]^2
E2 <- magnetic[K] - yhat2
e2[k]<-sum(abs(E2))
J3 <- lm(log(y) ~ x)
logyhat3 <- J3$coef[1] + J3$coef[2] * chemical[K]
yhat3 <- exp(logyhat3)
E3<- magnetic[K] - yhat3
e3[k]<-sum(abs(E3))
J4 <- lm(log(y) ~ log(x))
logyhat4 <- J4$coef[1] + J4$coef[2] * log(chemical[K])
yhat4 <- exp(logyhat4)
E4<- magnetic[K] - yhat4
e4[k]<-sum(abs(E4))
}
c(mean(e1^2), mean(e2^2), mean(e3^2), mean(e4^2))
##According to the prediction error criterion, Model 2, the quadratic model,
##would be the best fit for the data.

## ----setup, include=FALSE------------------------------------------------
library(latticeExtra)
library(RANN)
library(energy)
library(Ball)
library(boot)
library(ggplot2)

## ----Code chunk15, echo = FALSE------------------------------------------
##function:two-sample Cramer-von Mises test for equal distributions
cvm <- function(x,y,data){
  r <- 1000 #permutation samples
  reps <- numeric(r)
  n <- length(x)
  m <- length(y)
  v.n <- numeric(n)
  v1.n <- numeric(n)
  v.m <- numeric(m)
  v1.m <- numeric(m)
  z <- c(x,y)
  N <- length(z)
  Ix <- seq(1:n)
  Iy <- seq(1:m)
  v.n <- (x-Ix)**2
  v.m <- (y-Iy)**2
  #test statistic
  reps_0 <- ((n * sum(v.n)+m * sum(v.m))/(m * n * N))-(4 * m * n - 1)/(6 * N)
  for (k in 1:r){#permutation samples
    w <- sample(N,size=n,replace=FALSE)
    x1 <- sort(z[w])
    y1 <- sort(z[-w])
    v1.n <- (x1-Ix)**2
    v1.m <- (y1-Iy)**2
    reps[k] <- ((n * sum(v1.n)+m * sum(v1.m))/(m * n * N))-(4 * m * n - 1)/(6 * N)
  }
  p <- mean(c(reps_0,reps) >= reps_0)
  return(
    histogram(c(reps_0,reps),
              type="density",
              col="#0080ff",
              xlab="Replicates of Cramer-Von Mises statistic",
              ylab=list(rot=0),
              main=paste0("Data:",data),
              sub=list(substitute(paste(hat(p),"=",pvalue),list(pvalue=p)),col=2),
              panel=function(...){
                panel.histogram(...)
                panel.abline(v=reps_0,col=2,lwd=2)
              })
  )
}

##Data: Example 8.1
attach(chickwts)
x <- sort(as.vector(weight[feed == "soybean"]))
y <- sort(as.vector(weight[feed == "linseed"]))

cvm1 <- cvm(x,y,"Example 8.1")

##Data: Example 8.2
x <- sort(as.vector(weight[feed == "sunflower"]))
y <- sort(as.vector(weight[feed == "linseed"]))
detach(chickwts)

cvm2 <- cvm(x,y,"Example 8.2")

##Results
print(cvm1)
print(cvm2)

## ------------------------------------------------------------------------
## variable definition
m <- 500 #permutation samples
p<-2 # dimension of data
n1 <- n2 <- 50 #the sample size of x and y
R<-999 #boot parameter
k<-3 #boot parameter
n <- n1 + n2
N = c(n1,n2)
# the function of NN method
Tn <- function(z, ix, sizes,k){
  n1 <- sizes[1]; n2 <- sizes[2]; n <- n1 + n2
  if(is.vector(z)) z <- data.frame(z,0);
  z <- z[ix, ];
  NN <- nn2(data=z, k=k+1) 
  block1 <- NN$nn.idx[1:n1,-1]
  block2 <- NN$nn.idx[(n1+1):n,-1]
  i1 <- sum(block1 < n1 + .5)
  i2 <- sum(block2 > n1+.5)
  (i1 + i2) / (k * n)
}

eqdist.nn <- function(z,sizes,k){
  boot.obj <- boot(data=z,statistic=Tn,R=R,sim = "permutation", 
                 sizes = sizes,k=k)
  ts <- c(boot.obj$t0,boot.obj$t)
  p.value <- mean(ts>=ts[1])
  list(statistic=ts[1],p.value=p.value)
}
p.values <- matrix(NA,m,3) #p<U+05B5>

## ------------------------------------------------------------------------
##(1)Unequal variances and equal expectations
set.seed(1)
sd <- 1.5
for(i in 1:m){
  x <- matrix(rnorm(n1*p),ncol=p)
  y <- matrix(rnorm(n2*p,sd=sd),ncol=p)
  z <- rbind(x,y)
  p.values[i,1] <- eqdist.nn(z,N,k)$p.value#NN method
  p.values[i,2] <- eqdist.etest(z,sizes=N,R=R)$p.value#energy methods
  p.values[i,3] <- bd.test(x=x,y=y,R=999,seed=i*12345)$p.value# ball method
}
alpha <- 0.05;
pow <- colMeans(p.values<alpha)
power <- data.frame(methods = c('NN','energy','Ball'),pow)
ggplot(power,aes(methods,pow))+#plot
  geom_col(fill = 'palegreen3')+
  coord_flip()

## ------------------------------------------------------------------------
##(2)Unequal variances and unequal expectations
set.seed(1)
mu <- 0.5
sd <- 1.5
for(i in 1:m){
  x <- matrix(rnorm(n1*p),ncol=p)
  y <- matrix(rnorm(n2*p,mean=mu,sd=sd),ncol=p)
  z <- rbind(x,y)
  p.values[i,1] <- eqdist.nn(z,N,k)$p.value#NN method
  p.values[i,2] <- eqdist.etest(z,sizes=N,R=R)$p.value#energy methods
  p.values[i,3] <- bd.test(x=x,y=y,R=999,seed=i*12345)$p.value# ball method
}
alpha <- 0.05;
pow <- colMeans(p.values<alpha)
pow
power <- data.frame(methods = c('NN','energy','Ball'),pow)
ggplot(power,aes(methods,pow))+#plot
  geom_col(fill = 'palegreen3')+
  coord_flip()

## ------------------------------------------------------------------------
##Non-normal distributions: t distribution with 1 df (heavy-tailed
##distribution), bimodal distribution (mixture of two normal
##distributions)
set.seed(1)
mu <- 0.5
sd <- 2
for(i in 1:m){
  x <- matrix(rt(n1*p,df=1),ncol=p)
  y1 = rnorm(n2*p);  y2 = rnorm(n2*p,mean=mu,sd=sd)
  w = rbinom(n, 1, .5) # 50:50 random choice
  y <- matrix(w*y1 + (1-w)*y2,ncol=p)# normal mixture
  z <- rbind(x,y)
  p.values[i,1] <- eqdist.nn(z,N,k)$p.value#NN method
  p.values[i,2] <- eqdist.etest(z,sizes=N,R=R)$p.value#energy methods
  p.values[i,3] <- bd.test(x=x,y=y,R=999,seed=i*12345)$p.value# ball method
}
alpha <- 0.05;
pow <- colMeans(p.values<alpha)
pow
power <- data.frame(methods = c('NN','energy','Ball'),pow)
ggplot(power,aes(methods,pow))+#plot
  geom_col(fill = 'palegreen3')+
  coord_flip()

## ------------------------------------------------------------------------
##Unbalanced samples 
set.seed(1)
mu <- 0.5
N = c(n1,n2*2)
for(i in 1:m){
  x <- matrix(rnorm(n1*p),ncol=p);
  y <- cbind(rnorm(n2*2),rnorm(n2*2,mean=mu));
  z <- rbind(x,y)
  p.values[i,1] <- eqdist.nn(z,N,k)$p.value#NN method
  p.values[i,2] <- eqdist.etest(z,sizes=N,R=R)$p.value#energy methods
  p.values[i,3] <- bd.test(x=x,y=y,R=999,seed=i*12345)$p.value# ball method
}
alpha <- 0.05;
pow <- colMeans(p.values<alpha)
pow
power <- data.frame(methods = c('NN','energy','Ball'),pow)
ggplot(power,aes(methods,pow))+#plot
  geom_col(fill = 'palegreen3')+
  coord_flip()

## ----Code chunk16, echo = FALSE------------------------------------------
set.seed(1)
n <- 10000 #Sample size
x <- numeric(n)
u <- runif(n)
theta=1
eta=0

x[1] <- rnorm(1)
k <- 0

# cauchy functions
f <- function(x, theta=1, eta=0){
  out <- 1/(pi * theta * (1+((x-eta)/theta)^2))
  return(out)
}

for(i in 2:n){
  xt <- x[i-1]
  y <- rnorm(1,mean=xt)
  R <- f(y)*dnorm(xt,mean=y)/(f(xt)*dnorm(y,mean=xt))
  if(u[i] <= R){
    x[i] <- y
  }else{
    x[i] <- xt
    k <- k+1
  }
}

is <- 1001:n
par(mfrow=c(1,2))
plot(is,x[is],type="l")
hist(x[is], probability=TRUE,breaks=100)
plot.x <- seq(min(x[is]),max(x[is]),0.01)
lines(plot.x,f(plot.x))
par(mfrow=c(1,1))
#compare the deciles
observations <- quantile(x[is],seq(0,1,0.1))
expectations <- qcauchy(seq(0,1,0.1))
decile <- data.frame(observations,expectations)
decile

## ----Code chunk19_1,echo = FALSE-----------------------------------------
Sk_1 <- function(a,k){
  q <- sqrt(a^2*(k-1)/(k-a^2)) 
  return (1-pt(q,df=k-1))
}
Sk <- function(a,k){
  q <- sqrt(a^2*k/(k+1-a^2)) 
  return (1-pt(q,df=k))
}
difSK <- function(x,k) { 
  Sk_1(x,k)-Sk(x,k)
}
kset <- c(4:25,100,500,1000)
out <- 1:length(kset)
for (i in 1:length(kset)){
  out[i] <- uniroot( difSK
            , lower = 0+1e-5, upper = sqrt(kset[i])-1e-5,k=kset[i]) $root
}
out


kset[ abs(out-sqrt(kset)) < sqrt(kset)*0.01]
 
n <- 1:length(kset)
Kwrongnum <- n[abs(out-sqrt(kset)) < sqrt(kset)*0.01]

 
#Example : k=23
k=23
xx <- seq(0.01,sqrt(k)-1e-5,length=1000)
y <- difSK(xx,k)
plot(xx,y,type="l")

#Example : k=1000
k=1000
xx <- seq(0.01,sqrt(k)-1e-5,length=1000)
y <- difSK(xx,k)
plot(xx,y,type="l")

#change upper to 3

for (i in Kwrongnum){
  out[i] <- uniroot( difSK
                     , lower = 0+1e-5, upper =3,k=kset[i]) $root
}
names(out) <- kset

out


## ----Code chunk20, echo = FALSE------------------------------------------
f<-function(y,theta,eta){

1/(theta*3.141592653*(1+((y-eta)/theta)^2))
}


pdf<-function(x,theta,eta,lower.tail=TRUE){
 if(lower.tail) res<-integrate(f,lower = -Inf,upper = x,rel.tol=.Machine$double.eps^0.25,theta=theta,eta=eta)
 else res<-integrate(f,lower = x,upper = Inf,rel.tol=.Machine$double.eps^0.25,theta=theta,eta=eta)
  return(res$value)
}
pdf(x=0,theta = 1,eta = 0)
pcauchy(0,location = 0,scale = 1)

pdf(x=2,theta = 2,eta =1,lower.tail = F )
pcauchy(2,location = 1,scale = 2,lower.tail = F)

## ----echo=FALSE----------------------------------------------------------
        dat <- rbind(Genotype=c('AA','BB','OO','AO','BO','AB'),
                     Frequency=c('p2','q2','r2','2pr','2qr','2pq',1),
                     Count=c('nAA','nBB','nOO','nAO','nBO','nAB','n'))
    knitr::kable(dat,format='markdown',caption = "Comparation of them",align = "c")

## ----Code chunk21_1,echo = FALSE-----------------------------------------
library(nloptr)
# Mle 

eval_f0 <- function(x,x1,n.A=28,n.B=24,nOO=41,nAB=70) {
  
  r1<-1-sum(x1)
  nAA<-n.A*x1[1]^2/(x1[1]^2+2*x1[1]*r1)
  nBB<-n.B*x1[2]^2/(x1[2]^2+2*x1[2]*r1)
  r<-1-sum(x)
  return(-2*nAA*log(x[1])-2*nBB*log(x[2])-2*nOO*log(r)-
           (n.A-nAA)*log(2*x[1]*r)-(n.B-nBB)*log(2*x[2]*r)-nAB*log(2*x[1]*x[2]))
}


# constraint
eval_g0 <- function(x,x1,n.A=28,n.B=24,nOO=41,nAB=70) {
  return(sum(x)-0.999999)
}

opts <- list("algorithm"="NLOPT_LN_COBYLA",
             "xtol_rel"=1.0e-8)
mle<-NULL
r<-matrix(0,1,2)
r<-rbind(r,c(0.2,0.35))# the beginning value of p0 and q0
j<-2
while (sum(abs(r[j,]-r[j-1,]))>1e-8) {
res <- nloptr( x0=c(0.3,0.25),
               eval_f=eval_f0,
               lb = c(0,0), ub = c(1,1), 
               eval_g_ineq = eval_g0, 
               opts = opts, x1=r[j,],n.A=28,n.B=24,nOO=41,nAB=70 )
j<-j+1
r<-rbind(r,res$solution)
mle<-c(mle,eval_f0(x=r[j,],x1=r[j-1,]))
}
r  #the result of EM algorithm
mle #the max likelihood values


## ----Code chunk23, echo = FALSE------------------------------------------
attach(mtcars)
formulas <- list(
  mpg ~ disp,
  mpg ~ I(1 / disp),
  mpg ~ disp + wt,
  mpg ~ I(1 / disp) + wt
)
#1 for loops
f3<- vector("list", length(formulas))
for (i in seq_along(formulas)){
  f3[[i]] <- lm(formulas[[i]], data = mtcars)
}
f3
#2 lapply
la3<-lapply(formulas, function(x) lm(formula = x, data = mtcars))
la3

## ----Code chunk24,echo = FALSE-------------------------------------------

set.seed(123)
bootstraps <- lapply(1:10, function(i) {
  rows <- sample(1:nrow(mtcars), rep = TRUE)
  mtcars[rows, ]
})
# for loops
f4<- vector("list", length(bootstraps))
for (i in seq_along(bootstraps)){
  f4[[i]] <- lm(mpg ~ disp, data = bootstraps[[i]])
}
f4
# lapply without anonymous function
la4<- lapply(bootstraps, lm, formula = mpg ~ disp)
la4

## ----Code chunk25,echo = FALSE-------------------------------------------
rsq <- function(mod) summary(mod)$r.squared
#3
sapply(la3, rsq)
sapply(f3, rsq)
#4
sapply(la4,rsq)
sapply(f4,rsq)

## ----Code chunk26,echo = FALSE-------------------------------------------
set.seed(123)
trials <- replicate(
100,
t.test(rpois(10, 10), rpois(7, 10)),
simplify = FALSE
)
# anonymous function:
sapply(trials, function(x) x[["p.value"]])

## ----Code chunk27,echo = FALSE-------------------------------------------
#example
options(warn = -1)
testlist <- list(iris, mtcars, cars)
lapply(testlist, function(x) vapply(x, mean, numeric(1)))
#a more specialized function:
lmapply <- function(X, FUN, FUN.VALUE, simplify = FALSE){
  out <- Map(function(x) vapply(x, FUN, FUN.VALUE), X)
  if(simplify == TRUE){return(simplify2array(out))}
  out
}
lmapply(testlist, mean, numeric(1))

## ------------------------------------------------------------------------
chisq.test2 <- function(x, y){
  
  # Input
  if (!is.numeric(x)) {
    stop("x must be numeric")}
  if (!is.numeric(y)) {
    stop("y must be numeric")}
  if (length(x) != length(y)) {
    stop("x and y must have the same length")}
  if (length(x) <= 1) {
    stop("length of x must be greater one")}
  if (any(c(x, y) < 0)) {
    stop("all entries of x and y must be greater or equal zero")}
  if (sum(complete.cases(x, y)) != length(x)) {
    stop("there must be no missing values in x and y")}
  if (any(is.null(c(x, y)))) {
    stop("entries of x and y must not be NULL")}
  
  # compute the theoretical value
  m <- rbind(x, y)#the actual value
  margin1 <- rowSums(m)
  margin2 <- colSums(m)
  n <- sum(m)
  me <- tcrossprod(margin1, margin2) / n #the theoretical value
  
  # Output
  STATISTIC = sum((m - me)^2 / me)
  dof <- (length(margin1) - 1) * (length(margin2) - 1)#degree of freedom
  p <- pchisq(STATISTIC, df = dof, lower.tail = FALSE)
  return(list(X_squared = STATISTIC, df = dof, `p-value` = p))
}


## ------------------------------------------------------------------------
a <- 11:15
b <- c(11,12.5,13.5,14.5,15.5)
m_test <- cbind(a,b)

identical(chisq.test(m_test),chisq.test2(a, b))

## ------------------------------------------------------------------------
chisq.test(m_test)

chisq.test2(a, b)


## ------------------------------------------------------------------------
chisq.test2c <- compiler::cmpfun(chisq.test2)

microbenchmark::microbenchmark(
  chisq.test(m_test),
  chisq.test2(a,b),
  chisq.test2c(a,b)
)

## ------------------------------------------------------------------------
table2 <- function(x,y){
  
  x_val <- unique(x)
  y_val <- unique(y)
  mat <- matrix(0L, length(x_val), length(y_val))
  for (i in seq_along(x)) {
    mat[which(x_val == x[[i]]), which(y_val == y[[i]])] <-
      mat[which(x_val == x[[i]]),  which(y_val == y[[i]])] + 1L
  }
  dimnames <- list(x_val, y_val)
  names(dimnames) <- as.character(as.list(match.call())[-1])  # R has names for dimnames... :/
  tab <- array(mat, dim = dim(mat), dimnames = dimnames)
  class(tab) <- "table"
  tab
}

## ------------------------------------------------------------------------
x <- c(1, 2, 3, 1, 2, 3)
y <- c(2, 3, 4, 2, 3, 4)
identical(table(x,y), table2(x,y))

## ------------------------------------------------------------------------
microbenchmark::microbenchmark(table(x,y), table2(x,y))


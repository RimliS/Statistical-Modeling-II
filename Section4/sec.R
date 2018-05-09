library(MASS)

x <- seq(0,100,len=200)

k <- matrix(0,200,200)
alpha <- 0.1
l <- 1

for (i in 1:200){
  for (j in 1:200) {
    k[i,j] <- (alpha^2)*exp(-(1/l^2)*(abs(x[i]-x[j]))^2)
  }
}

f <- rep(0,200)
f <- mvrnorm(1,rep(0,200),as.matrix(k))
plot(x,f,ty="l",main="length scale=1 and alpha=1")

nsample <- 5
fx <- matrix(0,200,nsample)
fnew <- rep(0,200)
count=0
for (i in 1:nsample) {
  fnew <- mvrnorm(1,rep(0,200),as.matrix(k))
  fx[,i] <- fnew
  count=count+1
}


plot(x,fx[,1],ty="l",ylim=c(-4,4),ylab="Functions")
lines(x,fx[,2],ty="l",col="blue")
lines(x,fx[,3],ty="l",col="green")
lines(x,fx[,4],ty="l",col="red")
lines(x,fx[,5],ty="l",col="purple")
legend('bottomleft',c("f1","f2","f3","f4","f5"),
       col=c("black","blue","green","red","purple"),lty=1:2,lwd=1,cex=0.5)
mtext("length scale = 0.1", side=3, line=2, cex.lab=1,las=1, col="blue")


#Ex-4.5
faith <- read.table("/Users/rimlisengupta/Dropbox/UT Austin/Spring 2018/SDS383D-StatMod2/Section4/faithful.txt", fill = TRUE, 
                    header=TRUE)

library(MASS)

x <- faith$waiting

k <- matrix(0,272,272)
alpha <- 0.1
l <- 1

for (i in 1:272){
  for (j in 1:272) {
    k[i,j] <- (alpha^2)*exp(-(1/l^2)*(abs(x[i]-x[j]))^2)
  }
}

y <- faith$eruptions
sigma <- 1
x.star <- 80

k.star <- rep(0,272)

for (i in 1:272) {
  k.star[i] <- (alpha^2)*exp(-(1/l^2)*(abs(x.star-x[i]))^2)
}

kstar <- (alpha^2)*exp(-(1/l^2)*(abs(x.star-x.star))^2)

mu.star <- t(k.star)%*%(solve(k+sigma^2*diag(272)))%*%y
Sig.star <- kstar+sigma^2-t(k.star)%*%(solve(k+sigma^2*diag(272)))%*%k.star

y.star <- rnorm(272,mu.star,sqrt(Sig.star))

plot(x,y.star,pch=16)
par(new=TRUE)
plot(x,y,col="red",pch=16)

n<-272
ynew <- matrix(0,1000,n)
count=0
for(m in 1:1000) {
  ynew[m,] <- rnorm(n,mu.star,sqrt(Sig.star))
  count <- count+1
}

post.mean <- apply(ynew,1,mean)
std <- apply(ynew,1,sd)
int.l <- post.mean-1.96*std
int.u <- post.mean+1.96*std

tab <- cbind(post.mean, int.l, int.u)
colnames(tab) <- c("posterior mean","lower inetrval","upper interval")
tab <- data.frame(tab)

library("reshape2")
library("ggplot2")

ggplot(data=tab,
       aes(x=int.l, y=int.u, colour=post.mean)) +
       geom_point() +
       geom_smooth()
  
  

#4.6
smp_size <- floor(0.5 * nrow(faith))

## set the seed to make your partition reproductible
set.seed(123)
train_ind <- sample(seq_len(nrow(faith)), size = smp_size)

train <- faith[train_ind, ]
test <- faith[-train_ind, ]

x <- train$waiting
y <- train$eruptions

library(GPfit)
library(lhs)
n <- 136
GPmodel <- GP_fit(x,y)
print(GPmodel)

Number Of Observations: n = 272
Input Dimensions: d = 1

Correlation: Exponential (power = 1.95)
Correlation Parameters: 
  beta_hat
[1] -7.491063

sigma^2_hat: [1] 520660.5

delta_lb(beta_hat): [1] 2.803141e-07

nugget threshold parameter: 20

xnew <- test$waiting
pred <- predict(GPmodel,xnew)
yhat <- pred$Y_hat
ytest <- test$eruptions

plot(xnew,ytest,col="red",ylab="fitted values",xlab="x")
points(xnew,yhat,col="blue",pch=16)

newdata <- faith[1:10,]

x <- newdata$waiting
y <- newdata$eruptions

GPmodel <- GP_fit(x,y)
print(GPmodel)

Number Of Observations: n = 272
Input Dimensions: d = 1

Correlation: Exponential (power = 1.95)
Correlation Parameters: 
  beta_hat
[1] -7.491063

sigma^2_hat: [1] 520660.5

delta_lb(beta_hat): [1] 2.803141e-07

nugget threshold parameter: 20

xnew <- test$waiting
pred <- predict(GPmodel,xnew)
yhat <- pred$Y_hat
ytest <- test$eruptions

plot(xnew,ytest,col="red",ylab="fitted values",xlab="x")
points(xnew,yhat,col="blue",pch=16)



#Ex:4.12
iris <- read.csv("/Users/rimlisengupta/Dropbox/UT Austin/Spring 2018/SDS383D-StatMod2/Section4/iris.csv",
                  sep = ",", header=TRUE)
iris <- iris[,-1]

ind1 <- which(iris$Species=="setosa")
ind2 <- which(iris$Species=="virginica")

temp1 <- iris[ind1,]
temp2 <- iris[ind2,]
temp <- rbind(temp1,temp2)
temp$Sepal.Length <- scale(temp$Sepal.Length)
temp$Sepal.Width <- scale(temp$Sepal.Width)
temp$Petal.Length <- scale(temp$Petal.Length)
temp$Petal.Width <- scale(temp$Petal.Width)


X <- as.matrix(cbind(rep(1, nrow(temp)),temp[,-5]))
colnames(X) <- c("intercept","Sepal.Length","Sepal.Width","Petal.Length","Petal.Width")
y <- as.numeric(temp$Species == "virginica")

#Define sigmoid
g <- function (z) {
  return (1 / (1 + exp(-z) ))
} 

#Define hypothesis 
h <- function (x,beta) {
  return(g(x %*% beta))
}

#Define log-posterior
log.p <- function (x,y,beta,lambda) {
  return(t(y)%*%log(h(x,beta))+t((1-y))%*%log(1-h(x,beta))-(lambda/2)*t(beta)%*%beta)
}

#Gradient descent
grad <- function (x,y,beta,lambda) {
  return(t(x)%*%(y-h(x,beta))-lambda*beta)
}

#Define Hessian
H <- function (x,y,beta,lambda) {
  return (-t(x)%*%x*diag(h(x,beta))*diag(1-h(x,beta))-lambda)
}

map <- optim(c(0, 0), function(beta) -log.p(X,y,beta,1), method = "L-BFGS", gr = function(beta) -grad(X,y,beta,1))

Mean <- map$par
Sigma <- solve(-H(X,y,Mean,1))
w <- 1.96 * sqrt(diag(Sigma))
Mean
#-0.3468521 -0.1246658
l <- Mean - w
#-0.4910317 -0.2689405
u <- Mean + w
#-0.20267238  0.01960878

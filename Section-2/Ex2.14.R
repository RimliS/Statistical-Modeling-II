dental <- read.csv("/Users/rimlisengupta/Dropbox/UT/Spring-2018/SDS383D-StatMod2/Section2/dental.csv",
                   sep = ",", header=TRUE)

loc <- c(1,4)
temp1 <- dental[,-loc]
temp1 <- data.frame(temp1)


library(MASS)

# length of chain
nit <- 5000
n <- 108
p <- 5

inds0 <- which(temp1$Sex=="Male")
inds1 <- which(temp1$Sex=="Female")

indsy10 <- which(temp1$age==10)
indsy12 <- which(temp1$age==12)
indsy14 <- which(temp1$age==14)

gender <- rep(0,n)
gender[inds0] <- 1
gender[inds1] <- 0
temp1$gender <- gender

y10 <- rep(0,n)
y12 <- rep(0,n)
y14 <- rep(0,n)

y10[indsy10] <- 1
y12[indsy12] <- 1
y14[indsy14] <- 1

temp1$y10 <- y10
temp1$y12 <- y12
temp1$y14 <- y14

#reorder by column index
temp1 <- temp1[c(1,2,3,5,6,7,4)]

X <- cbind(beta0=rep(1,n),temp1[,4:7])
X <- as.matrix(X)
K <- diag(p)

y <- temp1$distance

a <- 2
b <- 1
a.star <- a+(n+p)/2

# start chain
count = 0
beta <- matrix(0,p,nit)
w <- rep(0,nit)
lambda <- matrix(0,n,nit)
mu <- c(0,0,0,0,0)
beta.new <- rep(0,p)
w.new <- rep(0,nit)
w[1] <- 1
w.new[1] <- w[1]
lambda.new <- rep(0,n)
tau <- 1

L <- matrix(0,n,n)
D <- matrix(0,p,p)
S <- matrix(0,p,p)
XX.inv <- matrix(0,p,p)
M <- matrix(0,n,n)
mu.beta <- NULL

for(m in 2:nit) {
  for (i in 1:n) {
  lambda.new[i] <- rgamma(1,tau+0.5,tau+0.5*w[m]*(y[i]-as.numeric(X[i,]%*%beta[,m-1])))
  }
  L <- diag(lambda.new)
  D <- t(X)%*%L%*%X
  S <- solve(D+K)
  XX.inv <- solve(crossprod(X,X))
  M <- X%*%XX.inv%*%t(X)
  mu.beta <- as.vector(crossprod(S,(t(X) %*% M %*% y+K%*%mu)))
  
  beta.new <- mvrnorm(1,mu.beta, S*(1/w[m-1]))
  w.new <- rgamma(1,a.star,b+0.5*t((y-X%*%beta.new))%*%L%*%(y-X%*%beta.new)
                  +0.5*crossprod((beta.new-mu),(beta.new-mu)))
  
  beta[,m] <- beta.new
  w[m] <- w.new
  lambda[,m] <- lambda.new
  count = count + 1
}

dental <- read.csv("/Users/rimlisengupta/Dropbox/UT/Spring-2018/SDS383D-StatMod2/Section2/dental.csv",
                    sep = ",", header=TRUE)

loc <- c(1,4)
temp1 <- dental[,-loc]
temp1 <- data.frame(temp1)

lmfit <- lm(distance ~ factor(age)+factor(Sex),data=temp1)

library(MASS)
lmridge <- lm.ridge(distance ~ factor(age)+factor(Sex),lambda=0,data=temp1)

#lmridge
#age    SexMale 
#15.3856902  0.6601852  2.3210227


# length of chain
nit <- 10000
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

#temp2 <- as.matrix(temp1)[,-3]
X <- cbind(beta0=rep(1,n),temp1[,4:7])
X <- as.matrix(X)
K <- diag(p)
L <- diag(n)
mu <- c(0,0,0,0,0)

D <- t(X)%*%L%*%X
S <- solve(D+K)
y <- temp1$distance
XX.inv <- solve(crossprod(X,X))
M <- X%*%XX.inv%*%t(X)
mu.beta <- as.vector(crossprod(S,(t(X) %*% M %*% y+K%*%mu)))

a <- 2
b <- 1
a.star <- a+(n+p)/2

beta <- matrix(0,p,nit)
w <- rep(0,nit)

# start chain
count = 0
beta.new <- rep(0,p)
w.new <- rep(0,nit)
w[1] <- 1
w.new[1] <- w[1]

for(m in 2:nit) {
  beta.new <- mvrnorm(1,mu.beta, S*(1/w[m-1]))
  w.new <- rgamma(1,a.star,b+0.5*crossprod((y-X%*%beta.new),(y-X%*%beta.new))
                  +0.5*crossprod((beta.new-mu),(beta.new-mu)))
  
  beta[,m] <- beta.new
  w[m] <- w.new
  count = count + 1
}

col <- 1:5000

beta.hat <- c(mean(beta[1,-col]),mean(beta[2,-col]),mean(beta[3,-col]),
              mean(beta[4,-col]),mean(beta[5,-col]))
names(beta.hat) <- c("Intercept","beta1","beta2","beta3","beta4")

par(mfrow=c(3,2))
plot(beta[1,-col],ty="l",ylab="Intercept",main="traceplot of intercept")
plot(beta[2,-col],ty="l",ylab="beta1",main="traceplot of beta.1")
plot(beta[3,-col],ty="l",ylab="beta2",main="traceplot of beta.2")
plot(beta[4,-col],ty="l",ylab="beta3",main="traceplot of beta.3")
plot(beta[5,-col],ty="l",ylab="beta4",main="traceplot of beta.4")
hist(w[-col],breaks=30,xlab="w")


data <- read.table("/Users/rimlisengupta/Dropbox/UT/Spring-2018/SDS383D-StatMod2/Section3/pima.txt", fill = TRUE, 
                    header=TRUE)

df <- subset(data, select = c(3,6,8,9))
X.mat <- subset(data, select=c(3,6,8))
X.mat <- scale(X.mat)
final <- cbind(df[,4],X.mat)
final <- data.frame(final)

fit <- glm(final$V1~X.mat-1,family=binomial)
beta.hat <- coef(fit)

# starting point
#lambda0 <- 1
p <- 3
beta <- matrix(0,p,nit)

# length of chain
nit <- 5000
n <- 768

z <- matrix(1,n,nit)
inds0 <- which(final$V1==0)
z[inds0,1] <- -1

inds1 <- which(final$V1==1)
z[inds1,1] <- 1

# initialize
#x[1] <- x0
#u[1] <- u0


XX.inv <- solve(t(X.mat)%*%X.mat)

library(MASS)
library(truncnorm)

# start chain
count = 0
beta.new <- rep(0,p)
z.new <- rep(0,n)
for(m in 1:nit) {
  beta.new <- mvrnorm(1,XX.inv%*%t(X.mat)%*%z[,m], XX.inv)
  
 for (j in inds1) {
  z.new[j] <- rtruncnorm(n, a=0, b=Inf, mean=X.mat[j,]%*%beta.new, sd=1)
    }
    
  for (j in inds0) {
    z.new[j] <- rtruncnorm(n, a=-Inf, b=0, mean=X.mat[j,]%*%beta.new, sd=1)
    }
  
  beta[,m+1] <- beta.new
  z[,m+1] <- z.new
  count = count + 1
}


plot(beta[1,],ty="l")
plot(beta[2,],ty="l")
plot(beta[3,],ty="l")


plot(z[1,],ty="l")
hist(z[1,],freq=F,xlab="Z1",main="Histogram for Z1")
mean(z[1,])

plot(z[2,],ty="l")
hist(z[2,],freq=F,xlab="Z2",main="Histogram for Z2")
mean(z[2,])

plot(z[10,],ty="l")
hist(z[10,],freq=F,xlab="Z10",main="Histogram for Z10")
mean(z[10,])

plot(z[600,],ty="l")
hist(z[600,],freq=F,xlab="Z600",main="Histogram for Z600")
mean(z[600,])

mean(beta[1,]==0)
mean(beta[2,]==0)
mean(beta[3,]==0)

loc=1:1000
betaf <- beta[,-loc]

mean(betaf[1,])

#beta1: -0.1053866
#beta2: 0.4280609
#beta3: 0.3179886


#3.2
titanic <- read.csv("/Users/rimlisengupta/Dropbox/UT Austin/Spring 2018/SDS383D-StatMod2/Section3/titanic.csv",
                    sep = ",", header=TRUE)

df <- data.frame(titanic)
temp1 <- df[!is.na(df$Age),]

##Design matrix
x <- temp1[,3]
x.sc <- scale(x) ##Standardize

lookup <- c("No" = 0,"Yes"=1)
temp1$new_y <- lookup[temp1$Survived]
y <- temp1[,6] ##Response

##GLM Estimates:
fit <- glm(y~x.sc-1,family=binomial)
theta_glm <- coef(fit)
-0.1211981

#x <- as.matrix(x)
x.sc <- as.matrix(x.sc)

##Number of iterations:
num.iterations <- 5000

y <- matrix(y)
m <- nrow(x.sc)
#x.sc <- cbind(rep(1, m), x.sc) ##Including 1 in design matrix
num.features <- ncol(x.sc)

# Gradient descent function
grad <- function(x, y, theta) {
  gradient <- (1 / nrow(y)) * t(x) %*% ((1 / (1 + exp(-x %*% t(theta)))) - y)-theta
  return(t(gradient))
}

  num.iterations <- 10000
  alpha <- 0.04
  theta <- matrix(rep(0, num.features), nrow=1)
  count <- 0
  for (i in 2:num.iterations) {
    theta.new <- theta[i-1] - alpha * grad(x,y,theta[i-1])
    theta[i] <- theta.new
    count <- count+1
    }
  
  mean(theta)
  # -0.1286303 it=1000
  # -0.1286033 it=3000
  # -0.1284868 it=5000
  # -0.1284008 it=6000
  # -0.1282017 it=7000
  # -0.128117 it=10000
  
  l <- c(1000,3000,5000,6000,7000,10000)
  k <- c(-0.1286303,-0.1286033,-0.1284868,-0.1284008,-0.1282017,-0.128117)
  plot(l,k,ty="l")
  
  
#Prob 3.3-3.5
X <- as.matrix(cbind(rep(1, nrow(temp1)), scale(temp1$Age)))
y <- as.numeric(temp1$Survived == "Yes")

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
l <- Mean - 1.96 * sqrt(diag(Sigma))
u <- Mean + 1.96 * sqrt(diag(Sigma))


#3.7
tea <- tea_discipline_oss
tea <- tea[tea$ACTIONS > 0, ]

X <- as.matrix(cbind(rep(1, nrow(tea)), scale(as.numeric(tea$GRADE))))
y <- tea$ACTIONS

map <- optim(c(0, 0), function(beta) -log.p(X,y,beta,1), 
             gr = function(beta) -grad(X,y,beta,1))

mu <- map$par
cov <- solve(-H(X,y,mu,1))
l <- mu - 1.96 * sqrt(diag(cov))
u <- mu + 1.96 * sqrt(diag(cov))



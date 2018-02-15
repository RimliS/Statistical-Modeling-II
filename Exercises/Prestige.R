
#Reading data
adata <- read.table("/Users/rimlisengupta/Dropbox/SDS383D-StatMod2/prestige.txt", fill = TRUE, header=TRUE)

#Getting rid of categorical variable
prestige = adata[,1:5]
prestige<-na.omit(prestige)

#Getting rid of the first column
prestige <- prestige[, -1]

#Standardize
X.mat<-prestige[,-2]
X.mat<-scale(X.mat)

final<-cbind(prestige[,2],X.mat)
final<-data.frame(final)
final$income <- final$V1
final <- final[,-1]

mod <- lm(income ~ education+women+prestige, data=final)
summary(mod)

#Coefficients:
#            Estimate  Std. Error  t value    Pr(>|t|)    
#(Intercept)   6797.9    254.9     26.665     < 2e-16 ***
#education     483.5     511.9      0.944       0.347    
#women        -1614.7    271.4     -5.948       4.19e-08 ***
#prestige      2433.3    514.6      4.729       7.58e-06 ***

#sigma^2 unknown
beta0 <- rep(1,102) #Define intercept column
Xmat <- cbind(beta0,X.mat) #Design matrix
y <- final$income #response

qr(Xmat)$rank #check rank, full rank (X'X)^(-1) is invertible
#4

A=t(Xmat)%*%Xmat
B=solve(A)
beta=B%*%t(Xmat)%*%y #estimates of coefficients

M=Xmat%*%B%*%t(Xmat) #Hat matrix 

#Using {y'(I-M)y}/(n-p) where p is the rank of X
C=diag(102)-M 
e.sq=as.vector(t(y)%*%C%*%y)

sigma.hat.sq=e.sq/(102-4)

#variance of least square:
var.beta=sigma.hat.sq*B

sqrt(var.beta[1,1])
# 254.9342
sqrt(var.beta[2,2])
# 511.9442
sqrt(var.beta[3,3])
# 271.4444
sqrt(var.beta[4,4])
# 514.5795

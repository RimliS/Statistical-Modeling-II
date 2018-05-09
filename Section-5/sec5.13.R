library(MCMCpack)

vect <- rep(10,100)
pi <- rdirichlet(1,vect)
z <- rmultinom(10, 1, pi)

k <- rep(0,100)

for (i in 1:100) {
  k[i] <- length(which(z[i,]>0))
}

inds1 <- which(k > 0)
inds2 <- which(k > 0)
inds3 <- which(k > 0)
inds4 <- which(k > 0)
inds5 <- which(k > 0)

# 3  6 24 29 30 37 43 54 95
# 8 14 33 35 37 45 56 82 88
# 3 12 26 31 38 53 76 88 93
# 2  12  16  20  35  36  55  89 100
# 13 15 24 30 32 41 47 61 85 97



vect <- rep(1,100)
pi <- rdirichlet(1,vect)
z <- rmultinom(10, 1, pi)

k <- rep(0,100)

for (i in 1:100) {
  k[i] <- length(which(z[i,]>0))
}

inds1 <- which(k > 0)
inds2 <- which(k > 0)
inds3 <- which(k > 0)
inds4 <- which(k > 0)
inds5 <- which(k > 0)

# 10 27 36 40 43 46 53 64 85
# 8 14 33 35 37 45 56 82 88
# 3 12 26 31 38 53 76 88 93
# 2  12  16  20  35  36  55  89 100
# 13 15 24 30 32 41 47 61 85 97

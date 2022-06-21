rrr.sim3<-function (n = n, 
                    p = p, rho.x = 0.0,
                    q.mix = c(5, 20, 5), 
                    nrank = 2, 
                    intercept = rep(0.5,q), 
                    mis.prop = 0.2) 
{
  q <- sum(q.mix)
  q1 <- q.mix[1]
  S = matrix(rho.x, ncol = p, nrow = p); diag(S) = 1
  out = eigen(S, symmetric = TRUE)
  S.sqrt = out$vectors %*% diag(out$values^0.5) %*% t(out$vectors)
  
  X <- matrix(rnorm(n * p), n, p)%*% S.sqrt
  u0 <- matrix(nrow = p, ncol = nrank, rnorm(p * nrank, 0,1))
  u0 <- qr.Q(qr(u0))
  v0 <- matrix(nrow = q, ncol = nrank, rnorm(q * nrank, 0,1))
  C <- u0 %*% t(v0) 
  C0 <- rbind(intercept, C)
  X0 <- cbind(1, X)
  MU <- X0 %*% C0
  family <- list(gaussian(), binomial(), poisson())
  familygroup <- c(rep(1, q1))
  cfamily <- unique(familygroup)
  nfamily <- length(cfamily)
  Y <- matrix(nrow = n, ncol = q, 0)
  sigma <- 1
  if (sum(familygroup == 1) > 0) {
    Y[, familygroup == 1] <- MU[, familygroup == 1] + 
      matrix(nrow = n,ncol = q1, rnorm(n * q1, 0, sigma)) # rt(n*q1, df=3)
  }
  
  N <- n * q
  if (is.numeric(mis.prop) & mis.prop < 1 & mis.prop > 0) {
    N_mis <- round(N * (1-mis.prop) )
    index_mis <- sample(1:N, size = N_mis, replace = TRUE) # FALSE
    Y_mis <- matrix(nrow = n, ncol = q, NA)#Y
    Y_mis[index_mis] <- Y[index_mis]#NA
  }
  else {
    index_mis <- NULL
    Y_mis <- Y
  }
  list(Y = Y, Y.mis = Y_mis, index.miss = index_mis, X = X, 
       C = C, family = family[q.mix != 0], familygroup = familygroup)
}

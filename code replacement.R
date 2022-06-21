library(rrpack)
source('/Users/thetm/Dropbox/ongoing_works/B RRR/BRRR with missing Y/Rcodes/rrr_sim3.R')
library(Rcpp)
sourceCpp("/Users/thetm/Dropbox/ongoing_works/B RRR/BRRR with missing Y/Rcodes/test.cpp")

q = 8
p = 12
n = 100 # ell
r = 2
missfrac = 0.8
rho.x = 0

u0 <- matrix(nrow = p, ncol = r, rnorm(p * r, 0,1))
u0 <- qr.Q(qr(u0))
v0 <- matrix(nrow = q, ncol = r, rnorm(q * r, 0,1))
C <- u0 %*% t(v0) *2 + rnorm(p*q,0, 0.1)
X <- matrix(rnorm(n*p), n, p)
tX = t(X)
Y = matrix(0, nrow = n,ncol = q)

C0 <- rbind(rep(1,q), C)
X0 <- cbind(1, X)
MU <- X %*% C
N <- n * q
Nobs <- round(N * (1-missfrac) )
index_obs <- sample(1:N, size = Nobs, replace = TRUE) # FALSE
Yobs <- rep(NA, Nobs)
for (i in 1:Nobs) Yobs[i] = MU[ index_obs[i] ] + rnorm(1)

Y[unique(index_obs)] <- sapply(unique(index_obs) , function(i)mean(Yobs[index_obs==i]))

Yna = matrix(NA, nrow = n,ncol = q)
Yna[unique(index_obs)]=Y[unique(index_obs)]

simdata <- rrr.sim3(n = n, p = p, q.mix = c(q, 0, 0),intercept = rep(1,q),rho.x = rho.x,
                    nrank = r, mis.prop = missfrac)
family <- simdata$family ; familygroup <- simdata$familygroup
control = list(epsilon = 1e-4, sv.tol = 1e-2, maxit = 2000, trace = FALSE, gammaC0 = 1.1, plot.cv = F,  conv.obj = TRUE)
fit.cv.mrrr <- cv.mrrr(Yna , X, family = family, familygroup = familygroup,
                       maxrank = 8,control = control, 
                       penstr = list(penaltySVD = "rankCon"),
                       nfold = 5, nlam = 50)
hatc = coef(fit.cv.mrrr$fit)[-1,]

P_i = function(Yobs,XM,entri){ 
  A = matrix(0, nrow = n,ncol = q)
  newid = index_obs[entri]
  A[newid] = Yobs[entri] -XM[newid]
  return(A)
}
P_i <- compiler::cmpfun(P_i)

Iters = 3000
burnin = 1000
h =  1/(p*q)^1

M = hatc
tau = 10
ystar = diag(p)*tau^2
Bm = matrix(0,nr=p,nc=q)
a = 0; sig2 = 1
h_tX = h*tX 
for(s in 1:Iters){
  ystar_Mt_M = ystar + tcrossprod(M)
  tam1 = solve(ystar_Mt_M , M )
  XM = X%*%M
  ss1 = 0; for (i in 1:Nobs)ss1 = ss1 + eigenMapMatMult(h_tX , P_i(Yobs,XM,i))
  tam = M + ss1 /sig2 + h*(p+q+2)*tam1 + sqrt(2*h)*rnorm(p*q)
  #
  ystar_tamt_tam = ystar + tcrossprod(tam)
  logdet = determinant( ystar_tamt_tam )
  
  Xtam = X%*%tam
  pro.tam = -sum((Yobs - Xtam[index_obs])^2 )/(2*sig2) - 0.5*(p+q+2)*logdet$modulus*logdet$sign
  logdet = determinant( ystar_Mt_M )
  pro.M = -sum((Yobs - XM[index_obs])^2 )/(2*sig2) - 0.5*(p+q+2)*logdet$modulus*logdet$sign
  
  ss= 0; for (i in 1:Nobs)ss = ss + eigenMapMatMult(h_tX , P_i(Yobs,Xtam,i))
  tran.m = -sum(( M-tam- ss/sig2 - h*(p+q+2)*solve(a = ystar_tamt_tam,tam) )^2 )/(4*h)
  
  tran.tam = -sum((tam-M-ss1/sig2 -h*(p+q+2)*tam1 )^2 )/(4*h)
  pro.trans = pro.tam+tran.m-pro.M-tran.tam
  if(log(runif(1)) < pro.trans){
    M = tam
    a = a+1
  } 
  if (s>burnin) Bm = Bm + M/(Iters-burnin)
}
a/Iters



### LMC
M = hatc
h =  h/2
ystar = diag(p)*tau^2
M.lmc = matrix(data=0,nr=p,nc=q)
start_time <- Sys.time()
for(s in 1:Iters){
  XM = eigenMapMatMult(X,M)
  ss1 =  0; for (i in 1:Nobs)ss1 = ss1 + eigenMapMatMult(h_tX , P_i(Yobs,XM,i))
  M = M + ss1/sig2 + h*(p+q+2)*solve(ystar + tcrossprod(M),M) + sqrt(2*h)*rnorm(p*q)
  if (s>burnin) M.lmc = M.lmc + M/(Iters-burnin)
}
end_time <- Sys.time();end_time - start_time
a/Iters
mean((X%*%Bm - MU )^2);
mean((X%*%hatc - MU )^2);
mean((X%*%M.lmc - MU )^2)

(sum((X%*%Bm - MU )^2) -sum(( (X%*%Bm)[unique(index_obs)] - MU[unique(index_obs)] )^2) )/(prod(dim(Y)) - length(unique(index_obs)))
(sum((X%*%hatc - MU )^2) -sum(( (X%*%hatc)[unique(index_obs)] - MU[unique(index_obs)] )^2) )/(prod(dim(Y)) - length(unique(index_obs)))
(sum((X%*%M.lmc - MU )^2) -sum(( (X%*%M.lmc)[unique(index_obs)] - MU[unique(index_obs)] )^2) )/(prod(dim(Y)) - length(unique(index_obs)))



library(rrr) # to get data
setwd("~/Dropbox/ongoing_works/B RRR/BRRR with missing Y/Rcodes")
load("~/Dropbox/ongoing_works/B RRR/BRRR with missing Y/Rcodes/galaxy_full.rda")
library(Rcpp)
sourceCpp("~/Documents/testtest/test.cpp")
P_Omega = function(a,entri){ a[entri] = 0; return(a) }
P_Omega <- compiler::cmpfun(P_Omega)
library(dplyr)
data(COMBO17)
galaxy <- as_data_frame(COMBO17); rm(COMBO17)
galaxy <- select(galaxy, -starts_with("e."), -Nr, -UFS:-IFD)
galaxy <- na.omit(galaxy)
X <- scale( as.matrix(select(galaxy, -Rmag:-chi2red) ))
Y <- scale( as.matrix(select(galaxy, Rmag:chi2red) )); rm(galaxy)
n <- nrow(Y)
p <- ncol(X)
q <-  ncol(Y)
tX =t(X)

B.ols <- solve(tX %*%X, tX %*%Y) ; dimnames(B.ols) <- NULL
missfrac = .8

library(rrpack)
control = list(epsilon = 1e-4, sv.tol = 1e-2, maxit = 1000,
               trace = FALSE, gammaC0 = 1.1, plot.cv = F, conv.obj = TRUE)

mrrr = mala = lmc = list()
ac = c()

for (ss in 87:100) {

  Yna = Y 
  imiss = sample(seq(n*q),n*q*missfrac,replace=FALSE)
  Yna[imiss] = NA
  fit.cv.mrrr <- cv.mrrr(Yna, X, family = list(gaussian()),control = control,
                         penstr = list(penaltySVD = "rankCon",lambdaSVD = c(1:min(p,q)) ),
                         maxrank = min(p,q) , nfold = 5)$fit
  hatc = fit.cv.mrrr$coef[- 1 ,]
  Iters = 15000
  burnin = 5000
  h =  1/(p*q)^1.9
  h_tX = h*tX ; 
  M =  hatc
  lam = 10
  ystar = diag(p)*lam^2
  Bm = matrix(data=0,nr=p,nc=q)
  a = 0; sig2 = 1
  for(s in 1:Iters){
    ystar_Mt_M = ystar + tcrossprod(M)
    tam1 = solve(a = ystar_Mt_M,M)
    P_Ome_YnaXM = P_Omega(Yna - eigenMapMatMult(X,M) ,imiss)
    tam = M + eigenMapMatMult(h_tX,P_Ome_YnaXM )/sig2 + h*(p+q+2)*tam1 + sqrt(2*h)*rnorm(p*q)
    #
    ystar_tamt_tam = ystar + tcrossprod(tam)
    logdet = determinant( ystar_tamt_tam )
    Xtam = eigenMapMatMult(X,tam)
    pro.tam = -sum((Yna - Xtam)^2,na.rm = TRUE)/(2*sig2) - 0.5*(p+q+2)*logdet$modulus*logdet$sign
    logdet = determinant( ystar_Mt_M )
    pro.M = -sum(P_Ome_YnaXM^2, na.rm = TRUE)/(2*sig2) - 0.5*(p+q+2)*logdet$modulus*logdet$sign
    
    tam2 = solve(a = ystar_tamt_tam,tam)
    tran.m = -sum((M-tam-eigenMapMatMult(h_tX,P_Omega(Yna - Xtam,imiss) )/sig2 - h*(p+q+2)*tam2 )^2)/(4*h)
    tran.tam = -sum((tam-M-eigenMapMatMult(h_tX, P_Ome_YnaXM )/sig2 - h*(p+q+2)*tam1 )^2)/(4*h)
    
    pro.trans = pro.tam+tran.m-pro.M-tran.tam
    if(log(runif(1)) < pro.trans){
      M = tam
      a = a+1
    } 
    if (s>burnin){
      Bm = Bm + M/(Iters-burnin)
    } 
  }

  M = hatc #
  h =  h/2
  M.lmc = matrix(data=0,nr=p,nc=q)
  for(s in 1:Iters){
    tam = solve(a = ystar + tcrossprod(M),M)
    M = M + eigenMapMatMult(h_tX,P_Omega(Yna - X%*%M,imiss) )/sig2 + h*(p+q+2)*tam + sqrt(2*h)*rnorm(p*q)
    if (s>burnin){
      M.lmc = M.lmc + M/(Iters-burnin)
    } 
  }
  mala[[ss]] = c(mean((X%*%Bm - X%*%B.ols)^2) ,  mean(( (X%*%Bm)[imiss] - Y[imiss])^2)  )
  mrrr[[ss]] =  c(mean((X%*%hatc - X%*%B.ols)^2) , mean(( (X%*%hatc)[imiss] - Y[imiss])^2) )
  lmc[[ss]] = c(mean((X%*%M.lmc - X%*%B.ols)^2), mean(( (X%*%M.lmc)[imiss] - Y[imiss])^2) )
  print(paste(ss, ac[ss]<- a/Iters) )
}
save.image(file = "~/Dropbox/ongoing_works/B RRR/BRRR with missing Y/Rcodes/outsimulations/galaxy_out_08.rda")


mean(sapply(lmc, function(x)x[1]) ) ; sd(sapply(lmc, function(x)x[1]) )
mean(sapply(lmc, function(x)x[2]) ) ; sd(sapply(lmc, function(x)x[2]) )

mean(sapply(mala, function(x)x[1]) ) ; sd(sapply(mala, function(x)x[1]) )
mean(sapply(mala, function(x)x[2]) ) ; sd(sapply(mala, function(x)x[2]) )

mean(sapply(mrrr, function(x)x[1]) ) ; sd(sapply(mrrr, function(x)x[1]) )
mean(sapply(mrrr, function(x)x[2]) ) ; sd(sapply(mrrr, function(x)x[2]) )




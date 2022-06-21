library(rrr) # to get data
setwd("~/Dropbox/ongoing_works/B RRR/BRRR with missing Y/Rcodes")
library(Rcpp)
sourceCpp("~/Documents/testtest/test.cpp")
P_Omega = function(a,entri){ a[entri] = 0; return(a) }
P_Omega <- compiler::cmpfun(P_Omega)
library(rrpack)

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
B.ols <- solve(t(X) %*%X, t(X) %*%Y) ; dimnames(B.ols) <- NULL

missfrac = .001
Yna = Y
imiss = sample(seq(n*q),n*q*missfrac,replace=FALSE)
Yna[imiss] = NA

control = list(epsilon = 1e-4, sv.tol = 1e-2, maxit = 2000,
               trace = FALSE, gammaC0 = 1.1, plot.cv = F,  conv.obj = TRUE)
fit.cv.mrrr <- cv.mrrr(Yna, X, family = list(gaussian()),control = control,
                       penstr = list(penaltySVD = "rankCon",lambdaSVD = c(1:min(p,q)) ),
                       maxrank = min(p,q) , nfold = 5,nlam = 50)$fit
hatc = fit.cv.mrrr$coef[- 1 ,]

mean((X%*%hatc - X%*%B.ols)^2)
mean(( (X%*%hatc)[imiss] - Y[imiss])^2) 


######################################
#####  MALA  full  ############
Iters = 5000
burnin = 2000
h =  1/(p*q)^2.25
M = hatc 
tau = 10
ystar = diag(p)*tau^2
Bm = matrix(data=0,nr=p,nc=q)
a = 0; sig2 = 1
tX =t(X); h_tX = h*tX ; 
start_time <- Sys.time()
entri.12 = entri.102 = c()
for(s in 1:Iters){
  ystar_Mt_M = ystar + tcrossprod(M)
  tam1 = solve(ystar_Mt_M, M)
  YnaXM = Yna - eigenMapMatMult(X,M)
  tam = M + h_tX%*%P_Omega(YnaXM,imiss)/sig2 + h*(p+q+2)*tam1 + sqrt(2*h)*rnorm(p*q)
  #
  ystar_tamt_tam = ystar + tcrossprod(tam)
  logdet = determinant( ystar_tamt_tam )
  Yna_Xtam = Yna -eigenMapMatMult(X,tam)
  pro.tam = -sum(Yna_Xtam^2,na.rm = TRUE)/(2*sig2) - 0.5*(p+q+2)*logdet$modulus*logdet$sign
  logdet = determinant( ystar_Mt_M )
  pro.M = -sum(YnaXM^2, na.rm = TRUE)/(2*sig2) - 0.5*(p+q+2)*logdet$modulus*logdet$sign
  
  tam2 = solve(a = ystar_tamt_tam,tam)
  tran.m = -sum((M-tam-h_tX%*%P_Omega(Yna_Xtam,imiss)/sig2 - h*(p+q+2)*tam2 )^2)/(4*h)
  tran.tam = -sum((tam-M-h_tX%*%P_Omega(YnaXM,imiss)/sig2 - h*(p+q+2)*tam1 )^2)/(4*h)
  
  pro.trans = pro.tam+tran.m-pro.M-tran.tam
  if(log(runif(1)) < pro.trans){
    M = tam
    a = a+1
  } 
  if (s>burnin){   Bm = Bm + M/(Iters-burnin) } 
 # entri.12[s] = M[1,2]
 # entri.102[s] = M[10,2]
}
end_time <- Sys.time();end_time - start_time
a/Iters


### Langevin MC full
M = hatc #
h =  h/2
M.lmc = matrix(data=0,nr=p,nc=q)
for(s in 1:Iters){
  tam = solve(a = ystar + tcrossprod(M),M)
  M = M + h_tX%*%P_Omega(Yna - X%*%M,imiss)/sig2 + h*(p+q+2)*tam + sqrt(2*h)*rnorm(p*q)
  if (s>burnin)M.lmc = M.lmc + M/(Iters-burnin)
}

mean((X%*%hatc - X%*%B.ols)^2); mean(( (X%*%hatc)[imiss] - Y[imiss])^2) 
mean((X%*%Bm - X%*%B.ols)^2); mean(( (X%*%Bm)[imiss] - Y[imiss])^2) 
mean((X%*%M.lmc - X%*%B.ols)^2); mean(( (X%*%M.lmc)[imiss] - Y[imiss])^2) 

a/Iters

hatc.rr <- cv.rrr(Y,X,nfold = 5)
hatc.rr = coef(hatc.rr)
mean((X%*%hatc.rr - X%*%B.ols)^2); mean(( (X%*%hatc.rr)[imiss] - Y[imiss])^2) 


library(rrpack)
library(Rcpp)
sourceCpp("/Users/thetm/Dropbox/ongoing_works/B RRR/BRRR with missing Y/Rcodes/test.cpp")
P_Omega = function(a,entri){ a[entri] = 0; return(a) }; P_Omega <- compiler::cmpfun(P_Omega)
source('/Users/thetm/Dropbox/ongoing_works/B RRR/BRRR with missing Y/Rcodes/rrr_sim3.R')

q = 8
p = 12
n = 100
r = 2
missfrac = 0.2
rho.x = 0.5

lmc = mala = mrrpack = list()
ac = c()
for (ss in 1:200) {
  simdata <- rrr.sim3(n = n, p = p, q.mix = c(q, 0, 0),intercept = rep(1,q),rho.x = rho.x,
                      nrank = r, mis.prop = missfrac)
  Y <- simdata$Y
  Yna <- simdata$Y.mis
  X <- simdata$X
  X0 <- cbind(1,X)
  imiss <- simdata$index.miss
  C <- simdata$C
  family <- simdata$family
  familygroup <- simdata$familygroup
  svdX0d1 <- svd(X0)$d[1]
  init1 = list(kappaC0 = svdX0d1 * 5)
  offset = NULL
  control = list(epsilon = 1e-4, sv.tol = 1e-2, maxit = 2000,
                 trace = FALSE, gammaC0 = 1.1, plot.cv = F, conv.obj = TRUE)
  fit.cv.mrrr <- cv.mrrr(Yna , X, family = family,
                         familygroup = familygroup,
                         maxrank = 8,
                         penstr = list(penaltySVD = "rankCon",
                                       lambdaSVD = c(1 : 6)),
                         control = control, init = init1,
                         nfold = 5, nlam = 50)
  hatc = coef(fit.cv.mrrr$fit)[-1,]
  
  Iters = 5000
  burnin = 2000
  h =  1/(p*q)^1.41
  M =  hatc
  tau = 10
  ystar = diag(p)*tau^2
  Bm = matrix(data=0,nr=p,nc=q)
  a = 0; sig2 = 1
  tX =t(X); h_tX = h*tX 
  Yna = simdata$Y.mis 
  for(s in 1:Iters){
    ystar_Mt_M = ystar + tcrossprod(M)
    tam1 = solve(a = ystar_Mt_M,M)
    YnaXM = Yna - eigenMapMatMult(X,M)
    tam = M + eigenMapMatMult(h_tX,P_Omega(YnaXM,imiss) )/sig2 + h*(p+q+2)*tam1 + sqrt(2*h)*rnorm(p*q)
    #
    ystar_tamt_tam = ystar + tcrossprod(tam)
    logdet = determinant( ystar_tamt_tam )
    Xtam = eigenMapMatMult(X,tam)
    pro.tam = -sum((Yna - Xtam)^2,na.rm = TRUE)/(2*sig2) - 0.5*(p+q+2)*logdet$modulus*logdet$sign
    logdet = determinant( ystar_Mt_M )
    pro.M = -sum(YnaXM^2, na.rm = TRUE)/(2*sig2) - 0.5*(p+q+2)*logdet$modulus*logdet$sign
    
    tam2 = solve(a = ystar_tamt_tam,tam)
    tran.m = -sum((M-tam- eigenMapMatMult(h_tX,P_Omega(Yna - Xtam,imiss))/sig2 - h*(p+q+2)*tam2 )^2)/(4*h)
    tran.tam = -sum((tam-M-eigenMapMatMult(h_tX,P_Omega(YnaXM,imiss) )/sig2 - h*(p+q+2)*tam1 )^2)/(4*h)
    
    pro.trans = pro.tam+tran.m-pro.M-tran.tam
    if(log(runif(1)) < pro.trans){
      M = tam
      a = a+1
    } 
    if (s>burnin){
      Bm = Bm + M/(Iters-burnin)
    } 
  }
  
  ### LMC
  M = hatc #
  h =  h/2
  ystar = diag(p)*tau^2
  M.lmc = matrix(data=0,nr=p,nc=q)
  for(s in 1:Iters){
    tam = solve(a = ystar + tcrossprod(M),M)
    M = M + eigenMapMatMult(h_tX,P_Omega(Yna - X%*%M,imiss) )/sig2 + h*(p+q+2)*tam + sqrt(2*h)*rnorm(p*q)
    if (s>burnin) M.lmc = M.lmc + M/(Iters-burnin)
  }
  print(ac[ss]<- a/Iters)
  
  lmc[[ss]] =c(mean((X%*%M.lmc - X%*%C )^2), mean(( (X%*%M.lmc)[imiss] - Y[imiss])^2))
  mrrpack[[ss]] = c(mean((X%*%C - X%*%hatc)^2), mean(( (X%*%hatc)[imiss] - Y[imiss])^2) )
  mala[[ss]] = c(mean((X%*%Bm - X%*%C )^2), mean(( (X%*%Bm)[imiss] - Y[imiss])^2) )
  print(ss)
}
save.image(file ='~/Dropbox/ongoing_works/B RRR/BRRR with missing Y/Rcodes/outsimulations/simu4_rx05_mf02.rda')


######################################



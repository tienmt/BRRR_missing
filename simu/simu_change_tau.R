library(rrpack)
library(Rcpp)
sourceCpp("~/Documents/testtest/test.cpp")
P_Omega = function(a,entri){ a[entri] = 0; return(a) }
P_Omega <- compiler::cmpfun(P_Omega)
n = 100  # samples
q=l = 8  # response
p = 12  # predictors
r = 3   # true rank
missfrac = 0.5
rho.x = 0.0
S = matrix(rho.x, ncol = p, nrow = p); diag(S) = 1
out = eigen(S, symmetric = TRUE)
S.sqrt = out$vectors %*% diag(out$values^0.5) %*% t(out$vectors)
control = list(epsilon = 1e-4, sv.tol = 1e-2, maxit = 1000, trace = FALSE, gammaC0 = 1.1, plot.cv = F, conv.obj = TRUE)

ac = c() ; lmc_tau10 = mala_tau10  = mrrr =
  lmc_tau100 = mala_tau100 = 
  lmc_tau1 = mala_tau1=
  lmc_tau01 = mala_tau01 = list() 
for (ss in 1:100) {
  X = matrix(rnorm(n*p),nr=n) %*% S.sqrt
  B = matrix(rnorm(p*r),nr=p)%*%t(matrix(rnorm(l*r),nr=l))# + rnorm(p*l)
  Yna = Y = X%*%B + rnorm(n*l) # rt(n*l,df=2)#
  imiss = sample(seq(n*l),n*l*missfrac,replace=FALSE)
  Yna[imiss] = NA
  Y_mis = Yna
  # mRRR
  fit.cv.mrrr <- cv.mrrr(Yna, X, family = list(gaussian()),control = control,
                         penstr = list(penaltySVD = "rankCon",lambdaSVD = c(1:min(p,q)) ),
                         maxrank = min(p,q) , nfold = 5)$fit
  hatc = fit.cv.mrrr$coef[-1 ,]
  mrrr[[ss]] = c(mean((X%*%hatc - X%*%B)^2), mean((hatc - B )^2)/mean(B^2),mean(( (X%*%hatc)[imiss] - Y[imiss])^2) )
  #####  MALA    ############
  Iters = 5000
  burnin = 2000
  h =  1/(p*q)^1
  M = hatc
  lam = 10
  ystar = diag(p)*lam^2
  M.mala = matrix(data=0,nr=p,nc=q)
  a = 0; sig2 = 2
  tX =t(X); h_tX = h*tX
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
      M.mala = M.mala + M/(Iters-burnin)
    } 
  }
  mala_tau10[[ss]] = c(mean((X%*%M.mala - X%*%B )^2), mean((M.mala - B )^2)/mean(B^2) , mean(( (X%*%M.mala)[imiss] - Y[imiss])^2))
  
  print(ac[ss] <- a/Iters)
  #Langevin MC for BRRR
  M = hatc #
  Iters = 5000
  burnin = 2000
  h =  h/1.5
  M.lmc = matrix(data=0,nr=p,nc=q)
  for(s in 1:Iters){
    tam = solve(a = ystar + tcrossprod(M),M)
    M = M + eigenMapMatMult(h_tX,P_Omega(Yna - X%*%M,imiss) )/sig2 + h*(p+q+2)*tam + sqrt(2*h)*rnorm(p*q)
    if (s>burnin){
      M.lmc = M.lmc + M/(Iters-burnin)
    } 
  }
  lmc_tau10[[ss]] =c(mean((M.lmc -B)^2),  mean((M.lmc -B )^2)/mean(B^2), mean(( (X%*%M.lmc)[imiss] - Y[imiss])^2) )
  print(ss)
  
  
  h =  1/(p*q)^.98
  M = hatc
  lam = 100
  ystar = diag(p)*lam^2
  M.mala = matrix(data=0,nr=p,nc=q)
  a = 0; sig2 = 2
  tX =t(X); h_tX = h*tX
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
      M.mala = M.mala + M/(Iters-burnin)
    } 
  }
  mala_tau100[[ss]] = c(mean((X%*%M.mala - X%*%B )^2), mean((M.mala - B )^2)/mean(B^2) , mean(( (X%*%M.mala)[imiss] - Y[imiss])^2))
  
  print(ac[ss] <- a/Iters)
  #Langevin MC for BRRR
  M = hatc #
  Iters = 5000
  burnin = 2000
  h =  h/1.5
  M.lmc = matrix(data=0,nr=p,nc=q)
  for(s in 1:Iters){
    tam = solve(a = ystar + tcrossprod(M),M)
    M = M + eigenMapMatMult(h_tX,P_Omega(Yna - X%*%M,imiss) )/sig2 + h*(p+q+2)*tam + sqrt(2*h)*rnorm(p*q)
    if (s>burnin){
      M.lmc = M.lmc + M/(Iters-burnin)
    } 
  }
  lmc_tau100[[ss]] =c(mean((M.lmc -B)^2),  mean((M.lmc -B )^2)/mean(B^2), mean(( (X%*%M.lmc)[imiss] - Y[imiss])^2) )
  print(ss)
  
  
  h =  1/(p*q)^1.6
  M = hatc
  lam = 1
  ystar = diag(p)*lam^2
  M.mala = matrix(data=0,nr=p,nc=q)
  a = 0; sig2 = 2
  tX =t(X); h_tX = h*tX
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
      M.mala = M.mala + M/(Iters-burnin)
    } 
  }
  mala_tau1[[ss]] = c(mean((X%*%M.mala - X%*%B )^2), mean((M.mala - B )^2)/mean(B^2) , mean(( (X%*%M.mala)[imiss] - Y[imiss])^2))
  
  print(ac[ss] <- a/Iters)
  #Langevin MC for BRRR
  M = hatc #
  Iters = 5000
  burnin = 2000
  h =  h/1.5
  M.lmc = matrix(data=0,nr=p,nc=q)
  for(s in 1:Iters){
    tam = solve(a = ystar + tcrossprod(M),M)
    M = M + eigenMapMatMult(h_tX,P_Omega(Yna - X%*%M,imiss) )/sig2 + h*(p+q+2)*tam + sqrt(2*h)*rnorm(p*q)
    if (s>burnin){
      M.lmc = M.lmc + M/(Iters-burnin)
    } 
  }
  lmc_tau1[[ss]] =c(mean((M.lmc -B)^2),  mean((M.lmc -B )^2)/mean(B^2), mean(( (X%*%M.lmc)[imiss] - Y[imiss])^2) )
  print(ss)
  
  
  h =  1/(p*q)^2.62
  M = hatc
  lam = .1
  ystar = diag(p)*lam^2
  M.mala = matrix(data=0,nr=p,nc=q)
  a = 0; sig2 = 2
  tX =t(X); h_tX = h*tX
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
      M.mala = M.mala + M/(Iters-burnin)
    } 
  }
  mala_tau01[[ss]] = c(mean((X%*%M.mala - X%*%B )^2), mean((M.mala - B )^2)/mean(B^2) , mean(( (X%*%M.mala)[imiss] - Y[imiss])^2))
  
  print(ac[ss] <- a/Iters)
  #Langevin MC for BRRR
  M = hatc #
  Iters = 5000
  burnin = 2000
  h =  h/1.5
  M.lmc = matrix(data=0,nr=p,nc=q)
  for(s in 1:Iters){
    tam = solve(a = ystar + tcrossprod(M),M)
    M = M + eigenMapMatMult(h_tX,P_Omega(Yna - X%*%M,imiss) )/sig2 + h*(p+q+2)*tam + sqrt(2*h)*rnorm(p*q)
    if (s>burnin){
      M.lmc = M.lmc + M/(Iters-burnin)
    } 
  }
  lmc_tau01[[ss]] =c(mean((M.lmc -B)^2),  mean((M.lmc -B )^2)/mean(B^2), mean(( (X%*%M.lmc)[imiss] - Y[imiss])^2) )
  print(ss)
  
}


save.image(file ='~/Dropbox/ongoing_works/B RRR/BRRR with missing Y/Rcodes/outsimulations/simu_change_tau.rda')


######################################



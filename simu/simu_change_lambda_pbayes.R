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
  control = list(epsilon = 1e-4, sv.tol = 1e-2, maxit = 1500, trace = FALSE, gammaC0 = 1.1, plot.cv = F, conv.obj = TRUE)
  
  ac = ac2 = ac8 = ac32 = c()  
  lmc_lam1 = mala_lam1 =
    lmc_lam2 = mala_lam2  = 
    lmc_lam8 = mala_lam8 = 
    lmc_lam32 = mala_lam32 =list() 
  for (ss in 1:100) {
    X = matrix(rnorm(n*p),nr=n) %*% S.sqrt
    B = matrix(rnorm(p*r),nr=p)%*%t(matrix(rnorm(l*r),nr=l))*2 + rnorm(p*l,sd = .1)
    Yna = Y = X%*%B + rnorm(n*l) # rt(n*l,df=2)#
    imiss = sample(seq(n*l),n*l*missfrac,replace=FALSE)
    Yna[imiss] = NA
    Y_mis = Yna
    fit.cv.mrrr <- cv.mrrr(Yna, X, family = list(gaussian()),control = control,
                           penstr = list(penaltySVD = "rankCon",lambdaSVD = c(1:min(p,q)) ),
                           maxrank = min(p,q) , nfold = 5)$fit
    hatc = fit.cv.mrrr$coef[- 1 ,]
  
    ######################################
    #####  MALA    ############
    Iters = 10000
    burnin = 2000
    h =  1/(p*q)^1.28
    M =  hatc
    tau = 10
    ystar = diag(p)*tau^2
    Bm = matrix(data=0,nr=p,nc=q)
    a = 0; sig2 = .5
    tX =t(X); h_tX = h*tX ; 
    entri.mala = list()
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
      entri.mala[[s]] = M
    }
    print(ac[ss] <- a/Iters )
    mala_lam1$mse[ss] = mean((X%*%Bm - X%*%B )^2)
    aa = 0
    for (i in 1:(p*q))aa = aa+data.table::between(B[i], quantile(sapply(entri.mala, function(x) x[i]),probs = 0.025 ),quantile(sapply(entri.mala, function(x) x[i]),probs = 0.975))
    mala_lam1$credible[ss] = aa/(p*q)
    
    #Langevin MC for BRRR
    M = hatc #
    h =  h/1.5
    M.lmc = matrix(data=0,nr=p,nc=q)
    entri.lmc = list()
    for(s in 1:Iters){
      tam = solve(a = ystar + tcrossprod(M),M)
      M = M + eigenMapMatMult(h_tX,P_Omega(Yna - X%*%M,imiss) )/sig2 + h*(p+q+2)*tam + sqrt(2*h)*rnorm(p*q)
      if (s>burnin){
        M.lmc = M.lmc + M/(Iters-burnin)
      } 
      entri.lmc[[s]] = M
    }
    lmc_lam1$mse[ss] = mean((X%*%M.lmc -X%*%B)^2)
    aa = 0
    for (i in 1:(p*q))aa = aa+data.table::between(B[i], quantile(sapply(entri.lmc, function(x) x[i]),probs = 0.025 ),quantile(sapply(entri.lmc, function(x) x[i]),probs = 0.975))
    lmc_lam1$credible[ss] = aa/(p*q)
    ######################################
    #####  MALA    ############
    Iters = 10000
    burnin = 2000
    h =  1/(p*q)^1.14
    M =  hatc
    tau = 10
    ystar = diag(p)*tau^2
    Bm = matrix(data=0,nr=p,nc=q)
    a = 0; sig2 = 1
    tX =t(X); h_tX = h*tX ; 
    entri.mala = list()
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
      entri.mala[[s]] = M
    }
    print(ac2[ss] <- a/Iters )
    mala_lam2$mse[ss] = mean((X%*%Bm - X%*%B )^2)
    aa = 0
    for (i in 1:(p*q))aa = aa+data.table::between(B[i], quantile(sapply(entri.mala, function(x) x[i]),probs = 0.025 ),quantile(sapply(entri.mala, function(x) x[i]),probs = 0.975))
    mala_lam2$credible[ss] = aa/(p*q)
    
    #Langevin MC for BRRR
    M = hatc #
    h =  h/1.5
    M.lmc = matrix(data=0,nr=p,nc=q)
    entri.lmc = list()
    for(s in 1:Iters){
      tam = solve(a = ystar + tcrossprod(M),M)
      M = M + eigenMapMatMult(h_tX,P_Omega(Yna - X%*%M,imiss) )/sig2 + h*(p+q+2)*tam + sqrt(2*h)*rnorm(p*q)
      if (s>burnin){
        M.lmc = M.lmc + M/(Iters-burnin)
      } 
      entri.lmc[[s]] = M
    }
    lmc_lam2$mse[ss] = mean((X%*%M.lmc -X%*%B)^2)
    aa = 0
    for (i in 1:(p*q))aa = aa+data.table::between(B[i], quantile(sapply(entri.lmc, function(x) x[i]),probs = 0.025 ),quantile(sapply(entri.lmc, function(x) x[i]),probs = 0.975))
    lmc_lam2$credible[ss] = aa/(p*q)
    
    ######################################
    #####  MALA    ############
    Iters = 10000
    burnin = 2000
    h =  1/(p*q)^.85
    M =  hatc
    tau = 10
    ystar = diag(p)*tau^2
    Bm = matrix(data=0,nr=p,nc=q)
    a = 0; sig2 = 2^2
    tX =t(X); h_tX = h*tX ; 
    entri.mala = list()
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
      entri.mala[[s]] = M
    }
    print(ac8[ss] <- a/Iters )
    mala_lam8$mse[ss] = mean((X%*%Bm - X%*%B )^2)
    aa = 0
    for (i in 1:(p*q))aa = aa+data.table::between(B[i], quantile(sapply(entri.mala, function(x) x[i]),probs = 0.025 ),quantile(sapply(entri.mala, function(x) x[i]),probs = 0.975))
    mala_lam8$credible[ss] = aa/(p*q)
    
    #Langevin MC for BRRR
    M = hatc #
    h =  h/1.5
    M.lmc = matrix(data=0,nr=p,nc=q)
    entri.lmc = list()
    for(s in 1:Iters){
      tam = solve(a = ystar + tcrossprod(M),M)
      M = M + eigenMapMatMult(h_tX,P_Omega(Yna - X%*%M,imiss) )/sig2 + h*(p+q+2)*tam + sqrt(2*h)*rnorm(p*q)
      if (s>burnin){
        M.lmc = M.lmc + M/(Iters-burnin)
      } 
      entri.lmc[[s]] = M
    }
    lmc_lam8$mse[ss] = mean((X%*%M.lmc -X%*%B)^2)
    aa = 0
    for (i in 1:(p*q))aa = aa+data.table::between(B[i], quantile(sapply(entri.lmc, function(x) x[i]),probs = 0.025 ),quantile(sapply(entri.lmc, function(x) x[i]),probs = 0.975))
    lmc_lam8$credible[ss] = aa/(p*q)
    
    ######################################
    #####  MALA    ############
    Iters = 10000
    burnin = 2000
    h =  1/(p*q)^.65
    M =  hatc
    tau = 10
    ystar = diag(p)*tau^2
    Bm = matrix(data=0,nr=p,nc=q)
    a = 0; sig2 = 2^4
    tX =t(X); h_tX = h*tX ; 
    entri.mala = list()
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
      entri.mala[[s]] = M
    }
    print(ac32[ss] <- a/Iters )
    mala_lam32$mse[ss] = mean((X%*%Bm - X%*%B )^2)
    aa = 0
    for (i in 1:(p*q))aa = aa+data.table::between(B[i], quantile(sapply(entri.mala, function(x) x[i]),probs = 0.025 ),quantile(sapply(entri.mala, function(x) x[i]),probs = 0.975))
    mala_lam32$credible[ss] = aa/(p*q)
    
    #Langevin MC for BRRR
    M = hatc #
    h =  h/1.5
    M.lmc = matrix(data=0,nr=p,nc=q)
    entri.lmc = list()
    for(s in 1:Iters){
      tam = solve(a = ystar + tcrossprod(M),M)
      M = M + eigenMapMatMult(h_tX,P_Omega(Yna - X%*%M,imiss) )/sig2 + h*(p+q+2)*tam + sqrt(2*h)*rnorm(p*q)
      if (s>burnin){
        M.lmc = M.lmc + M/(Iters-burnin)
      } 
      entri.lmc[[s]] = M
    }
    lmc_lam32$mse[ss] = mean((X%*%M.lmc -X%*%B)^2)
    aa = 0
    for (i in 1:(p*q))aa = aa+data.table::between(B[i], quantile(sapply(entri.lmc, function(x) x[i]),probs = 0.025 ),quantile(sapply(entri.lmc, function(x) x[i]),probs = 0.975))
    lmc_lam32$credible[ss] = aa/(p*q)
    
    print(ss)
  }
  
  
  save.image(file ='~/Dropbox/ongoing_works/B RRR/BRRR with missing Y/Rcodes/outsimulations/simu_change_lambda_pcbayes.rda')
  
  
  ######################################
  
  

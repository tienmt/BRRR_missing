library(rrr) # to get data
setwd("~/Dropbox/ongoing_works/B RRR/BRRR with missing Y/Rcodes")
library(Rcpp)
sourceCpp("~/Documents/testtest/test.cpp")
P_Omega = function(a,entri){ a[entri] = 0; return(a) }
P_Omega <- compiler::cmpfun(P_Omega)
library(rrpack)
control = list(epsilon = 1e-4, sv.tol = 1e-2, maxit = 1000,
               trace = FALSE, gammaC0 = 1.1, plot.cv = F, conv.obj = TRUE)
library(secure) # to get data
x = scale( as.matrix(CellCycle$X) )
tX =t(x)
dimnames(x) <- NULL
y <- scale( as.matrix(CellCycle$Y) ) ; sum(is.na(y))/prod(dim(y))*100  # missing rate in Y
n <- nrow(y)
p <- ncol(x)
q <-  ncol(y)

mrrr = mala = lmc = ac = c()
for (ss in 2:100) {
  test = sample(1:n, size = n*0.2)
  X = x[-test,]
  Yna = y[-test,]
  imiss = is.na(Yna)

  fit.cv.mrrr <- cv.mrrr(Yna, X, family = list(gaussian()),control = control,
                         penstr = list(penaltySVD = "rankCon",lambdaSVD = c(1:min(p,q)) ),
                         maxrank = min(p,q) , nfold = 5)$fit
  hatc = fit.cv.mrrr$coef[- 1 ,]
  mrrr[ss] = mean((y[test,]- x[test,]%*%hatc)^2 ,na.rm = TRUE )
  

  ######################################
  #####  MALA  full  ############
  Iters = 30000
  burnin = 5000
  h =  1/(p*q)^1.36 #hatc#
  M =  hatc
  lam = 10
  ystar = diag(p)*lam^2
  Bm = matrix(data=0,nr=p,nc=q)
  a = 0; sig2 = 1
  tX =t(X);h_tX = h*tX 
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
  print(paste(ss, ac[ss] <- a/Iters))
  mala[ss] = mean((y[test,]- x[test,]%*%Bm)^2 ,na.rm = TRUE )
  
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

  lmc[ss] = mean((y[test,]- x[test,]%*%M.lmc)^2 ,na.rm = TRUE ) 
}

save.image(file = '~/Dropbox/ongoing_works/B RRR/BRRR with missing Y/Rcodes/outsimulations/cellcycle_out_1.rda')


mean(lmc); sd(lmc)
mean(mala); sd(mala)
mean(mrrr); sd(mrrr)
par(mar = c(2, 2, .7, 0.1))
boxplot(cbind('LMC' = lmc,
              'MALA' = mala,
              'mRRR' = mrrr))
grid(lwd=1.5)



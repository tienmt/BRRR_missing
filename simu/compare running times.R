library(rrpack)
library(Rcpp)
sourceCpp("~/Documents/testtest/test.cpp")
P_Omega = function(a,entri){ a[entri] = 0; return(a) }
P_Omega <- compiler::cmpfun(P_Omega)

n = 2500  # samples
q=l = 10  # response
p = 10  # predictors
r = 3   # true rank
missfrac = 0.5
rho.x = 0.0
S = matrix(rho.x, ncol = p, nrow = p); diag(S) = 1
out = eigen(S, symmetric = TRUE)
S.sqrt = out$vectors %*% diag(out$values^0.5) %*% t(out$vectors)
control = list(epsilon = 1e-4, sv.tol = 1e-2, maxit = 1000,
               trace = FALSE, gammaC0 = 1.1, plot.cv = F, conv.obj = TRUE)

X = matrix(rnorm(n*p),nr=n) %*% S.sqrt
B = matrix(rnorm(p*r),nr=p)%*%t(matrix(rnorm(l*r),nr=l))  #*2 + rnorm(p*l,sd = .1)
Yna = Y = X%*%B + rt(n*l,df=2) # rnorm(n*l) # 
imiss = sample(seq(n*l),n*l*missfrac,replace=FALSE)
Yna[imiss] = NA

start_time <- Sys.time()
fit.cv.mrrr <- cv.mrrr(Yna, X, family = list(gaussian()),control = control,
                       penstr = list(penaltySVD = "rankCon",lambdaSVD = c(1:min(p,q)) ),
                       maxrank = min(p,q) , nfold = 5)$fit
hatc = fit.cv.mrrr$coef[- 1 ,]
end_time <- Sys.time();end_time - start_time
######################################
#####  MALA    ############
Iters = 5000
burnin = 0
h =  1/(p*q)^1.7
M =  hatc
tau = 1
ystar = diag(p)*tau^2
Bm = matrix(data=0,nr=p,nc=q)
a = 0; sig2 = 2
tX =t(X); h_tX = h*tX ; 
start_time <- Sys.time()
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
end_time <- Sys.time();end_time - start_time
a/Iters


#Langevin MC for BRRR
M = hatc #
h =  h/2
M.lmc = matrix(data=0,nr=p,nc=q)
start_time <- Sys.time()
for(s in 1:Iters){
  tam = solve(a = ystar + tcrossprod(M),M)
  M = M + eigenMapMatMult(h_tX,P_Omega(Yna - X%*%M,imiss) )/sig2 + h*(p+q+2)*tam + sqrt(2*h)*rnorm(p*q)
  if (s>burnin){
    M.lmc = M.lmc + M/(Iters-burnin)
  } 
}
end_time <- Sys.time();end_time - start_time
c(mean((Bm - B )^2),Matrix::rankMatrix(Bm,.1)[1],  
  mean((Bm - B )^2)/mean(B^2) , mean(( (X%*%Bm)[imiss] - Y[imiss])^2))
c(mean((hatc - B)^2), Matrix::rankMatrix(hatc,.1)[1],
  mean((hatc - B )^2)/mean(B^2),mean(( (X%*%hatc)[imiss] - Y[imiss])^2) )
c(mean((M.lmc -B)^2), Matrix::rankMatrix(M.lmc,.1)[1],
  mean((M.lmc -B )^2)/mean(B^2), mean(( (X%*%M.lmc)[imiss] - Y[imiss])^2) )
a/Iters


#fixed n,p and change q
lmc = c(0.1415231, 0.4033508,1.596942)
mala = c(0.5353441,1.568539, 4.909745)
mrrr = c(0.6074359, 1.537456,6.328265 )
par(mar = c(2.5,5,3,1) , mfrow = c(1,3) )
plot(log(lmc),type = 'b',ylim = c(-2,2),lty=3, lwd = 2,
     xaxt = "n",xlab = '',ylab = 'time (s) in log-scale ',
     main = '(a)',cex = 3,cex.lab= 2,cex.axis=1.5)
grid(lwd = 2)
axis(1, at=1:3,labels=c('q = 10','q = 50','q = 250'),cex.axis=1.5)
lines(log(mala),lty=3,lwd = 2,type = 'b',pch=2,cex = 3)
lines(log(mrrr),lty=3,lwd = 2,type = 'b',pch=3,cex = 3)

# fixed p,q change n
lmc = c(0.1808319,0.3324721,1.072639 )
mala = c(0.4276309,1.306748,3.526047 )
mrrr = c(0.4568031,0.5426159,2.112822 )
plot(log(lmc),type = 'b',ylim = c(-2,2),lty=3, lwd = 2,
     xaxt = "n",xlab = '',ylab = '',
     main = '(b)',cex = 3,cex.axis=1.5)
grid(lwd = 2)
axis(1, at=1:3,labels=c('n = 100','n = 500','n = 2500'),cex.axis=1.5)
lines(log(mala),lty=3,lwd = 2,type = 'b',pch=2,cex = 3)
lines(log(mrrr),lty=3,lwd = 2,type = 'b',pch=3,cex = 3)


#fixed n,q and change p
lmc = c(0.1773829,0.6151769,15.56312)
mala = c(0.4482291,1.517594,51.1786)
mrrr = c(0.4392419,0.9671741,1.611686 )
plot(log(lmc),type = 'b',ylim = c(-2,4),lty=3, lwd = 2,
     xaxt = "n",xlab = '',ylab = '',
     main = '(c)',cex = 3,cex.axis=1.5)
legend("topleft",inset=.12,
           x.intersp = 0,
           c(" LMC "," MALA",' mRRR'),lwd=c(1,1,1),
           pch = c(1,2,3), horiz=F, bty='n', cex=2)
grid(lwd = 2)
axis(1, at=1:3,labels=c('p = 10','p = 50','p = 250'),cex.axis=1.5)
lines(log(mala),lty=3,lwd = 2,type = 'b',pch=2,cex = 3)
lines(log(mrrr),lty=3,lwd = 2,type = 'b',pch=3,cex = 3)



  
  

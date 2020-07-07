### It is suggested to have X standardized before model fitting
library(RSpectra)

net.gen.from.P <- function(P,mode="undirected"){
  n <- nrow(P)
  if(mode=="undirected"){
    upper.index <- which(upper.tri(P))
    upper.p <- P[upper.index]
    upper.u <- runif(n=length(upper.p))
    upper.A <- rep(0,length(upper.p))
    upper.A[upper.u < upper.p] <- 1
    
    
    A <- matrix(0,n,n)
    A[upper.index] <- upper.A
    A <- A + t(A)
    
  }else{
    A <- matrix(0,n,n)
    R <- matrix(runif(n^2),n,n)
    A[R<P] <- 1
    diag(A) <- 0
  }
  return(A)
}


SP.Inf <- function(X,Y,A,K,r=NULL,sigma2=NULL,thr=NULL,alpha.CI=0.05,boot.thr=TRUE,boot.n= 50){
  n <- nrow(X)
  p <- ncol(X)
  eigen.A <- eigs_sym(A=A,k=K)
  X.svd <- svd(X)
  x.proj <- X.svd$v%*%(t(X.svd$u)/X.svd$d)
  R <- X.svd$u
  Uhat <- matrix(eigen.A$vectors[,1:K],ncol=K)
  hat.SVD <- svd(t(R)%*%Uhat,nv=K,nu=p)
  Vhat <- hat.SVD$v
  U.curl <- Uhat%*%Vhat
  R.curl <- R%*%hat.SVD$u
  if(is.null(r)){
  if(is.null(thr)){
    if(boot.thr){
      P0 <- eigen.A$vectors %*%t(eigen.A$vectors*eigen.A$values)
      P0[P0>1] <- 1
      P0[P0<0] <- 0
      P0 <- P0/mean(rowSums(P0))*mean(rowSums(A))
      eig.P0 <- eigs_sym(P0,k=K)
      fake.svd <- svd(t(R)%*%matrix(eig.P0$vectors[,1:K],ncol=K))
      print("Use bootstrapping to find r-threshold....")
      sim.s <- matrix(0,boot.n,min(p,K))
      for(II in 1:boot.n){
        A.sim <- net.gen.from.P(P0)
        eig.A.sim <- eigs_sym(A.sim,k=K)
        sim.svd <- svd(t(R)%*%matrix(eig.A.sim$vectors[,1:K],ncol=K))
        sim.s[II,] <- sim.svd$d
      }
      
      thr <- 1-max(abs(t(sim.s)-fake.svd$d))
      print(paste("Select r by threshold",thr))
      
    }else{
      dhat <- max(rowSums(A))
      thr <- 1- 4*sqrt(p*K*log(n))/dhat
      if(thr < 0.05){
        thr <- 0.05
      }
      print(paste("Select r by asymptotic threshold",thr))
    }
  }
    print("Observed singular values are ")
    print(hat.SVD$d)
    r <- sum(hat.SVD$d>=thr)
    print(paste("Select r =",r))
  }
  if(r >0){
  Z <- cbind(R.curl[,-(1:r)],U.curl[,-(1:r)])
  Z.proj <- Z%*%solve(t(Z)%*%Z,t(Z))
  PW.hat <- matrix(R.curl[,1:r],ncol=r)%*%t(matrix(R.curl[,1:r],ncol=r))
  S <- cbind(R.curl[,-(1:r)],matrix(0,nrow=n,ncol=K-r))%*%solve(t(Z)%*%Z,t(Z))
  S.perp <- Z.proj - S
  H <- Z.proj + PW.hat
  theta.hat <- x.proj%*%PW.hat%*%Y
  PU.hat <- t(U.curl[,-(1:r)])
  }else{
    Z <- cbind(R.curl,U.curl)
    Z.proj <- Z%*%solve(t(Z)%*%Z,t(Z))
    S <- cbind(R.curl,matrix(0,nrow=n,ncol=K))%*%solve(t(Z)%*%Z,t(Z))
    S.perp <- Z.proj - S
    H <- Z.proj
    theta.hat <- matrix(0,nrow=p,ncol=1)
    PU.hat <- t(U.curl)
  }
  beta.z <- S%*%Y
  beta.hat <- x.proj%*%beta.z
  sigma2.hat <- sum((Y-H%*%Y)^2)/(n-p-K+r)
  if(!is.null(sigma2)){
      sigma2.hat <- sigma2
  }
  tmp <- x.proj%*%S
  cov.hat <- sigma2.hat*tmp%*%t(tmp)
  fitted.Y <- H%*%Y
  alpha.hat <- fitted.Y - X%*%beta.hat - X%*%theta.hat
  diag.sd <- sqrt(diag(cov.hat))
  abs.t.val <- abs(beta.hat/diag.sd)
  t.pval <- 1-pnorm(abs.t.val)
  CI.lower <- beta.hat - qnorm(1-alpha.CI/2)*diag.sd
  CI.upper <- beta.hat + qnorm(1-alpha.CI/2)*diag.sd
  coef.mat <- cbind(beta.hat,CI.lower,CI.upper,t.pval)
  colnames(coef.mat) <- c("coefficient","CI-lower","CI-upper","p-val")
  cov.gamma <- sigma2.hat*PU.hat%*%S.perp%*%t(S.perp)%*%t(PU.hat)
  gamma <- PU.hat%*%S.perp%*%Y
  if((K-r)>1){
  eig.cov.gamma <- eigen(cov.gamma,symmetric=TRUE)
  inv.sqrt.cov.gamma <- t(eig.cov.gamma$vectors)*sqrt(1/eig.cov.gamma$values)
  adjusted.gamma <- inv.sqrt.cov.gamma%*%gamma
  }else{
    adjusted.gamma <- 0
  }
  chisq.val <- sum(adjusted.gamma^2)
  chisq.pval <- pchisq(chisq.val,df=K-r,lower.tail=FALSE)
  if(is.null(colnames(X))){
    rownames(coef.mat) <- paste("V",1:p,sep="")
  }else{
    rownames(coef.mat) <- colnames(X)
  }
  return(list(beta=beta.hat,alpha=alpha.hat,theta=theta.hat,r=r,
              sigma2=sigma2.hat,cov.hat=cov.hat,coef.mat=coef.mat,
              fitted=fitted.Y,chisq.val=chisq.val,chisq.p=chisq.pval))
}




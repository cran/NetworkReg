### It is suggested to have X standardized before model fitting
library(RSpectra)

#' @title generates a network from the given connection probability
#' @description Generates an adjacency matrix from a given probability matrix,
#' according independent Bernoulli -- the so-called inhomogeneous
#' Erdos-Renyi model. It is used to generate new networks from a given model.
#'
#' @usage net.gen.from.P(P, mode = "undirected")
#' @param P connection probability between nodes
#' @param mode "undirected" (default) if the network is undirected, so the
#' adjacency matrix will be symmetric with only upper diagonal entries
#' being generated as independent Bernoulli. Otherwise, the adjacency matrix gives independent Bernoulli everywhere.
#'
#' @return An adjacency matrix
#' @export
#'
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

SP.Inf.Linear <- function(X,Y,A,K,r=NULL,thr=NULL,alpha.CI=0.05,boot.thr=TRUE,boot.n= 50,sigma2=NULL){
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
        message("Use bootstrapping to find r-threshold....")
        sim.s <- matrix(0,boot.n,min(p,K))
        for(II in 1:boot.n){
          A.sim <- net.gen.from.P(P0)
          eig.A.sim <- eigs_sym(A.sim,k=K)
          sim.svd <- svd(t(R)%*%matrix(eig.A.sim$vectors[,1:K],ncol=K))
          sim.s[II,] <- sim.svd$d
        }

        thr <- 1-max(abs(t(sim.s)-fake.svd$d))
        message(paste("Select r by threshold",thr))

      }else{
        dhat <- max(rowSums(A))
        thr <- 1- 4*sqrt(p*K*log(n))/dhat
        if(thr < 0.05){
          thr <- 0.05
        }
        message(paste("Select r by asymptotic threshold",thr))
      }
    }
    r <- sum(hat.SVD$d>=thr)
    message(paste("Select r =",r))
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
  chisq.p <- pchisq(chisq.val,df=K-r,lower.tail=FALSE)
  if(is.null(colnames(X))){
    rownames(coef.mat) <- paste("V",1:p,sep="")
  }else{
    rownames(coef.mat) <- colnames(X)
  }
  return(list(beta=beta.hat,alpha=alpha.hat,theta=theta.hat,r=r,
              sigma2=sigma2.hat,cov.hat=cov.hat,coef.mat=coef.mat,
              fitted=fitted.Y,chisq.val=chisq.val,chisq.p=chisq.p,thr=thr))
}



SP.Inf.Logistic <- function(X,Y,A,K,r=NULL,thr=NULL,boot.thr=TRUE,alpha.CI=0.05,boot.n= 50){
  n <- nrow(X)
  p <- ncol(X)
  eigen.A <- eigs_sym(A=A,k=K,which = "LA")
  X.svd <- svd(X)
  x.proj <- X.svd$v%*%(t(X.svd$u)/X.svd$d)
  R <- X.svd$u
  Uhat <- matrix(eigen.A$vectors[,1:K],ncol=K)
  hat.SVD <- svd(t(R)%*%Uhat,nv=K,nu=p)
  Vhat <- hat.SVD$v
  U.curl <- sqrt(n)*Uhat%*%Vhat
  R.curl <- sqrt(n)*R%*%hat.SVD$u
  if(is.null(r)){
    if(is.null(thr)){
      if(boot.thr){
        P0 <- eigen.A$vectors %*%t(eigen.A$vectors*eigen.A$values)
        P0[P0>1] <- 1
        P0[P0<0] <- 0
        P0 <- P0/mean(rowSums(P0))*mean(rowSums(A))
        eig.P0 <- eigs_sym(P0,k=K)
        fake.svd <- svd(t(R)%*%matrix(eig.P0$vectors[,1:K],ncol=K))
        message("Use bootstrapping to find r-threshold....")
        sim.s <- matrix(0,boot.n,min(p,K))
        for(II in 1:boot.n){
          A.sim <- net.gen.from.P(P0)
          eig.A.sim <- eigs_sym(A.sim,k=K)
          sim.svd <- svd(t(R)%*%matrix(eig.A.sim$vectors[,1:K],ncol=K))
          sim.s[II,] <- sim.svd$d
        }

        thr <- 1-max(abs(t(sim.s)-fake.svd$d))
        message(paste("Select r by threshold",thr))

      }else{
        dhat <- max(rowSums(A))
        thr <- 1- 4*sqrt(p*K*log(n))/dhat
        if(thr < 0.05){
          thr <- 0.05
        }
        message(paste("Select r by asymptotic threshold",thr))
      }
    }
    r <- sum(hat.SVD$d>=thr)
    message(paste("Select r =",r))
  }
  if(r >0){
    Z <- cbind(R.curl,U.curl[,-(1:r)])
    fit<-glm(Y~0+Z,family = "binomial")
    theta<-as.matrix(fit$coefficients[1:r])
    beta<-as.matrix(fit$coefficients[(r+1):p])
    VCOV<-vcov(fit)
    cov_theta<-VCOV[1:r,1:r]
    cov_beta<-VCOV[(r+1):p,(r+1):p]
    Xtheta<-R.curl[,(1:r)]%*%theta
    Xbeta<-R.curl[,(r+1):p]%*%beta
    theta.hat<-x.proj%*%Xtheta
    cov_theta<-x.proj%*%R.curl[,(1:r)]%*%cov_theta%*%t(x.proj%*%R.curl[,(1:r)])
    alpha.hat<-U.curl[,-(1:r)]%*%fit$coefficients[(p+1):(p+K-r)]
  }else{
    Z <- cbind(R.curl,U.curl)
    fit<-glm(Y~0+Z,family = "binomial")
    beta<-as.matrix(fit$coefficients[1:p])
    VCOV<-vcov(fit)
    cov_beta<-VCOV[1:p,1:p]
    Xbeta<-R.curl[,1:p]%*%beta
    theta.hat<-NA
    alpha.hat<-U.curl%*%fit$coefficients[(p+1):(p+K)]
  }
  beta.hat <- x.proj%*%Xbeta
  cov_beta<-x.proj%*%R.curl[,(r+1):p]%*%cov_beta%*%t(x.proj%*%R.curl[,(r+1):p])
  diag.sd<-sqrt(diag(cov_beta))
  CI.lower <- beta.hat - qnorm(1-alpha.CI/2)*diag.sd
  CI.upper <- beta.hat + qnorm(1-alpha.CI/2)*diag.sd
  SEE<-diag(cov_beta)
  cov_network<-VCOV[(p+1):(p+K-r),(p+1):(p+K-r)]
  network_effect_test_statistics<-t(as.matrix(fit$coefficients[(p+1):(p+K-r)]))%*%solve(cov_network)%*%as.matrix(fit$coefficients[(p+1):(p+K-r)])
  network<-(network_effect_test_statistics>qchisq(1-alpha.CI,K-r))
  chisq.p<-pchisq(network_effect_test_statistics,df=K-r,lower.tail=FALSE)
  abs.t.val <- abs(beta.hat/diag.sd)
  t.pval <- 1-pnorm(abs.t.val)
  coef.mat <- cbind(beta.hat,CI.lower,CI.upper,t.pval)
  colnames(coef.mat) <- c("coefficient","CI-lower","CI-upper","p-val")
  return(list(theta=theta.hat,alpha=alpha.hat,r=r,
              beta=beta.hat,cov.hat=cov_beta,coef.mat=coef.mat,
              chisq.val=network_effect_test_statistics,chisq.p= chisq.p,
              fitted=fit$fitted.values))
}

SP.Inf.Poisson <- function(X,Y,A,K,r=NULL,thr=NULL,boot.thr=TRUE,alpha.CI=0.05,boot.n= 50){
  n <- nrow(X)
  p <- ncol(X)
  eigen.A <- eigs_sym(A=A,k=K,which = "LA")
  X.svd <- svd(X)
  x.proj <- X.svd$v%*%(t(X.svd$u)/X.svd$d)
  R <- X.svd$u
  Uhat <- matrix(eigen.A$vectors[,1:K],ncol=K)
  hat.SVD <- svd(t(R)%*%Uhat,nv=K,nu=p)
  Vhat <- hat.SVD$v
  U.curl <- sqrt(n)*Uhat%*%Vhat
  R.curl <- sqrt(n)*R%*%hat.SVD$u
  if(is.null(r)){
    if(is.null(thr)){
      if(boot.thr){
        P0 <- eigen.A$vectors %*%t(eigen.A$vectors*eigen.A$values)
        P0[P0>1] <- 1
        P0[P0<0] <- 0
        P0 <- P0/mean(rowSums(P0))*mean(rowSums(A))
        eig.P0 <- eigs_sym(P0,k=K)
        fake.svd <- svd(t(R)%*%matrix(eig.P0$vectors[,1:K],ncol=K))
        message("Use bootstrapping to find r-threshold....")
        sim.s <- matrix(0,boot.n,min(p,K))
        for(II in 1:boot.n){
          A.sim <- net.gen.from.P(P0)
          eig.A.sim <- eigs_sym(A.sim,k=K)
          sim.svd <- svd(t(R)%*%matrix(eig.A.sim$vectors[,1:K],ncol=K))
          sim.s[II,] <- sim.svd$d
        }

        thr <- 1-max(abs(t(sim.s)-fake.svd$d))
        message(paste("Select r by threshold",thr))

      }else{
        dhat <- max(rowSums(A))
        thr <- 1- 4*sqrt(p*K*log(n))/dhat
        if(thr < 0.05){
          thr <- 0.05
        }
        message(paste("Select r by asymptotic threshold",thr))
      }
    }
    r <- sum(hat.SVD$d>=thr)
    message(paste("Select r =",r))
  }
  if(r >0){
    Z <- cbind(R.curl,U.curl[,-(1:r)])
    fit<-glm(Y~0+Z,family = poisson(link = "log"))
    theta<-as.matrix(fit$coefficients[1:r])
    beta<-as.matrix(fit$coefficients[(r+1):p])
    VCOV<-vcov(fit)
    cov_theta<-VCOV[1:r,1:r]
    cov_beta<-VCOV[(r+1):p,(r+1):p]
    Xtheta<-R.curl[,(1:r)]%*%theta
    Xbeta<-R.curl[,(r+1):p]%*%beta
    theta.hat<-x.proj%*%Xtheta
    cov_theta<-x.proj%*%R.curl[,(1:r)]%*%cov_theta%*%t(x.proj%*%R.curl[,(1:r)])
    alpha.hat<-U.curl[,-(1:r)]%*%fit$coefficients[(p+1):(p+K-r)]
  }else{
    Z <- cbind(R.curl,U.curl)
    fit<-glm(Y~0+Z,family = poisson(link = "log"))
    beta<-as.matrix(fit$coefficients[1:p])
    VCOV<-vcov(fit)
    cov_beta<-VCOV[1:p,1:p]
    Xbeta<-R.curl[,1:p]%*%beta
    theta.hat<-NA
    alpha.hat<-U.curl%*%fit$coefficients[(p+1):(p+K)]
  }
  beta.hat <- x.proj%*%Xbeta
  cov_beta<-x.proj%*%R.curl[,(r+1):p]%*%cov_beta%*%t(x.proj%*%R.curl[,(r+1):p])
  diag.sd<-sqrt(diag(cov_beta))
  CI.lower <- beta.hat - qnorm(1-alpha.CI/2)*diag.sd
  CI.upper <- beta.hat + qnorm(1-alpha.CI/2)*diag.sd
  SEE<-diag(cov_beta)
  cov_network<-VCOV[(p+1):(p+K-r),(p+1):(p+K-r)]
  network_effect_test_statistics<-t(as.matrix(fit$coefficients[(p+1):(p+K-r)]))%*%solve(cov_network)%*%as.matrix(fit$coefficients[(p+1):(p+K-r)])
  network<-(network_effect_test_statistics>qchisq(1-alpha.CI,K-r))
  chisq.p<-pchisq(network_effect_test_statistics,df=K-r,lower.tail=FALSE)
  abs.t.val <- abs(beta.hat/diag.sd)
  t.pval <- 1-pnorm(abs.t.val)
  coef.mat <- cbind(beta.hat,CI.lower,CI.upper,t.pval)
  colnames(coef.mat) <- c("coefficient","CI-lower","CI-upper","p-val")
  return(list(theta=theta.hat,alpha=alpha.hat,r=r,
              beta=beta.hat,cov.hat=cov_beta,coef.mat=coef.mat,
              chisq.val=network_effect_test_statistics,chisq.p= chisq.p,
              fitted=fit$fitted.values))
}

#' @title Fitting Generalized Linear Model on Network-Linked Data
#' @description SP.Inf is used to the regression model on network-linked data by subspace project and produce the inference result.
#' @usage SP.Inf(X, Y, A, K, model = "linear", r = NULL, sigma2 = NULL, thr = NULL,
#' alpha.CI = 0.05, boot.thr = TRUE, boot.n = 50)
#' @param X the covariate matrix where each row is an observation and each column is a covariate. If an intercept is to be included in the model, the column of ones should be in the matrix.
#' @param Y the column vector of response.
#' @param A the network information. The most natural choice is the adjacency matrix of the network. However, if the network is assumed to be noisy and a better estimate of the structural connection strength, it can also be used. This corresponds to the Phat matrix in the original paper. A Laplacian matrix can also be used, but it should be flipped. See 'Details'.
#' @param K the dimension of the network eigenspace for network effect.
#' @param model the type of Generalized Linear Regression. The "linear" ,"logistic" and "poisson" represents linear regression, logistic regression and poisson regression. The default is linear regression.
#' @param r the covariate-network cofounding space dimension. This is typically unknown and can be unspecified by using the default value 'NULL'. If so, the user should provide a threshold or resort to a tuning procedure by either the theoretical rule or a bootstrapping method, as described in the paper.
#' @param sigma2 the variance of random noise for linear regression. Typically unknown.
#' @param thr threshold for r estimation. If r is unspecified, we will use the thereshold to select r. If this is also 'NULL', aa theoretical threshold or a bootsrapping method can be evoked to estimate it.
#' @param alpha.CI the 1-alpha.CI confidence level will be produced for the parameters.
#' @param boot.thr logical. Only effective if both r and thr are NULLs. If FALSE, the theoretical threshold will be used to select r. Otherwise, the bootstrapping procedure will be used to find the threshold.
#' @param boot.n the number of bootstrapping samples used when boot.thr is TRUE.
#' @details The model fitting procedure is following the paper exactly, so please check the procedure and theory in the paper. If the Laplacian matrix L=D-A is the network quantity to use, notice that typically we treat the smallest values and their corresponding eigenvectors as network cohesive space. Therefore, one should consider flip the Laplacian matrix by using cI - L as the value for A, where c is sufficiently large to ensure PSD of cI-L.
#'
#' @return A list object with
#' \item{\code{beta}}{estimate of beta, the covariate effects}
#' \item{\code{alpha}}{individual effects}
#' \item{\code{theta}}{coefficients of confounding effects with respect to the covariates}
#' \item{\code{r}}{confounding dimension}
#' \item{\code{sigma}}{estimated random noise variance for linear regression}
#' \item{\code{cov.hat}}{covariance matrix of beta}
#' \item{\code{coef.mat}}{beta and the confidence intervals according to alpha.CI and the p-values of the significance test}
#' \item{\code{fitted}}{fitted value of response}
#' \item{\code{chisq.val}}{the value of the chi-square statistic for the significance test for network effect}
#' \item{\code{chisq.p}}{the p-value of the significance test for network effect}
#' @import RSpectra randnet stats
#' @references Le, C. M., & Li, T. (2022). Linear regression and its inference on noisy network-linked data. Journal of the Royal Statistical Society Series B: Statistical Methodology, 84(5), 1851-1885.
#'
#' Wang J, Le C M, Li T. Perturbation-Robust Predictive Modeling of Social Effects by Network Subspace Generalized Linear Models. arXiv preprint arXiv:2410.01163, 2024.
#'
#' @author
#' Jianxiang Wang, Can M. Le, and Tianxi Li.
#' Maintainer: Jianxiang Wang <jw1881@scarletmail.rutgers.edu>
#'
#' @examples
#' set.seed(1)
#' library(randnet)
#' library(RSpectra)
#' ### Example data generation procedure from Section 5.3 of the paper with logistic regression
#' n <- 1000
#' big.model <- BlockModel.Gen(lambda=n^(2/3), n=n, beta=0.2, K=4)
#' P <- big.model$P
#' big.X <- cbind(rnorm(n), runif(n), rexp(n))
#' eigen.P <- eigs_sym(A=P, k=4)
#' X.true <- big.X
#' X.true <- scale(X.true, center=TRUE, scale=TRUE) * sqrt(n/(n-1))
#' beta <- matrix(c(1,1,1), ncol=1)
#' Xbeta <- X.true %*% beta
#' U <- eigen.P$vectors[,1:4]
#' alpha.coef <- matrix(sqrt(n) * c(1, 1, 1, 1), ncol=1)
#' alpha <- U %*% alpha.coef
#' EY <- (1 + exp(-Xbeta - alpha))^(-1)
#' ## Model fitting
#' A <- net.gen.from.P(P)
#' Y <- rbinom(n, 1, EY)
#' fit <- SP.Inf(X.true, Y, A, K=4, model=c("logistic"), alpha=0.05, boot.thr=FALSE)
#' fit$coef.mat
#'
#' @keywords models
#' @keywords regression
#' @export
#'
SP.Inf <- function(X,Y,A,K,model="linear",r=NULL,sigma2=NULL,thr=NULL,alpha.CI=0.05,boot.thr=TRUE,boot.n= 50){
  if(length(model)>1) model <- model[1]
  if(!any(c("linear","logistic","poisson")==model)) stop("model must be one of linear, logistic and poisson.")
  if(model=="linear"){
    result <- SP.Inf.Linear(X,Y,A,K,r,sigma2=sigma2,thr=thr,alpha.CI=alpha.CI,boot.thr=boot.thr,boot.n= boot.n)
    return(result)
  }
  if(model=="logistic"){
    result <- SP.Inf.Logistic(X,Y,A,K,r,thr=thr,alpha.CI=alpha.CI,boot.thr=boot.thr,boot.n= boot.n)
    return(result)
  }
  if(model=="poisson"){
    result <- SP.Inf.Poisson(X,Y,A,K,r,thr=thr,alpha.CI=alpha.CI,boot.thr=boot.thr,boot.n= boot.n)
    return(result)
  }
}




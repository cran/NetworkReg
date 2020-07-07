\name{SP.Inf}
\alias{SP.Inf}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Fitting Linear Regression Models on Network-Linked Data
}
\description{
SP.Inf is used to the regression model on network-linked data by subspace project and produce the inference result.
}
\usage{
SP.Inf(X, Y, A, K, r = NULL, sigma2 = NULL, thr = NULL, alpha.CI = 0.05,
boot.thr = TRUE, boot.n = 50)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{X}{
the covariate matrix where each row is an observation and each column is a covariate. If an intercept is to be included in the model, the column of ones should be in the matrix.
}
  \item{Y}{
the column vector of response.
}
  \item{A}{
the network information. The most natural choice is the adjacency matrix of the network. However, if the network is assumed to be noisy and a better estimate of the structural connection strength, it can also be used. This corresponds to the Phat matrix in the original paper. A Laplacian matrix can also be used, but it should be flipped. See 'Details'. 
}
  \item{K}{
the dimension of the network eigenspace for network effect. 
}
  \item{r}{
the covariate-network cofounding space dimension. This is typically unknown and can be unspecified by using the default value 'NULL'. If so, the user should provide a threshold or resort to a tuning procedure by either the theoretical rule or a bootstrapping method, as described in the paper.
}
  \item{sigma2}{
the variance of random noise. Typically unknown.
}
  \item{thr}{
threshold for r estimation. If r is unspecified, we will use the thereshold to select r. If this is also 'NULL', aa theoretical threshold or a bootsrapping method can be evoked to estimate it.
}
  \item{alpha.CI}{
the 1-alpha.CI confidence level will be produced for the parameters.
}
  \item{boot.thr}{
logical. Only effective if both r and thr are NULLs. If FALSE, the theoretical threshold will be used to select r. Otherwise, the bootstrapping procedure will be used to find the threshold.
}
  \item{boot.n}{
the number of bootstrapping samples used when boot.thr is TRUE.
}
}
\details{
The model fitting procedure is following the paper exactly, so please check the procedure and theory in the paper. If the Laplacian matrix L=D-A is the network quantity to use, notice that typically we treat the smallest values and their corresponding eigenvectors as network cohesive space. Therefore, one should consider flip the Laplacian matrix by using cI - L as the value for A, where c is sufficiently large to ensure PSD of cI-L.
}
\value{
A list object with
%%  If it is a LIST, use
  \item{beta }{estimate of beta, the covariate effects}
  \item{alpha }{individual effects}
  \item{theta}{coefficients of confounding effects with respect to the covariates}
  \item{r}{confounding dimension}
  \item{sigma}{estimated random noise variance}
  \item{cov.hat}{covariance matrix of beta}
  \item{coef.mat}{beta and the confidence intervals according to alpha.CI and the p-values of the significance test}
  \item{fitted}{fitted value of response}
  \item{chisq.val}{the value of the chi-square statistic for the significance test for network effect}
  \item{chisq.p}{the p-value of the significance test for network effect}
%% ...
}
\references{
Can M. Le and Tianxi Li.  Linear regression and its inference on noisy network-linked data. arXiv:2007.00803
}
\author{
Can M. Le and  Tianxi Li.\cr

Maintainer: Tianxi Li <tianxili@virginia.edu>
}

\examples{
library(randnet)
library(RSpectra)
### data generating procedure in Section 5.3 of the paper

n <- 1000
big.model <- BlockModel.Gen(lambda=n^(1/2),n=n,beta=0.2,K=4)
P <- big.model$P
big.X <- cbind(rnorm(n),runif(n),rexp(n))

eigen.P <- eigs_sym(A=P,k=4)
X.true <- big.X
X.true <- scale(X.true,center=TRUE,scale=TRUE)*sqrt(n/(n-1))
X.true <- cbind(sqrt(n)*eigen.P$vectors[,1],X.true)
X.svd <- svd(X.true)
x.proj <- X.svd$v\%*\%(t(X.svd$u)/X.svd$d)
Theta <- X.svd$v\%*\%(t(X.svd$v)/(X.svd$d^2))*n
R <- X.svd$u
U <- eigen.P$vectors[,1:4]
true.SVD <- svd(t(R)\%*\%U,nu=4,nv=4)
V <- true.SVD$v
r <- 1
U.tilde <- U\%*\%V
R.tilde <- R\%*\%true.SVD$u
theta.tilde <- matrix(c(sqrt(n),0,0,0),ncol=1)
beta.tilde <- matrix(sqrt(n)*c(0,1,1,1),ncol=1)
Xtheta <- R.tilde\%*\%theta.tilde
Xbeta <- R.tilde\%*\%beta.tilde

theta <- solve(t(X.true)\%*\%X.true,t(X.true)\%*\%Xtheta)
beta <- solve(t(X.true)\%*\%X.true,t(X.true)\%*\%Xbeta)
alpha.coef <- matrix(sqrt(n)*c(0,1,1,1),ncol=1)
alpha <- U.tilde\%*\%alpha.coef

EY <- Xtheta+Xbeta + alpha


#### model fitting


A <- net.gen.from.P(P)
Khat <- BHMC.estimate(A, K.max = 15)$K  ### estimate K to use

## model fitting
Y <- EY + rnorm(n)
fit <- SP.Inf(X.true,Y,A,K=Khat,alpha=0.05,boot.thr=FALSE)
### In general, boot.thr = T works better for small sample but is slower.
### It was used in the paper.
fit$coef.mat  
### notice that beta1 inference is meaningful here. Check the paper.
beta
fit$chisq.p

## find a parametric estimation of the network. This is generally not available.
rsc <- reg.SP(A,K=Khat,tau=0.1)
est <- SBM.estimate(A,rsc$cluster)
Phat <- est$Phat
fit2 <- SP.Inf(X.true,Y,Phat,K=Khat,alpha=0.05,boot.thr=FALSE)
fit2$coef.mat  
### notice that beta1 inference is meaningful here. Check the paper.

}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ models }% use one of  RShowDoc("KEYWORDS")
\keyword{ regression }% __ONLY ONE__ keyword per line

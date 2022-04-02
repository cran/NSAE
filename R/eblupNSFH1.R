#' EBLUP under nonstationary Fay-Herriot model for sample area
#'
#' @description This function gives the EBLUP and the estimate of mean squared error (mse)
#'     based on a nonstationary Fay-Herriot model for sample area.
#'
#' @param formula an object of class list of formula, describe the model to be fitted
#' @param vardir a vector of sampling variances of direct estimators for each small area
#' @param lat a vector of latitude for each small area
#' @param long a vector of longitude for each small area
#' @param method type of fitting method, default is "REML" methods
#' @param data a data frame comprising the variables named in formula, vardir, lat and long
#'
#' @return The function returns a list with the following objects:
#' \describe{
#'   \item{eblup}{a vector with the values of the estimators for each small area}
#'   \item{mse}{a vector of the mean squared error estimates for each small area}
#'   \item{sample}{a matrix consist of area code, eblup, mse, SE and CV}
#'   \item{fit}{a list containing the following objects:}
#'   \itemize{
#'     \item estcoef : a data frame with the estimated model coefficients in the first column (beta),their asymptotic standard errors in the second column (std.error), the t statistics in  the third column (tvalue) and the p-values of the significance of each coefficient in last column (pvalue)
#'     \item refvar : estimated random effects variance
#'     \item spatialcorr : spatial correlation parameter
#'     \item randomeffect : a data frame with the values of the random effect estimators
#'   }
#'  }
#'
#' @export eblupNSFH1
#'
#' @examples
#' # Load data set
#' data(paddysample)
#' # Fit nonstationary Fay-Herriot model using sample part of paddy data
#' result <- eblupNSFH1(y ~ x1+x2, var, latitude, longitude, method="REML", data = paddysample)
#' result
eblupNSFH1 <- function (formula, vardir, lat,long, method = "REML", data) {
  namevar <- deparse(substitute(vardir))
  namelat <- deparse(substitute(lat))
  namelong <- deparse(substitute(long))
  if (!missing(data)) {
    formuladata <- model.frame(formula, na.action = na.omit,data)
    X <- model.matrix(formula, data)
    vardir <- data[, namevar]
    lat <- data[, namelat]
    long <- data[, namelong]
  }
  else {
    formuladata <- model.frame(formula, na.action = na.omit)
    X <- model.matrix(formula)
  }
  if (attr(attributes(formuladata)$terms, "response") == 1)
    textformula <- paste(formula[2], formula[1], formula[3])
  else textformula <- paste(formula[1], formula[2])
  if (length(na.action(formuladata)) > 0)
    stop("Argument formula=", textformula, " contains NA values.")
  if (any(is.na(vardir)))
    stop("Argument vardir=", namevar, " contains NA values.")
  y <- formuladata[, 1]
  m <- length(y)
  p <- dim(X)[2]
  direct <- y
  I<-diag(1,m)
  distance<-matrix(0,m,m)
  distance<-as.matrix(dist(cbind(as.vector(lat),as.vector(long))))
  W <- 1/(1+distance)
  W.sam <- W[1:m,1:m]
  logl=function(delta)  {
    area=m
    psi=matrix(c(vardir),area,1)
    Y=matrix(c(direct),area,1)
    Z.area=diag(1,area)
    lambda<-delta[1]
    sigma.u<-delta[2]
    C<-lambda*diag(1,p)
    Cov<-(X%*%C%*%t(X))*W.sam+sigma.u*Z.area%*%t(Z.area)
    V<-Cov+I*psi[,1]
    Vi<-solve(V)
    Xt=t(X)
    XVi<-Xt%*%Vi
    Q<-solve(XVi%*%X)
    P<-Vi-(Vi%*%X%*%Q%*%XVi)
    b.s<-Q%*%XVi%*%Y
    ee=eigen(V)
    -(area/2)*log(2*pi)-0.5*sum(log(ee$value))-(0.5)*log(det(t(X)%*%Vi%*%X))-(0.5)*t(Y)%*%P%*%Y
  }
  grr=function(delta) {
    lambda<-delta[1]
    sigma.u<-delta[2]
    area=m
    psi=matrix(c(vardir),area,1)
    Y=matrix(c(direct),area,1)
    Z.area=diag(1,area)
    I<-diag(1,area)
    C<-lambda*diag(1,p)
    Cov<-(X%*%C%*%t(X))*W.sam+sigma.u*Z.area%*%t(Z.area)
    V<-Cov+I*psi[,1]
    Vi<-solve(V)
    Xt=t(X)
    XVi<-Xt%*%Vi
    Q<-solve(XVi%*%X)
    P<-Vi-(Vi%*%X%*%Q%*%XVi)
    derLambda<-(X%*%t(X))*W
    derSigmau<-Z.area%*%t(Z.area)
    s<-matrix(0,2,1)
    PG<-P%*%derLambda
    PU<-P%*%derSigmau
    Pdir<-P%*%Y
    s[1,1]<-((-0.5)*sum(diag(PG)))+((0.5)*(t(Y)%*%PG%*%Pdir))
    s[2,1]<-((-0.5)*sum(diag(PU)))+((0.5)*(t(Y)%*%PU%*%Pdir))
    c(s[1,1],s[2,1])
  }
  opt=constrOptim(c(0.1,0.2),logl,grr,method="Nelder-Mead",ui=rbind(c(1,0),c(0,1)),
                  ci=c(0.001,0),control=list(fnscale=-1))
  lambda.stim.S=opt$par[1]
  sigma2.u.stim.S=opt$par[2]
  C.est<-lambda.stim.S*diag(1,p)
  Sigma.l<-kronecker(C.est,W.sam)
  z.mat = list()
  for (i in 1:p) {
    z.mat[[i]] <- diag(X[,i])
  }
  Z <- rlist::list.cbind(z.mat)
  Cov.est<-Z%*%Sigma.l%*%t(Z)+sigma2.u.stim.S*I%*%t(I)
  V<-Cov.est+I*vardir
  Vi<-solve(V)
  Q<-solve(t(X)%*%Vi%*%X)
  Beta.hat<-Q%*%t(X)%*%Vi%*%direct
  P<-Vi-Vi%*%X%*%solve(t(X)%*%Vi%*%X)%*%t(X)%*%Vi
  res<-direct-c(X%*%Beta.hat)
  Sigma.u=sigma2.u.stim.S*I
  spatial.hat=Sigma.l%*%t(Z)%*%Vi%*%res
  u.hat=Sigma.u%*%t(I)%*%Vi%*%res
  EBLUP.Mean<-X%*%Beta.hat+Z%*%spatial.hat+I%*%u.hat
  zvalue <- Beta.hat/sqrt(diag(Q))
  pvalue <- 2 * pnorm(abs(zvalue), lower.tail = FALSE)
  coef <- data.frame(beta = Beta.hat, std.error = sqrt(diag(Q)),tvalue = zvalue, pvalue)
  Sigma.w<-matrix(0,((p+1)*m),((p+1)*m))
  Sigma.w[1:(p*m),1:(p*m)]<-Sigma.l
  Sigma.w[(p*m+1):((p+1)*m),(p*m+1):((p+1)*m)]<-Sigma.u
  w.i<-cbind(Z,I)
  c.i<-X-w.i%*%Sigma.w%*%t(w.i)%*%Vi%*%X
  g1<-matrix(0,m,1)
  for (i in 1:m)  {
    g1[i,1]<-c.i[i,]%*%solve(t(X)%*%Vi%*%X)%*%cbind(c.i[i,])
  }
  g2<-matrix(0,m,1)
  for(i in 1:m) {
    g2[i,1]<-w.i[i,]%*%Sigma.w%*%(diag((p+1)*m)-t(w.i)%*%Vi%*%w.i%*%Sigma.w)%*%cbind(w.i[i,])
  }
  g3<-matrix(0,m,1)
  Ds.1<-matrix(0,((p+1)*m),((p+1)*m))
  Ds.1[1:(p*m),1:(p*m)]<-kronecker(diag(1,p),W.sam)
  Ds.1[(p*m+1):((p+1)*m),(p*m+1):((p+1)*m)]<-0
  Ds.2<-diag(c(rep(0,p*m),rep(1,m)))
  B.1<-Z%*%(kronecker(diag(1,p),W.sam))%*%t(Z)
  B.2<-I%*%t(I)
  B<-list(B.1,B.2)
  Dv.1<--Vi%*%B.1%*%Vi
  Dv.2<--Vi%*%B.2%*%Vi
  II<-matrix(0,2,2)
  P<-Vi-Vi%*%X%*%solve(t(X)%*%Vi%*%X)%*%t(X)%*%Vi
  for(i in 1:2) {
    for(j in 1:2) {
      II[i,j]<--0.5*sum(diag(P%*%B[[i]]%*%P%*%B[[j]]))
    }
  }
  II<--II
  ESS<-matrix(0,2,m)
  for (i in 1:m)  {
    ESS[1,]<-w.i[i,]%*%(Ds.1%*%t(w.i)%*%Vi+Sigma.w%*%t(w.i)%*%Dv.1)
    ESS[2,]<-w.i[i,]%*%(Ds.2%*%t(w.i)%*%Vi+Sigma.w%*%t(w.i)%*%Dv.2)
    g3[i,1]<-2*t(direct-X%*%Beta.hat)%*%t(ESS)%*%solve(II)%*%ESS%*%(direct-X%*%Beta.hat)
  }
  EBLUP.MSE.PR<-c(g1+g2+g3)
  areacode=1:m
  FH.SE=round(sqrt(EBLUP.MSE.PR),2)
  FH.CV=round(100*(sqrt(EBLUP.MSE.PR)/EBLUP.Mean),2)
  result1= cbind(areacode,EBLUP.Mean, EBLUP.MSE.PR,FH.SE,FH.CV)
  colnames(result1)=c("area","EBLUP","EBLUP.MSE","EBLUP.SE","EBLUP.CV")
  result <- list(eblup = NA, mse = NA, sample = NA,
                 fit = list(estcoef = NA, refvar = NA, spatialcorr = NA, randomeffect = NA))
  result$fit$estcoef <- coef
  result$fit$refvar <- sigma2.u.stim.S
  result$fit$spatialcorr <- lambda.stim.S
  result$fit$randomeffect <- u.hat
  result$eblup <- EBLUP.Mean
  result$mse <- EBLUP.MSE.PR
  result$sample <- result1
  return(result)
}

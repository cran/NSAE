
#' EBLUP under stationary Fay-Herriot model for sample area
#'
#'@description This function gives the EBLUP and the estimate of mean squared error (mse)
#'    based on a stationary Fay-Herriot model for sample area.
#'
#' @param formula an object of class list of formula, describe the model to be fitted
#' @param vardir a vector of sampling variances of direct estimators for each small area
#' @param method type of fitting method, default is "REML" method
#' @param data a data frame comprising the variables named in formula and vardir
#' @return The function returns a list with the following objects:
#' \describe{
#'   \item{eblup}{a vector with the values of the estimators for each small area}
#'   \item{mse}{a vector of the mean squared error estimates for each small area}
#'   \item{sample}{a matrix consist of area code, eblup, mse, standard error (SE) and coefficient of variation (CV)}
#'   \item{fit}{a list containing the following objects:}
#'   \itemize{
#'     \item estcoef : a data frame with the estimated model coefficients in the first column (beta),
#'         their asymptotic standard errors in the second column (std.error), the t statistics in
#'         the third column (tvalue) and the p-values of the significance of each coefficient in
#'         last column (pvalue)
#'     \item refvar : estimated random effects variance
#'     \item randomeffect : a data frame with the values of the random effect estimators
#'   }
#'  }
#' @export eblupFH1
#'
#' @examples
#' # Load data set
#' data(paddysample)
#' # Fit Fay-Herriot model using sample part of paddy data
#' result <- eblupFH1(y ~ x1+x2, var, method="REML", data = paddysample)
#' result
eblupFH1 <- function (formula, vardir, method = "REML", data) {
  namevar <- deparse(substitute(vardir))
  if (!missing(data)) {
    formuladata <- model.frame(formula, na.action = na.omit, data)
    X <- model.matrix(formula, data)
    vardir <- data[, namevar]
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
  direct <- y
  I<-diag(1,m)
  p<-dim(X)[2]
  logl=function(delta){
    area=m
    psi=matrix(c(vardir),area,1)
    Y=matrix(c(direct),area,1)
    Z.area=diag(1,area)
    sigma.u<-delta[1]
    V<-sigma.u*Z.area%*%t(Z.area)+I*psi[,1]
    Vi<-solve(V)
    Xt<-t(X)
    XVi<-Xt%*%Vi
    Q<-solve(XVi%*%X)
    P<-Vi-(Vi%*%X%*%Q%*%XVi)
    b.s<-Q%*%XVi%*%Y
    ee<-eigen(V)
    -(area/2)*log(2*pi)-0.5*sum(log(ee$value))-(0.5)*log(det(t(X)%*%Vi%*%X))-(0.5)*t(Y)%*%P%*%Y
  }
  opt<-optimize(logl,c(0.001,100),maximum = TRUE)
  estsigma2u<-opt$maximum
  Xt <-t(X)
  D=diag(1,m)
  V<-estsigma2u*D%*%t(D)+I*vardir
  Vi<-solve(V)
  Q<-solve(Xt%*%Vi%*%X)
  Beta.hat<-Q%*%Xt%*%Vi%*%direct
  res<-direct-X%*%Beta.hat
  Sigma.u=estsigma2u*I
  u.hat=Sigma.u%*%t(D)%*%Vi%*%res
  EBLUP.Mean<-X%*%Beta.hat+D%*%u.hat
  zvalue <- Beta.hat/sqrt(diag(Q))
  pvalue <- 2 * pnorm(abs(zvalue), lower.tail = FALSE)
  coef <- data.frame(beta = Beta.hat, std.error = sqrt(diag(Q)),tvalue = zvalue, pvalue)
  g1<-NULL
  Ga<-Sigma.u-Sigma.u%*%Vi%*%Sigma.u
  for (i in 1:m) {
    g1[i]<-Ga[i,i]
  }
  g2<-NULL
  Gb<-Sigma.u%*%Vi%*%X
  Xa<-matrix(0,1,p)
  for (i in 1:m) {
    Xa[1,]<-X[i,]-Gb[i,]
    g2[i]<-Xa%*%Q%*%t(Xa)
  }
  g3<-NULL
  Deriv1=solve((estsigma2u*I)+ diag(c(vardir),m))
  II<-((1/2)*sum(diag(Deriv1%*%Deriv1)))^(-1)
  for (i in 1:m) {
    g3[i]=(vardir[i]^2)*(vardir[i]+estsigma2u)^(-3)*II
  }
  EBLUP.MSE.PR<-c(g1+g2+g3)
  areacode=1:m
  FH.SE=round(sqrt(EBLUP.MSE.PR),2)
  FH.CV=round(100*(sqrt(EBLUP.MSE.PR)/EBLUP.Mean),2)
  result1= cbind(areacode,EBLUP.Mean, EBLUP.MSE.PR,FH.SE,FH.CV)
  colnames(result1)=c("area","EBLUP","EBLUP.MSE","EBLUP.SE","EBLUP.CV")
  result <- list(eblup = NA, mse = NA, sample = NA,fit = list(estcoef = NA, refvar = NA))
  result$fit$estcoef <- coef
  result$fit$refvar <- estsigma2u
  result$fit$randomeffect <- u.hat
  result$eblup <- EBLUP.Mean
  result$mse <- EBLUP.MSE.PR
  result$sample <- result1
  return(result)
}

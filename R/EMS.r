##' Compute the projection onto the column space of the matrix X.
##' @title Column Space of a Matrix
##' @param X a matrix.
##' @return A projection matrix.
qr.projection <- function(X) {
  QR <- qr(X)
  Q <- qr.Q(QR)[,seq_len(QR$rank),drop=F]
  attr(Q,"pivot") <- QR$pivot[seq_len(QR$rank)]
  Q
}



##' Compute the projection onto the orthogonal complement of the
##' column space of the matrix X.
##' @title Complement of Column Space of a Matrix
##' @param X a matrix.
##' @return A projection matrix.
qr.complement <- function(X) {
  ## orthonormal basis for column space
  Q <- qr.projection(X)
  ## Form P = I - QQ^T
  P <- -tcrossprod(Q)
  diag(P) <- diag(P)+1
  P
}


##' Compute the projection onto the model space.
##'
##' @title Projection onto Model Space
##' @param formula A formula describing the model.
##' @param data A data frame in which the variables specified in the
##'   model formula will be found.
##' @param type The "type" of test.
##' @importFrom stats model.matrix contrasts<- contr.sum delete.response
##' @return a list with elements
##' \item{Q}{the projection matrix.}
##' \item{assign}{a vector assigning the columns of Q to terms in the model.}
project.effects <- function(formula,data,type=c("I","II","III")) {

  ## Extract model frame and terms
  mf <- match.call(expand.dots = FALSE)
  mf <- mf[c(1L, match(c("formula", "data"), names(mf), 0L))]
  mf[[1L]] <- as.name("model.frame")
  mf <- eval(mf, parent.frame())
  mt <- delete.response(attr(mf, "terms"))
  ## Enforce true contrasts
  for(k in seq_len(ncol(mf)))
    if(is.factor(mf[,k])) contrasts(mf[,k]) <- contr.sum
  ## Model matrix and assignment of columns to terms
  X <- model.matrix(formula,mf)
  X.assign <- attr(X,"assign")

  ## Select SS type
  type <- match.arg(type)

  ## Type I sums of squares
  if(type=="I") {
    ## Compute the projection onto the model space, and the
    ## assignment of terms to columns
    Q <- qr.projection(X)
    assign <- X.assign[attr(Q,"pivot")]
  }

  ## Type II sums of squares
  if(type=="II") {
    ## Intercept term
    Q <- qr.projection(X[,X.assign==0,drop=F])
    assign <- rep(0,length=ncol(Q))
    ## For each remaining term
    for(k in seq_len(max(X.assign))) {
      ## Determine terms marginal to this one
      f <- ifelse(attr(mt,"factors") > 0,1,0)
      f <- apply(f[,k] <= f,2,all)
      f <- !(X.assign %in% which(f))
      ## Compute the projection of the model space for this term onto
      ## the orthogonal complement of the space of the marginal terms
      Xk <- qr.complement(X[,f,drop=F])%*%X[,X.assign==k,drop=F]
      Qk <- qr.projection(Xk)
      Q <- cbind(Q,Qk)
      assign <- c(assign,rep(k,length=ncol(Qk)))
    }
  }

  ## Type III sums of squares
  if(type=="III") {
    ## Intercept term
    Q <- qr.projection(X[,X.assign==0,drop=F])
    assign <- rep(0,length=ncol(Q))
    ## For each remaining term
    for(k in seq_len(max(X.assign))) {
      Xk <- qr.complement(X[,X.assign!=k,drop=F])%*%X[,X.assign==k,drop=F]
      Qk <- qr.projection(Xk)
      Q <- cbind(Q,Qk)
      assign <- c(assign,rep(k,length=ncol(Qk)))
    }
  }
  ## Return projection and the assignment of columns to terms
  list(Q=Q,assign=assign)
}


##' Construct the submatrix of the design matrix corresponding to a
##' single term of the model.
##'
##' @title Design Matrix for a Model Term
##' @param frame A model frame with columns corresponding to the
##'   factors in the term.
##' @param random A logical vector indicating which factors are
##'   random.
##' @param restricted A logical indicating whether to assume a
##'   restricted or unrestricted model.
##' @return A matrix
term.matrix <- function(frame,random,restricted) {
  if(!restricted && any(random)) {
    g <- interaction(frame)
    ifelse(outer(g,levels(g),"=="),1,0)
  } else {
    X <- 1
    for(k in 1:ncol(frame)) {
      m <- nlevels(frame[[k]])
      if(random[k])
        X <- kronecker(diag(1,m,m),X)
      else
        X <- kronecker(sqrt(m/(m-1))*(diag(1,m,m-1)+matrix(-1/m,m,m-1)),X)
    }
    X[unclass(interaction(frame)),]
  }
}


##' Compute expressions for the expeted means squares in a mixed model.
##'
##' Computes the coefficients in the expressions for the expected mean
##' squares for a given dataset.
##'
##' This approach to analysis of variance is not recommended.
##' @title Expected Mean Squares for Anova
##' @param formula A formula describing the model.
##' @param data A data frame in which the variables specified in the
##'   modle formula will be found.
##' @param random A character vector listing the predictors to be
##'   considered random.
##' @param type The "type" of test.
##' @param restricted Should EMS be computed assuming a restricted or
##'   unrestricted model.
##' @importFrom stats delete.response setNames
##' @return Returns an object of class "ems" with components
##' \item{\code{call}}{the matched call.}
##' \item{\code{restricted}}{logical indicating whether the EMS were computed for the restricted or unrestricted model.}
##' \item{variables}{a logical vector indicating which model variables are random.}
##' \item{terms}{a logical vector indicating which model terms are random.}
##' \item{EMS}{a table of the coefficients in the expansions of the EMS.}
##' @export
ems <- function(formula,data,random=NULL,type=c("III","I","II"),restricted=FALSE) {

  ## Model frame and terms
  cl <- match.call()
  mf <- match.call(expand.dots = FALSE)
  mf <- mf[c(1L, match(c("formula", "data"), names(mf), 0L))]
  mf[[1L]] <- as.name("model.frame")
  mf <- eval(mf, parent.frame())
  mt <- delete.response(attr(mf, "terms"))

  ## Projection onto effects and degrees of freedom
  type <- match.arg(type)
  eff <- project.effects(formula,mf,type=type)
  df <- table(eff$assign)
  rdf <- nrow(mf)-sum(df)
  if(attr(mt,"intercept")==1) {
    df <- df[-1]
    eff$Q <- eff$Q[,eff$assign>0]
    eff$assign <- eff$assign[eff$assign>0]
  }

  ## Term factors, variable names and term names
  fac <- attr(mt,"factors")
  vrs <- rownames(fac)
  tms <- colnames(fac)

  ## Which variables and terms are considered random
  rvrs <- setNames(vrs %in% random,vrs)
  rtms <- setNames(sapply(seq_along(tms),function(k) any(vrs[fac[,k]!=0] %in% random)),tms)

  ## EMS are the sums of squares of the projected columns of the
  ## (overparameterized) design matrix
  EMS <- matrix(0,length(tms),length(tms),dimnames=list(tms,tms))
  for(k in 1:ncol(fac)) {
    fs <- fac[,k]!=0
    Zk <- term.matrix(mf[,vrs[fs],drop=F],rvrs[fs],restricted)
    EMS[,k] <- tapply(rowSums(crossprod(eff$Q,Zk)^2),eff$assign,mean)
  }
  EMS <- cbind(Df=as.vector(c(df,rdf)),
               rbind(EMS,Error=0),
               Error=if(rdf==0) 0 else 1)

  r <- list(call=cl,restricted=restricted,
            variables=rvrs,terms=rtms,EMS=EMS)
  class(r) <- "ems"
  r
}


##' Print method for ems objects
##'
##' @title Print Method for ems Objects
##' @param x an ems object.
##' @param ... currently ignored.
##' @importFrom stats printCoefmat
##' @export
print.ems <-function(x,...) {
  cat("\nCall:\n", deparse(x$call), "\n", sep = "")
  cat("\nTerms:\n",
      paste(encodeString(names(x$terms),width=NULL),
            ifelse(x$terms,"Random","Fixed"),
            sep=": ",collapse="\n"),"\n\n",sep="")
  if(x$restricted)
    cat("EMS (Restricted model):\n")
  else
    cat("EMS (Unrestricted model):\n")
  printCoefmat(x$EMS,zap.ind=1:ncol(x$EMS))
}


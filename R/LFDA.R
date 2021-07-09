

LFDA<- function(x, y, r,prior = proportions, CV=FALSE,usekernel = TRUE, fL = 0,tol,kernel="gaussian",metric = c("orthonormalized","plain","weighted"),knn = 5,...) {
 requireNamespace("lfda")
   Y=as.factor(y)
  cl <- match.call()
  tol=tol
  metric <- match.arg(metric) # the type of the transforming matrix (metric)
  x=as.matrix(x)
  p <- ncol(x)
  n <- nrow(x) # number of samples
  if(n != length(y))
    stop("nrow(x) and length(y) are different")
  g <- as.factor(y)
  lev <- lev1 <- levels(g)
  counts <- as.vector(table(g))
  if(!missing(prior)) {
    if(any(prior < 0) || round(sum(prior), 5) != 1) stop("invalid 'prior'")
    if(length(prior) != nlevels(g)) stop("'prior' is of incorrect length")
    prior <- prior[counts > 0L]
    
  }
  if(any(counts == 0L)) {
    empty <- lev[counts == 0L]
    warning(sprintf(ngettext(length(empty),
                             "group %s is empty",
                             "groups %s are empty"),
                    paste(empty, collapse = " ")), domain = NA)
    lev1 <- lev[counts > 0L]
    g <- factor(g, levels = lev1)
    counts <- as.vector(table(g))
  }
  
  proportions <- counts/n
  ng <- length(proportions)
  names(prior) <- names(counts) <- lev1
  group.means <- tapply(c(x), list(rep(g, p), col(x)), mean)
  ### new section  from local da
  x <- t(as.matrix(x)) # transpose of original samples
  y <- t(as.matrix(y)) # transpose of original class labels
  d <- nrow(x) # number of predictors
  
  if(is.null(r)) r <- d # if no dimension reduction requested, set r to d
  
  tSb <- mat.or.vec(d, d) # initialize between-class scatter matrix (to be maximized)
  tSw <- mat.or.vec(d, d) # initialize within-class scatter matrix (to be minimized)
  
  # compute the optimal scatter matrices in a classwise manner
  for (i in unique(as.vector(t(y)))) {
    
    Xc <- x[, y == i] # data for this class
    nc <- ncol(Xc)
    
    # determine local scaling for locality-preserving projection
    Xc2 <- t(as.matrix(colSums(Xc^2)))  ## square of each class for each variable (predictor)
    # calculate the distance, using a self-defined repmat function that's the same
    # as repmat() in Matlab ,repmat is basically a X %x% Y, nc X nc matrix
    distance2 <- repmat(Xc2, nc, 1) + repmat(t(Xc2), 1, nc) - 2 * t(Xc) %*% Xc
    
    # Get affinity matrix, using knn, which is nearest neibger distancethe larger the element in the matrix, the closer two data points are
    A <- getAffinityMatrix(distance2, knn, nc)
    
    Xc1 <- as.matrix(rowSums(Xc))
    G <- Xc %*% (repmat(as.matrix(colSums(A)), 1, d) * t(Xc)) - Xc %*% A %*% t(Xc)
    tSb <- tSb + (G/n) + Xc %*% t(Xc) * (1 - nc/n) + Xc1 %*% (t(Xc1)/n)
    tSw <- tSw + G/nc
  }
  
  X1 <- as.matrix(rowSums(x))
  tSb <- tSb - X1 %*% t(X1)/n - tSw
  
  tSb <- (tSb + t(tSb))/2 # final between-class cluster matrix
  tSw <- (tSw + t(tSw))/2 # final within-class cluster matrix
  #A "scatter matrix" is just a covariance matrix without devidedness by sample_size-1.
  # find generalized eigenvalues and normalized eigenvectors of the problem
  if (r == d) {
    # without dimensionality reduction
    eigTmp <- eigen(ginv(tSw) %*% tSb)  # eigenvectors here are normalized, here we used ginv from MASS package to subsitute solve because solve can not solve singular values
  } else {
    # dimensionality reduction (select only the r largest eigenvalues of the problem)
    eigTmp <- suppressWarnings(rARPACK::eigs(A = ginv(tSw) %*% tSb, k = r, which = 'LM')) # r largest magnitude eigenvalues
  }
  eigVec <- Re(eigTmp$vectors) # the raw transforming matrix
  eigVal <- as.matrix(Re(eigTmp$values))
  
  # options to require a particular type of returned transform matrix
  # transforming matrix (do not change the "=" in the switch statement)
  Tr <- getMetricOfType(metric, eigVec, eigVal, d)
  # scaling is equal to Tr
  Z <- t(t(Tr) %*% x)
 # require(klaR) ## FOR Nativebaye function
 bayes_judgement=Mabayes(Z,Y,var.equal = FALSE,tol=tol)
  bayes=klaR::NaiveBayes(as.data.frame(Z), Y, prior=prior, usekernel=usekernel, fL = 0,kernel=kernel,...)
 
   bayes_assigment=predict.NaiveBayes(bayes)
  
  #  group.means1 <- tapply(c(Z), list(rep(g, p), col(Z)), mean)
  if(CV) {
    x <- t(t(Tr) %*% x) # transformed data Z
    
    dm <- group.means %*% Tr
    ### K is deternmined from n.pc, ng is number of gropus
    K = ng 
    dist <- matrix(0, n, ng) # n row and ng col
    
    ## dev matrix
    for(i in 1L:ng) {
      dev <- x - matrix(dm[i,  ], n, r, byrow = TRUE)
      dist[, i] <- rowSums(dev^2)### dis or devation of each class
    }
    ind <- cbind(1L:n, g) ### give the order number (sequence) and group(class )
    nc <- counts[g] ## number of class
    cc <- nc/((nc-1)*(n-K)) ### proportation of eac
    dist2 <- dist
    for(i in 1L:ng) {
      dev <- x - matrix(dm[i,  ], n, r, byrow = TRUE)
      dev2 <- x - dm[g, ]
      tmp <- rowSums(dev*dev2) 
      dist[, i] <- (n-1L-K)/(n-K) * (dist2[, i] +  cc*tmp^2/(1 - cc*dist2[ind]))
    }
    ### dist should be discriminat function 
    dist[ind] <- dist2[ind] * (n-1L-K)/(n-K) * (nc/(nc-1))^2 /(1 - cc*dist2[ind])
    dist <- 0.5 * dist - matrix(log(prior), n, ng, byrow = TRUE) # proboloty of discrim density function 
    
    dist <- exp(-(dist - min(dist, na.rm = TRUE))) #### proporbility of pi this is distribution Normal distribution function
    cla <- factor(lev1[max.col(dist)], levels = lev)
    ## calculate distance is   
    #    deltaTrain[,i]=deltaTrain[,i]-0.5*t(Mu_k[,i])%*%sigmaMinus1%*%Mu_k[,i]+log(Pi_k[i])
    ##  convert to posterior probabilities
    
    posterior <- dist/drop(dist %*% rep(1, length(prior)))
    dimnames(posterior) <- list(rownames(x), lev1)
    return(CV=list(prior=prior,tol=tol,usekernel=usekernel,fL=fL,kernel=kernel,counts = counts, means = group.means, Y=Y,post_class = cla,lev=lev,bayes=bayes,bayes_judgement=bayes_judgement,bayes_assigment=bayes_assigment,lev=lev, posterior = posterior,N=n,call=cl,Z=Z,T=Tr))
  }
  xbar <- colSums(prior %*% group.means)
  fac <-1/(ng - 1)
  X <- sqrt((n * prior)*fac) * scale(group.means, center = xbar, scale = FALSE) %*% Tr
  X.s <- svd(X, nu = 0L)
  # rank <- sum(X.s$d > tol * X.s$d[1L])
  #rank=r-1
  #T_lda <- Tr %*% X.s$v[, 1L:rank]
 # if(is.null(dimnames(x)))
  #  dimnames(T_lda) <- list(NULL, paste("LD", 1L:rank, sep = ""))
  #else {
   # dimnames(T_lda) <- list(colnames(x), paste("LD", 1L:rank, sep = ""))
    #dimnames(group.means)[[2L]] <- colnames(x)
#  }
  cl <- match.call()
  #svd = X.s$d[1L:rank]
  #PLD1=svd[1]^2/sum(svd^2)
  #PLD2=svd[2]^2/sum(svd^2)
  cl[[1L]] <- as.name("LFAD")
  structure(list(prior = prior,tol=tol,usekernel=usekernel,fL=fL,kernel=kernel,bayes=bayes,bayes_judgement=bayes_judgement, bayes_assigment=bayes_assigment,counts = counts, means = group.means, Z=Z, Tr=Tr,lev = lev, N = n,Y=Y, call = cl),
            class = "LFDA")
}






dkernel<-function(x, kernel=stats::density(x), interpolate=FALSE, ...)
{
  foo<-function(x,kernel,n)
  {
    if (x <= kernel$x[1]) return(kernel$y[1])
    else if (x >= kernel$x[n]) return(kernel$y[n])
    else 
    {
      pos<-which.min((kernel$x-x)^2)
      if (kernel$x[pos]>x) return(mean(c(kernel$y[pos],kernel$y[(pos-1)])))
      else return(mean(c(kernel$y[pos],kernel$y[(pos+1)])))
    }
  }
  n<-length(kernel$x)
  if (interpolate) y<-sapply(x,foo,kernel=kernel,n=n)
  else y<-sapply(x,FUN=function(y){kernel$y[(which.min((kernel$x-y)^2))]})
  return(y)
}

predict.NaiveBayes <- function (object, newdata, threshold = 0.001, ...)
{
  if (missing(newdata))
    newdata <- object$x
  if (sum(is.element(colnames(newdata), object$varnames)) < length(object$varnames))
    stop("Not all variable names used in object found in newdata") 
  ## (both colnames & varnames are given) & (varnames is a subset of colnames):
  newdata <- data.frame(newdata[, object$varnames])
  nattribs <- ncol(newdata)
  islogical <- sapply(newdata, is.logical)
  isnumeric <- sapply(newdata, is.numeric)
  #   L <- sapply(
  #      1:nrow(newdata),
  #      function(i)
  #      {
  #         ndata <- as.numeric(newdata[i, ])
  #         L <-  sapply(
  #               1:nattribs,
  #               function(v)
  #               {
  #                  nd <- ndata[v]
  #                  if (is.na(nd))
  #                  {
  #                     rep(1, length(object$apriori))
  #                  } else {
  #                     prob <- if (isnumeric[v])
  #                     {
  #                        msd <- object$tables[[v]]
  #                        if (object$usekernel) sapply(
  #                           msd,
  #                           FUN = function(y)
  #                           {
  #                              dkernel(x = nd, kernel = y, ...)
  #                           })
  #                        else dnorm(nd, msd[, 1], msd[, 2])
  #                     }
  #                     else object$tables[[v]][, nd]
  #
  #                     prob
  #                  }
  #               }
  #               )
  #
  #         L <- ifelse(L < threshold, threshold, L)
  #         # normalize by p(x) = p(x_1|y) + ... + p(x_p|y)
  #         Lnorm <- apply(L, 2, function(x, y) x/sum(x * y), y = as.vector(object$apriori))
  #         # get product
  #         Lprod <- apply(Lnorm, 1, prod)
  #         # normalize by posterior
  #         Lpost <- object$apriori * Lprod
  #         Lpost <- Lpost/sum(Lpost)
  #         Lpost
  #      }
  #   )
  newdata <- data.matrix(newdata)
  Lfoo <- function(i) {
    tempfoo <- function(v) {
      nd <- ndata[v]
      if (is.na(nd))
        return(rep(1, length(object$apriori)))
      prob <-
        if (isnumeric[v]) {
          msd <- object$tables[[v]]
          if (object$usekernel)
            sapply(msd, FUN = function(y)
              dkernel(x = nd, kernel = y, ...))
          else stats::dnorm(nd, msd[, 1], msd[, 2])
        } else if (islogical[v]) {
          object$tables[[v]][, nd + 1]
        } else {
          object$tables[[v]][, nd]
        }
      prob[prob == 0] <- threshold
      return(prob)
    }
    
    ndata <- newdata[i, ]
    tempres <- log(sapply(1:nattribs, tempfoo))
    L <- log(object$apriori) + rowSums(tempres)
    
    #        L <- exp(L)
    #        L/sum(L)
    
    if(isTRUE(all.equal(sum(exp(L)), 0)))
      warning("Numerical 0 probability for all classes with observation ", i)
    L
  }
  L <- sapply(1:nrow(newdata), Lfoo)
  
  classdach <- factor(object$levels[apply(L, 2, which.max)],
                      levels = object$levels)
  posterior <- t(apply(exp(L), 2, function(x) x/sum(x)))
  
  #                  print(str(posterior))
  colnames(posterior) <- object$levels
  rownames(posterior) <- names(classdach) <- rownames(newdata)
  return(list(class = classdach, posterior = posterior))
}

#' @export
predict.LFDA=function(object,prior=NULL,testData,...){
  tol=object$tol
  testData=as.matrix(testData)
  Z=object$Z
  Y=object$Y
  Trans=as.matrix(object$T)
  Z2=testData %*% Trans 
  if (is.null(prior)==TRUE){
    prior=object$prior
  }
  usekernel=object$usekernel
  fL=object$fL
  kernel=object$kernel
  # require(klaR) ## FOR Nativebaye function
  bayes_jud_pred=Mabayes(TrnX=Z,TrnG=Y,TstX = Z2,var.equal = FALSE,tol=tol)
  bayes=klaR::NaiveBayes(as.data.frame(Z), Y, prior, usekernel, fL,kernel,bw = "nrd0", adjust = 1,weights = NULL, window = kernel, give.Rkern = FALSE,...)
  Nbayes_assig_pred=predict.NaiveBayes(bayes,newdata=as.data.frame(Z2))
  ng=length(object$lev)
  
  means <- colSums(prior*object$means)
  scaling <- object$scaling
  x <-Z2
  dm <- scale(object$means, center = means, scale = FALSE) %*% Trans
  
  dimen <- length(object$svd) 
  N <- object$N
  
  dm <- dm[, 1L:dimen, drop = FALSE]
  dist <- matrix(0.5 * rowSums(dm^2) - log(prior), nrow(x),
                 length(prior), byrow = TRUE) - x[, 1L:dimen, drop=FALSE] %*% t(dm)
  dist <- exp( -(dist - apply(dist, 1L, min, na.rm=TRUE)))
  
  posterior <- dist / drop(dist %*% rep(1, length(prior)))
  nm <- names(object$prior)
  cl <- factor(nm[max.col(posterior)], levels = object$lev)
  dimnames(posterior) <- list(rownames(x), nm)
  list(class = cl, posterior = posterior, x = x[, 1L:dimen, drop = FALSE],bayes_jud_pred=bayes_jud_pred,Nbayes_assig_pred=Nbayes_assig_pred)
}

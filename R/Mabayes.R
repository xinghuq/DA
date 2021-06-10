### citation :	Bingpei Wu, 2012, WMDB 1.0: Discriminant Analysis Methods by Weight Mahalanobis Distance and bayes
### I made some modification from Bingpei Wu, 2012, WMDB, which some of the functions do not work properly when used in discriminant analysis

Mabayes<- function (TrnX, TrnG, p = rep(1, length(levels(TrnG))),TstX = NULL, var.equal = FALSE,tol){
  tol=tol
 
  mahalanobis=function (x, center, cov, inverted = FALSE,tol, ...) 
  {
    x <- if (is.vector(x)) 
      matrix(x, ncol = length(x))
    else as.matrix(x)
    if (!isFALSE(center)) 
      x <- sweep(x, 2L, center)
    if (!inverted) 
      cov <- solve(cov, tol,...)
    setNames(rowSums(x %*% cov * x), rownames(x))
  }
  
  mx <- nrow(TrnX)
  mg <- nrow(TrnG)
  TrnG <- as.factor(TrnG)
  G=as.factor(as.numeric(TrnG))
  counts <- as.vector(table(TrnG))
  lev1=levels(as.factor(TrnG))
  if (is.null(TstX) == TRUE) TstX <- TrnX
  if (is.vector(TstX) == TRUE) TstX <- t(as.matrix(TstX))
  else if 
  (is.matrix(TstX) != TRUE)
    TstX <- as.matrix(TstX)
  if (is.matrix(TrnX) != TRUE) TrnX <- as.matrix(TrnX)
  nx <- nrow(TstX)
  
  blong <- matrix(rep(0, nx), nrow=1,dimnames=list("blong", 1:nx))
  g <- length(levels(TrnG))
  mu <- matrix(0, nrow=g, ncol=ncol(TrnX))
  TrnX <- as.matrix(TrnX)
  
  # mu=aggregate(TrnX~TrnG, data=TrnX,mean)
  #mu=as.matrix(mu[,-1])
  for (i in 1:g) {
    mu[i, ] <- colMeans(TrnX[TrnG == lev1[i], ])
  }
  D <- matrix(0, nrow=g, ncol=nx)
  p = counts/mx
  if (var.equal == TRUE || var.equal == TRUE){
    for (i in 1:g){
      d2 <- mahalanobis(TstX, mu[i,], var(TrnX),tol)
      D[i,] <- d2 - 2*log(p[i])
    }
  }
  else{
    
    for (i in 1:g){
      S <- var(TrnX[TrnG == lev1[i], ])
      d2 <- mahalanobis(TstX, mu[i,], S,tol)
      D[i,] <- d2 - 2*log(p[i])-log(det(S))
    }
  }
  
  for (j in 1:nx){
    dmin <- Inf
    for (i in 1:g) {
      if ((!is.na(D[i,j]) < dmin)==TRUE){
      dmin <- D[i,j] 
      blong[j] <- i
    }
    }
  }
  
  dist <- exp(-(D/2 - min(D/2, na.rm = TRUE)))
  D1=t(D)
  posterior <- D1/drop(D1 %*% rep(1, g))
  posterior2 <- t(dist)/drop(t(dist) %*% rep(1, g))
  # if (is.null(TstX)!=TRUE){
  dimnames(posterior) <- list(rownames(TstX), lev1)
  dimnames(posterior2) <- list(rownames(TstX), lev1)
  
  cl <- factor(lev1[max.col(D1)], levels = lev1)
  cl2 <- factor(lev1[max.col(t(dist))], levels = lev1)
  #Assigments=factor(blong, labels = lev1)
  if (is.null(TstX) == TRUE){
  wrong_j=(cl2!=TrnG)
  sample_assign_F=(cl2[which(cl2!=TrnG)])
  sample_act_T=(TrnG[which(cl2!=TrnG)])
  success=(1-length(which(cl2!=TrnG))/length(cl2))
  
  
  
  print(blong)
  print("num of wrong discrimiant")
  print(which(blong!=G))
  print("samples wrongly assigned to")
  print(blong[which(blong!=G)])
  print("samples actually belong to")
  print(G[which(blong!=G)])
  print("percent of right judgement")
  print(1-length(which(blong!=G))/length(blong))    
  return(list(posterior2=posterior2,cl2=cl2,blong=blong,posterior=posterior,post_class=cl,wrong_j=wrong_j,sample_assign_F=sample_assign_F,sample_act_T=sample_act_T,success=success))
  }
 else
   return(list(blong=blong,D=D,posterior2=posterior2,cl2=cl2,posterior=posterior,post_class=cl))
  
}


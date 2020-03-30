

#### now here x can be matrix, data.frame, or kernel matrix, here we can use the kernel version of discriminant， these function also correct the existing function at R CRAN

LDAKPC <- function(x,y, n.pc,usekernel = FALSE, fL = 0,kernel.name = "rbfdot", kpar=list(0.001), kernel="gaussian",threshold = 1e-5,...){
  
  LDAKPC <- list()
  class(LDAKPC) <- "Linear Discriminant Analysis of Kernel principle components"
  
  # kpca
  require(kernlab)
  LDAKPC.train <- kernlab::kpca(~.,data=x,kernel = kernel.name,
                       kpar = kpar,
                       th = threshold,...)
  if (is.null(n.pc)){
    LDAKPC.rotation.train <- as.data.frame(LDAKPC.train@rotated)
  } else {
    LDAKPC.rotation.train <- as.data.frame(LDAKPC.train@rotated[,1:n.pc])}
  # KPC + lda 
  lda.rotation.train <- MASS::lda(LDAKPC.rotation.train,y,...)
  
  LDs <- as.matrix(LDAKPC.rotation.train) %*% as.matrix(lda.rotation.train$scaling)
  labels <- as.factor(y)
  
  LDAKPC$kpca<- LDAKPC.train
  LDAKPC$kpc=LDAKPC.rotation.train
  LDAKPC$LDAKPC<- lda.rotation.train
  LDAKPC$LDs <- LDs
  LDAKPC$label <- labels
  LDAKPC$n.pc=n.pc
  return(LDAKPC)
}

### once predict, r or n.pc should be the same with the input data, or the transformation will not work
predict.LDAKPC <- function(object = obj,prior=NULL, testData = data){
  n.pc=object$n.pc
  # kpca
  if(is.null(prior)==TRUE){
    prior=object$LDAKPC$prior
  }
  require(kernlab)
  predict.kpca <- kernlab::predict(object = object$kpca,
                          testData)[,1:n.pc]
  
  # kpca + lfda = lfdakpc
 
  predicted_LDs <- predict.kpca %*% as.matrix(object$LDAKPC$scaling)
  predict.LDAKPC <- predict(object$LDAKPC,prior,
                                 newdata = predict.kpca)
  
  return(list(predicted_LDs=predicted_LDs,predict.LDAKPC=predict.LDAKPC))
  
}

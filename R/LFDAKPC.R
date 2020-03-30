
####  this is fucntion for locial fisher DA of kernel PCA

#### now here x can be matrix, data.frame, or kernel matrix, here we can use the kernel version of discriminant

LFDAKPC <- function(x,y, n.pc,usekernel = FALSE, fL = 0,kernel.name = "rbfdot", kpar=list(0.001), kernel="gaussian",threshold = 1e-5,...){
  
  LFDAKPC  <- list()
  class(LFDAKPC ) <- "Local Fisher Discriminant Analysis of Kernel principle components"
  
  # kpca
  require(kernlab)
  LFDAKPC.train <- kpca(~.,data=x,kernel = kernel.name,
                       kpar = kpar,
                       th = threshold,...)
  if (is.null(n.pc)){
  LFDAKPC.rotation.train <- as.data.frame(LFDAKPC.train@rotated)
  } else {
    LFDAKPC.rotation.train <- as.data.frame(LFDAKPC.train@rotated[,1:n.pc])}
  # KPC + klda 
  klda.rotation.train <- LFDA(LFDAKPC.rotation.train,y,r=n.pc,usekernel =usekernel, 
                                fL = fL,kernel=kernel,...)
  
  LDs <- as.matrix(LFDAKPC.rotation.train) %*% as.matrix(klda.rotation.train$T)
  labels <- as.factor(y)
  
  LFDAKPC $kpca<- LFDAKPC.train
  LFDAKPC $kpc=LFDAKPC.rotation.train
  LFDAKPC $LFDAKPC<- klda.rotation.train
  LFDAKPC $LDs <- LDs
  LFDAKPC $label <- labels
  LFDAKPC $n.pc=n.pc
  return(LFDAKPC )
}

### once predict, r or n.pc should be the same with the input data, or the transformation will not work
predict.LFDAKPC  <- function(object = obj,prior, testData = data){
  n.pc=object$n.pc
  # kpca
  if(is.null(prior)==TRUE){
  prior=object$LFDAKPC$prior
  }
  predict.kpca <- kernlab::predict(object = object$kpca,
                          testData)[,1:n.pc]
  
  # kpca + lfda = klfdapc
  predicted_LDs <- predict.kpca %*% object$LFDAKPC$Tr
  predict.LFDAKPC <- predict.LFDA(object$LFDAKPC,prior,
                 newdata = predict.kpca)
  
  return(list(predicted_LDs=predicted_LDs,predict.LFDAKPC=predict.LFDAKPC))
  
}


multinomial_kernel=function(data1,data2,order, beta, gamma){
  ## Multinomial Kernel
  # Arguments:
  # - data1: data matrix where rows are observations and columns are features
  # - data1: same as data1
  # - order: the order of the multinomial expansion
  # - beta and gamma: no idea how to call them... coefficients? Default value
  # removes multiple instances of the same polynome
  
  if (ncol(data1)!=ncol(data2)) stop("data1 and data2 should have the same number of variables")
  
  nargin <- length(as.list(match.call())) -1 
  if (nargin<3) {
    if (order ==1)
      beta = order^(1/(order-1))
    gamma = order^(order/(1-order))
  }
  else beta = 1 
  gamma = 1
  
  
  K = gamma * ((1 + beta * (data2 %*% t(data1)) ^ order) - 1)
  
  return(K)
}

###### This is a traditional implementation of KLFDA using self-defined classifier, see Pierre Enel (2019). Kernel Fisher Discriminant Analysis (https://www.github.com/p-enel/MatlabKFDA), GitHub. Retrieved March 29, 2019.

KLFDA_mk=function(X, Y, r,order, regParam, usekernel=TRUE,fL=0.5,priors,tol,reg,metric,plotFigures=FALSE,
              verbose,...){
  ## Kernel Local Fisher Discriminant Analysis
  # (with polynomial kernel)
  # Arguments
  # - X: the predictors data, a matrix where rows are
  # observations and columns are features
  # - Y: the class of each observation, a cell of strings whose
  # length being equal to the number of rows in X
  # - order: the order of the polynomial expansion, integer
  # - regParam: optional, the regularization parameter. More
  # specifically, the regularization term is equal to the mean
  # value of the kernel matrix times regParam. Default value is
  # 0.25
  # - priors: optional, a structure, each field corresponding to
  # a class and each value coresponding to the prior probability
  # associated with this class. Default value is the empirical
  # probability: number of obs in a class / total number of obs
  requireNamespace("lfda")
  requireNamespace("MASS")
   n=ncol(X)
  if (nrow(X)!=nrow(Y))
    stop(" Number of observations for training data ,X,Y should be the same")
  
  
  obj.trainData = X
  obj.trainClass = Y
  obj.order = order;
  plotFigures = plotFigures;
  verbose = verbose;
  obj.classes = sort(unique(obj.trainClass));
  obj.nClasses = length(obj.classes);
  
  nargin <- length(as.list(match.call())) -1
  if (nargin < 4) {
  regParam = 0.25}
  if (nargin == 5) {
   if (length(priors)!= obj.nClasses)
    stop ("Number of prior should be equal to number of classes") 
}
#  in matlabe size is dim in R
  # Assessing the regularity of the number of features
  obj.nObservations=dim(obj.trainData)[1]
  obj.nFeatures = dim(obj.trainData)[2]
  
  # Concatenate training activity and generate class vectors (vectors where
                                                              # an entry is one if the pattern belong to this class)
  classVecsTrain = matrix(nrow=obj.nObservations, ncol=obj.nClasses);
  obj.nObsPerClas = matrix(nrow=1, ncol=obj.nClasses);
  for (i in seq_along(obj.classes)) {
  clas = obj.classes[i]
  classVecsTrain[, i] = match(obj.trainClass, clas,nomatch = 0)
  
  obj.nObsPerClas[i] = sum(classVecsTrain[,i])
  }
  

 
  ## II.1) Compute 'gram' (kernel, K) Matrix i.e. dot product of any two high-dim vectors.

if (verbose==TRUE) {
  print('Computing K...')
  
 K = multinomial_kernel(obj.trainData,obj.trainData, obj.order)
 K = K / obj.nObservations
}
if (verbose==FALSE){
  print('done!\n')
  }

## II.2) Compute the 'dual' to the means difference between classes ('M' matrix) and
## of pooled covariance ('N' matrix) on a high-dimensional space.

#First of all we need a 'dual' version of the averages per class.
#For that means, an auxiliary column of ones and zeros has to be
#firstly constructed per class.
MuQ = matrix(data=0,nrow=obj.nObservations,ncol= obj.nClasses)# The high dimensional mean of each class
Mu = matrix(data=0,nrow=obj.nObservations, ncol=1) # The mean of the means
for (i in 1 : obj.nClasses){
#See e.g. Scholkopf & Smola 02 ch 14 for justification of the next step.
MuQ[,i] = K %*% classVecsTrain[,i] / obj.nObsPerClas[i]
Mu = Mu + MuQ[,i]
}

Mu = Mu/ obj.nClasses; #Just the mean of the class-means (or centroids)

# Again, see e.g. Scholkopf & Smola 02 ch 14 for justification of the next
# step. M will represent the "dual" mean-differences (numerator in the FD ration)
# and N the "dual" pooled covariance.

M = matrix(data=0,obj.nObservations, obj.nObservations);
N = K %*% t(K)

for (i in 1 : obj.nClasses){
M = M + (MuQ[,i] - Mu) %*% t((MuQ[,i] - Mu))
N = N - (obj.nObsPerClas[i] * (MuQ[,i] %*% t(MuQ[,i])))
}
M = M * (obj.nClasses - 1) # across-class unbiased covariance
                                     
 # Regularizing
mK = abs(mean(K))
                                     
if (verbose==TRUE) {
  print(paste('Mean K is', mK))
}
C = regParam * mK
# The value of 'C' is such that the maximum eigenvector
                                    
 #         is proportional to alpha. This is taken as a stable solution.
                                    
 #         This value seems ok. C cannot be much smaller than
                                    
 #         0.01*mean(mean(K))in MSUA-like data. The bigger, the worse
                                   
  #         in-sample classification.
                                     
                                     
N = N + C * K; # Emulates SVM-like penalization (i.e. complexity).
                                     
                                     
# Extracting eigenvalues and eigenvectors
# find generalized eigenvalues and normalized eigenvectors of the problem
eigTmp <- suppressWarnings(rARPACK::eigs(A = ginv(N + reg * diag(1, nrow(N), ncol(N))) %*% M,
                                         k = r,which ='LM')) # r largest magnitude eigenvalues
eigVec <- Re(eigTmp$vectors) # the raw transforming matrix
eigVal <- as.matrix(Re(eigTmp$values))
                                     
#[Vtmp, lambda] = eig(M, N);
                                    
 #lambda = real(diag(lambda));
                                     
 # Warning: eigenvalues have to be ordered.
#[~, index] = sort(abs(lambda), 'descend');
#obj.V = Vtmp(:, index);
         
 projData = K %*% eigVec
 Tr <- getMetricOfType(metric, eigVec, eigVal, n)                              
  Z <- t(t(Tr) %*% K)
  
  for (i in 1 : obj.nClasses){
   projDataC = projData[as.logical(classVecsTrain[,i]),]
  }
 if (plotFigures==TRUE){
   
  
   
  # library(plotly)
   cols=grDevices::rainbow(obj.nClasses)
   p1 <- plotly::plot_ly(as.data.frame(Z), x =Z[,1], y =Z[,2], z =Z[,3], color = obj.trainClass,colors=cols,...) %>% 
     add_markers() %>%
     layout(scene = list(xaxis = list(title = 'LDA1'),
                         yaxis = list(title = 'LDA2'),
                         zaxis = list(title = 'LDA3')))
   #cols=rainbow(obj.nClasses)
   
   ## change the color by adding colors = c("#7570B3",'#BF382A', "#0C4B8E"  color=cols[obj.trainClass]
  # p2 <- plot_ly(as.data.frame(projDataC), x =projDataC[,1], y =projDataC[,2], z =projDataC[,3],color=obj.trainClass) %>% 
   #  add_markers() %>%
   #  layout(scene = list(xaxis = list(title = 'LDA1'),
     #                    yaxis = list(title = 'LDA2'),
     #                    zaxis = list(title = 'LDA3')))
   
   print(p1)
   #print(p2)
 }            
                                     
                                     # Reduced projected data to a dimensionality of obj.nClasses - 1,
                                     # sufficient for separation of the points in obj.nClasses
  #redProjData = matrix(data=0,nrow=obj.nObservations,ncol=obj.nClasses);
  redProjData=list()
  obj.means = matrix(data=0,nrow=1,ncol=obj.nClasses);
  obj.means0 = list()
  # obj.covariances = matrix(data=0,nrow=1,ncol=obj.nClasses);
  obj.covariances=list()
  for (i in 1 : obj.nClasses){
  clas = obj.classes[i];
  redProjData[[i]]=list()
   #redProjData[[i]] = projData[as.logical(classVecsTrain[,i]), (1 : (obj.nClasses - 1))]### here obj,nClasses -1 can be subsetuted by d, the reduced number of dim
  redProjData[[i]] = projData[as.logical(classVecsTrain[,i]), (1 : r)]
  obj.means[,i] = mean(redProjData[[i]])
  obj.means0[[i]]= colMeans(redProjData[[i]])
   obj.covariances[[i]]=list()
  obj.covariances[[i]] = cov(redProjData[[i]])
  }

  #### here we test whether the covariance is colinear and set tolerance, in some cases like our Shannon information data, the colinearity will abort the next step   
  
 # for (i in 1:obj.nClasses){
  #calculate singular value
 # sX <- svd(obj.covariances[[i]], nu = 0L)
 # rank <- sum(sX$d > tol^2)
 # if(rank == 0L) stop("rank = 0: variables are numerically constant")
 # if(rank < obj.nFeatures) warning("variables are collinear")  
 # }
  
  
  ###  aternational
  redZ=list()
  
  obj.means1 = list()
  obj.covariances1=list()
  for (i in 1 : obj.nClasses){
   
    redZ[[i]]=list()
    #redProjData[[i]] = projData[as.logical(classVecsTrain[,i]), (1 : (obj.nClasses - 1))]### here obj,nClasses -1 can be subsetuted by d, the reduced number of dim
    redZ[[i]] = Z[as.logical(classVecsTrain[,i]), (1 : r)]
    
    obj.means1[[i]]= colMeans(redZ[[i]])
    obj.covariances1[[i]]=list()
    obj.covariances1[[i]] = cov(redZ[[i]])
  }
## Generate priors and likelihood functions
  # Priors
  nargin <- length(as.list(match.call())) -1 
  if (is.null(priors)==TRUE){
   obj.priors = matrix(data=0,nrow=1, ncol=obj.nClasses);
  for (i in  1 : obj.nClasses){
  obj.priors[,i] = obj.nObsPerClas[,i] / obj.nObservations;
  }
  }
  if(is.null(priors)==FALSE) {
   priorsTmp = priors;
      obj.priors = matrix(data=0,nrow=1, ncol=obj.nClasses)
    for (i in 1 : obj.nClasses){
     clas = obj.classes[i]
    obj.priors[,i]= priorsTmp[i]
    }
  }
                                    
 # Likelihood
obj.likelihoodFuns = list()
# warning('error', 'nearlySingularMatrix')
obj.covInv=list()
obj.covDet=list()
factor=list()
#obj.covInv=solve(obj.covariances[[i]])
 for (i in 1 : obj.nClasses){
clas = obj.classes[i]
obj.covInv[[i]]=list()
obj.covDet[[i]]=list()
obj.covInv[[i]] = ginv(obj.covariances[[i]]);# ginv(obj.co) get inverse matrix 
obj.covDet[[i]] = det(obj.covariances[[i]]);

#pvaCov = vpa(obj.covariances[[i]]);
#obj.covInv[[i]] = double((pvaCov));
#obj.covDet.(clas) = double(det(pvaCov));
factor[[i]]=list()
factor[[i]] = (2 * pi) ^ (-(obj.nClasses - 1) / 2) * (obj.covDet[[i]] ^ -0.5)
 }

### the methos of Z
obj.covInv1=list()
obj.covDet1=list()
factor1=list()
#obj.covInv=solve(obj.covariances[[i]])
for (i in 1 : obj.nClasses){
  
  obj.covInv1[[i]]=list()
  obj.covDet1[[i]]=list()
  obj.covInv1[[i]] = ginv(obj.covariances1[[i]]);# ginv(obj.co)
  obj.covDet1[[i]] = det(obj.covariances1[[i]]);
  
  #pvaCov = vpa(obj.covariances[[i]]);
  #obj.covInv[[i]] = double((pvaCov));
  #obj.covDet.(clas) = double(det(pvaCov));
  factor1[[i]]=list()
  factor1[[i]] = (2 * pi) ^ (-(obj.nClasses - 1) / 2) * (obj.covDet1[[i]] ^ -0.5)
}

likelihoods = matrix(data=0,nrow = obj.nObservations, ncol=obj.nClasses);
likelihoods1 = matrix(data=0,nrow=obj.nObservations, ncol=obj.nClasses);
posteriors = matrix(data=0,nrow=obj.nObservations, ncol=obj.nClasses);
posteriors1 = matrix(data=0,nrow=obj.nObservations, ncol=obj.nClasses)
dist=matrix(data=0,nrow=obj.nObservations, ncol=obj.nClasses)
dist1=matrix(data=0,nrow=obj.nObservations, ncol=obj.nClasses)
## here from the  scripts  of Perr, each row of new data has nclass of 
# (rpdTest[j,] - obj.means[[i]]) %*% obj.covInv[[i]] * t(rpdTest[j,] -obj.means[[i]])
for (j in 1 : obj.nObservations){
  
  for (i in 1 : obj.nClasses){
    # clas = obj.classes[i]
    
    
    #mean((rpdTest[j,] - obj.means0[[i]]) %*% obj.covInv[[i]] * t(rpdTest[j,] -obj.means0[[i]]))
    dist[j,i]=mean((projData[j,] - obj.means0[[i]]) %*% obj.covInv[[i]] * t(projData[j,] -obj.means0[[i]]))
    dist1[j,i]=mean((Z[j,] - obj.means1[[i]]) %*% obj.covInv1[[i]] * t(Z[j,] -obj.means1[[i]]))
    # likelihoods[[j]][[i]] =factor[[i]] * exp(-0.5 * (rpdTest[j,] - obj.means[i]) %*% obj.covInv[[i]] * t(rpdTest[j,] -obj.means[i]))
    #likelihoods1 =factor[[i]] * exp(-0.5 * (rpdTest[j,] - obj.means0[[i]]) %*% obj.covInv[[i]] * t(rpdTest[j,] -obj.means0[[i]]))
    likelihoods[j,i] =factor[[i]] * exp(-0.5 * dist[j,i])
    likelihoods1[j,i] =factor1[[i]] * exp(-0.5 * dist1[j,i])
  
    
  }
  posteriors[j,] = likelihoods[j,] * obj.priors / sum(likelihoods[j,] * obj.priors)
  posteriors1[j,] = likelihoods1[j,] * obj.priors / sum(likelihoods1[j,] * obj.priors)
}
#}
### posteriors.class is the normal prodiction using projData posteriors.classZ using Z

## Predicting the class of each data point

posteriors.class=factor(obj.classes[max.col(posteriors)], levels = obj.classes)
posteriors.classZ=factor(obj.classes[max.col(posteriors1)], levels = obj.classes)

#requireNamespace("klaR")
bayes_proj=klaR::NaiveBayes(as.data.frame(projData), as.factor(obj.trainClass), prior=obj.priors, usekernel=usekernel, fL = fL,...)
kern_bayes_assigment_proj=predict.NaiveBayes(bayes_proj)
bayes_Z=klaR::NaiveBayes(as.data.frame(Z), as.factor(obj.trainClass), prior=obj.priors, usekernel=usekernel, fL = fL,...)
kern_bayes_assigment_Z=predict.NaiveBayes(bayes_Z)

return(list(bayes_proj=bayes_proj,bayes_Z=bayes_Z,kern_bayes_assigment_proj=kern_bayes_assigment_proj,kern_bayes_assigment_Z=kern_bayes_assigment_Z,posteriors=posteriors,posteriorsZ=posteriors1,posteriors.class=posteriors.class,posteriors.classZ=posteriors.classZ, obj.order=obj.order, obj.priors=obj.priors,obj.trainData=obj.trainData, obj.nClasses=obj.nClasses, obj.nFeatures=obj.nFeatures,classVecsTrain=classVecsTrain,obj.classes=obj.classes,obj.nObsPerClas=obj.nObsPerClas,obj.means=obj.means,obj.means0=obj.means0,obj.means1=obj.means1,obj.covInv=obj.covInv,obj.covInv1=obj.covInv1,factor=factor,factor1=factor1,eigVec=eigVec,eigVal=eigVal,Z=Z,Tr=Tr,LD=projDataC))
#obj.likelihoodFuns[[i]] =  factor * exp(-0.5 * (x - obj.means[[i]] * obj.covInv[[i]] * (x - t(obj.means[[i]]))

}
  

#' @export
predict.KLFDA_mk=function(object,prior=NULL, testData,...){
  ## Predict the class of new data
  # Arguments are the same as the X and Y in the constructor
  # function
  ## X is test data
  # Checking arguments
  
  #if (nrow(X)!= length(Y)) stop("number of class and number of observations are not equal")
  nObsTest=dim(testData)[1]
  nFeaturesTest = dim(testData)[2]
  obj.nClasses=object$obj.nClasses
  obj.nFeatures=object$obj.nFeatures
  if (obj.nFeatures != nFeaturesTest)
    stop('The number of features must be constant across classes and train/test data')
  
  
  #   classVecsTest = matrix(data=NA,nrow=nObsTest, ncol=obj.nClasses);
  #  Yid = matrix(data=0,nrow=nObsTest, ncol=1);
  # for (i in 1:obj.nClasses){
  #clas = obj.classes[i];
  #classVecsTest[,i] = match(Y, clas,nomatch=0)
  #Yid[as.logical(classVecsTest[,i])] = i;
  #}
  obj.trainData=object$obj.trainData
  obj.order=object$obj.order
  ## Computing the test kernel matrix
  K2 = multinomial_kernel(obj.trainData, testData, obj.order);
  # K2 = K2/ obj.nObservations;
  eigVec=object$eigVec
  Tr=object$Tr
  ## Projecting data onto the discriminant axes
  # rpdTest = K2 %*% eigVec[, (1 : (obj.nClasses - 1))] # Reduced projected data test
  rpdTest = K2 %*% eigVec
  Z2 <- K2 %*% Tr
  
  
  # requireNamespace("klaR") ## FOR Nativebaye function
  
  #bayes=NaiveBayes(Z, Y, prior, usekernel, fL,kernel,bw = "nrd0", adjust = 1,weights = NULL, window = kernel, give.Rkern = FALSE,...)
  Kern_Nbayes_assig_pred_proj=predict.NaiveBayes(object$bayes_proj,as.data.frame(rpdTest),...)
  Kern_Nbayes_assig_pred_Z=predict.NaiveBayes(object$bayes_Z,as.data.frame(Z2),...)
  ## Retrieving the likelihood of each test point
  likelihoods0 = matrix(data=0,nrow=nObsTest, ncol=obj.nClasses);
  likelihoods = matrix(data=0,nrow=nObsTest, ncol=obj.nClasses);
  likelihoods1 = matrix(data=0,nrow=nObsTest, ncol=obj.nClasses);
  posteriors0 = matrix(data=0,nrow=nObsTest, ncol=obj.nClasses);
  posteriors = matrix(data=0,nrow=nObsTest, ncol=obj.nClasses);
  posteriors1 = matrix(data=0,nrow=nObsTest, ncol=obj.nClasses)
  dist0=matrix(data=0,nrow=nObsTest, ncol=obj.nClasses)
  dist=matrix(data=0,nrow=nObsTest, ncol=obj.nClasses)
  dist1=matrix(data=0,nrow=nObsTest, ncol=obj.nClasses)
  ## here from the  scripts  of Perr, each row of new data has nclass of 
  # (rpdTest[j,] - obj.means[[i]]) %*% obj.covInv[[i]] * t(rpdTest[j,] -obj.means[[i]])
  obj.means=object$obj.means
  obj.means0=object$obj.means0
  obj.means1=object$obj.means1
  obj.covInv=object$obj.covInv
  
  obj.covInv1=object$obj.covInv1
  factor=object$factor
  factor1=object$factor1
  obj.classes=object$obj.classes
  if(is.null(prior)) 
    obj.priors=object$obj.priors
  else obj.priors=prior
  
  for (j in 1 : nObsTest){
    
    for (i in 1 : obj.nClasses){
      # clas = obj.classes[i]
      
      dist[j,i]=mean((rpdTest[j,] - obj.means[i]) %*% obj.covInv[[i]] * t(rpdTest[j,] -obj.means[i]))
      #mean((rpdTest[j,] - obj.means0[[i]]) %*% obj.covInv[[i]] * t(rpdTest[j,] -obj.means0[[i]]))
      dist0[j,i]=mean((rpdTest[j,] - obj.means0[[i]]) %*% obj.covInv[[i]] * t(rpdTest[j,] -obj.means0[[i]]))
      dist1[j,i]=mean((Z2[j,] - obj.means1[[i]]) %*% obj.covInv1[[i]] * t(Z2[j,] -obj.means1[[i]]))
      # likelihoods[[j]][[i]] =factor[[i]] * exp(-0.5 * (rpdTest[j,] - obj.means[i]) %*% obj.covInv[[i]] * t(rpdTest[j,] -obj.means[i]))
      #likelihoods1 =factor[[i]] * exp(-0.5 * (rpdTest[j,] - obj.means0[[i]]) %*% obj.covInv[[i]] * t(rpdTest[j,] -obj.means0[[i]]))
      likelihoods[j,i] =factor[[i]] * exp(-0.5 * dist[j,i])
      likelihoods0[j,i] =factor[[i]] * exp(-0.5 * dist0[j,i])
      likelihoods1[j,i] =factor1[[i]] * exp(-0.5 * dist1[j,i])
      
      
    }
    posteriors[j,] = likelihoods[j,] * obj.priors / sum(likelihoods[j,] * obj.priors)## rpdtest means using 
    posteriors0[j,] = likelihoods0[j,] * obj.priors / sum(likelihoods0[j,] * obj.priors)##rpdtest
    posteriors1[j,] = likelihoods1[j,] * obj.priors / sum(likelihoods1[j,] * obj.priors)## Z
  }
  #}
  
  
  ## Predicting the class of each data point
  
  posteriors.class=factor(obj.classes[max.col(posteriors)], levels = obj.classes)
  posteriors.class0=factor(obj.classes[max.col(posteriors0)], levels = obj.classes)
  posteriors.class1=factor(obj.classes[max.col(posteriors1)], levels = obj.classes)
  
  
  ## Creating the output variables
  
  results = list();
  results$Kern_Nbayes_assig_pred_proj=Kern_Nbayes_assig_pred_proj
  results$Kern_Nbayes_assig_pred_Z=Kern_Nbayes_assig_pred_Z
  results$likelihood = likelihoods1;
  results$posteriors = posteriors;
  results$posteriors0 = posteriors0;
  results$posteriors1 = posteriors1;
  results$posteriors.class = posteriors.class
  results$posteriors.class0 = posteriors.class0
  results$posteriors.class1 = posteriors.class1
  results$LD=rpdTest
  results$Z=Z2
  return(results)
}


      
############## KLFDA function that adopted from Tang et al 2016


kmatrixGauss=function (x, sigma = 1)
{
  x <- t(as.matrix(x))
  d <- nrow(x)
  n <- ncol(x)
  X2 <- t(as.matrix(colSums(x^2)))
  distance2 <- repmat(X2, n, 1) + repmat(t(X2), 1, n) - 2 *
    t(x) %*% x
  K <- exp(-distance2/(2 * sigma^2))
  return(K)
}



KLFDA=function(kdata, y, r,  metric = c("weighted", "orthonormalized",
                                    "plain"),tol=1e-5, knn = 6, reg = 0.001)
{
  requireNamespace("lfda")
  k=kdata
  
  obj.classes = sort(unique(y));
  obj.nClasses = length(obj.classes);
 
 
  
  metric <- match.arg(metric)
  y <- t(as.matrix(y))
  n <- nrow(k)
  if (is.null(r))
    r <- n
  tSb <- mat.or.vec(n, n)
  tSw <- mat.or.vec(n, n)
  for (i in unique(as.vector(t(y)))) {
    Kcc <- k[y == i, y == i]
    Kc <- k[, y == i]
    nc <- nrow(Kcc)
    Kccdiag <- diag(Kcc)
    distance2 <- repmat(Kccdiag, 1, nc) + repmat(t(Kccdiag),
                                                 nc, 1) - 2 * Kcc
    A <- getAffinityMatrix(distance2, knn, nc)
    Kc1 <- as.matrix(rowSums(Kc))
    Z <- Kc %*% (repmat(as.matrix(colSums(A)), 1, n) * t(Kc)) -
      Kc %*% A %*% t(Kc)
    tSb <- tSb + (Z/n) + Kc %*% t(Kc) * (1 - nc/n) + Kc1 %*%
      (t(Kc1)/n)
    tSw <- tSw + Z/nc
  }
  K1 <- as.matrix(rowSums(k))
  tSb <- tSb - K1 %*% t(K1)/n - tSw
  tSb <- (tSb + t(tSb))/2
  tSw <- (tSw + t(tSw))/2
  F=tSb/tSw

  eigTmp <- suppressWarnings(rARPACK::eigs(A = solve(tSw + reg * diag(1, nrow(tSw), ncol(tSw)),tol=tol) %*% tSb, k = r,
                                        which = "LM"))
  eigVec <- Re(eigTmp$vectors)
  eigVal <- as.matrix(Re(eigTmp$values))

  Tr <- getMetricOfType(metric, eigVec, eigVal, n)
  Z <- t(t(Tr) %*% k)
  out <- list(T = Tr, Z = Z,obj.classes = obj.classes,
  obj.nClasses = obj.nClasses, kmat=kdata)
  
  class(out) <- "KLFDA"
  return(out)
}



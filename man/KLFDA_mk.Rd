\name{KLFDA_mk}
\alias{KLFDA_mk}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{ Kernel Local Fisher Discriminant Analysis (KLFDA) with Multinomial kernel
%%  ~~function to do ... ~~
}
\description{Kernel Local Fisher Discriminant Analysis (KLFDA). This function implements the Kernel Local Fisher Discriminant Analysis with a Multinomial kernel. 
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
KLFDA_mk(X, Y, r, order, regParam, 
usekernel = TRUE, fL = 0.5, 
priors, tol, reg, metric, 
plotFigures = FALSE, verbose, ...)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{X}{The input training data
%%     ~~Describe \code{X} here~~
}
  \item{Y}{The training labels
%%     ~~Describe \code{Y} here~~
}
  \item{r}{ The number of reduced features
%%     ~~Describe \code{r} here~~
}
  \item{order}{ The order passing to Multinomial Kernel
%%     ~~Describe \code{order} here~~
}
  \item{regParam}{The regularization parameter for kernel matrix
%%     ~~Describe \code{regParam} here~~
}
  \item{usekernel}{ Whether to used kernel classifier
%%     ~~Describe \code{usekernel} here~~
}
  \item{fL}{ pass to kernel classifier if usekenel is TRUE
%%     ~~Describe \code{fL} here~~
}
  \item{priors}{The weight of each class
%%     ~~Describe \code{priors} here~~
}
  \item{tol}{ The tolerance for rejecting uni-variance
%%     ~~Describe \code{tol} here~~
}
  \item{reg}{ The regularization parameter 
%%     ~~Describe \code{reg} here~~
}
  \item{metric}{	Type of metric in the embedding space (default: 'weighted') 'weighted' - weighted eigenvectors 'orthonormalized' - orthonormalized 'plain' - raw eigenvectors
%%     ~~Describe \code{metric} here~~
}
  \item{plotFigures}{ whether to plot the reduced features, 3D plot
%%     ~~Describe \code{plotFigures} here~~
}
  \item{verbose}{ silence the processing
%%     ~~Describe \code{plotFigures} here~~
}
  \item{...}{ additional arguments for the classifier
%%     ~~Describe \code{plotFigures} here~~
}

}
\details{This function uses Multinomial Kernel, users can replace the Multinomial Kernel based on your own purpose. The final discrimination employs three classifiers, the basic linear classifier, the Mabayes (Bayes rule and the Mahalanobis distance), and Niave Bayes classifier.
%%  ~~ If necessary, more details than the description above ~~
}
\value{
%%  ~Describe the value returned
%%  If it is a LIST, use
\item{class}{The class labels from linear classifier}
\item{posterior}{The posterior possibility of each class from linear classifier}
\item{bayes_judgement}{Discrimintion results using the Mabayes classifier}
\item{bayes_assigment}{Discrimintion results using the Naive bayes classifier}
\item{Z }{The reduced features}
%% ...
}
\references{
Sugiyama, M (2007). Dimensionality reduction of multimodal labeled data by local Fisher discriminant analysis. Journal of Machine Learning Research, vol.8, 1027-1061.

Sugiyama, M (2006). Local Fisher discriminant analysis for supervised dimensionality reduction. In W. W. Cohen and A. Moore (Eds.), Proceedings of 23rd International Conference on Machine Learning (ICML2006), 905-912.

Original Matlab Implementation: http://www.ms.k.u-tokyo.ac.jp/software.html#LFDA

Tang, Y., & Li, W. (2019). lfda: Local Fisher Discriminant Analysis inR. Journal of Open Source Software, 4(39), 1572.

Moore, A. W. (2004). Naive Bayes Classifiers. In School of Computer Science. Carnegie Mellon University.

Pierre Enel (2020). Kernel Fisher Discriminant Analysis (https://www.github.com/p-enel/MatlabKFDA), GitHub. Retrieved March 30, 2020.

Karatzoglou, A., Smola, A., Hornik, K., & Zeileis, A. (2004). kernlab-an S4 package for kernel methods in R. Journal of statistical software, 11(9), 1-20.

%% ~put references to the literature/web site here ~
}
\author{qinxinghu@gmail.com
%%  ~~who you are~~
}


\seealso{predict.KLFDA_mk, klfda_1
%% ~~objects to See Also as \code{\link{help}}, ~~~
}
\examples{
btest=KLFDA_mk(X=as.matrix(iris[,1:4]),
Y=as.matrix(as.data.frame(iris[,5])),r=3,order=2,regParam=0.25, 
usekernel=TRUE,fL=0.5,
priors=NULL,tol=1e-90,reg=0.01,metric =  'plain',plotFigures=FALSE,
verbose=TRUE)
#pred=predict.KLFDA_mk(btest,as.matrix(iris[1:10,1:4]))
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{KLFDA_mk}% use one of  RShowDoc("KEYWORDS")


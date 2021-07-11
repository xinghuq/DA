\name{KLFDA}
\alias{KLFDA}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{ Kernel Local Fisher Discriminant Analysis (KLFDA)
%%  ~~function to do ... ~~
}
\description{Kernel Local Fisher Discriminant Analysis (KLFDA). This function implements the Kernel Local Fisher Discriminant Analysis with an unified Kernel function. Different from KLFDA function, which adopts the Multinomial Kernel as an example, this function empolys the kernel function that allows you to choose various types of kernels. See the kernel function from "kernelMatriax" (kernlab).
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
KLFDA(x, y, kernel = kernlab::polydot(degree = 1, scale = 1, offset = 1), 
r = 20, tol, prior, CV = FALSE, usekernel = TRUE, 
fL = 0.5, metric = c("weighted", "orthonormalized", "plain"), 
knn = 6, reg = 0.001, ...)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{x}{ The input training data
%%     ~~Describe \code{x} here~~
}
  \item{y}{ The training labels
%%     ~~Describe \code{y} here~~
}
  \item{kernel}{ The kernel function used to calculate kernel matrix. Choose the corresponding kernel you want, see details.
%%     ~~Describe \code{kernel} here~~
}
  \item{r}{ The number of reduced features you want to keep.
%%     ~~Describe \code{r} here~~
}
  \item{tol}{ The tolerance used to reject the uni-variance. This is important when the variance between classes is small, and setting the large tolerance will avoid the data distortion.
%%     ~~Describe \code{tol} here~~
}
  \item{prior}{ The weight of each class, or the proportion of each class.
%%     ~~Describe \code{prior} here~~
}
  \item{CV}{ Whether to do cross validation.
%%     ~~Describe \code{CV} here~~
}
  \item{usekernel}{whether to use kernel classifier, if TRUE, pass to Naive Bayes classifier.
%%     ~~Describe \code{usekernel} here~~
}
  \item{fL}{ If usekernel is TRUE, pass to the kernel function.
%%     ~~Describe \code{fL} here~~
}
  \item{metric}{type of metric in the embedding space (default: 'weighted') 'weighted' - weighted eigenvectors 'orthonormalized' - orthonormalized 'plain' - raw eigenvectors
%%     ~~Describe \code{metric} here~~
}
  \item{knn}{ The number of nearest neighbours
%%     ~~Describe \code{knn} here~~
}
  \item{reg}{ The regularization parameter
%%     ~~Describe \code{reg} here~~
}
  \item{...}{ additional arguments for the classifier
%%     ~~Describe \code{plotFigures} here~~
}
}
\details{
This function empolys three different classifiers, the basic linear classifier, the Mabayes (Bayes rule and the Mahalanobis distance), and Niave Bayes classifier.
The argeument "kernel" in the klfda function is the kernel function used to calculate the kernel matrix. If usekernel is TRUE, the corresponding kernel parameters will pass the the Naive Bayes kernel classifier.
The kernel parameter can be set to any function, of class kernel, which computes the inner product in feature space between two vector arguments. kernlab provides the most popular kernel functions which can be initialized by using the following functions:

rbfdot Radial Basis kernel function

polydot Polynomial kernel function

vanilladot Linear kernel function

tanhdot Hyperbolic tangent kernel function

laplacedot Laplacian kernel function

besseldot Bessel kernel function

anovadot ANOVA RBF kernel function

splinedot the Spline kernel

(see example.)

kernelFast is mainly used in situations where columns of the kernel matrix are computed per invocation. In these cases, evaluating the norm of each row-entry over and over again would cause significant computational overhead.
%%  ~~ If necessary, more details than the description above ~~
}
\value{ The results give the classified classes and the posterior possibility of each class using different classifier. 
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

Sugiyama, M (2007).Dimensionality reduction of multimodal labeled data by local Fisher discriminant analysis. Journal of Machine Learning Research, vol.8, 1027-1061.

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

\seealso{predict.KLFDA, KLFDAM
%% ~~objects to See Also as \code{\link{help}}, ~~~
}
\examples{
require(kernlab)
btest=KLFDA(as.matrix(iris[,1:4]),as.matrix(as.data.frame(iris[,5])),
kernel=kernlab::rbfdot(sigma = 0.1),
r=3,prior=NULL,tol=1e-90,
reg=0.01,metric =  'plain')
pred=predict.KLFDA(btest,testData=as.matrix(iris[1:10,1:4]),prior=NULL)
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{KLFDA}% use one of  RShowDoc("KEYWORDS")


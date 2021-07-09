\name{LFDAKPC}
\alias{LFDAKPC}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{Local Fisher Discriminant Analysis of Kernel principle components (LFDAKPC)
%%  ~~function to do ... ~~
}
\description{Local Fisher Discriminant Analysis of Kernel principle components
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
LFDAKPC(x, y, n.pc, 
usekernel = FALSE, fL = 0, 
kernel.name = "rbfdot", 
kpar = list(0.001), kernel = "gaussian", 
threshold = 1e-05, ...)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{x}{ Input traing data
%%     ~~Describe \code{x} here~~
}
  \item{y}{ Input labels
%%     ~~Describe \code{y} here~~
}
  \item{n.pc}{ number of pcs that will be kept in analysis
%%     ~~Describe \code{n.pc} here~~
}
  \item{usekernel}{ Whether to use kernel function, if TRUE, it will pass to the kernel.names
%%     ~~Describe \code{usekernel} here~~
}
  \item{fL}{ if using kernel, pass to kernel function
%%     ~~Describe \code{fL} here~~
}
  \item{kernel.name}{ if usekernel is TURE, this will take the kernel name and use the parameters set as you defined
%%     ~~Describe \code{kernel.name} here~~
}
  \item{kpar}{	
the list of hyper-parameters (kernel parameters). This is a list which contains the parameters to be used with the kernel function. Valid parameters for existing kernels are :

sigma inverse kernel width for the Radial Basis kernel function "rbfdot" and the Laplacian kernel "laplacedot".

degree, scale, offset for the Polynomial kernel "polydot"

scale, offset for the Hyperbolic tangent kernel function "tanhdot"

sigma, order, degree for the Bessel kernel "besseldot".

sigma, degree for the ANOVA kernel "anovadot".

Hyper-parameters for user defined kernels can be passed through the kpar parameter as well.
%%     ~~Describe \code{kpar} here~~
}
  \item{kernel}{ kernel name if all the above are not used
%%     ~~Describe \code{kernel} here~~
}
  \item{threshold}{	
the threshold for kpc: value of the eigenvalue under which principal components are ignored (only valid when features = 0). (default : 0.0001)
}
  \item{...}{ additional arguments for the classifier
%%     ~~Describe \code{plotFigures} here~~
}
}

\value{
\item{kpca}{Results of kernel principal component analysis. Kernel Principal Components Analysis is a nonlinear form of principal component analysis}
\item{kpc}{Kernel principal components. The scores of the components}
\item{LFDAKPC}{LOcal linear discriminant anslysis of kernel principal components}
\item{LDs}{The discriminant function. The scores of the components}
\item{label}{The corresponding class of the data}
\item{n.pc}{Number of Pcs kept in analysis}

}
\references{
Sugiyama, M (2007). Dimensionality reduction of multimodal labeled data by local Fisher discriminant analysis. Journal of Machine Learning Research, vol.8, 1027-1061.

Sugiyama, M (2006). Local Fisher discriminant analysis for supervised dimensionality reduction. In W. W. Cohen and A. Moore (Eds.), Proceedings of 23rd International Conference on Machine Learning (ICML2006), 905-912.

Tang, Y., & Li, W. (2019). lfda: Local Fisher Discriminant Analysis inR. Journal of Open Source Software, 4(39), 1572.


Karatzoglou, A., Smola, A., Hornik, K., & Zeileis, A. (2004). kernlab-an S4 package for kernel methods in R. Journal of statistical software, 11(9), 1-20.

%% ~put references to the literature/web site here ~
}
\author{qinxinghu@gmail.com
%%  ~~who you are~~
}

\examples{

train=LFDAKPC(iris[,1:4],y=iris[,5],tol=1,n.pc=3,kernel.name = "rbfdot")
pred=predict.LFDAKPC(train,prior=NULL,testData = iris[1:10,1:4])
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{LFDAKPC}% __ONLY ONE__ keyword per line

\name{LDAKPC}
\alias{LDAKPC}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{Linear Fisher discriminant analysis of kernel principal components (DAKPC)
%%  ~~function to do ... ~~
}
\description{ Linear Fisher discriminant analysis of kernel principal components (DAKPC). This function empolies the LDA and kpca. This function is called Kernel Fisher Discriminant Analysis (KFDA) in other package (kfda). "KFDA" is the misleading name and "KFDA" has crucial error in package kfda. This function rectifies the current existing error for kfda.
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
LDAKPC(x, y, n.pc, usekernel = FALSE, 
fL = 0, kernel.name = "rbfdot", 
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
%%     ~~Describe \code{threshold} here~~
}
  \item{...}{ additional arguments for the classifier
%%     ~~Describe \code{plotFigures} here~~
}
}
\value{
%%  ~Describe the value returned
%%  If it is a LIST, use
\item{kpca}{Results of kernel principal component analysis. Kernel Principal Components Analysis is a nonlinear form of principal component analysis}
\item{kpc}{Kernel principal components. The scores of the components}
\item{LDAKPC}{Linear discriminant anslysis of kernel principal components}
\item{LDs}{The discriminant function. The scores of the components}
\item{label}{The corresponding class of the data}
\item{n.pc}{Number of Pcs kept in analysis}

%% ...
}
\references{
Karatzoglou, A., Smola, A., Hornik, K., & Zeileis, A. (2004). kernlab-an S4 package for kernel methods in R. Journal of statistical software, 11(9), 1-20.

Mika, S., Ratsch, G., Weston, J., Scholkopf, B., & Mullers, K. R. (1999, August). Fisher discriminant analysis with kernels. In Neural networks for signal processing IX: Proceedings of the 1999 IEEE signal processing society workshop (cat. no. 98th8468) (pp. 41-48). Ieee.
%% ~put references to the literature/web site here ~
}
\author{
qinxinghu@gmail.com
%%  ~~who you are~~
}

\examples{
data(iris)
train=LDAKPC(iris[,1:4],y=iris[,5],n.pc=3,kernel.name = "rbfdot")
pred=predict.LDAKPC(train,testData = iris[1:10,1:4])
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{LDAKPC}% use one of  RShowDoc("KEYWORDS")


\name{LFDA}
\alias{LFDA}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{ Local Fisher Discriminant Analysis (LFDA)
%%  ~~function to do ... ~~
}
\description{ This function implements local Fisher discriminant analysis. It gives the discriminant function with the posterior possibility of each class.
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
LFDA(x, y, r, prior = proportions, CV = FALSE, usekernel = TRUE, fL = 0, tol, kernel = "gaussian", metric = c("orthonormalized", "plain", "weighted"), knn = 5, ...)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{x}{ Input training data
%%     ~~Describe \code{x} here~~
}
  \item{y}{ Training labels
%%     ~~Describe \code{y} here~~
}
  \item{r}{ Number of reduced features that will be kept
%%     ~~Describe \code{r} here~~
}
  \item{prior}{ Prior possibility of each class
%%     ~~Describe \code{prior} here~~
}
  \item{CV}{ Whether to do cross validation
%%     ~~Describe \code{CV} here~~
}
  \item{usekernel}{ Whether to use the kernel discrimination in native bayes classifier
%%     ~~Describe \code{usekernel} here~~
}
  \item{fL}{ Feed to native bayes classifier. Factor for Laplace correction, default factor is 0, i.e. no correction.
%%     ~~Describe \code{fL} here~~
}
  \item{tol}{ The tolerance used in Mabayes discrimination, see Mabayes
%%     ~~Describe \code{tol} here~~
}
  \item{kernel}{ If usekernel is TRUE, specifying the kernel names, see NaiveBaye.
%%     ~~Describe \code{kernel} here~~
}
  \item{metric}{ The type of metric in the embedding space (no default), e.g., 'weighted',  weighted eigenvectors; 'orthonormalized' , orthonormalized; 'plain',  raw eigenvectors.
%%     ~~Describe \code{metric} here~~
}
  \item{knn}{ Number of nearest neighbors
%%     ~~Describe \code{knn} here~~
}
 

}
\details{The results give the classified classes and the posterior possibility of each class using different classifier.
%%  ~~ If necessary, more details than the description above ~~
}
\value{
%%  ~Describe the value returned
%%  If it is a LIST, use
\item{class}{The class labels}
\item{posterior}{The posterior possibility of each class}
\item{bayes_judgement}{Discrimintion results using the Mabayes classifier}
\item{bayes_assigment}{Discrimintion results using the Naive bayes classifier}
\item{Z }{The reduced features}

%% ...
}
\references{
Sugiyama, M (2007). Dimensionality reduction of multimodal labeled data by local Fisher discriminant analysis. Journal of Machine Learning Research, vol.8, 1027-1061.

Sugiyama, M (2006). Local Fisher discriminant analysis for supervised dimensionality reduction. In W. W. Cohen and A. Moore (Eds.), Proceedings of 23rd International Conference on Machine Learning (ICML2006), 905-912.

Tang, Y., & Li, W. (2019). lfda: Local Fisher Discriminant Analysis inR. Journal of Open Source Software, 4(39), 1572.

Moore, A. W. (2004). Naive Bayes Classifiers. In School of Computer Science. Carnegie Mellon University.

Pierre Enel (2020). Kernel Fisher Discriminant Analysis (https://www.github.com/p-enel/MatlabKFDA), GitHub. Retrieved March 30, 2020.
%% ~put references to the literature/web site here ~
}
\author{qinxinghu@gmail.com
%%  ~~who you are~~
}
\note{
%%  ~~further notes~~
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
%% ~~objects to See Also as \code{\link{help}}, ~~~
}
\examples{
LFDAtest=LFDA(iris[,1:4],y=iris[,5],r=3, CV=FALSE,usekernel = TRUE, fL = 0,kernel="gaussian",metric = "plain",knn = 6,tol = 1)
LFDApred=predict.LFDA(LFDAtest,newdata=iris[1:10,1:4],prior=NULL)
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{LFDA}% use one of  RShowDoc("KEYWORDS")


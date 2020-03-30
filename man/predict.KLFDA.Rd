\name{predict.KLFDA}
\alias{predict.KLFDA}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{Prediction function for Kernel Local Fisher Discriminant Analysis (KLFDA)
%%  ~~function to do ... ~~
}
\description{Prediction function for Kernel Local Fisher Discriminant Analysis (KLFDA) with a Multinomial kernel. Predictions are based on three classifiers. See KLFDA for detail. Results give the class labels and posterior possibility of the tested data.
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
predict.KLFDA(obj, newdata)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{obj}{The object from KLFDA model
%%     ~~Describe \code{obj} here~~
}
  \item{newdata}{ The tested data
%%     ~~Describe \code{newdata} here~~
}
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


%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{KLFDA, predict.klfda_1
%% ~~objects to See Also as \code{\link{help}}, ~~~
}
\examples{
btest=KLFDA(X=as.matrix(iris[,1:4]),Y=as.matrix(as.data.frame(iris[,5])),r=3,order=2,regParam=0.25, usekernel=TRUE,fL=0.5,priors=NULL,tol=1e-90,reg=0.01,metric =  'plain',plotFigures=FALSE,
           verbose=TRUE)
pred=predict.KLFDA(btest,as.matrix(iris[1:10,1:4]))
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{predict.KLFDA}% use one of  RShowDoc("KEYWORDS")


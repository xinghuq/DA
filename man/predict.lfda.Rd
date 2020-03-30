\name{predict.LFDA}
\alias{predict.LFDA}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{ Prediction for Local Fisher Discriminant Analysis (LFDA)
%%  ~~function to do ... ~~
}
\description{Prediction function for Local Fisher Discriminant Analysis (LFDA)
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
predict.LFDA(object, newdata, prior, dimen, ...)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{object}{ The LFDA object
%%     ~~Describe \code{object} here~~
}
  \item{newdata}{ The newdata you want to predict
%%     ~~Describe \code{newdata} here~~
}
  \item{prior}{ The prior of the new data
%%     ~~Describe \code{prior} here~~
}
  \item{dimen}{ The predicted dimen based on training model
%%     ~~Describe \code{dimen} here~~
}

}
\details{
%%  ~~ If necessary, more details than the description above ~~
}
\value{
%%  ~Describe the value returned
%%  If it is a LIST, use
\item{class}{The class labels from liner classifier}
\item{posterior}{The posterior possibility from linear classifier}
\item{bayes_jud_pred}{The predicted results using Mabayes classifier}
\item{Nbayes_assig_pred}{The predicted results using Naive Bayes classifier}
%% ...
}
\references{
Sugiyama, M (2007). Dimensionality reduction of multimodal labeled data by local Fisher discriminant analysis. Journal of Machine Learning Research, vol.8, 1027-1061.

Sugiyama, M (2006). Local Fisher discriminant analysis for supervised dimensionality reduction. In W. W. Cohen and A. Moore (Eds.), Proceedings of 23rd International Conference on Machine Learning (ICML2006), 905-912.

Tang, Y., & Li, W. (2019). lfda: Local Fisher Discriminant Analysis inR. Journal of Open Source Software, 4(39), 1572.

Moore, A. W. (2004). Naive Bayes Classifiers. In School of Computer Science. Carnegie Mellon University.
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
\keyword{predict.LFDA}% use one of  RShowDoc("KEYWORDS")


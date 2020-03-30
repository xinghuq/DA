\name{predict.LFDAKPC}
\alias{predict.LFDAKPC}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{ Prediction of Local Fisher Discriminant Analysis of Kernel principle components (LFDAKPC)
%%  ~~function to do ... ~~
}
\description{This is the prediction function of Local Fisher Discriminant Analysis of Kernel principle components (LFDAKPC)
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
predict.LFDAKPC(object = obj, prior, testData = data)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{object}{ LFDAKPC object
%%     ~~Describe \code{object} here~~
}
  \item{prior}{ The prior proportion or weight
%%     ~~Describe \code{prior} here~~
}
  \item{testData}{ The test data
%%     ~~Describe \code{testData} here~~
}
}
\details{
%%  ~~ If necessary, more details than the description above ~~
}
\value{
%%  ~Describe the value returned
%%  If it is a LIST, use
\item{predicted_LDs}{predicted discriminant function}
\item{predict.LFDAKPC}{predicted results of the test data}
%% ...
}
\references{
Karatzoglou, A., Smola, A., Hornik, K., & Zeileis, A. (2004). kernlab-an S4 package for kernel methods in R. Journal of statistical software, 11(9), 1-20.

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
\keyword{ predict.LFDAKPC }% use one of  RShowDoc("KEYWORDS")


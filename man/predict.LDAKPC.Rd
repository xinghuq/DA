\name{predict.LDAKPC}
\alias{predict.LDAKPC}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{Predict function of Linear Fisher discriminant analysis of kernel principal components (DAKPC)
%%  ~~function to do ... ~~
}
\description{ Predict function of Linear Fisher discriminant analysis of kernel principal components (DAKPC)
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
predict.LDAKPC(object = obj, prior, testData = data)
}
%- maybe also 'usage' for other objects documented here.
\arguments{

\item{object}{ The trained LDAKPC object}
\item{prior}{The prior weight of the predicted data}
\item{testData}{The data you want to test/predict}
}
\details{
%%  ~~ If necessary, more details than the description above ~~
}
\value{
%%  ~Describe the value returned
%%  If it is a LIST, use
\item{posterior}{The posterior probabilities of the predicted data.}
\item{x }{The discriminant function of the predicted data.}
\item{predicted_LDs}{The predicted scores of discriminant function, is always the same with x if there is no transformation.}
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
\keyword{predict.LDAKPC}% use one of  RShowDoc("KEYWORDS")


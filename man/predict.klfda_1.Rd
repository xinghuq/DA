\name{predict.klfda_1}
\alias{predict.klfda_1}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{ The predict method for KLFDA(klfda_1)
%%  ~~function to do ... ~~
}
\description{Prediction function for Kernel Local Fisher Discriminant Analysis KLFDA(klfda_1). The prediction uses three classifiers, as explained in the manual of klfda_1.
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
predict.klfda_1(object, newdata, dimen, ...)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{object}{The object from Kernel Local Fisher Discriminant Analysis KLFDA(klfda_1).
%%     ~~Describe \code{object} here~~
}
  \item{newdata}{The test data. The data you want to predict.
%%     ~~Describe \code{newdata} here~~
}
  \item{dimen}{ The data dimen, if NULL, the function will automatically learns from the data.
%%     ~~Describe \code{dimen} here~~
}

}
\details{
%%  ~~ If necessary, more details than the description above ~~
}
\value{The values include the class and posterior possibility of each class predicted by three classifiers.
%%  ~Describe the value returned
%%  If it is a LIST, use
\item{class}{The class labels from liner classifier}
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

\seealso{klfda_1, KLFDA, predict.KLFDA
%% ~~objects to See Also as \code{\link{help}}, ~~~
}
\examples{
btest=klfda_1(as.matrix(iris[,1:4]),as.matrix(as.data.frame(iris[,5])),kernel=kernlab::rbfdot(sigma = 0.1),r=3,prior=NULL,tol=1e-90,reg=0.01,metric =  'plain')
pred=predict.klfda_1(btest,as.matrix(iris[1:10,1:4]))
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{predict.klfda_1}% use one of  RShowDoc("KEYWORDS")


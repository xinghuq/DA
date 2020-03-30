\name{Mabayes}
\alias{Mabayes}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{ Membership assignment by weighted Mahalanobis distance and bayes rule
%%  ~~function to do ... ~~
}
\description{ The function gives the discrimintion of the potential classes based on Bayes rule and the Mahalanobis distance. This function adopts the function from  Bingpei Wu, 2012, WMDB 1.0 with some corrections of the judement rule.
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
Mabayes(TrnX, TrnG, p = rep(1, length(levels(TrnG))), TstX = NULL, var.equal = FALSE, tol)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{TrnX}{ Training data
%%     ~~Describe \code{TrnX} here~~
}
  \item{TrnG}{ Training label
%%     ~~Describe \code{TrnG} here~~
}
  \item{p}{ prior or proportion of each class
%%     ~~Describe \code{p} here~~
}
  \item{TstX}{ Test data
%%     ~~Describe \code{TstX} here~~
}
  \item{var.equal}{ whether the variance or the weight is equal between classes
%%     ~~Describe \code{var.equal} here~~
}
  \item{tol}{ The threshold or tolerance value for the covariance and distance
%%     ~~Describe \code{tol} here~~
}
}
\details{
%%  ~~ If necessary, more details than the description above ~~
}
\value{
%%  ~Describe the value returned
%%  If it is a LIST, use
\item{posterior and class}{The posterior possibility and class labels}
%%  \item{comp2 }{Description of 'comp2'}
%% ...
}
\references{

Bingpei Wu, 2012, WMDB 1.0: Discriminant Analysis Methods by Weight Mahalanobis Distance and bayes.

Ito, Y., Srinivasan, C., & Izumi, H. (2006, September). Discriminant analysis by a neural network with Mahalanobis distance. In International Conference on Artificial Neural Networks (pp. 350-360). Springer, Berlin, Heidelberg.

Wölfel, M., & Ekenel, H. K. (2005, September). Feature weighted Mahalanobis distance: improved robustness for Gaussian classifiers. In 2005 13th European signal processing conference (pp. 1-4). IEEE.

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

data(iris)
train=Mabayes(iris[,1:4],iris[,5],TstX= iris[1:10,1:4],tol = 1)

}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.


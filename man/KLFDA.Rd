\name{KLFDA}
\alias{KLFDA}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{Kernel local Fisher discriminant analysis
%%  ~~function to do ... ~~
}
\description{This function performs Kernel Local Fisher Discriminant Analysis. The function provided here allows users to carry out the KLFDA using other input datasets. We used the gaussan matrix as example. Users can input the
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
KLFDA(kdata, y, r,
metric = c("weighted", "orthonormalized", "plain"),
tol=1e-5,knn = 6, reg = 0.001)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{kdata}{The input dataset (kernel matrix). The input data can be a genotype matrix, dataframe, species occurence matrix, or principal components. The dataset have to convert to a kernel matrix before feed into this function.
%%     ~~Describe \code{pcdata} here~~
}
  \item{y}{ The group lables
%%     ~~Describe \code{y} here~~
}
  \item{r}{ Number of reduced features
%%     ~~Describe \code{r} here~~
}
  \item{metric}{Type of metric in the embedding space (default: 'weighted') 'weighted' - weighted eigenvectors 'orthonormalized' - orthonormalized 'plain' - raw eigenvectors
%%     ~~Describe \code{metric} here~~
}
  \item{knn}{The number of nearest neighbours
%%     ~~Describe \code{knn} here~~
}
  \item{tol}{Tolerance to avoid singular values
%%     ~~Describe \code{reg} here~~
}
  \item{reg}{The regularization parameter
%%     ~~Describe \code{reg} here~~
}
}
\details{Kernel Local Fisher Discriminant Analysis for any kernel matrix. It was proposed in Sugiyama, M (2006, 2007) as a non-linear improvement for discriminant analysis.
%%  ~~ If necessary, more details than the description above ~~
}
\value{
 \item{Z}{The reduced features}
\item{Tr}{The transformation matrix}
%% ...
}
\references{
Tang, Y., & Li, W. (2019). lfda: Local Fisher Discriminant Analysis inR. Journal of Open Source Software, 4(39), 1572.

Sugiyama, M (2007). Dimensionality reduction of multimodal labeled data by local Fisher discriminant analysis. Journal of Machine Learning Research, vol.8, 1027-1061.

Sugiyama, M (2006). Local Fisher discriminant analysis for supervised dimensionality reduction. In W. W. Cohen and A. Moore (Eds.), Proceedings of 23rd International Conference on Machine Learning (ICML2006), 905-912.

%% ~put references to the literature/web site here ~
}

\seealso{
KLFDAPC
}
\examples{
kmat <- kmatrixGauss(iris[, -5],sigma=1)
zklfda=KLFDA(kmat, iris[, 5], r=3,metric = "plain",tol=1e-5 )
print(zklfda$Z)
}


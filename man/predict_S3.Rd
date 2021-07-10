\name{predict}
\alias{predict}
\alias{predict.klfda_1}
\alias{predict.KLFDA_mk}
\alias{predict.LDAKPC}
\alias{predict.LFDA}
\alias{predict.LFDAKPC}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{ Predict method in DA for discriminant analysis
%%  ~~function to do ... ~~
}
\description{ Predict method for DA. 
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
\method{predict}{KLFDA_mk}(object,prior,testData, ...)
\method{predict}{klfda_1}(object,prior,testData, ...)
\method{predict}{LDAKPC}(object,prior,testData, ...)
\method{predict}{LFDA}(object,prior,testData, ...)
\method{predict}{LFDAKPC}(object,prior,testData, ...)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{object}{One of the trained object from discriminant analysis
}
\item{prior}{ The weights of the groups.
}
 \item{testData}{The test data or new data 
%%     ~~Describe \code{\dots} here~~
}
  \item{\dots}{Arguments passed to the classifiers
}

}

\value{ The predict function will output the predicted points and their predicted possibility
}
\keyword{ predict }% use one of  RShowDoc("KEYWORDS")

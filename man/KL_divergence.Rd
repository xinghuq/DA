\name{KL_divergence}
\alias{KL_divergence}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{Symmetrised Kullback - Leibler divergence (KL-Divergence)
%%  ~~function to do ... ~~
}
\description{ This function calculates Symmetrised Kullback - Leibler divergence (KL-Divergence) between each class. Designed for KLFDA. 
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
KL_divergence(obj)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{obj}{ The KLFDA object. Users can mdify it to adapt your own purpose.
%%     ~~Describe \code{obj} here~~
}
}
\details{
%%  ~~ If necessary, more details than the description above ~~
}
\value{
  Returns a symmetrised version of the KL divergence between each pair of class
%%  ~Describe the value returned
%%  If it is a LIST, use
%%  \item{comp1 }{Description of 'comp1'}
%%  \item{comp2 }{Description of 'comp2'}
%% ...
}
\references{
Van Erven, T., & Harremos, P. (2014). Renyi divergence and Kullback-Leibler divergence. IEEE Transactions on Information Theory, 60(7), 3797-3820.

Pierre Enel (2019). Kernel Fisher Discriminant Analysis (https://www.github.com/p-enel/MatlabKFDA), GitHub.
%% ~put references to the literature/web site here ~
}
\author{qinxinghu@gmail.com

%%  ~~who you are~~
}
\note{ This function is useful for extimating the loss between reduced features and the original features. It has been adopted in TSNE to determine its projection performance.
%%  ~~further notes~~
}


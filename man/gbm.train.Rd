\name{gbm.train}
\alias{gbm.train}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Train GBM Model
}
\description{
Train GBM Model
}
\usage{
gbm.train(X.train, Surv.train, train.par)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{X.train}{
%%     ~~Describe \code{X.train} here~~
}
  \item{Surv.train}{
%%     ~~Describe \code{Surv.train} here~~
}
  \item{train.par}{
%%     ~~Describe \code{train.par} here~~
}
}
\details{
%%  ~~ If necessary, more details than the description above ~~
}
\value{
%%  ~Describe the value returned
%%  If it is a LIST, use
%%  \item{comp1 }{Description of 'comp1'}
%%  \item{comp2 }{Description of 'comp2'}
%% ...
}
\references{
%% ~put references to the literature/web site here ~
}
\author{
Tai-Hsien Ou Yang
}
\note{
%%  ~~further notes~~
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
%% ~~objects to See Also as \code{\link{help}}, ~~~
}
\examples{
##---- Should be DIRECTLY executable !! ----
##-- ==>  Define data, use random,
##--	or do  help(data=index)  for the standard data sets.

## The function is currently defined as
function (X.train, Surv.train, train.par) 
{
    gbmmodel = gbm.cvrun(Surv.train ~ ., data = X.train, distribution = "coxph", 
        shrinkage = 0.01, n.trees = 1000, interaction.depth = 6, 
        cv.folds = 5, verbose = F, seed = 2013)
    return(gbmmodel)
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line

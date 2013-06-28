\name{cox.predict}
\alias{cox.predict}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Predict Cox Model
}
\description{
Predict Cox Model
}
\usage{
cox.predict(coxmodel, X.test, predict.par)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{coxmodel}{
%%     ~~Describe \code{coxmodel} here~~
}
  \item{X.test}{
%%     ~~Describe \code{X.test} here~~
}
  \item{predict.par}{
%%     ~~Describe \code{predict.par} here~~
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
function (coxmodel, X.test, predict.par) 
{
    predicted = predict(coxmodel, data.frame(X.test))
    return(predicted)
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line

\name{load.clnc}
\alias{load.clnc}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Load Clinical Data
}
\description{
Load Clinical Data
}
\usage{
load.clnc(file.path)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{file.path}{
%%     ~~Describe \code{file.path} here~~
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
function (file.path) 
{
    clinical.raw = read.table(file.path, header = TRUE, sep = "\t")
    rownames(clinical.raw) = clinical.raw[, 1]
    clnc = clinical.raw[, -1]
    cat("Converting Clinical data...\n")
    clnc[, "grade"] = as.numeric(gsub("G", "", clnc[, "grade"]))
    clnc[, "stage"] = gsub("IIIC", "3", clnc[, "stage"])
    clnc[, "stage"] = gsub("IIIB", "3", clnc[, "stage"])
    clnc[, "stage"] = gsub("IIIA", "3", clnc[, "stage"])
    clnc[, "stage"] = gsub("IV", "4", clnc[, "stage"])
    clnc[, "stage"] = gsub("IIB", "2", clnc[, "stage"])
    clnc[, "stage"] = gsub("IIC", "2", clnc[, "stage"])
    clnc[, "stage"] = gsub("IC", "1", clnc[, "stage"])
    clnc[, "stage"] = gsub("IB", "1", clnc[, "stage"])
    clnc[, "stage"] = gsub("IA", "1", clnc[, "stage"])
    clnc[, "stage"] = gsub("[Not Available]", 0, clnc[, "stage"])
    clnc[, "stage"] = as.numeric(clnc[, "stage"])
    return(clnc)
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line

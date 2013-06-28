\name{progModelCV}
\alias{progModelCV}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Cross Validate the performance of Integrated Prognosis Prediction Models
}
\description{
Cross Validate the performance of Integrated Prognosis Prediction Models
}
\usage{
progModelCV(mRNA.ID, cli.ID, sur.ID, model.train, model.predict, ci.list.size, ci.sample.size, xlist = c("mitotic", "mt", "ls"), clist, cv.fold = 100, train.par = NULL, predict.par = NULL)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{mRNA.ID}{
%%     ~~Describe \code{mRNA.ID} here~~
}
  \item{cli.ID}{
%%     ~~Describe \code{cli.ID} here~~
}
  \item{sur.ID}{
%%     ~~Describe \code{sur.ID} here~~
}
  \item{model.train}{
%%     ~~Describe \code{model.train} here~~
}
  \item{model.predict}{
%%     ~~Describe \code{model.predict} here~~
}
  \item{ci.list.size}{
%%     ~~Describe \code{ci.list.size} here~~
}
  \item{ci.sample.size}{
%%     ~~Describe \code{ci.sample.size} here~~
}
  \item{xlist}{
%%     ~~Describe \code{xlist} here~~
}
  \item{clist}{
%%     ~~Describe \code{clist} here~~
}
  \item{cv.fold}{
%%     ~~Describe \code{cv.fold} here~~
}
  \item{train.par}{
%%     ~~Describe \code{train.par} here~~
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

synapseLogin()

lusc.cv=progModelCV("syn1710382", "syn1715826", "syn1710384", gbm.train, gbm.predict, 100, 100, xlist = c("mitotic"), clist=c("age","stage"), cv.fold = 100)
lusc.cv=progModelCV("syn1710382", "syn1715826", "syn1710384", gbm.train, gbm.predict, 100, 100, xlist = c("mt","ls"), clist=c("age","stage"), cv.fold = 100)
lusc.cv=progModelCV("syn1710382", "syn1715826", "syn1710384", gbm.train, gbm.predict, 100, 100, xlist = NULL, clist=c("age","stage"), cv.fold = 100)
lusc.cv=progModelCV("syn1710382", "syn1715826", "syn1710384", single.train, single.predict, 100, 100, xlist = c("mitotic"), clist=NULL, cv.fold = 100)


ci.submit<-list("gisl[modelname]"="Mitotic Only","gisl[ci]"=round(lusc.cv["mean"],6),"gisl[std]"=round(lusc.cv["std"],6),"gisl[ctype]"="LUSC","gisl[date]"="04032013","gisl[desc]"="Mitotic Only Tai-Hsien"  )
postToHost("128.59.65.84","/gisls",ci.submit,port=2013)


## The function is currently defined as
function (mRNA.ID, cli.ID, sur.ID, model.train, model.predict, 
    ci.list.size, ci.sample.size, xlist = c("mitotic", "mt", 
        "ls"), clist, cv.fold = 100, train.par = NULL, predict.par = NULL) 
{
    cat("Downloading Synapse Entities...\n")
    syn = loadEntity(mRNA.ID)
    cli = loadEntity(cli.ID)
    sur = loadEntity(sur.ID)
    cat("Loading mRNA Datasets...\n")
    ge = load.exp(file.path(syn$cacheDir, syn$files[[1]][1]))
    ge = t(ge)
    cat("Imputing Expression data...\n")
    for (i in 1:nrow(ge)) {
        ge[i, is.na(ge[i, ])] <- mean(ge[i, ], na.rm = TRUE)
    }
    rn.ge = rownames(ge)
    rn.ge = substr(rn.ge, 6, nchar(rn.ge))
    rownames(ge) = rn.ge
    data(attractome.minimalist)
    rn.ge = as.matrix(rn.ge)
    dim(rn.ge) = c(length(rn.ge), 1)
    colnames(rn.ge) = "Gene.Symbol"
    rownames(rn.ge) = rn.ge
    metagene = CreateMetageneSpace(ge, attractome.minimalist, 
        rn.ge)$metaSpace
    cat("Loading Clinical Data...\n")
    clnc = load.clnc(file.path(cli$cacheDir, cli$files[[1]][1]))
    clnc.imputed = lazyImputeDFClncOslo(clnc)
    cat("Loading Survival Data...\n")
    survival.ge = load.exp(file.path(sur$cacheDir, sur$files[[1]][1]))
    ci.list = matrix(0, 1:ci.list.size)
    cat("Generating Cross Validation Sample Sets...\n")
    sample.id = replicate(ci.list.size, sample(1:ncol(ge), ncol(ge), 
        replace = TRUE))
    test.id = replicate(ci.list.size, sample(1:ncol(ge), ncol(ge), 
        replace = TRUE))
    cat("Cross validating...\n")
    for (i in 1:cv.fold) {
        sample.list = colnames(ge)[sample.id[1:ci.sample.size, 
            i]]
        Surv.train = Surv(survival.ge[sample.list, "OS_OS"], 
            survival.ge[sample.list, "OS_vital_status"])
        X.train = cbind(t(metagene[xlist, sample.list]), clnc.imputed[sample.list, 
            clist])
        trainmodel <- match.fun(model.train)
        trainedmodel = trainmodel(X.train, Surv.train, train.par)
        X.test = cbind(t(metagene[xlist, test.id[, i]]), clnc.imputed[test.id[, 
            i], clist])
        predictmodel <- match.fun(model.predict)
        p1 = predictmodel(trainedmodel, X.test, predict.par)
        ci.list[i] = concordance.index(p1, survival.ge[test.id[, 
            i], "OS_OS"], survival.ge[test.id[, i], "OS_vital_status"])$c.index
        if (i\%\%100 == 0) {
            setTxtProgressBar(txtProgressBar(style = 3), i/cv.fold)
        }
    }
    ci.mean = mean(ci.list, na.rm = T)
    ci.std = sd(ci.list, na.rm = T)
    cat("Done!\n")
    cat("MEAN=", ci.mean, "STD=", ci.std, "\n")
    return(c(mean = ci.mean, std = ci.std))
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line

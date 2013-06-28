\name{progEval-package}
\alias{progEval-package}
\alias{progEval}
\docType{package}
\title{ Integrated Prognosis Prediction Model Cross Validation Tool
}
\description{
Prognosis Prediction Model Cross Validation Tool
}
\details{
\tabular{ll}{
Package: \tab progEval\cr
Type: \tab Package\cr
Version: \tab 1.0\cr
Date: \tab 2013-03-29\cr
License: \tab GPL \cr
}



}
\author{
Tai-Hsien Ou Yang
Maintainer: Tai-Hsien Ou Yang <to2232@columbia.edu>

}
\references{
DREAMBOX7
Sage Bionetworks SynapseClient
}
\keyword{ prognosis }
\seealso{
DREAMBOX7 by Weiyi Cheng
\code{\link[DreamBox7]{https://github.com/weiyi-bitw/DreamBox7}}
}
\examples{
aiccox.ov<-progModelCV( "syn1710361", "syn1715828", "syn1710363", aiccox.train, aiccox.predict, ci.list.size=100, ci.sample.size=100, xlist=c( "mitotic","mt","ls" ), clist=c("age", "grade", "stage"), cv.fold=100)
gbm.ov<-progModelCV( "syn1710361", "syn1715828", "syn1710363", gbm.train, gbm.predict, ci.list.size=100, ci.sample.size=100, xlist=c( "mitotic","mt","ls" ), clist=c("age", "grade", "stage"), cv.fold=100)
}

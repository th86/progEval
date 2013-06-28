#GBM MODEL
#Tai-Hsien Ou Yang 
#03292013

gbm.train<-function ( X.train, Surv.train,  train.par ) {
	gbmmodel = gbm.cvrun(Surv.train ~ ., data = X.train, distribution = "coxph", 
        shrinkage = 0.01, n.trees =1000, interaction.depth = 6, 
        cv.folds = 5, verbose = F, seed = 2013)
return(gbmmodel)
}

gbm.predict<-function ( gbmmodel, X.test, predict.par ) {
	best.iter = gbm.perf(gbmmodel, method = "cv", plot.it = FALSE)
	predicted=predict.gbm(gbmmodel, X.test, best.iter)
return(predicted)
}








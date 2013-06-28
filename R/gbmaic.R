#GBM-AIC MODEL
#Tai-Hsien Ou Yang 
#03302013

gbmaic.train<-function ( X.train, Surv.train,  train.par ) {

	upper = terms(Surv.train ~ (.), data = X.train)
    	coxmodel = step(coxph(Surv.train ~ ., data = X.train), scope = upper, direction = "both", k = 2, trace = FALSE)

	gbmmodel = gbm.cvrun(Surv.train ~ ., data = X.train, distribution = "coxph", 
        shrinkage = 0.01, n.trees =1000, interaction.depth = 6, 
        cv.folds = 5, verbose = F, seed = 2013)
	
	model=list(aic=coxmodel, gbm=gbmmodel )
return(model)
}

gbmaic.predict<-function ( model, X.test, predict.par ) {

	predicted.aic = predict(model$aic, data.frame(X.test))

	best.iter = gbm.perf(model$gbm, method = "cv", plot.it = FALSE)
	predicted.gbm =predict.gbm(model$gbm, X.test, best.iter)
	predicted= predicted.aic + predicted.gbm
return(predicted)
}








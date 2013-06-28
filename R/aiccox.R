#AIC-COX MODEL
#Tai-Hsien Ou Yang 
#03292013

aiccox.train<-function ( X.train, Surv.train,  train.par ) {
	upper = terms(Surv.train ~ (.), data = X.train)
    	coxmodel = step(coxph(Surv.train ~ ., data = X.train), scope = upper, direction = "both", k = 2, trace = FALSE)
return(coxmodel)
}

aiccox.predict<-function ( coxmodel, X.test, predict.par ) {
	predicted = predict(coxmodel, data.frame(X.test))
return(predicted)
}


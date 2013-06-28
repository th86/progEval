#SIMPLE-COX MODEL
#Tai-Hsien Ou Yang 
#03292013

cox.train<-function ( X.train, Surv.train,  train.par ) {
	coxmodel=coxph(Surv.train ~ ., data = X.train)
return(coxmodel)
}

cox.predict<-function ( coxmodel, X.test, predict.par ) {
	predicted = predict(coxmodel, data.frame(X.test))
return(predicted)
}


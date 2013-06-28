#Prognosis Model Prediction Cross Validation Package
#Tai-Hsien Ou Yang
#03242013

library(synapseClient)
library(DreamBox7)
library(survival)
library(survcomp)  

source(".\\R\\loadclnc.R")
source(".\\R\\direct.R")
source(".\\R\\simplecox.R")
source(".\\R\\aiccox.R")
source(".\\R\\gbm.R")
source(".\\R\\gbmaic.R")
#Login Synapse
synapseLogin()


progModelCV <-function( mRNA.ID, cli.ID, sur.ID,  train.ID, test.ID   ,  model.train, model.predict,  xlist=c( "mitotic","mt","ls" ), clist=NULL, glist=NULL , cv.fold=100, train.par=NULL, predict.par=NULL ) {


	#Download entities
	cat("Downloading Synapse Entities...\n")
	syn = loadEntity(mRNA.ID) 
	cli = loadEntity(cli.ID) 
	sur = loadEntity(sur.ID) 
	trainset = loadEntity(train.ID)
	testset = loadEntity(test.ID)


	#load the matrix and transpose 
	cat("Loading mRNA Datasets...\n")
	ge = load.exp(file.path(syn$cacheDir, syn$files[[1]][1]))
	ge=t(ge)

	cat("Imputing Expression data...\n")
	#mean imputation
	for(i in 1:nrow(ge)){
 		ge[i, is.na(ge[i,] )] <- mean(ge[i,], na.rm = TRUE) 
 		#ge[i, is.nan(ge[i,] )] <- mean(ge[i,], na.rm = TRUE) 
	}

	#remove prefix
	rn.ge=rownames(ge)
	rn.ge=substr(rn.ge,6,nchar(rn.ge))



	rn.ge = sapply(rn.ge, 
		function(x){	
		if(regexpr("\\?", x) > 0){
			o = strsplit(x, "\\|")[[1]][2]
		}else{
			o = strsplit(x, "\\|")[[1]][1]
		}
		return (o)
		}
	)

	rownames(ge)=rn.ge 
	
	#build gene symbol map and generate metagene
	data(attractome.minimalist)
	rn.ge=as.matrix(rn.ge)
	dim(rn.ge)=c(length(rn.ge),1)
	colnames(rn.ge)="Gene.Symbol"
	rownames(rn.ge)=rn.ge
	metagene=CreateMetageneSpace(ge,attractome.minimalist, rn.ge)$metaSpace
	
	#load clinical
	cat("Loading Clinical Data...\n")
	
	clnc = load.clnc(file.path(cli$cacheDir, cli$files[[1]][1])) 

	clnc.imputed = lazyImputeDFClncOslo(clnc)


	#condition
	if( "mtStage3" %in% xlist   ){
	 stagec=clnc.imputed[ ,"stage" ]
	 mt=metagene["mt",]
	 mt[which(stagec>3)]=0
	 mt.3.med=median(mt[which(stagec>0)]	)
	 cat("MT<4 Median=",mt.3.med,"\n")
	 mt[stagec>0]=mt[stagec>0]-mt.3.med

	 metagene=rbind(metagene,mtStage3=mt)
	}

	if( "lsxStage3" %in% xlist   ){
	 stagec=clnc.imputed[ , "stage" ]
	 ls=metagene["ls",]
         lsxStage3=ls*(3-stagec)
	 lsxStage3= lsxStage3- median(lsxStage3)

	 metagene=rbind(metagene, lsxStage3=lsxStage3)
	}




	#load survival
	cat("Loading Survival Data...\n")
	survival.ge = load.exp(file.path(sur$cacheDir, sur$files[[1]][1]))

	#cross validation CI

	ci.list=matrix(0,1:cv.fold)
	cat("Generating Cross Validation Sample Sets...\n")
	#sample.id=replicate(ci.list.size,sample(1:ncol(ge),ncol(ge),replace = TRUE))
	#test.id=replicate(ci.list.size,sample(1:ncol(ge),ncol(ge),replace = TRUE))



	sample.id=as.matrix(read.table( file.path(trainset$cacheDir, trainset$files[[1]][1]) ,header=F))
	test.id=as.matrix(read.table( file.path(testset$cacheDir, testset$files[[1]][1]) ,header=F))

	rownames( survival.ge)= colnames(ge) 

	#xlist=c( "mitotic","mt","ls" )
	#clist=c("age", "grade", "stage")

	cat("Cross validating...\n")


	for(i in 1:cv.fold){
		sample.list=sample.id[,i]

		Surv.train=Surv(survival.ge[sample.list,"OS_OS"], survival.ge[sample.list,"OS_vital_status"])

		#deal with the singular feature
		metagenelist=t(metagene[xlist,sample.list] )
		if(length(xlist)==1){
			names(metagenelist)=  sample.list
			metagenelist=metagenelist[1,]
		}


		X.train=cbind(metagenelist,t(ge[glist,sample.list]) , clnc.imputed[sample.list,clist])
	

		trainmodel <- match.fun(model.train) 
		trainedmodel=trainmodel(X.train, Surv.train, train.par  ) 		
		
		metagenelist=t(metagene[xlist,test.id[,i]] )
		
		if(length(xlist)==1){
			names(metagenelist)=test.id[,i]
			metagenelist=metagenelist[1,]
		}
		X.test=cbind(metagenelist,t(ge[glist,test.id[,i]]) , clnc.imputed[test.id[,i],clist] )

		predictmodel <- match.fun(model.predict) 
		p1=predictmodel( trainedmodel,  X.test, predict.par  )

		#ci.list[i]=concordance.index(p1,survival.ge[test.id[,i],"OS_OS"], survival.ge[test.id[,i],"OS_vital_status"])$c.index
		ci.list[i]=getCCDIdx(p1,cbind(survival.ge[test.id[,i],"OS_OS"], survival.ge[test.id[,i],"OS_vital_status"]))
	
 		if (i%%100 == 0) {
            	setTxtProgressBar(txtProgressBar(style = 3), i/cv.fold)
        	}
		 
	}


ci.mean=mean(ci.list, na.rm=T)
ci.std=sd(ci.list, na.rm=T)

cat("Done!\n")
cat("MEAN=",ci.mean,"STD=",ci.std,"\n")

return(c(mean=ci.mean, std=ci.std)    )

}


evalMatrixCI<- function (ge, surv ) {
	
	surv=surv[colnames(ge),]
	cilist=matrix(0,nrow(ge),1)

	for(i in 1:nrow(ge) ){
		cilist[i,1]= getCCDIdx(ge[i, ],surv   ) 
	}
	rownames(cilist)=rownames(ge)

	return(cilist)
}



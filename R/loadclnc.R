#Load Clinical Data
#Tai-Hsien Ou Yang
#04012013

require('IlluminaHumanMethylation27k.db')

if(!require(IlluminaHumanMethylation27k.db)){
	source("http://bioconductor.org/biocLite.R")
	biocLite("IlluminaHumanMethylation27k.db")
#  library(IlluminaHumanMethylation27k.db)
}



load.clnc<-function( file.path ) {

 	clinical.raw =read.table( file.path ,header = TRUE,sep="\t")
	rownames(clinical.raw)=clinical.raw[,1]
	clnc=clinical.raw[,-1]
	cat("Converting Clinical data...\n")
	if ( "grade" %in% colnames(clnc) ){
		clnc[,"grade"]=as.numeric(gsub("G", "", clnc[,"grade"]) )
	}
	
	if ( "stage" %in% colnames(clnc) ){
		clnc[,"stage"]=gsub("Stage", "", clnc[,"stage"]) 
		clnc[,"stage"]=gsub("IIIC", "3.6", clnc[,"stage"]) 
		clnc[,"stage"]=gsub("IIIB", "3.3", clnc[,"stage"]) 
		clnc[,"stage"]=gsub("IIIA", "3.0", clnc[,"stage"])
		clnc[,"stage"]=gsub("III", "3.1", clnc[,"stage"]) 
		clnc[,"stage"]=gsub("IV", "4.1", clnc[,"stage"])
		clnc[,"stage"]=gsub("IVC", "4.6", clnc[,"stage"])  
		clnc[,"stage"]=gsub("IVB", "4.3", clnc[,"stage"]) 
		clnc[,"stage"]=gsub("IVA", "4.0", clnc[,"stage"])   
		clnc[,"stage"]=gsub("IIC", "2.6", clnc[,"stage"]) 
		clnc[,"stage"]=gsub("IIB", "2.3", clnc[,"stage"]) 
		clnc[,"stage"]=gsub("IIA", "2.0", clnc[,"stage"]) 
		clnc[,"stage"]=gsub("II", "2.1", clnc[,"stage"]) 
		clnc[,"stage"]=gsub("IC", "1.6", clnc[,"stage"]) 
		clnc[,"stage"]=gsub("IB", "1.3", clnc[,"stage"]) 
		clnc[,"stage"]=gsub("IA", "1.0", clnc[,"stage"]) 
		clnc[,"stage"]=gsub("I", "1.1", clnc[,"stage"]) 
		clnc[,"stage"]=gsub("[Not Available]", 0, clnc[,"stage"]) 
		clnc[,"stage"]=as.numeric( clnc[,"stage"] )

	}


return(clnc)
}



preproc.exp<-function(file.path, tr=T){

ge = load.exp(file.path)

if(tr==T)
{
ge=t(ge)
 }

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

return(ge)
}


loadAllHTD<-function (ctype, methfile=T, rppafile=T){

meth=NULL
rppa=NULL

ge = preproc.exp(   paste(ctype, "_mRNA_core.txt",sep="")   )
clnc = load.clnc( paste(ctype, "_clinical_core.txt",sep="") )
surv = load.exp( paste(ctype, "_OS_core.txt",sep="") )
mi = preproc.exp(   paste(ctype, "_miRNA_core.txt",sep="")   )
cnv = preproc.exp(   paste(ctype, "_CNV_core.txt",sep="")   )

rownames(cnv) =colnames(load.exp(   paste(ctype, "_CNV_core.txt",sep="")   ))

shared=intersect(colnames(ge),rownames(clnc))

ge=ge[, shared]
ge.rn=rownames(ge)
clnc=clnc[shared,]
surv=surv[shared,]

mi=mi[, shared]
cnv=cnv[, shared]

if( rppafile==T ){
	rppa = preproc.exp(   paste(ctype, "_RPPA_core.txt",sep="")   )
	rppa=rppa[, shared]
}

if( methfile==T ){
	meth = preproc.exp(   paste(ctype, "_methylation_core.txt",sep="")   )
	meth=meth[, shared]
	x <- IlluminaHumanMethylation27kSYMBOL
	mapped_probes <- mappedkeys(x)
	xx <- as.list(x[mapped_probes])
	xxx=unlist(xx)
	rn.meth=substr(rownames(meth),8,17)
	rn=xxx[rn.meth]
	rn[is.na(rn)]= rownames(meth)[is.na(rn)]
	rownames(meth) = rn
}

o=list( mrna=ge, cnv=cnv, meth=meth, mirna=mi, rppa=rppa, clinical=clnc, surv=surv   )

return(o)
}



#From WY Cheng's DreamBox7
load.exp<-function (file, sep = "\t") 
{
    line = readLines(file)
    tokens = strsplit(line[1], "\t")[[1]]
    n = length(tokens) - 1
    m = length(line) - 1
    mset = matrix(NA, m, n)
    cname = tokens[2:(n + 1)]
    rname = rep(NA, m)
    b = txtProgressBar(style = 3)
    for (i in 1:m) {
        tokens = strsplit(line[i + 1], "\t")[[1]]
        tokens[tokens == "NA"] = NA
        rname[i] = tokens[1]
        mset[i, ] = as.numeric(tokens[2:(n + 1)])
        if (i%%100 == 0) {
            setTxtProgressBar(b, i/m)
        }
    }
    cat("\n")
    colnames(mset) <- cname
    rownames(mset) <- rname
    mset
}


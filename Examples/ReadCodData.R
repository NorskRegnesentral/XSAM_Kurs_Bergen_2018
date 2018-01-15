#Specialized code for reading NSS herring data. Structure inherited from WGWIDE history...
#

#---------------------------------------------------------#
#First some functions specific for the XSA files
#---------------------------------------------------------#
readXSAfiles0<-function(file){
	ncol <- max(count.fields(file, sep = ""))
	tab<-read.table(file,fill=T,flush=T,as.is=T,header=F,col.names = paste0("V", seq_len(ncol)))
	#k0<<-tab
	options(warn=-1)
	FileType<-CheckFileType(tab)
	options(warn=0)
	#print(FileType)
	if(FileType=="NumTab"){
		out<-readXSAtab1(tab)
	}
	if(FileType=="FleetFile"){
		out<-readXSAtab2(tab)
	}
	if(FileType=="OtherFile"){
		cat("This function cannot read this file!!\n")
		out<-NULL
	}
	out
}

CheckFileType<-function(tab){
	nrows<-nrow(tab)
	ntextrows<-sum(is.na(as.numeric(tab[,1])))
	nnumrows<-nrows-ntextrows
	if(ntextrows==1)FileType<-"NumTab"
	if(ntextrows>1 & ntextrows<nnumrows)FileType<-"FleetFile"
	if(ntextrows>nnumrows)FileType<-"OtherFile"
	FileType
}

readXSAtab1<-function(tab,na.rm=T){
	desc<-paste(tab[1,][!is.na(tab[1,])],collapse=" ")
	print(paste("Reads:", desc,collapse=" "))
	options(warn=-1)
	years<-as.numeric(tab[3,])
	ages<-as.numeric(tab[4,])
	code<-as.numeric(tab[5,])
	options(warn=0)
	years<-years[!is.na(years)]
	ages<-ages[!is.na(ages)]
	code<-code[!is.na(code)]
	years<-years[1]:years[2]
	ages<-ages[1]:ages[2]
	nlines<-length(years)
	if(code==5){
		ncols<-1
		ages<-1
	}
	else ncols<-length(ages)
	#print(years)
	#print(ages)
	out<-tab[6:(6+nlines-1),1:ncols,drop=F]
	#Convert to numeric (by column since as.numeric ikke funker på data.frame...)
	out<-apply(out,MAR=2,FUN=as.numeric)
	dimnames(out)[[2]]<-ages
	dimnames(out)[[1]]<-years
	if(na.rm)out[out<0]<-NA
	attributes(out)$ages<-ages
	attributes(out)$years<-years
	out
}

readXSAtab2<-function(tab,na.rm=T){
	desc<-paste(tab[1,][!is.na(tab[1,])],collapse=" ")
	print(paste("Reads:", desc,collapse=" "))
	options(warn=-1)
	nfleets<-as.numeric(tab[2,])
	options(warn=0)
	nfleets<-nfleets[!is.na(nfleets)]
	nfleets<-nfleets-100
	print(paste("File contains:",nfleets,"fleets.",collapse=" "))
	options(warn=-1)
	pointer<-3
	names<-c()
	out<-list()
	for(i in 1:nfleets){
		names[i]<-paste(tab[pointer,][!is.na(tab[pointer,])],collapse=" ")
		pointer<-pointer+1
		years<-as.numeric(tab[pointer,])
		pointer<-pointer+1
		sampleTime<-as.numeric(tab[pointer,])
		pointer<-pointer+1
		ages<-as.numeric(tab[pointer,])
		pointer<-pointer+1
		#
		years<-years[!is.na(years)]
		sampleTime<-sampleTime[!is.na(sampleTime)]
		ages<-ages[!is.na(ages)]
		years<-years[1]:years[2]
		ages<-ages[1]:ages[2]
		#print(years);print(sampleTime);print(ages)
		out[[i]]<-tab[pointer:(pointer+length(years)-1),2:(length(ages)+1),drop=F]
		#Convert to numeric (by column since as.numeric ikke funker på data.frame...)
		out[[i]]<-apply(out[[i]],MAR=2,FUN=as.numeric)
		pointer<-pointer+length(years)

		dimnames(out[[i]])[[2]]<-ages
		dimnames(out[[i]])[[1]]<-years
		if(na.rm)out[[i]][out[[i]]<0]<-NA
		attributes(out[[i]])$ages<-ages
		attributes(out[[i]])$years<-years
		attributes(out[[i]])$sampleTime<-sampleTime[3]

	}
	names(out)<-names
	out
}


AddData<-function(dobj,years,val=NULL,selectyr=max(dimnames(dobj)[[1]])){
	#dobj: dataobject as returned by readXSAfiles0
	#years: Which year(s) to add to dobj
	#val: which value to add. OBS: If selectyear>1, then val must correpsond to a matrix with correct dims.
	#If is.null(.) Then
	#select: average from years given by select. Default is nrow(dobj), i.e. last year
	dimnames01<-dimnames(dobj)[[1]]
	dimnames02<-dimnames(dobj)[[2]]
	if(is.null(val)){
		tmp<-dobj[paste(selectyr),,drop=F]
		tmp<-colMeans(tmp)
		for(i in 1:length(years)){
			dobj<-rbind(dobj,tmp)
		}
	}
	else dobj<-rbind(dobj,val)
	dimnames(dobj)[[1]]<-paste(c(dimnames01,years))
	attributes(dobj)$years<-as.numeric(dimnames(dobj)[[1]])
	dimnames(dobj)[[2]]<-dimnames02
	attributes(dobj)$ages<-as.numeric(dimnames02)

	dobj
}

#---------------------------------------------------------#
#End XSA reading functions
#---------------------------------------------------------#
#Creates a datalist on a format that can be used by XSAM utility functions that creates data objsects that can be passed on to the XSAM template.


#Set up paths and create datalist
#path<-"C:/Users/aanes/Documents/XSAMcourse/
mypath<-getwd()
datapath<-paste(mypath,"/Data/cod/xsa/",sep="")

files<-list.files(datapath)	
files

datalist<-list()

datalist$caa<-readXSAfiles0(paste(datapath,files[2],sep=""))
attributes(datalist$caa)$k<-1

datalist$SurveyIndices<-readXSAfiles0(paste(datapath,files[4],sep=""))
for(i in 1:length(datalist$SurveyIndices))attributes(datalist$SurveyIndices[[i]])$k<-1

datalist$matprop<-readXSAfiles0(paste(datapath,files[8],sep=""))
datalist$west<-readXSAfiles0(paste(datapath,files[12],sep=""))
datalist$weca<-readXSAfiles0(paste(datapath,files[11],sep=""))
datalist$caton<-readXSAfiles0(paste(datapath,files[3],sep=""))
attributes(datalist$caton)$k<-1

datalist$natMor<-readXSAfiles0(paste(datapath,files[10],sep=""))
datalist$propF<-readXSAfiles0(paste(datapath,files[5],sep=""))
datalist$propM<-readXSAfiles0(paste(datapath,files[9],sep=""))


CForecast<-data.frame(y=as.numeric(rownames(datalist$caton)),p=datalist$caton,a=datalist$caton)
names(CForecast)<-c("y","p","a")
CForecast<-rbind(CForecast,c(2016,NA,NA))
CForecast<-rbind(CForecast,c(2017,NA,NA))

CForecast$p
CatchPrediction<-data.frame(Prediction=CForecast$p,caton=CForecast$a,se=NA,rse=0.05)
rownames(CatchPrediction)<-CForecast$y

###
CatchPrediction["2016","Prediction"]<-894000
CatchPrediction["2016","rse"]<-0.05
#CatchPrediction[,"se"]<-sqrt(CatchPrediction[,"rse"]^2+1)
CatchPrediction[,"se"]<-sqrt(log(CatchPrediction[,"rse"]^2+1))

CatchPrediction["2017","Prediction"]<-805000
CatchPrediction["2017","rse"]<-0.05
CatchPrediction[,"se"]<-sqrt(log(CatchPrediction[,"rse"]^2+1))

datalist$CatchPrediction<-CatchPrediction
attributes(datalist$CatchPrediction)$k<-1


######-----------------------------------------####
#OOOOOOOOOOBBBBBBBBBBBBBBBBSSSSSSSSSSS!!!!!!
datalist$SurveyIndices[[1]][is.na(datalist$SurveyIndices[[1]])]<-0

datalist$SurveyIndices[[2]][is.na(datalist$SurveyIndices[[2]])]<-0

datalist$SurveyIndices[[3]]

datalist$SurveyIndices[[4]]
###########################################################################################

#ADD DATA FOR 2016 BY USING AVERAGES OF LAST 3 YEARS
datalist$matprop<-AddData(dobj=datalist$matprop,years=2016,val=NULL,selectyr=2013:2015)
datalist$west<-AddData(dobj=datalist$west,years=2016,val=NULL,selectyr=2013:2015)
datalist$natMor<-AddData(dobj=datalist$natMor,years=2016,val=NULL,selectyr=2013:2015)
datalist$propF<-AddData(dobj=datalist$propF,years=2016,val=NULL,selectyr=2013:2015)
datalist$propM<-AddData(dobj=datalist$propM,years=2016,val=NULL,selectyr=2013:2015)

########################################




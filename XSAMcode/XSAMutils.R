ReadXSAMsettings<-function(file="C:/Users/aanes/Documents/WGWIDE/XSAM/XSAMsetup.txt",CheckFile=T){
	variables<-c(
	"agerange",
	"yearrange",
	"maxAgePlusGroup",
	"stockRecruitmentModelCode",
	"CatchConstraint",
	"corFlaglogF",
	"LatentEffort",
	"TSsel",
	"corFlagU",
	"EstimateARInterceptEffort",
	"EstimateUInterceptEffort",
	"ymeantype",
	"am",
	"Am",
	"AT",
	"RecruitmentProcess",
	"keyLogFconst",
	"keyVarF",
	"fleetIndex",
	"Modelq",
	"keyLogqpar",
	"keyVarObs",
	"Fbarrange",
	"agerangeC",
	"agerangeI",
	"UseCatchPred"
	)
	tt<-scan(file,what="",sep="\n",comment.char = "#",flush=T,quiet=T)
	#------------------------------------------#
	#Check that required headers are present
	#------------------------------------------#
	fileheaders<-tt[tt%in%variables]
	missingheaders<-variables[!variables%in%fileheaders]
	if(length(missingheaders)>0)stop(paste("Required headers:<",paste(missingheaders,collapse=", "),"> are missing from input file",collapse=" "))
	#------------------------------------------#
	#Check for redundant headers
	#------------------------------------------#
	fileheaders<-lapply(tt,FUN=function(x){x<-strsplit(x,split=" ");x[[1]][[1]][1]})
	fileheaders<-lapply(fileheaders,FUN=function(x){if(x=="NA"){x<--999};x})
	options(warn=-1)
	fileheaders<-unlist(lapply(fileheaders,FUN=function(x)as.numeric(x)))
	options(warn=0)
	fileheaders<-tt[is.na(fileheaders)]
	redunantheaders<-fileheaders[!fileheaders%in%variables]
	if(length(redunantheaders)>0)stop(paste("Headers:<",paste(redunantheaders,collapse=", "),"> found in input file are redundant for XSAM and results in erroneous interpretation of model setup. \nNote that this function does not accept leading blanks for numerical values as they will be interpreted as a header!\nUse leading # to insert comments",collapse=" "))
	#Must sort variables according to ordering in file
	#h0<<-variables
	#h1<<-fileheaders
	#variables<-variables[match(variables,fileheaders)]
	variables<-variables[match(fileheaders,variables)]
	#h2<<-variables
	headerpositions<-match(variables,tt)
	#print(headerpositions)
	#k0<<-tt
	nlines<-diff(c(headerpositions,length(tt)+1))-1
	#print(nlines)
	out<-list()
	for(i in 1:length(variables)){
		#print(variables[i])
		#print(nlines[i])
		out[[i]]<-NULL
		ncols<-c()
		for(k in 1:nlines[i]){
			tmp<-unlist(strsplit(tt[headerpositions[i]+k],split=" "))
			tmp[tmp=="NA"]<-"-9999"
			tmp<-as.numeric(tmp)
			tmp<-tmp[!is.na(tmp)]
			ncols[k]<-length(tmp)
			if(k>1){
				if(ncols[k]!=ncols[k-1])stop(paste("Different number of columns for",variables[i],sep=""))
				out[[i]]<-rbind(out[[i]],tmp)
			}
			else out[[i]]<-tmp
		}
		if(k>1)dimnames(out[[i]])[[1]]<-NULL
		out[[i]][out[[i]]==-9999]<-NA
		names(out)[i]<-variables[i]
	}
	if(is.null(dim(out$keyLogqpar)))out$keyLogqpar<-matrix(out$keyLogqpar,nrow=1)
	if(is.null(dim(out$agerangeI)))out$agerangeI<-matrix(out$agerangeI,nrow=1)
	if(CheckFile)CheckXSAMsettings(out)
	out
}
#tmp<-ReadXSAMsettings()

#CheckOrdering(tmp$keyLogFconst)

#CheckOrdering(c(1:9,8))
#diff(c(1:9,8))

CheckOrdering<-function(vec,rmvals=T){
	#Check whether first value is 0 or 1 and then ordering of vector
	#If rmvals NAs and negative values are removed before checking
	n<-length(vec)
	error<-0
	if(rmvals){
		vec<-vec[!is.na(vec) & vec>=0]
	}
	x0<<-vec
	vec<-unique(vec)
	#Test 1
	if(!vec[1]%in%c(0,1))error<-1
	#Test 2
	if(sum(diff(vec))!=(length(vec)-1))error<-2
	start<-vec[1]
	stopp<-max(vec)#vec[length(vec)]
	#Test 3
	if(stopp>n)error<-3
	#Test4; the unique operation remove bits of vec that may be descending after reached maximum
	if(sum(diff(x0)<0)>0)error<-4
	
	list(error=error,start=start,stopp=stopp)
}


CheckXSAMsettings<-function(settings){
	n<-0
	#Check correspondance in age ranges
	nages<-diff(settings$agerange)+1
	dimF<-length(settings$keyLogFconst)
	dimq<-dim(settings$keyLogqpar)[2]
	dimvarobs<-dim(settings$keyVarObs)[2]
	if(dimF!=nages){warning("Age span in keyLogFconst is different than age span: May result in crash or erroneous results!");n<-n+1}
	if(dimq!=nages){warning("Age span in keyLogqpar is different than age span: May result in crash or erroneous results!");n<-n+1}
	if(dimvarobs!=nages){warning("Age span in keyVarObs is different than age span: May result in crash or erroneous results!");n<-n+1}
	#
	dimF<-length(settings$keyVarF)
	if(dimF!=nages){warning("Age span in keyVarF is different than age span: May result in crash or erroneous results!");n<-n+1}
	#
	nIndices<-length(settings$fleetIndex)-1
	dimq<-dim(settings$keyLogqpar)[1]
	if(nIndices>dimq){warning("Number of indices is larger than number of q vectors: Will result in crash!");n<-n+1}
	if(nIndices<dimq){warning("Number of indices is smaller than number of q vectors: They should correspond!");n<-n+1}
	dimq<-dim(settings$agerangeI)[1]
	if(nIndices>dimq){warning("Number of indices is larger than corresponding number of age spans provided: Will result in crash!");n<-n+1}	
	if(nIndices<dimq){warning("Number of indices is smaller than corresponding number of age spans provided: They should correspond!");n<-n+1}
	uniqueF<-unique(settings$keyLogFconst)
	#Check consistency in age settings
	minA<-settings$agerange[1]
	maxA<-settings$agerange[2]
	if(settings$Am>maxA){warning("Am must be <=maximum A. Will result in crash!");n<-n+1}
	if(settings$am>settings$Am){warning("am must be <=Am. Will result in crash!");n<-n+1}
	if(settings$AT<maxA){warning("AT must be >=maximum A. Will result in crash!");n<-n+1}
	tmptest<-CheckOrdering(settings$keyLogFconst)
	if(tmptest$error==1){warning("Start value of keyLogFconst should be 0 or 1. Will result in crash!");n<-n+1}
	if(tmptest$error%in%c(2,4)){warning("Ordering of keyLogFconst is not chronological. Is this correct? May cause undesireable results if not deliberately set!");n<-n+1}
	if(tmptest$error==3){warning("Largest value of index is larger than the dimension of keyLogFconst. Will result in crash!");n<-n+1}
	start1<-tmptest$start
	tmptest<-CheckOrdering(settings$keyVarF)
	if(tmptest$error==1){warning("Start value of keyVarF should be 0 or 1. Will result in crash!");n<-n+1}
	if(tmptest$error%in%c(2,4)){warning("Ordering of keyVarF is not chronological. Is this correct? May cause undesireable results if not deliberately set!");n<-n+1}
	if(tmptest$error==3){warning("Largest value of index is larger than the dimension of keyVarF. Will result in crash!");n<-n+1}
	start2<-tmptest$start
	tmptest<-CheckOrdering(t(settings$keyLogqpar))
	if(tmptest$error==1){warning("Start value of keyLogqpar should be 0 or 1. Will result in crash!");n<-n+1}
	if(tmptest$error%in%c(2,4)){warning("Ordering of keyLogqpar is not chronological. Will result in crash!");n<-n+1}
	if(tmptest$error==3){warning("Largest value of index is larger than the dimension of keyLogqpar. Will result in crash!");n<-n+1}
	start3<-tmptest$start
	tmptest<-CheckOrdering(t(settings$keyVarObs))
	if(tmptest$error==1){warning("Start value of keyVarObs should be 0 or 1. Will result in crash!");n<-n+1}
	if(tmptest$error%in%c(2,4)){warning("Ordering of keyVarObs is not chronological. Will result in crash!");n<-n+1}
	if(tmptest$error==3){warning("Largest value of index is larger than the dimension of keyVarObs. Will result in crash!");n<-n+1}
	start4<-tmptest$start
	if(!sum(c(start1,start2,start3,start4))%in%c(0,4)){warning("Start value of indexing of keyLogFconst, keyVarF, keyLogqpar and keyVarObs must be consistent and either 0 or 1. Will result in crash!");n<-n+1}
	#Check correspondance between agerangeI and keyLogqpar
	if(settings$Modelq==0){
		agerangeI<-settings$agerangeI
		#print(agerangeI)
		keyLogqpar<-settings$keyLogqpar
		for(i in 1:nrow(keyLogqpar)){
			idx<-1:length(minA:maxA)
			nonapos<-match(agerangeI[i,1]:agerangeI[i,2],minA:maxA)
			napos<-(1:length(idx))[is.na(match(minA:maxA,agerangeI[i,1]:agerangeI[i,2]))]
			#Check correspondance with keyLogqpar
			keyidx<-keyLogqpar[i,]
			#First test: If keyLogqpar states NA values for q for ages that is included by agerange
			test1<-keyidx[nonapos]
			test1[test1<0]<-NA
			if(sum(!is.na(test1))!=length(nonapos)){warning(paste("keyLogqpar dictates NA for ages contained in agerange for index:",i,"Include the correct index or correct age range!",sep=" "));n<-n+1}
			#Second test: If keyLogqpar  values for q outside ages defined by agerange
			test2<-keyidx[napos]
			test2[test2<0]<-NA
			if(sum(is.na(test2))!=length(napos)){warning(paste("keyLogqpar is specified for ages outside agerange for index:",i,"Include the correct indices or correct age range!",sep=" "));n<-n+1}
		}
	}
	if(n>0){
		warning(paste("'CheckXSAMsettings' found",n,"potential problems with the setting file that may need attention!",collapse=" "))
	}
	else{
		cat("Found the setting file OK!\n")
	}
}

ReadTASACSFleetFile<-function(file="C:/Users/aanes/Documents/SILD/Benchmark/WD/fleet.txt",na.rm=T,sep=" "){
	#Reads and interprets the TASACS format of a fleet file, i.e. a file containint multiple indices on a specific format.
	tt<-scan(file,what="",sep="\n",comment.char = "#",flush=T,quiet=T)
	tmp<-lapply(tt,FUN=function(x)strsplit(x,split=sep))
	#Remove empty places
	tt1<-lapply(tmp,FUN=function(x)unlist(x)[unlist(x)!=""])
	#Find identyfier
	id<-lapply(tt1,FUN=function(x)all.equal(c("-1","-1","-1","-1"),x))
	#Extract first element of each vector in list (since all.equal returns a vector of text if mismatch)
	id<-lapply(id,FUN=function(x)x[1])
	#Extract line positions of identyfier
	id<-(1:length(id))[unlist(id)==T]
	tt1<-lapply(tt1,FUN=function(x)unlist(x))
	print(id)
	out<-list()
	for(i in 1:length(id)){
		out[[i]]<-NULL
		name<-paste(unlist(tt1[id[i]-2]),collapse=" ")
		print(name)
		yrrange<-as.numeric(unlist(tt1[id[i]-1]))
		agerange<-as.numeric(unlist(tt1[id[i]+1]))
		nlines<-diff(yrrange)+1
		print(yrrange)
		ncols<-c()
		for(j in 1:nlines){
			tmp<-unlist(tt1[id[i]+1+j])
			tmp[tmp=="NA"]<-"-9999"
			tmp<-as.numeric(tmp)
			tmp<-tmp[!is.na(tmp)]
			ncols[j]<-length(tmp)
			if(j>1){
				if(ncols[j]!=ncols[j-1])stop(paste("Different number of columns for",name,sep=""))
				out[[i]]<-rbind(out[[i]],tmp)
			}
			else out[[i]]<-tmp
		}
		if(j>1)dimnames(out[[i]])[[1]]<-NULL
		out[[i]][out[[i]]==-9999]<-NA
		names(out)[i]<-name
		out[[i]]<-out[[i]][,2:ncol(out[[i]]),drop=F]
		dimnames(out[[i]])[[2]]<-agerange[1]:agerange[2]
		dimnames(out[[i]])[[1]]<-yrrange[1]:yrrange[2]
		if(na.rm)out[[i]][out[[i]]<0]<-NA
		attributes(out[[i]])$ages<-agerange[1]:agerange[2]
		attributes(out[[i]])$years<-yrrange[1]:yrrange[2]
	}
	out
}

ReadDataFiles0<-function(file,na.rm=T,sep=" "){
	#Reads: canum, weca, west, natmor, matprop, etc, i.e. a file which contains one data table with equal numbers of columns
	#where a descriptive name is found in the first line, the years on line 3 and ages on line 4, and the table itself starts at line 6.
	tt<-scan(file,what="",sep="\n",comment.char = "#",flush=T,quiet=T)
	#tmp<-lapply(tt,FUN=function(x)strsplit(x,split=" "))
	tmp<-lapply(tt,FUN=function(x)strsplit(x,split=sep))
	#Remove empty places
	tt1<-lapply(tmp,FUN=function(x)unlist(x)[unlist(x)!=""])
	tt1<-lapply(tt1,FUN=function(x)unlist(x))
	name<-paste(tt1[[1]],collapse=" ")
	years<-as.numeric(tt1[[3]])
	years<-years[!is.na(years)]
	years<-years[1]:years[2]
	ages<-as.numeric(tt1[[4]])
	ages<-ages[!is.na(ages)]
	ages<-ages[1]:ages[2]
	#print(ages);print(years)
	nlines<-length(years)
	out<-NULL
	ncols<-c()
	for(j in 1:nlines){
		tmp<-unlist(tt1[5+j])
		tmp[tmp=="NA"]<-"-9999"
		tmp<-as.numeric(tmp)
		tmp<-tmp[!is.na(tmp)]
		ncols[j]<-length(tmp)
		if(j>1){
			if(ncols[j]!=ncols[j-1])stop(paste("Different number of columns for",name,sep=""))
			out<-rbind(out,tmp)
			}
		else out<-tmp
	}
	if(j>1)dimnames(out)[[1]]<-NULL
	out[out==-9999]<-NA
	dimnames(out)[[2]]<-ages
	dimnames(out)[[1]]<-years
	if(na.rm)out[out<0]<-NA
	attributes(out)$ages<-ages
	attributes(out)$years<-years
	out
}


CheckDataList<-function(datalist){
	reqobj<-c(
	"caa",
	"SurveyIndices",
	"matprop",
	"west",
	"weca",
	"natMor",
	"propF",
	"propM",
	"caton",
	"CatchPrediction"
	)
	listobj<-names(datalist)
	missingobj<-reqobj[!reqobj%in%listobj]
	if(length(missingobj)>0)warning(paste("Required objects:<",paste(missingobj,collapse=", "),"> are missing from data list and will cause the function to crash or produce erroneous results.",collapse=" "))
	redundantobj<-listobj[!listobj%in%reqobj]
	if(length(redundantobj)>0)warning(paste("Objects:<",paste(redundantobj,collapse=", "),"> found in data list./nThese will not be used by XSAM!",collapse=" "))	
	if(length(missingobj)==0 & length(redundantobj)==0)errorcode<-0
	if(length(missingobj)>0 & length(redundantobj)==0)errorcode<-1
	if(length(missingobj)>0 & length(redundantobj)>0)errorcode<-2
	if(length(missingobj)==0 & length(redundantobj)>0)errorcode<-3
	list(errorcode=errorcode,missingobj=missingobj,redundantobj=redundantobj)
}

CheckDataandSettings<-function(settings,datalist){
	errorcode<-0
	if(settings$UseCatchPred==1){
		stopyear<-max(settings$yearrange)
		useprediction<-datalist$CatchPrediction[paste(stopyear),c("Prediction","se")]
		if(is.na(useprediction[1]) | useprediction[1]<0 | is.na(useprediction[2]) | useprediction[2]<0){
			errorcode<-1
		}
		
	}
	#Insert test for data age range vs setup...
	list(errorcode=errorcode)
}
#####################################
#CreateDataObject<-function(ages,years,caa,deltaC,SurveyIndices,indices,deltaI,kcaa=1000,ki=1/1000,
#	IsampleTimes,matprop,west,weca,natMor,propF,caton,kC=1/1000,NACval=NULL,NAIval=NULL,settings=NULL){
CreateDataObject<-function(settings=XSAMsettings,datalist,deltaC=0,deltaI=0,NAval=NULL,actualAge=T){
	#----------------------------------------------------------------------------------------------------------------------#
	#Combines XSAM settings (model configuration) with data into an object that can be passed onto the XSAM template

	#TODO: Define settings here
	#TODO: Define datalist here
	#----------------------------------------------------------------------------------------------------------------------#
	tmp<-CheckDataList(datalist)
	if(tmp$errorcode%in%c(1,2))stop("Datalist misses necessary objects!")
	tmp<-CheckDataandSettings(settings,datalist)
	if(tmp$errorcode%in%c(1))stop("Catchprediction and corresponding SE must be available and positive for assessmentyear if to be used when fitting the model!")

	#unlist the objects in datalist
	caa<-datalist$caa
	SurveyIndices<-datalist$SurveyIndices
	matprop<-datalist$matprop
	west<-datalist$west
	weca<-datalist$weca
	natMor<-datalist$natMor
	propF<-datalist$propF
	propM<-datalist$propM
	caton<-datalist$caton
	CatchPrediction<-datalist$CatchPrediction
	#Takes input data plus some additional information and creates a data object that can be fed into the SAM model or the GENERALIZED GUDMUNDSSON
	#if !is.null(settings); it apply settings provided by the list
	#-----------------------------------------------------------------#
	#ages: vector of ages to keep
	#years: vector years to keep
	#caa: matrix of catch at age in numbers rows-years, columns-ages
	#deltaC: constant to add to 0 observations of caa
	#SurveyIndices: list of sruveyindices, each a matrix with numerbs at age in numbers rows-years, columns-ages
	#indices: vector of indices with length equal to length of SurveyIndices with value corresponding to SurveyIndices.
	#deltaI: vector of constants to add to 0 observations of SurveyIndices
	#kcaa: relative unit conversion of caa
	#ki: vector of relative unit conversion of SurveyIndices
	#IsampleTimes: vector of length equal to SurveyIndices with elements in corresponding position providing timing of the survey in fraction of the year
	#matprop: matrix of proportion mature at age (columns) and year (rows)
	#west: matrix of mean weight in stock at age (columns) and year (rows)
	#weca: matrix of mean weight in catch at age (columns) and year (rows)
	#natMor: matrix of natural mortality at age (columns) and year (rows)
	#propF: matrix of proportion of F before spawning at age (columns) and year (rows)
	#caton: 1 column matrix of total landings in year (rows)
	#kC: relative unit conversion of caton

	#All matrices must be data.frame with ages as column names and years as row names
	#-----------------------------------------------------------------#
	if(is.null(NAval))NAval<--999
	#PREPARE obs in data object

	#XSAMsettings
	minage<-min(settings$agerange)
	maxage<-max(settings$agerange)
	ages<-minage:maxage
	nages<-length(ages)
	startyear<-min(settings$yearrange)
	stopyear<-max(settings$yearrange)
	years<-startyear:stopyear
	nyears<-length(years)
	#---------------------------------#
	#Check if datarange is consistent with settings range
	#---------------------------------#
	#---------------------------------#
	#Reduce dataset to ages and years provided by settings
	#---------------------------------#
	agesC<-settings$agerangeC[1]:settings$agerangeC[2]
	kcaa<-attributes(caa)$k
	caa<-caa[as.numeric(dimnames(caa)[[1]])%in%years,as.numeric(dimnames(caa)[[2]])%in%agesC]
	#print(caa*kcaa)
	#Create dataobject
	agemat<-matrix(as.numeric(rep(colnames(caa),nrow(caa))),nrow=nrow(caa),byrow=T)
	agerangeC<-range(agemat)
	yearmat<-matrix(as.numeric(rep(rownames(caa),ncol(caa))),nrow=nrow(caa),byrow=F)
	#dataobj<-data.frame(year=c(yearmat),fleet=1,age=c(agemat),obs=unlist(c(caa))*kcaa)
	dataobj<-data.frame(year=c(yearmat),fleet=0,age=c(agemat),obs=unlist(c(caa))*kcaa)
	caa<-caa*kcaa
	#dataobj[dataobj$obs==0,"obs"]<-deltaC#sort(unique(dataobj$obs))[2]
	#dataobj[dataobj$obs==0,"obs"]<-sort(unique(dataobj$obs))[2]
	if(deltaC>0){
		dataobj[,"obs"]<-dataobj[,"obs"]+deltaC
		caa<-caa+deltaC
	}
	else{
		dataobj[dataobj$obs==0,"obs"]<-NAval
		caa[caa==0]<-NAval
	}
	#Must establish age range of survey indices since they may not include all ages requested by ages
	indices<-settings$fleetIndex
	agerangeI<-matrix(NA,nrow=length(indices)-1,ncol=2)
	sampleTimes<-0
	SurveyList<-list()
	#indices<-settings$fleetTypes
	indices<-indices[indices>0]#0 is reserved for catch@age
	for(i in 1:length(indices)){
		#print(i)
		#tmp<-SurveyIndices[[i]]
		#Only the fleetIndices that are specyfied in settings
		tmpSI<-SurveyIndices[indices[i]][[1]]
		indeks<-indices[i]
		#Reduce dataset to ages and years
		agesI<-settings$agerangeI[i,1]:settings$agerangeI[i,2]
		agesI<-agesI[agesI%in%ages]
		tmp<-tmpSI[as.numeric(dimnames(tmpSI)[[1]])%in%years,as.numeric(dimnames(tmpSI)[[2]])%in%agesI,drop=F]
		####
		colnms<-dimnames(tmp)[[2]]
		h0<<-tmp
		print("0")
		print(tmp)
		tmp<-tmp[!is.na(rowSums(tmp)),]#,drop=F]
		print("0.5")
		if(is.null(nrow(tmp))){
			rownms<-names(tmp)
			tmp<-as.data.frame(matrix(tmp),ncol=1)
			dimnames(tmp)[[2]]<-colnms
			dimnames(tmp)[[1]]<-rownms
			rownms<-dimnames(tmp)[[1]]
			rownms<-rownms[rowSums(tmp)>0]
			tmp<-tmp[rowSums(tmp)>0,]
			tmp<-as.data.frame(matrix(tmp),ncol=1)
			dimnames(tmp)[[2]]<-colnms
			dimnames(tmp)[[1]]<-rownms
		}
		else{
			tmp<-tmp[rowSums(tmp)>0,]
		}
		h1<<-tmp
		###
		print("1")
		#print(tmp)
		#Create dataobject
		agemat<-matrix(as.numeric(rep(colnames(tmp),nrow(tmp))),nrow=nrow(tmp),byrow=T)
		yearmat<-matrix(as.numeric(rep(rownames(tmp),ncol(tmp))),nrow=nrow(tmp),byrow=F)
		#h1<<-agemat
		#h2<<-yearmat
		#k<-ifelse(length(ki)>1,ki[i],ki)
		k<-attributes(tmpSI)$k
		#print(tmp*k)
		SurveyList[[i]]<-tmp*k
		tmp<-data.frame(year=c(yearmat),fleet=indices[i],age=c(agemat),obs=unlist(c(tmp))*k)
		delta<-ifelse(length(deltaI)>1,deltaI[i],deltaI)
		#tmp[tmp$obs==0,"obs"]<-delta
		if(delta>0){
			tmp[,"obs"]<-tmp[,"obs"]+delta
			SurveyList[[i]]<-SurveyList[[i]]+delta
		}
		else{
			tmp[tmp$obs==0,"obs"]<-NAval
			SurveyList[[i]][SurveyList[[i]]==0]<-NAval
		}
		dataobj<-rbind(dataobj,tmp)
		agerangeI[i,]<-range(agemat)
		sampleTimes[i+1]<-attributes(tmpSI)$sampleTime
	}
	tmp<-dataobj
	dataobj<-list()
	#---------------------------------------------------------#
	# INSERT DATA TO BE USED FOR DIAGNOSTICS
	#---------------------------------------------------------#
	dataobj$caa<-caa
	dataobj$SurveyList<-SurveyList
	#---------------------------------------------------------#
	# INSERT DATA TO BE USED IN THE C TEMPLATE
	#---------------------------------------------------------#
	#dataobj$fleetTypes<-settings$fleetTypes#c(0,indices)
	dataobj$fleetIndex<-settings$fleetIndex#c(0,indices)
	dataobj$sampleTimes<-sampleTimes#c(0,IsampleTimes)
	dataobj$noYears<-length(unique(tmp$year))#length(unique(data$obs[,1]))#data$noYears
	dataobj$years<-years
	dataobj$obs<-tmp
	dataobj$nobs<-nrow(dataobj$obs)
	dataobj$propMat<-matprop[as.numeric(rownames(matprop))>=startyear,as.numeric(colnames(matprop))>=minage]
	dataobj$stockMeanWeight<-west[as.numeric(rownames(west))>=startyear,as.numeric(colnames(west))>=minage]
	dataobj$catchMeanWeight<-weca[as.numeric(rownames(weca))>=startyear,as.numeric(colnames(weca))>=minage]
	dataobj$natMor<-natMor[as.numeric(rownames(natMor))>=startyear,as.numeric(colnames(natMor))>=minage]
	dataobj$propF<-propF[as.numeric(rownames(propF))>=startyear,as.numeric(colnames(propF))>=minage]
	dataobj$propM<-propM[as.numeric(rownames(propM))>=startyear,as.numeric(colnames(propM))>=minage]
	#print(caton);print(startyear);print(stopyear)
	dataobj$caton<-caton[paste(startyear:(stopyear-1)),]*attributes(caton)$k
	dataobj$minAge<-minage#data$minAge
	dataobj$maxAge<-maxage#data$maxAge
	dataobj$agerangeI<-agerangeI
	dataobj$agerangeC<-agerangeC
	dataobj$UseCatchPred<-settings$UseCatchPred
	CatchPrediction[,"Prediction"]<-CatchPrediction[,"Prediction"]*attributes(CatchPrediction)$k
	dataobj$CatchPrediction<-c(unlist(CatchPrediction[paste(stopyear),c("Prediction","se")]))
	#---------------------------------------------------------#
	# CREATE DEFAULT MODEL SETTINGS
	#---------------------------------------------------------#
	if(is.null(settings$maxAgePlusGroup))dataobj$maxAgePlusGroup<-1
	else dataobj$maxAgePlusGroup<-settings$maxAgePlusGroup
	
	if(is.null(settings$keyLogFconst)){
		#keyLogFconst<-matrix(-1,nrow=length(SurveyIndices),ncol=nages)
		keyLogFconst[1,]<-minage:maxage-minage
		keyLogFconst<-minage:maxage-minage
		#keyLogFconst[1,nages]<-nages-2
		keyLogFconst[nages]<-nages-2
	}
	else{
		keyLogFconst<-settings$keyLogFconst
	} 
	testindeks<-min(keyLogFconst,na.rm=T)
	if(testindeks==1){#Adjust to C indexing
		keyLogFconst[is.na(keyLogFconst)]<-0
		keyLogFconst<-keyLogFconst-1
	}
	else keyLogFconst[is.na(keyLogFconst)]<--1

	dataobj$keyLogFconst<-keyLogFconst#data$keyLogFsta
	#h3<<-agerangeI
	if(is.null(settings$keyLogqpar)){
		#keyLogqpar<-matrix(0,nrow=length(SurveyIndices)+1,ncol=nages)
		keyLogqpar<-matrix(0,nrow=length(SurveyIndices),ncol=nages)
		startid<-1
		#for(i in 2:nrow(keyLogqpar)){
		for(i in 1:nrow(keyLogqpar)){
			#sdata$keyLogqpar
			#tmpagerange<-agerangeI[i-1,1]:agerangeI[i-1,2]
			tmpagerange<-agerangeI[i,1]:agerangeI[i,2]
			keyLogqpar[i,ages%in%tmpagerange]<-startid:(length(tmpagerange)+startid-1)
			startid<-length(tmpagerange)+startid
		}
		#Adjust to C indexing
		keyLogqpar<-keyLogqpar-1
	}
	else keyLogqpar<-settings$keyLogqpar

	testindeks<-min(keyLogqpar,na.rm=T)
	if(testindeks==1){#Adjust to C indexing
		keyLogqpar[is.na(keyLogqpar)]<-0
		keyLogqpar<-keyLogqpar-1
	}
	else keyLogqpar[is.na(keyLogqpar)]<--1
	dataobj$keyLogqpar<-keyLogqpar

	if(is.null(settings$Modelq)){
		dataobj$Modelq<-0
	}
	else{
		dataobj$Modelq<-settings$Modelq
	}
	if(dataobj$Modelq==0)dataobj$nqs<-max(dataobj$keyLogqpar)+1
	else dataobj$nqs<-3*(nrow(keyLogqpar)-1)

	
	if(is.null(settings$keyVarObs)){
		keyVarObs<-matrix(-1,nrow=(length(SurveyIndices)+1),ncol=nages)
		for(i in 1:nrow(keyVarObs)){
			keyVarObs[i,]<-i-1
		}
	}
	else{
		keyVarObs<-settings$keyVarObs
	}
	testindeks<-min(keyVarObs,na.rm=T)
	if(testindeks==1){#Adjust to C indexing
		keyVarObs[is.na(keyVarObs)]<-0
		keyVarObs<-keyVarObs-1
	}
	else keyVarObs[is.na(keyVarObs)]<--1


	dataobj$keyVarObs<-keyVarObs#data$keyVarObs
	dataobj$nObsVar<-max(keyVarObs)+1

	
	if(is.null(settings$stockRecruitmentModelCode))dataobj$stockRecruitmentModelCode<-matrix(5)#data$stockRecruitmentModelCode
	else stockRecruitmentModelCode<-settings$stockRecruitmentModelCode
	dataobj$stockRecruitmentModelCode<-stockRecruitmentModelCode

	dataobj$nIndices<-length(indices)
	dataobj$SurveyIndex<-indices

	dataobj$dimC<-length(c(dataobj$minAge):c(dataobj$maxAge))
	
	dimI<-c()
	for(i in 1:nrow(agerangeI)){
		dimI[i]<-length(agerangeI[i,1]:agerangeI[i,2])
	}
	dataobj$dimI<-dimI

	#Create the following inputs:
	#DATA_ARRAY(sd_C);//matrix holding sd of catch by year (n,y)
	#DATA_ARRAY(R_C);//3 dim array holding correlation matrices of catch by year (n,n,y)
	#DATA_ARRAY(sd_I);//3 dim array holding sd of index by year and fleet (n,y,f)
	#DATA_ARRAY(R_I);//4 dim array holding correlations of index by year and fleet (n,n,y,f)
	dataobj$sd_C<-matrix(1,nrow=dataobj$dimC,ncol=dataobj$noYears)
	
	R_C<-array(NA,dim=c(dataobj$dimC,dataobj$dimC,dataobj$noYears))
	for(i in 1:dataobj$noYears)R_C[,,i]<-diag(dataobj$dimC)
	dataobj$R_C<-R_C

	#sd_I<-array(100,dim=c(max(dimI),dataobj$noYears,sum(dataobj$fleetTypes>0)))
	sd_I<-array(100,dim=c(max(dimI),dataobj$noYears,sum(dataobj$fleetIndex>0)))
	#for(i in 1:sum(dataobj$fleetTypes>0)){
	for(i in 1:sum(dataobj$fleetIndex>0)){
		for(j in 1:dataobj$noYears){
			sd_I[,j,i]<-rep(1,max(dimI))
		}
	}
	dataobj$sd_I<-sd_I

	#R_I<-array(0,dim=c(max(dimI),max(dimI),dataobj$noYears,sum(dataobj$fleetTypes>0)))
	R_I<-array(0,dim=c(max(dimI),max(dimI),dataobj$noYears,sum(dataobj$fleetIndex>0)))
	#for(i in 1:sum(dataobj$fleetTypes>0)){
	for(i in 1:sum(dataobj$fleetIndex>0)){
		for(j in 1:dataobj$noYears){
			R_I[1:max(dimI),1:max(dimI),j,i]<-diag(max(dimI))
		}
	}	
	dataobj$R_I<-R_I
	
	dataobj$NAval<-NAval

	#if(is.null(settings$Fbarminage))dataobj$Fbarminage<-5
	#else dataobj$Fbarminage<-settings$Fbarminage

	#if(is.null(settings$Fbarmaxage))dataobj$Fbarmaxage<-14
	#else dataobj$Fbarmaxage<-settings$Fbarmaxage

	if(is.null(settings$Fbarrange))dataobj$Fbarminage<-5
	else dataobj$Fbarminage<-settings$Fbarrange[1]

	if(is.null(settings$Fbarrange))dataobj$Fbarmaxage<-14
	else dataobj$Fbarmaxage<-settings$Fbarrange[2]

	#if(is.null(settings$biascorr))dataobj$biascorr<-0
	#else dataobj$biascorr<-settings$biascorr

	if(is.null(settings$CatchConstraint))dataobj$CatchConstraint<-0
	else dataobj$CatchConstraint<-settings$CatchConstraint


	if(is.null(settings$keyVarF)){
		keyVarF<-matrix(-1,nrow=2,ncol=nages)
		keyVarF[1,]<-0
	}
	else{
		keyVarF<-settings$keyVarF
	}
	testindeks<-min(keyVarF,na.rm=T)
	if(testindeks==1){#Adjust to C indexing
		keyVarF[is.na(keyVarF)]<-0
		keyVarF<-keyVarF-1
	}
	else keyVarF[is.na(keyVarF)]<--1


	dataobj$keyVarF<-keyVarF
	dataobj$nvf<-max(dataobj$keyVarF)+1

	#---------------------------------------------------------#
	# SETTINGS ONLY RELEVANT FOR SAM: OBS OBS NOT PROPERLY IMPLEMENTED AND TURNED OFF
	#---------------------------------------------------------#
	if(0){
	if(is.null(settings$keyVarLogN)){
		keyVarLogN<-matrix(c(0,rep(1,length((minage+1):maxage))),nrow=1)#data$keyVarLogN
	}
	else{
		keyVarLogN<-settings$keyVarLogN
	}
	dataobj$keyVarLogN<-keyVarLogN
	
	#dataobj$nlogF<-max(dataobj$keyLogFsta)+1
	dataobj$nlogF<-max(dataobj$keyLogFconst)+1

	
	dataobj$nlogN<-dataobj$maxAge-dataobj$minAge+1

	if(is.null(settings$corFlag))dataobj$corFlag<-matrix(0)#matrix(1)
	else corFlag<-settings$corFlag
	}
	#---------------------------------------------------------#
	# SETTINGS ONLY RELEVANT FOR XSAM
	#---------------------------------------------------------#
	if(is.null(settings$corFlaglogF))dataobj$corFlaglogF<-0#if(sdata$corFlaglogF>0): Estimate rhologF
	else dataobj$corFlaglogF<-settings$corFlaglogF

	if(is.null(settings$LatentEffort))dataobj$LatentEffort<-0#if(sdata$LatentEffort==1): Add latent V and estimate s2_3
	else dataobj$LatentEffort<-settings$LatentEffort

	if(is.null(settings$TSsel))dataobj$TSsel<-0#if(sdata$TSsel==1): Add latent for U (am x Y), estimate s2_2
	else dataobj$TSsel<-settings$TSsel

	if(is.null(settings$corFlagU))dataobj$corFlagU<-0#if(sdata$TSsel==1 & corFlagU>0): estimate rhoU
	else dataobj$corFlagU<-settings$corFlagU

	if(is.null(settings$EstimateARInterceptEffort))dataobj$EstimateARInterceptEffort<-1#if(sdata$EstimateARInterceptEffort==1): Estimate aY
	else dataobj$EstimateARInterceptEffort<-settings$EstimateARInterceptEffort

	if(is.null(settings$EstimateUInterceptEffort))dataobj$EstimateUInterceptEffort<-0#if(sdata$TSsel==1 & sdata$EstimateUInterceptEffort==1): Estimate aUpar(am) 
	else dataobj$EstimateUInterceptEffort<-settings$EstimateUInterceptEffort

	if(is.null(settings$ymeantype))dataobj$ymeantype<-0
	else dataobj$ymeantype<-settings$ymeantype

	if(is.null(settings$am))dataobj$am<-dataobj$maxAge-dataobj$minAge#sdata$nlogN-1
	else dataobj$am<-settings$am

	#if(is.null(settings$Am))dataobj$Am<-max(dataobj$keyLogFsta)+1
	if(is.null(settings$Am))dataobj$Am<-max(dataobj$keyLogFconst)+1
	else dataobj$Am<-settings$Am

	#if(is.null(settings$A))dataobj$A<-dataobj$nlogN
	#if(is.null(settings$A))dataobj$A<-dataobj$nlogN
	#else dataobj$A<-settings$A
	dataobj$A<-maxage

	dataobj$Fbarmaxage<-min(dataobj$Fbarmaxage,dataobj$A)
	#print(dataobj$Fbarmaxage)

	if(is.null(settings$AT))dataobj$AT<-dataobj$A
	else dataobj$AT<-settings$AT

	if(actualAge){
		dataobj$A<-dataobj$A-minage+1
		dataobj$Am<-dataobj$Am-minage+1
		dataobj$am<-dataobj$am-minage+1
		dataobj$AT<-dataobj$AT-minage+1
	}
	if(is.null(settings$RecruitmentProcess))dataobj$RecruitmentProcess<-1#if(sdata$RecruitmentProcess==0) Estimate initial values
	else dataobj$RecruitmentProcess<-settings$RecruitmentProcess

	#SKAL DET VÆRE -1 el ikke?
	if(is.null(settings$uam))dataobj$uam<-dataobj$am-1
	else dataobj$uam<-settings$uam

	dataobj$ReturnLL<-0
	
	#Type conversion to make it readable by C
	dataobj$obs<-as.matrix(dataobj$obs)
	dataobj$propMat<-as.matrix(dataobj$propMat)
	dataobj$stockMeanWeight<-as.matrix(dataobj$stockMeanWeight)
	dataobj$catchMeanWeight<-as.matrix(dataobj$catchMeanWeight)
	dataobj$natMor<-as.matrix(dataobj$natMor)
	dataobj$propF<-as.matrix(dataobj$propF)
	dataobj$propM<-as.matrix(dataobj$propM)


	dataobj
}


CreateParameterList<-function(dataobj,mapnames=NULL,sv=0){
	#-----------------------------------------------#
	#First determine dimension of latent variables
	#-----------------------------------------------#
	dimL<-dataobj$Am
	if(dataobj$TSsel==1)dimL<-dimL+(dataobj$am-1)
	if(dataobj$LatentEffort==1)dimL<-dimL+1
	dimL<-dimL+1#For Y
	if(dataobj$RecruitmentProcess==1)dimL<-dimL+1
	#-----------------------------------------------#
	#Then set up parameterlist
	#-----------------------------------------------#
	parameterlist<-list(
	initlogN=rep(0,dataobj$A-1),
	R=if(dataobj$RecruitmentProcess==0){rep(numeric(1),dataobj$noYears)}else{numeric(0)},
	logqpar=if(dataobj$Modelq==1){rep(-10,3*dataobj$nIndices)}else{rep(0,max(dataobj$keyLogqpar)+1)},
	#logs2_1=0,
	logs2_1=rep(0,max(dataobj$keyVarF)+1),
	logs2_2=if(dataobj$TSsel==1){numeric(1)}else{numeric(0)},
	logs2_3=if(dataobj$LatentEffort==1){numeric(1)}else{numeric(0)},
	logs2_4=0,
	logsR2=if(dataobj$RecruitmentProcess==1){numeric(1)}else{numeric(0)},
	loghvec=rep(0,max(dataobj$keyVarObs)+1),
	rec_loga=0,
 	rec_logb=0,
 	rhoU=if(dataobj$TSsel==1 & dataobj$corFlagU>0){numeric(1)}else{numeric(0)}, 
 	rhologF=if(dataobj$corFlaglogF>0){numeric(1)}else{numeric(0)}, 
 	aYpar=if(dataobj$EstimateARInterceptEffort==1){numeric(1)}else{numeric(0)},
 	bY=1,
	Upar=if(dataobj$TSsel==0){rep(numeric(1),dataobj$am-1)}else{numeric(0)},
	aUpar=if(dataobj$TSsel==1 & dataobj$EstimateUInterceptEffort==1){rep(numeric(1),dataobj$uam)}else{numeric(0)},
	bU=if(dataobj$TSsel==1){1}else{numeric(0)}#,
	#PredictCatchPars=dataobj$PredictCatchPars
	)
	#Convert to vector
	ParameterVector<-unlist(parameterlist)
	ParameterNames<-names(ParameterVector)
	#Identify parameters to map
	if(dataobj$RecruitmentProcess){
		if(dataobj$stockRecruitmentModelCode%in%0:1)mapnames<-c(mapnames,"rec_loga","rec_logb")
		if(dataobj$stockRecruitmentModelCode==5)mapnames<-c(mapnames,"rec_logb")
	}
	else mapnames<-c(mapnames,"rec_loga","rec_logb")
	#if(dataobj$PredictCatchPars[1]>0)mapnames<-c(mapnames,c("PredictCatchPars1","PredictCatchPars2"))
	#mapnames<-c(mapnames,c("PredictCatchPars1","PredictCatchPars2"))
	cat("Parameters mapped (also due to specified stock-recruitment model):\n")
	print(mapnames)
	#vector of names to map
	MapVector<-mapvec(ParameterVector,mapnames)
	#Parameterlist
	ParameterList <- list(
	FixedPars=ParameterVector,
	LP=matrix(sv, nrow=dimL ,ncol=dataobj$noYears)
	)
	list(ParameterList=ParameterList,MapVector=MapVector,ParameterNames=ParameterNames)
}

mapvec<-function(vec,fixname=NULL){
	map<-1:length(vec)
	if(!is.null(fixname)){
		map[names(vec)%in%fixname]<-NA
		#Reorder
		map[!is.na(map)]<-1:length(map[!is.na(map)])
	}	
	map	
}

CreateReportTable<-function(repobj,parobj){
	stab<-summary(repobj,"fixed")
	dimnames(stab)[[1]]<-names(parobj$ParameterList$FixedPars)[!is.na(parobj$MapVector)]
	stab
}

FitXSAMFunction<-function(dataobj,pars,control=list(eval.max=200,iter.max=150),lower=-Inf,upper=Inf,constraints=NULL,includeDiagnostics=T,returnObj=F){
	obj <- MakeADFun(dataobj,pars$ParameterList,random=c("LP"),DLL="XSAM",map=list(FixedPars=factor(pars$MapVector)))
	lower <- obj$par*0+lower
	upper <- obj$par*0+upper
	parnames<-pars$ParameterNames[!is.na(pars$MapVector)]
	print(parnames)
	print(lower)
	if(!is.null(constraints)){
		for(i in 1:length(constraints)){
			lower[parnames==constraints[[i]]$par]<-constraints[[i]]$lower
			upper[parnames==constraints[[i]]$par]<-constraints[[i]]$upper
		}
	}
	system.time(optobj<-nlminb(obj$par,obj$fn,obj$gr,control=control,lower=lower,upper=upper))
	#cbind(optXSAMobj$par,XSAMpars$ParameterNames[!is.na(XSAMpars$MapVector)])
	if(optobj$convergence==0){
		if(returnObj)rep<-sdreport(obj,getJointPrecision=T)
		else rep<-sdreport(obj)
		tab<-CreateReportTable(rep,pars)
		stats<-extractStatsG(report=summary(rep,"report"),dataobj=dataobj)
		if(includeDiagnostics)res<-ExtractResiduals(stats,tab,dataobj)
		else res<-NULL
	}
	else{
		rep<-tab<-stats<-res<-NULL
	}
	if(returnObj)obj<-obj
	else obj<-NULL
	convergence<-optobj$convergence
	list(convergence=convergence,optobj=optobj,rep=rep,tab=tab,stats=stats,res=res,obj=obj)
}

FitXSAMFunctionr<-function(dataobj,pars,control=list(eval.max=200,iter.max=150),lower=-Inf,upper=Inf,constraints=NULL,includeDiagnostics=T,returnObj=F){
	obj <- MakeADFun(dataobj,pars$ParameterList,random=c("LP"),DLL="XSAMr",map=list(FixedPars=factor(pars$MapVector)))
	lower <- obj$par*0+lower
	upper <- obj$par*0+upper
	parnames<-pars$ParameterNames[!is.na(pars$MapVector)]
	print(parnames)
	print(lower)
	if(!is.null(constraints)){
		for(i in 1:length(constraints)){
			lower[parnames==constraints[[i]]$par]<-constraints[[i]]$lower
			upper[parnames==constraints[[i]]$par]<-constraints[[i]]$upper
		}
	}
	system.time(optobj<-nlminb(obj$par,obj$fn,obj$gr,control=control,lower=lower,upper=upper))
	#cbind(optXSAMobj$par,XSAMpars$ParameterNames[!is.na(XSAMpars$MapVector)])
	if(optobj$convergence==0){
		if(returnObj)rep<-sdreport(obj,getJointPrecision=T)
		else rep<-sdreport(obj)
		tab<-CreateReportTable(rep,pars)
		stats<-extractStatsG(report=summary(rep,"report"),dataobj=dataobj)
		if(includeDiagnostics)res<-ExtractResiduals(stats,tab,dataobj)
		else res<-NULL
	}
	else{
		rep<-tab<-stats<-res<-NULL
	}
	if(returnObj)obj<-obj
	else obj<-NULL
	convergence<-optobj$convergence
	list(convergence=convergence,optobj=optobj,rep=rep,tab=tab,stats=stats,res=res,obj=obj)
}


extractStatsG<-function(report,dataobj){

	tmp<-report[dimnames(report)[[1]]%in%"N",1]
	N<-matrix(tmp,nrow=dataobj$AT,ncol=dataobj$noYears)
	tmp<-report[dimnames(report)[[1]]%in%"N",2]
	Nse<-matrix(tmp,nrow=dataobj$AT,ncol=dataobj$noYears)

	tmp<-report[dimnames(report)[[1]]%in%"logN",1]
	logN<-matrix(tmp,nrow=dataobj$AT,ncol=dataobj$noYears)
	tmp<-report[dimnames(report)[[1]]%in%"logN",2]
	logNse<-matrix(tmp,nrow=dataobj$AT,ncol=dataobj$noYears)

	tmp<-report[dimnames(report)[[1]]%in%"F",1]
	FF<-matrix(tmp,nrow=dataobj$A,ncol=dataobj$noYears)
	tmp<-report[dimnames(report)[[1]]%in%"F",2]
	FFse<-matrix(tmp,nrow=dataobj$A,ncol=dataobj$noYears)

	tmp<-report[dimnames(report)[[1]]%in%"logF",1]
	logF<-matrix(tmp,nrow=dataobj$Am,ncol=dataobj$noYears)
	tmp<-report[dimnames(report)[[1]]%in%"logF",2]
	logFse<-matrix(tmp,nrow=dataobj$Am,ncol=dataobj$noYears)

	if(dataobj$Am<dataobj$A){
	        for(a in dataobj$Am:dataobj$A){
			logF<-rbind(logF,logF[dataobj$Am,])
			logFse<-rbind(logFse,logFse[dataobj$Am,])
		}
	}


	Fbar<-report[dimnames(report)[[1]]%in%"Fbar",1]
	Fbarse<-report[dimnames(report)[[1]]%in%"Fbar",2]

	FbarW<-report[dimnames(report)[[1]]%in%"FbarW",1]
	FbarWse<-report[dimnames(report)[[1]]%in%"FbarW",2]

	ssb<-report[dimnames(report)[[1]]%in%"ssb",1]
	ssbse<-report[dimnames(report)[[1]]%in%"ssb",2]

	totb<-report[dimnames(report)[[1]]%in%"totb",1]
	totbse<-report[dimnames(report)[[1]]%in%"totb",2]

	list(N=N,Nse=Nse,FF=FF,FFse=FFse,Fbar=Fbar,Fbarse=Fbarse,FbarW=FbarW,FbarWse=FbarWse,ssb=ssb,ssbse=ssbse,totb=totb,totbse=totbse,logN=logN,logNse=logNse,logF=logF,logFse=logFse)
}
##############################################
#--------------------------------------------#
#	    FUNCTIONS FOR DIAGNOSTICS		  #
#--------------------------------------------#
##############################################

SummaryPlot0<-function(xp,xo){
	lims<-range(c(xp,xo),na.rm=T)
	plot(xo,xp,xlim=lims,ylim=lims,xlab="Observed",ylab="Predicted")
	x1<-lims[1]-lims[2]
	x2<-lims[2]+lims[2]
	lines(c(x1,x2),c(x1,x2),col=3)
}

SummaryPlot1<-function(a,xp,xo,seo,p=0.95,add=F,xlab="",ylab="",xlim=NULL,...){
	p1<-(1-p)/2
	p2<-p+(1-p)/2
	q1<-qnorm(p1)
	q2<-qnorm(p2)
	x<-c(xp,xo+q1*seo,xo+q2*seo)
	lims<-range(x,na.rm=T)
	if(!add){
		if(is.null(xlim))plot(a,xp,type="o",pch=16,ylim=lims,lwd=2,col=3,xlab=xlab,ylab=ylab)
		else plot(a,xp,type="o",pch=16,ylim=lims,xlim=xlim,lwd=2,col=3,xlab=xlab,ylab=ylab)
	}
	
	lines(a,xo,type="o",pch=16,...)
	lines(a,xo+q1*seo,lty=2,type="o",...)
	lines(a,xo+q2*seo,lty=2,type="o",...)
}

SummaryPlot2<-function(resobj,dataobj,summaryobj,by="year",i,caa=T,index=F,k=1,sdType=1,a=3:12,y=1988:2015,p=0.95,add=F,...){
	if(by=="age"){
		x<-a
		if(caa){
			if(sdType==1)seo<-dataobj$sd_C[,i]
			else seo<-sqrt(resobj$varmatC[,i])
			xp<-resobj$cp[,i]
			xo<-resobj$co[,i]
		}
		else{
			if(sdType==1)seo<-dataobj$sd_I[,i,k]
			else seo<-sqrt(resobj$varmatI[,i,k])
			xp<-resobj$ip[[k]][,i]
			xo<-resobj$io[[k]][,i]
		}
	}

	if(by=="year"){
		x<-y
		if(caa){
			if(sdType==1)seo<-dataobj$sd_C[i,]
			else seo<-sqrt(resobj$varmatC[i,])
			xp<-resobj$cp[i,]
			xo<-resobj$co[i,]
		}
		else{
			if(sdType==1)seo<-dataobj$sd_I[i,,k]
			else seo<-sqrt(resobj$varmatI[i,,k])
			xp<-resobj$ip[[k]][i,]
			xo<-resobj$io[[k]][i,]
		}
	}
	
	if(by=="cohort"){
		amat<-matrix(rep(dataobj$minAge:dataobj$maxAge,dataobj$noYears),nrow=dataobj$A)
		startyear<-min(y)
		ymat<-matrix(rep(1:dataobj$noYears,dataobj$A),nrow=dataobj$A,byrow=T)+startyear-1
		yclmat<-ymat-amat
		cohort<-cohorts[i]
		xlim<-range(a)
		if(caa){
			if(sdType==1)seo<-dataobj$sd_C[yclmat==cohort]
			else seo<-sqrt(resobj$varmatC[yclmat==cohort])
			xp<-resobj$cp[yclmat==cohort]
			xo<-resobj$co[yclmat==cohort]
			x<-amat[yclmat==cohort]
		}
		else{
			xp<-log(summaryobj$N[yclmat==cohort])
			xo<-resobj$logniaa[[k]][yclmat==cohort]
			if(sdType==1)tmpsd<-dataobj$sd_I[,,k]
			else tmpsd<-sqrt(resobj$varmatI[,,k])
			seo<-tmpsd[yclmat==cohort]
			x<-amat[yclmat==cohort]
		}
	}
	else{
		xlim<-NULL
	}

	SummaryPlot1(a=x,xp=xp,xo=xo,seo=seo,p=p,xlim=xlim,add=add,...)
}

RestoreFobj<-function(file){
	load(file,env=.GlobalEnv)

	dataobj<<-FFsummary$dataobj
	pars<<-FFsummary$pars
	optFF<<-FFsummary$optFF
	objFFtab<<-FFsummary$objFFtab
	objFFstats<<-FFsummary$objFFstats
	LLFF<<-FFsummary$LLFF
}


ExtractResiduals<-function(reslist,partab,dataobj){
        cres<-PC(reslist,dataobj)
        co<-log(cres$caa)
        cp<-log(cres$pcaa)
        logcaa<-log(cres$caa)
        cres<-(log(cres$caa)-log(cres$pcaa))
        varmatC<-matrix(NA,nrow=dataobj$dimC,ncol=dataobj$noYears)
        if(sum(dimnames(partab)[[1]]=="loghvec1")==1){
                h<-exp(partab[dimnames(partab)[[1]]=="loghvec1",1])
        }
        else{
                if(sum(dimnames(partab)[[1]]=="loghvec")==1){
                        h<-exp(partab[dimnames(partab)[[1]]=="loghvec",1])[1]
                }
                else h<-1
        }
        for(i in 1:dataobj$noYears)varmatC[,i]<-h*dataobj$sd_C[,i]^2
        varmatC[is.na(co)]<-NA
        wmatC<-1/varmatC
        cvC<-sqrt(exp(varmatC)-1)
        cstdres<-cres/sqrt(varmatC)
        logiaa<-ires<-istdres<-cvI<-logniaa<-list()
        loghnames<-grep("loghvec",dimnames(partab)[[1]])
        loghnames<-dimnames(partab)[[1]][loghnames]
        print(loghnames)
		
        #varmatI<-wmatI<-array(NA,c(dataobj$A,dataobj$noYears,dataobj$nIndices))
	  #varmatI<-wmatI<-array(NA,c(max(dataobj$dimI),dataobj$noYears,dataobj$nIndices))
	  varmatI<-wmatI<-array(NA,c(dataobj$A,dataobj$noYears,dataobj$nIndices))
        io<-ip<-list()
        for(i in 1:dataobj$nIndices){
		agerange<-dataobj$agerangeI[i,]
		ages<-agerange[1]:agerange[2]
		ageindeks<-ages-dataobj$minAge+1
		print(ageindeks)
                ires[[i]]<-PI(fitobj=reslist,partab=partab,dataobj=dataobj,indeks=i)
                io[[i]]<-log(ires[[i]]$iaa)
                ip[[i]]<-log(ires[[i]]$piaa)
                logiaa[[i]]<-log(ires[[i]]$iaa)
                logniaa[[i]]<-log(ires[[i]]$niaa)
                ires[[i]]<-(log(ires[[i]]$iaa)-log(ires[[i]]$piaa))
                #varmat<-matrix(NA,nrow=dataobj$dimI[1],ncol=dataobj$noYears)
                #varmat<-matrix(NA,nrow=dataobj$dimI[i],ncol=dataobj$noYears)
                #varmat<-matrix(NA,nrow=max(dataobj$A),ncol=dataobj$noYears)
                if(length(loghnames)>0){
                        if(length(loghnames)<(i+1))loghnames[i+1]<-loghnames[i]
                        if(sum(dimnames(partab)[[1]]==loghnames[i+1])==1){
                                h<-exp(partab[dimnames(partab)[[1]]==loghnames[i+1],1])
                        }
                        else h<-exp(partab[dimnames(partab)[[1]]=="loghvec",1])[i+1]
                }
                else{
                        h<-1
                }
                #for(y in 1:dataobj$noYears)varmat[,y]<-h*dataobj$sd_I[,y,i]^2
                #for(y in 1:dataobj$noYears)varmat[1:dataobj$dimI[i],y]<-h*dataobj$sd_I[1:dataobj$dimI[i],y,i]^2
                for(y in 1:dataobj$noYears){
				print(y)
				#print(varmatI[1:dataobj$dimI[i],y,i])
				#print(dataobj$sd_I[1:dataobj$dimI[i],y,i])
                        #varmatI[1:dataobj$dimI[i],y,i]<-h*dataobj$sd_I[1:dataobj$dimI[i],y,i]^2
				print(dataobj$sd_I[,y,i])
				#varmatI[ageindeks,y,i]<-h*dataobj$sd_I[ageindeks,y,i]^2
				varmatI[ageindeks,y,i]<-h*dataobj$sd_I[1:length(ageindeks),y,i]^2
				print(varmatI[,y,i])
				print(io[[i]][,y])
                        #varmatI[1:dataobj$dimI[i],y,i][is.na(io[[i]][,y])]<-NA
				#varmatI[1:dataobj$dimI[i],y,i][is.na(io[[i]][dataobj$agerangeI[i,1]:dataobj$agerangeI[i,2],y])]<-NA
				varmatI[ageindeks,y,i][is.na(io[[i]][ageindeks,y])]<-NA
                }
                wmatI[,,i]<-1/varmatI[,,i]
                #cvI[[i]]<-sqrt(exp(varmat)-1)
                cvI[[i]]<-sqrt(exp(varmatI[,,i])-1)
                print("i");print(i)
                print(ires[[i]]);print(varmatI[,,i])
                #istdres[[i]]<-ires[[i]]/sqrt(varmat)
                #istdres[[i]]<-ires[[i]]/sqrt(varmatI[,,i])
			print("blarh")
			#print(ires[[i]][dataobj$agerangeI[i,1]:dataobj$agerangeI[i,2],])
			print(ires[[i]][ageindeks,])
			print("blarh2")
			print(dim(varmatI));print(dataobj$dimI[i])
			#print(varmatI[1:dataobj$dimI[i],,i])
			print(varmatI[ageindeks,,i])
			#istdres[[i]]<-ires[[i]][dataobj$agerangeI[i,1]:dataobj$agerangeI[i,2],]/sqrt(varmatI[1:dataobj$dimI[i],,i])
			#istdres[[i]]<-ires[[i]][ageindeks,]/sqrt(varmatI[ageindeks,,i])
			istdres[[i]]<-ires[[i]]/sqrt(varmatI[,,i])
			print(istdres[[i]])
        }
        list(cres=cres,cstdres=cstdres,cvC=cvC,logcaa=logcaa,ires=ires,istdres=istdres,cvI=cvI,logiaa=logiaa,logniaa=logniaa,varmatC=varmatC,varmatI=varmatI,co=co,cp=cp,io=io,ip=ip,wmatC=wmatC,wmatI=wmatI)
}

PC<-function(fitobj,dataobj){
	N<-fitobj$N
	#print(N)
	FF<-fitobj$FF
	#print(FF)
	A<-dataobj$A
	Y<-dataobj$noYears
	M<-t(dataobj$natMor)
	am<-max(dataobj$keyLogFconst)+1
	if(am<A){
		for(a in (am+1):A){
			FF<-rbind(FF,FF[am,])
		}
	}
	#print(c(A,Y,am))
	caa<-matrix(NA,nrow=A,ncol=Y)
	years<-dataobj$years
	datayears<-as.numeric(dimnames(dataobj$caa)[[1]])
	pcaa<-matrix(NA,nrow=A,ncol=Y)
	for(y in 1:Y){
		pcaa[1:A,y]<-FF[1:A,y]/(FF[1:A,y]+M[1:A,y])*(1-exp(-(FF[1:A,y]+M[1:A,y])))*N[1:A,y]
	}
	caa[,match(datayears,years)]<-t(dataobj$caa)

	list(pcaa=pcaa,caa=caa)
}

PI<-function(fitobj,partab,dataobj,indeks=1){
	N<-fitobj$N
	FF<-fitobj$FF
	A<-dataobj$A
	Y<-dataobj$noYears
	M<-t(dataobj$natMor)
	am<-max(dataobj$keyLogFconst)+1
	if(am<A){
		for(a in (am+1):A){
			FF<-rbind(FF,FF[am,])
		}
	}

	#Get the correct qs
	agerange<-dataobj$agerangeI[indeks,]
	ages<-agerange[1]:agerange[2]
	ageindeks<-ages-dataobj$minAge+1
	#print(ageindeks)
	qpos<-grep("logqpar",dimnames(partab)[[1]])
	qs<-partab[qpos,1]
	agepos<-dataobj$keyLogqpar[indeks,]+1
	#print(agepos)
	#ageindeks<-ageindeks[agepos>0]
	#print(ageindeks)
	agepos<-agepos[agepos>0]
	logqs<-qs[agepos]
	print(logqs)
	deltaY<-dataobj$sampleTimes[indeks+1]
	piaa<-matrix(NA,nrow=A,ncol=Y)
	iaa<-matrix(NA,nrow=A,ncol=Y)
	niaa<-matrix(NA,nrow=A,ncol=Y)
	years<-dataobj$years
	datayears<-as.numeric(dimnames(dataobj$SurveyList[[indeks]])[[1]])
	for(y in 1:length(years)){
		piaa[ageindeks,y]<-exp(logqs)*N[ageindeks,y]*exp(-(FF[ageindeks,y]+M[ageindeks,y])*deltaY)
	}
	iaa[ageindeks,match(datayears,years)]<-t(dataobj$SurveyList[[indeks]])
	iaa[iaa==-999]<-NA
	for(y in 1:length(years)){
		niaa[ageindeks,y]<-iaa[ageindeks,y]/exp(logqs)*exp(+(FF[ageindeks,y]+M[ageindeks,y])*deltaY)
	}

	list(piaa=piaa,iaa=iaa,niaa=niaa)
}

###############
CreateSummaryTab<-function(sdrep,pars){
	stab<-summary(sdrep,"fixed")
	dimnames(stab)[[1]]<-names(pars$ParameterList$FixedPars)[!is.na(pars$MapVector)]
	stab
}

ExtractLikelihoodComponents<-function(sdrep){
	#First total likelihood contributions
	LLS<-summary(sdrep,"report")
	LLS<-LLS[dimnames(LLS)[[1]]%in%"LLS",]
	LLS<-as.numeric(LLS[,1])
	LLS<-cbind(LLS,1:length(LLS))
	LLS<-LLS[LLS[,1]!=0,]
	LLS<-data.frame(LLS)
	dimnames(LLS)[[2]]<-c("negloglik","id")
	LLS$comp<-NA
	LLS$comp[1]<-"Y"
	LLS$comp[2]<-"V"
	LLS$comp[3]<-"TSSEL"
	LLS$comp[4]<-"F"
#	LLS$comp[5]<-"M"
	LLS$comp[5]<-"R"
	LLS$comp[6]<-"CAA"
	LLS$comp[7]<-"CAA-CatchConst"
	LLS$comp[8]<-"CAA-CatchPred"
	LLS$comp[9]<-"IAA_1"
	if(nrow(LLS)>9)LLS$comp[10]<-"IAA_2"
	if(nrow(LLS)>10)LLS$comp[11]<-"IAA_3"
	if(nrow(LLS)>11)LLS$comp[12]<-"IAA_4"
	if(nrow(LLS)>12)LLS$comp[13]<-"IAA_5"
	if(nrow(LLS)>13)LLS$comp[14]<-"IAA_6"
	LLS[2:nrow(LLS),1]<-LLS[2:nrow(LLS),1]-LLS[1:(nrow(LLS)-1),1]
	
	LL1<-summary(sdrep,"report")
	LL1<-LL1[dimnames(LL1)[[1]]%in%"llcaa",1]
	LL1<-c(LL1[1]-sum(LLS[1:5,1]),LL1[2:length(LL1)]-LL1[1:(length(LL1)-1)])
	LLY<-data.frame(LLC=LL1)

	LL1<-summary(sdrep,"report")
	LL1<-LL1[dimnames(LL1)[[1]]%in%"lliaa1",1]
	LL1<-c(LL1[1]-sum(LLS[1:8,1]),LL1[2:length(LL1)]-LL1[1:(length(LL1)-1)])
	LLY$LLI1<-LL1

	if(nrow(LLS)>9){
		LL1<-summary(sdrep,"report")
		LL1<-LL1[dimnames(LL1)[[1]]%in%"lliaa2",1]
		#LL1<-c(LL1[1]-sum(LLS[1:8,1]),LL1[2:length(LL1)]-LL1[1:(length(LL1)-1)])
		LL1<-c(LL1[1]-sum(LLS[1:9,1]),LL1[2:length(LL1)]-LL1[1:(length(LL1)-1)])
		LLY$LLI2<-LL1
	}
	if(nrow(LLS)>10){
		LL1<-summary(sdrep,"report")
		LL1<-LL1[dimnames(LL1)[[1]]%in%"lliaa3",1]
		LL1<-c(LL1[1]-sum(LLS[1:10,1]),LL1[2:length(LL1)]-LL1[1:(length(LL1)-1)])
		LLY$LLI3<-LL1
	}
	if(nrow(LLS)>11){
		LL1<-summary(sdrep,"report")
		LL1<-LL1[dimnames(LL1)[[1]]%in%"lliaa4",1]
		LL1<-c(LL1[1]-sum(LLS[1:11,1]),LL1[2:length(LL1)]-LL1[1:(length(LL1)-1)])
		LLY$LLI4<-LL1
	}
	if(nrow(LLS)>12){
		LL1<-summary(sdrep,"report")
		LL1<-LL1[dimnames(LL1)[[1]]%in%"lliaa5",1]
		LL1<-c(LL1[1]-sum(LLS[1:12,1]),LL1[2:length(LL1)]-LL1[1:(length(LL1)-1)])
		LLY$LLI5<-LL1
	}
	if(nrow(LLS)>13){
		LL1<-summary(sdrep,"report")
		LL1<-LL1[dimnames(LL1)[[1]]%in%"lliaa6",1]
		LL1<-c(LL1[1]-sum(LLS[1:13,1]),LL1[2:length(LL1)]-LL1[1:(length(LL1)-1)])
		LLY$LLI6<-LL1
	}
	#LLY[LLY==0]<-NA
	list(LLS=LLS,LLY=LLY)
}

################################
RunXSAMRetro<-function(XSAMdatalist,XSAMsettings,n,constraints=NULL,control=NULL,sdlist=NULL,StartValFixed=NULL){
	startyear<-XSAMsettings$yearrange[2]
	fits<-list()
	for(i in 0:n){
		assessmentyear<-startyear-i
		print(assessmentyear)
		XSAMsettings$yearrange[2]<-assessmentyear
		if(i>0)XSAMdatalist$caa[paste(assessmentyear),]<--999#NA
		#print(XSAMdatalist$caa)
		#stop()
		dataobj<-CreateDataObject(settings=XSAMsettings,datalist=XSAMdatalist)
		if(!is.null(sdlist)){
			tmp<-sdlist[[1]]
			tmp<-tmp[,1:ncol(dataobj$sd_C)]
			tmp[,ncol(dataobj$sd_C)]<-tmp[,ncol(dataobj$sd_C)-1]
			dataobj$sd_C<-as.matrix(tmp)
			for(j in 2:length(XSAMsettings$fleetIndex)){
				tmp<-sdlist[[j]]
				tmp<-tmp[,1:ncol(dataobj$sd_I[,,j-1])]
				dataobj$sd_I[,,j-1]<-tmp
			}
		}
		dataobj$caa[dataobj$caa<0]<--999
		pars<-CreateParameterList(dataobj,mapnames=NULL)
		if(!is.null(StartValFixed)) pars$ParameterList$FixedPars[pars$ParameterNames%in%pars$ParameterNames[!is.na(pars$MapVector)]]<-StartValFixed
		if(is.null(constraints) & is.null(control))fits[[i+1]]<-FitXSAMFunction(dataobj,pars)
		if(!is.null(constraints) & is.null(control))fits[[i+1]]<-FitXSAMFunction(dataobj,pars,constraints=constraints)
		if(is.null(constraints) & !is.null(control))fits[[i+1]]<-FitXSAMFunction(dataobj,pars,control=control)
		if(!is.null(constraints) & !is.null(control))fits[[i+1]]<-FitXSAMFunction(dataobj,pars,constraints=constraints,control=control)
	}
	fits
}

PlotXSAMfitobj<-function(years,fitobj,what="ssb",p=0.95,add=F,ylim=NULL,xlim=NULL,col=1,lwd=1,plot=T,lty2=3,...){
        if(!what%in%c("ssb","totb","R","Fbar","FbarW"))stop("what must be one of ssb,R,Fbar or FbarW!")
        if(what=="ssb"){
                y<-fitobj$stats$ssb
                yse<-fitobj$stats$ssbse
        }
        if(what=="totb"){
                y<-fitobj$stats$totb
                yse<-fitobj$stats$totbse
        }
        if(what=="R"){
                y<-fitobj$stats$N[1,]
                yse<-fitobj$stats$Nse[1,]
        }
        if(what=="Fbar"){
                y<-fitobj$stats$Fbar
                yse<-fitobj$stats$Fbarse
        }
        if(what=="FbarW"){
                y<-fitobj$stats$FbarW
                yse<-fitobj$stats$FbarWse
        }
        pl<-(1-p)/2
        pu<-1-pl
        yl<-y+qnorm(pl)*yse
        yu<-y+qnorm(pu)*yse
	  yl[yl<0]<-0
        if(plot){
                if(!add & is.null(ylim))ylim<-c(0,max(yu))
                if(!add & is.null(xlim))xlim<-range(years)+c(-1,1)
                #if(!add)plot(years,y,ylim=ylim,xlim=xlim,type="l",col=col,lwd=lwd,xlab=xlab,ylab=ylab,...)
			if(!add)plot(years,y,ylim=ylim,xlim=xlim,type="l",...)
                #else lines(years,y,col=col,lwd=lwd)
			else lines(years,y,col=col)
                #lines(years,yl,col=col,lwd=lwd,lty=lty2)
                #lines(years,yu,col=col,lwd=lwd,lty=lty2)
			lines(years,yl,lty=lty2,col=col)
                lines(years,yu,lty=lty2,col=col)
        }
        invisible(list(y=y,yl=yl,yu=yu))
}

PlotXSAMfitobjOLD<-function(years,fitobj,what="ssb",p=0.95,add=F,ylim=NULL,xlim=NULL,col=1,lwd=1,plot=T,xlab="",ylab=""){
	if(!what%in%c("ssb","totb","R","Fbar","FbarW"))stop("what must be one of ssb,R,Fbar or FbarW!")
	if(what=="ssb"){
		y<-fitobj$stats$ssb
		yse<-fitobj$stats$ssbse
	}
	if(what=="totb"){
		y<-fitobj$stats$totb
		yse<-fitobj$stats$totbse
	}
	if(what=="R"){
		y<-fitobj$stats$N[1,]
		yse<-fitobj$stats$Nse[1,]
	}
	if(what=="Fbar"){
		y<-fitobj$stats$Fbar
		yse<-fitobj$stats$Fbarse
	}
	if(what=="FbarW"){
		y<-fitobj$stats$FbarW
		yse<-fitobj$stats$FbarWse
	}
	pl<-(1-p)/2
	pu<-1-pl
	yl<-y+qnorm(pl)*yse
	yu<-y+qnorm(pu)*yse
	if(plot){
		if(!add & is.null(ylim))ylim<-c(0,max(yu))
		if(!add & is.null(xlim))xlim<-range(years)+c(-1,1)
		if(!add)plot(years,y,ylim=ylim,xlim=xlim,type="l",col=col,lwd=lwd,xlab=xlab,ylab=ylab)
		else lines(years,y,col=col,lwd=lwd)
		lines(years,yl,col=col,lwd=lwd,lty=2)
		lines(years,yu,col=col,lwd=lwd,lty=2)
	}
	invisible(list(y=y,yl=yl,yu=yu))
}

myAIC<-function(ll,npar,k=2)-2*ll + k*npar



mybubble2<-function(res,inches=.25,xlabs=NULL,ylabs=NULL,addScale=T,xtxt=c(-3,-1,0.5),ytxt=1,axes=T){
ylim<-range(row(res))+c(-1,1)
xlim<-range(col(res))+c(-1,1)

if(axes)plot(0,0,xlim=xlim,ylim=ylim,type="n",xlab="Year",ylab="Age",axes=F)
else plot(0,0,xlim=xlim,ylim=ylim,type="n",xlab="",ylab="",axes=F)

if(axes){
	if(is.null(xlabs))axis(1)
	else axis(1,at=(xlim[1]+1):(xlim[2]-1),labels=xlabs)
}
if(axes){
	if(is.null(ylabs))axis(2)
	else axis(2,at=(ylim[1]+1):(ylim[2]-1),labels=ylabs)
}

if(axes)box()

colv<-c("blue","red")
colv<-c("red","blue")
colm<-res
colm<-colv[1]
colm[res<0]<-colv[2]
colm[res>=0]<-colv[1]

tmp<-res
tmp<-abs(res)
radius <- sqrt( tmp/ pi ) 
print(radius)
symbols(col(tmp), row(tmp), circles=c(radius),fg="white",bg=c(colm),add=T,inches=inches)

if(addScale){
rr <- sqrt( 1/ pi ) 
#text(xlim[2]-3,ylim[2]+1,"Scale:",xpd=T)
#symbols(xlim[2]-1, ylim[2]+1, circles=rr,fg="white",bg=colv[1],add=T,inches=F,xpd=NA)
#text(xlim[2]+.5,ylim[2]+1,"=+1",xpd=T)
text(xlim[2]+xtxt[1],ylim[2]+ytxt,"Scale:",xpd=T)
symbols(xlim[2]+xtxt[2], ylim[2]+ytxt, circles=rr,fg="white",bg=colv[1],add=T,inches=F,xpd=NA)
text(xlim[2]+xtxt[3],ylim[2]+ytxt,"=+1",xpd=T)
}
}


###################
error.bar<-function(x, y = NULL, lower, upper, incr = TRUE, bar.ends = TRUE, gap = TRUE, add = FALSE,
        horizontal = FALSE, ..., xlab = deparse(substitute(x)), xlim, ylim,size.bar=NULL)
{
        draw.null.warn <- function(draw, gap)
        {
                if(!any(draw)) {
                        warning("Not enough room for a gap.")
                        draw <- !draw
                        gap <- 0
                }
                invisible(list(draw = draw, gap = gap))
        }
        if(missing(x))
                stop("no data for x or y")
        if(missing(y)) {
                if(missing(xlab))
                        xlab <- "Index"
                y <- x
                x <- time(x)
        }
        n <- length(x)
        if(length(y) != n)
                stop("length of y must equal the length of x")
        center <- if(horizontal) x else y
        if(missing(lower))
                stop("you must provide lower")
        if(length(lower) > 1 && length(lower) != n)
                stop("length of lower must be 1 or equal to the length of x")
        #if incr=T lower is assumed >=0
        if(incr) lower <- center - abs(lower) else lower <- rep(lower, length
                         = n)
        if(any(lower >= center))
                warning(paste(
                        "There are values of 'lower' which are greater or equal to ",
                        if(horizontal) "x" else "y"))
        if(missing(upper))
                upper <- 2 * center - lower
        else {
                if(length(upper) > 1 && length(upper) != n)
                        stop("length of upper must be 1 or equal to the length of x"
                                )
                if(incr)
                        upper <- center + upper
                else upper <- rep(upper, length = n)
        }
        if(any(upper <= center))
                warning(paste(
                        "There are values of 'upper' which are smaller or\nequal to ",
                        if(horizontal) "x" else "y"))
        if(!add)
                if(horizontal) {
                        if(missing(ylim))
                                plot(x, y, xlim = if(missing(xlim)) range(
                                                c(lower, upper), na.rm = TRUE)
                                         else xlim, xlab = xlab, ...)
                        else plot(x, y, xlim = if(missing(xlim)) range(c(lower,
                                                upper), na.rm = TRUE) else xlim,
                                        ylim = ylim, xlab = xlab, ...)
                }
                else {
                        if(missing(xlim))
                                plot(x, y, ylim = if(missing(ylim)) range(
                                                c(lower, upper), na.rm = TRUE)
                                         else ylim, xlab = xlab, ...)
                        else plot(x, y, ylim = if(missing(ylim)) range(c(lower,
                                                upper), na.rm = TRUE) else ylim,
                                        xlim = xlim, xlab = xlab, ...)
                }
        if(horizontal) {
                if(gap)
                        gap <- 0.75 * par("cxy")[1]
                draw <- x - lower > gap
                z <- draw.null.warn(draw, gap)
                draw <- z$draw
                gap <- z$gap
                segments(lower[draw], y[draw], x[draw] - gap, y[draw],...)
                draw <- upper - x > gap
                z <- draw.null.warn(draw, gap)
                draw <- z$draw
                gap <- z$gap
                segments(x[draw] + gap, y[draw], upper[draw], y[draw],...)
                if(bar.ends) {
                        if(is.null(size.bar))size.bar <- par("cxy")[2]
                        segments(lower, y - size.bar, lower, y + size.bar,...)
                        segments(upper, y - size.bar, upper, y + size.bar,...)
                }
        }
        else {
                if(gap)
                        gap <- 0.75 * par("cxy")[2]
                draw <- upper - y > gap
                z <- draw.null.warn(draw, gap)
                draw <- z$draw
                gap <- z$gap
                segments(x[draw], y[draw] + gap, x[draw], upper[draw],...)
                draw <- y - lower > gap
                z <- draw.null.warn(draw, gap)
                draw <- z$draw
                gap <- z$gap
                segments(x[draw], y[draw] - gap, x[draw], lower[draw],...)
                if(bar.ends) {
                        if(is.null(size.bar))size.bar <- par("cxy")[1]
                        segments(x - size.bar, upper, x + size.bar, upper,...)
                        segments(x - size.bar, lower, x + size.bar, lower,...)
                }
        }
}


#######
#More utility functions

#-------------------------------------------------------------------------------#
# MakePlusGroup functions
#-------------------------------------------------------------------------------#
makeNPlusGroup<-function(datatab,pg,keepattr=T){
	k<-attributes(datatab)$k
	sampleTime<-attributes(datatab)$sampleTime
	maxA<-max(as.numeric(colnames(datatab)))
	minA<-min(as.numeric(colnames(datatab)))
	if(maxA>=pg){
		pgnames<-paste(pg:maxA)
		notpgnames<-paste(minA:(pg-1))
		usenames<-paste(minA:pg)
		tmppg<-datatab[,pgnames,drop=F]
		tmppg<-rowSums(tmppg)
		datatab<-cbind(datatab[,notpgnames],tmppg)
		colnames(datatab)<-usenames
		if(keepattr){
			attributes(datatab)$ages<-as.numeric(colnames(datatab))
			attributes(datatab)$years<-as.numeric(rownames(datatab))
			attributes(datatab)$k<-k
			attributes(datatab)$sampleTime<-sampleTime
		}
	}
	datatab
}

makeWPlusGroup<-function(datatab,wtab,pg,keepattr=T){
	maxA<-max(as.numeric(colnames(datatab)))
	minA<-min(as.numeric(colnames(datatab)))
	if(is.null(wtab)){
		wtab<-datatab
		wtab[]<-1
	}
	pgnames<-paste(pg:maxA)
	notpgnames<-paste(minA:(pg-1))
	usenames<-paste(minA:pg)
	tmppg<-datatab[,pgnames,drop=F]
	tmppgW<-wtab[,pgnames,drop=F]
	#tmppg<-rowSums(tmppg*tmppgW)/rowSums(tmppgW)
	#Are 0's placed at the same places?
	z1<-tmppg==0
	z2<-tmppgW==0
	use<-!z1 & !z2
	#print(use)
	#tmpsum1<-rowSums(tmppg*tmppgW)
	#tmpsum2<-rowSums(tmppgW)
	tmpcol<-c()
	for(i in 1:nrow(datatab)){
		tmppguse<-tmppg[i,use[i,]]
		tmppgWuse<-tmppgW[i,use[i,]]
		if(sum(use[i,])>0){
			tmpsum1<-sum(tmppguse*tmppgWuse)
			tmpsum2<-sum(tmppgWuse)
			tmpcol[i]<-tmpsum1/tmpsum2
		}
		else{
			tmpcol[i]<-0
		}
	}
	#tmppg[tmpsum1==0 & tmpsum2==0]<-0
	print(usenames)
	print(notpgnames)
	#datatab<-cbind(datatab[,notpgnames],tmppg)
	datatab<-cbind(datatab[,notpgnames],tmpcol)
	kkk<<-datatab
	colnames(datatab)<-usenames
	if(keepattr){
		attributes(datatab)$ages<-as.numeric(colnames(datatab))
		attributes(datatab)$years<-as.numeric(rownames(datatab))
	}

	datatab
}

CreateXSAMPlusGroups<-function(datalist,pg){
	datalist$weca<-makeWPlusGroup(datalist$weca,datalist$caa,pg)
	datalist$west<-makeWPlusGroup(datalist$west,NULL,pg)
	datalist$caa<-makeNPlusGroup(datalist$caa,pg)
	for(i in 1:length(datalist$SurveyIndices)){
		datalist$SurveyIndices[[i]]<-makeNPlusGroup(datalist$SurveyIndices[[i]],pg)
	}
	datalist
}

addZlines<-function(z,col="lightblue",lty=1){
	hss<-pretty(0:1000,n=500)
	for(i in 1:length(hss))abline(hss[i],z,col=col,lty=lty)
	#for(i in 1:length(hss))abline(hss[i],-.3)
}


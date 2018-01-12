#Specialized code for reading NSS herring data. Structure inherited from WGWIDE history...
#
#Creates a datalist on a format that can be used by XSAM utility functions that creates data objsects that can be passed on to the XSAM template.


#Set up paths and create datalist
path<-"C:/Users/aanes/Documents/XSAMcourse/"

codepath<-paste(path,"XSAMcode/",sep="")
#Sorurce utility functions
source(paste(codepath,"XSAMutils.R",sep=""))

#datapath<-"C:/Users/aanes/Documents/XSAMcourse/Data/herring/"
datapath<-(paste(path,"Data/herring/",sep=""))


#datalist

datalist<-list()

datalist$caa<-ReadDataFiles0(paste(datapath,"canum.txt",sep=""),sep="\t")
attributes(datalist$caa)$k<-1000

datalist$SurveyIndices<-ReadTASACSFleetFile(paste(datapath,"fleet.txt",sep=""),sep="\t")
for(i in 1:length(datalist$SurveyIndices))attributes(datalist$SurveyIndices[[i]])$k<-1/1000
# Manually insert SampleTimes
attributes(datalist$SurveyIndices[[1]])$sampleTime<-2/12
attributes(datalist$SurveyIndices[[2]])$sampleTime<-6/12
attributes(datalist$SurveyIndices[[3]])$sampleTime<-1/12
attributes(datalist$SurveyIndices[[4]])$sampleTime<-5/12
attributes(datalist$SurveyIndices[[5]])$sampleTime<-5/12
attributes(datalist$SurveyIndices[[6]])$sampleTime<-9/12
attributes(datalist$SurveyIndices[[7]])$sampleTime<-9/12

datalist$matprop<-ReadDataFiles0(paste(datapath,"matprop.txt",sep=""))
datalist$west<-ReadDataFiles0(paste(datapath,"west.txt",sep=""),sep="\t")
datalist$weca<-ReadDataFiles0(paste(datapath,"weca.txt",sep=""),sep="\t")
datalist$caton<-ReadDataFiles0(paste(datapath,"caton.txt",sep=""))
attributes(datalist$caton)$k<-1/1000
#-------------------------------------------------------------------------------#
# manually create natural mortality
#-------------------------------------------------------------------------------#
ages<-0:15
years<-1950:2017
natmor<-matrix(rep(c(rep(0.9,3),rep(0.15,13)),length(years)),nrow=length(years),ncol=length(ages),byrow=T)
natmor<-data.frame(natmor)
colnames(natmor)<-min(ages):max(ages)
rownames(natmor)<-min(years):max(years)

#-------------------------------------------------------------------------------#
# manually create  propF
#-------------------------------------------------------------------------------#
propF<-matrix(0,nrow=length(years),ncol=length(ages))
propF<-data.frame(propF)
colnames(propF)<-min(ages):max(ages)
rownames(propF)<-min(years):max(years)

#-------------------------------------------------------------------------------#
# manually create  propM
#-------------------------------------------------------------------------------#
propM<-matrix(0,nrow=length(years),ncol=length(ages))
propM<-data.frame(propM)
colnames(propM)<-min(ages):max(ages)
rownames(propM)<-min(years):max(years)

datalist$natMor<-natmor
datalist$propF<-propF
datalist$propM<-propM

#-------------------------------------------------------------------------------#
# CatchPrediction
#-------------------------------------------------------------------------------#
#
#Read WG's historical catchprediction vs actual recorded catch
CForecast<-read.table(paste(datapath,"catch_forecasts.csv",sep=""),skip=3,header=F,sep=";")
names(CForecast)<-c("y","p","a","pp")
#Manually insert most recent years 
CForecast<-rbind(CForecast,c(2015,283013,328740,NA))
CForecast<-rbind(CForecast,c(2016,376612,NA,NA))
CForecast<-rbind(CForecast,c(2017,805142.4,NA,NA))



if(0)plot(p~a,data=CForecast)
f0<-lm(log(p)~-1+log(a),data=CForecast)
se<-sqrt(var(f0$res))
rse<-sqrt(exp(se^2)-1)

CatchPrediction<-data.frame(Prediction=CForecast$p,caton=CForecast$a,se=se,rse=rse)
rownames(CatchPrediction)<-CForecast$y

datalist$CatchPrediction<-CatchPrediction
attributes(datalist$CatchPrediction)$k<-1/1000



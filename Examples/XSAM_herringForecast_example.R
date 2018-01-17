################################################################
#FORECASTING
################################################################
#Set up parameters in ManagementPlan to do forecasting according to MP
ManagementPlanPars<-list()
ManagementPlanPars$Fs<-NULL
ManagementPlanPars$Bmp<-5000
ManagementPlanPars$Fmp<-0.125
ManagementPlanPars$Fmin<-0.05
ManagementPlanPars$Btrigger<-5000
ManagementPlanPars$Blim<-2500
ManagementPlanPars$Bpa<-NULL
ManagementPlanPars$Fmsy<-0.15
ManagementPlanPars$Flim<-NULL
ManagementPlanPars$Fpa<-0.15
ManagementPlanPars$type<-"MP"
####################################################

source("C:/Users/aanes/Documents/XSAMcourse/XSAMcode/XSAMutils.R")
source("C:/Users/aanes/Documents/XSAMcourse/XSAMcode/XSAMProjectionProgs.R")


#Load the XSAMdataobj,XSAMpars and XSAMfit. They are necessary in the forecast..



#How many years ahead?
np<-20
#How many replicates
n<-1000

#Setup parameters (stock- and catch- weight@age, proportion mature at age) for the prediction
Bpars<-SetBpars(XSAMdataobj,n=np)

#Deterministic projection using TS model for F
fbarrange<-c(XSAMdataobj$Fbarminage,XSAMdataobj$Fbarmaxage)
fc0<-ProjectPoP(XSAMdataobj,XSAMfit$rep,XSAMpars,Bpars,Fbarrange=fbarrange,n=np,TACs=NULL,Ftargets=NULL,DeterministicProjection=T,ConstraintType=3,UseManagementPlan=F,ManagementPlanPars=NULL,forceC=T)
names(fc0)

par(mfrow=c(1,2))
plot(fc0$years,fc0$SSB,type="l")
PlotXSAMfitobj(years=XSAMdataobj$years,fitobj=XSAMfit,what="ssb",p=0.95,add=T,col=3)

plot(fc0$years[1:length(fc0$FbarW)],fc0$FbarW,type="l")
PlotXSAMfitobj(years=XSAMdataobj$years,fitobj=XSAMfit,what="FbarW",p=0.95,add=T,col=3)

plot(fc0$years[1:length(fc0$Ctot)],fc0$Ctot,type="l")
tmp<-summary(XSAMfit$rep,"report")
tmp<-tmp[grep("predW",rownames(tmp)),1]
lines(XSAMdataobj$years,tmp,col=2)
lines(XSAMdataobj$years[1:length(XSAMdataobj$caton)],XSAMdataobj$caton,col=3)
points(2017,XSAMdataobj$CatchPrediction[1],col=3,pch=16)



#Stochastic projection using TS model for F, but assume known parameters..
fc1<-list()
for(i in 1:1000){
	fc1[[i]]<-ProjectPoP(XSAMdataobj,XSAMfit$rep,XSAMpars,Bpars,Fbarrange=fbarrange,n=np,TACs=NULL,Ftargets=NULL,DeterministicProjection=F,ConstraintType=3,UseManagementPlan=F,ManagementPlanPars=NULL,forceC=T)
}

#Stochastic projection using TS model for F, but unknown parameters..

fc2<-list()
parrepl<-GenerateSimDist(XSAMfit$rep,XSAMdataobj,XSAMpars,n=n)
for(i in 1:1000){
	tmprepobj<-list()
	tmprepobj$par.fixed<-parrepl$par.fixed[i,]
	tmprepobj$par.random<-parrepl$par.random[i,]
	fc2[[i]]<-ProjectPoP(XSAMdataobj,tmprepobj,XSAMpars,Bpars,Fbarrange=fbarrange,n=np,TACs=NULL,Ftargets=NULL,DeterministicProjection=F,ConstraintType=3,UseManagementPlan=F,ManagementPlanPars=NULL,forceC=T)
}

ssbs1<-lapply(fc1,FUN=function(x)x$SSB)
ssbs1<-matrix(unlist(ssbs1),ncol=length(ssbs1[[1]]),byrow=T)
ssbs1q<-apply(ssbs1,MAR=2,FUN=quantile,probs=c(0.025,0.975))
ssbs1m<-colMeans(ssbs1)

ssbs2<-lapply(fc2,FUN=function(x)x$SSB)
ssbs2<-matrix(unlist(ssbs2),ncol=length(ssbs2[[1]]),byrow=T)
ssbs2q<-apply(ssbs2,MAR=2,FUN=quantile,probs=c(0.025,0.975))
ssbs2m<-colMeans(ssbs2)

fbars1<-lapply(fc1,FUN=function(x)x$FbarW)
fbars1<-matrix(unlist(fbars1),ncol=length(fbars1[[1]]),byrow=T)
fbars1q<-apply(fbars1,MAR=2,FUN=quantile,probs=c(0.025,0.975))
fbars1m<-colMeans(fbars1)

fbars2<-lapply(fc2,FUN=function(x)x$FbarW)
fbars2<-matrix(unlist(fbars2),ncol=length(fbars2[[1]]),byrow=T)
fbars2q<-apply(fbars2,MAR=2,FUN=quantile,probs=c(0.025,0.975))
fbars2m<-colMeans(fbars2)
#
ctots1<-lapply(fc1,FUN=function(x)x$Ctot)
ctots1<-matrix(unlist(ctots1),ncol=length(ctots1[[1]]),byrow=T)
ctots1q<-apply(ctots1,MAR=2,FUN=quantile,probs=c(0.025,0.975))
ctots1m<-colMeans(ctots1)

ctots2<-lapply(fc2,FUN=function(x)x$Ctot)
ctots2<-matrix(unlist(ctots2),ncol=length(ctots2[[1]]),byrow=T)
ctots2q<-apply(ctots2,MAR=2,FUN=quantile,probs=c(0.025,0.975))
ctots2m<-colMeans(ctots2)

par(mfrow=c(1,3))
plot(fc0$years,fc0$SSB,type="l")
PlotXSAMfitobj(years=XSAMdataobj$years,fitobj=XSAMfit,what="ssb",p=0.95,add=T,col=3)
lines(fc0$years,ssbs1m,col=2)
lines(fc0$years,ssbs1q[1,],col=2,lty=3)
lines(fc0$years,ssbs1q[2,],col=2,lty=3)
lines(fc0$years,ssbs2q[1,],col=4,lty=3)
lines(fc0$years,ssbs2q[2,],col=4,lty=3)


plot(fc0$years[1:length(fc0$FbarW)],fc0$FbarW,type="l")
PlotXSAMfitobj(years=XSAMdataobj$years,fitobj=XSAMfit,what="FbarW",p=0.95,add=T,col=3)
lines(fc0$years[1:length(fc0$FbarW)],fbars1m,col=2)
lines(fc0$years[1:length(fc0$FbarW)],fbars1q[1,],col=2,lty=3)
lines(fc0$years[1:length(fc0$FbarW)],fbars1q[2,],col=2,lty=3)
lines(fc0$years[1:length(fc0$FbarW)],fbars2q[1,],col=4,lty=3)
lines(fc0$years[1:length(fc0$FbarW)],fbars2q[2,],col=4,lty=3)


plot(fc0$years[1:length(fc0$Ctot)],fc0$Ctot,type="l")
tmp<-summary(XSAMfit$rep,"report")
tmp<-tmp[grep("predW",rownames(tmp)),1]
lines(XSAMdataobj$years,tmp,col=2)
lines(XSAMdataobj$years[1:length(XSAMdataobj$caton)],XSAMdataobj$caton,col=3)
points(2017,XSAMdataobj$CatchPrediction[1],col=3,pch=16)
lines(fc0$years[1:length(fc0$FbarW)],ctots1m,col=2)
lines(fc0$years[1:length(fc0$FbarW)],ctots1q[1,],col=2,lty=3)
lines(fc0$years[1:length(fc0$FbarW)],ctots1q[2,],col=2,lty=3)
lines(fc0$years[1:length(fc0$FbarW)],ctots2q[1,],col=4,lty=3)
lines(fc0$years[1:length(fc0$FbarW)],ctots2q[2,],col=4,lty=3)

####################################################################################################

#Projections using management plan!
#Function adapted to herring
XSAMForecastFMP<-DoXSAMForecast(dataobj=XSAMdataobj,fitobj=XSAMfit,parobj=XSAMpars,StochasticForcast=T,Bpars=NULL,Fbarrange=NULL,np=np,n=n,ManagementPlanPars,BparsType=2,weightedF=T,forceC=T)

ManagementPlanPars$type<-"MSY"
XSAMForecastFMSY<-DoXSAMForecast(dataobj=XSAMdataobj,fitobj=XSAMfit,parobj=XSAMpars,StochasticForcast=T,Bpars=NULL,Fbarrange=NULL,np=np,n=n,ManagementPlanPars,BparsType=2,weightedF=T,forceC=T)

ManagementPlanPars$type<-"setF"
ManagementPlanPars$Fset<-0.15
XSAMForecastFpa<-DoXSAMForecast(dataobj=XSAMdataobj,fitobj=XSAMfit,parobj=XSAMpars,StochasticForcast=T,Bpars=NULL,Fbarrange=NULL,np=np,n=n,ManagementPlanPars,BparsType=2,weightedF=T,forceC=T)

ManagementPlanPars$type<-"setF"
ManagementPlanPars$Fset<-0
XSAMForecastF0<-DoXSAMForecast(dataobj=XSAMdataobj,fitobj=XSAMfit,parobj=XSAMpars,StochasticForcast=T,Bpars=NULL,Fbarrange=NULL,np=np,n=n,ManagementPlanPars,BparsType=2,weightedF=T,forceC=T)

#Fish in 2018 such that SSB(2019)=Blim... Iterative procedure for finding the correct F: here manual search for F...
ManagementPlanPars$type<-"setF"
#ManagementPlanPars$Fset<-0.605735
ManagementPlanPars$Fset<-0.4656
XSAMForecastF25<-DoXSAMForecast(dataobj=XSAMdataobj,fitobj=XSAMfit,parobj=XSAMpars,StochasticForcast=T,Bpars=NULL,Fbarrange=NULL,np=np,n=n,ManagementPlanPars,BparsType=2,weightedF=T,forceC=T)
objF25<-ManagementOptionsTable(XSAMForecastF25,XSAMdataobj,weightedF=T)
printManagementOptionsTable(objF25,printBasis=T)

########-------------------------------------------#######

#Extract some stats on SSB
hmp<-XSAMForecastStats(FC=XSAMForecastFMP,fitobj=XSAMfit,dataobj=XSAMdataobj,what="ssb",p=0.95,StatsType=3)
hmsy<-XSAMForecastStats(FC=XSAMForecastFMSY,fitobj=XSAMfit,dataobj=XSAMdataobj,what="ssb",p=0.95,StatsType=3)
hfpa<-XSAMForecastStats(FC=XSAMForecastFpa,fitobj=XSAMfit,dataobj=XSAMdataobj,what="ssb",p=0.95,StatsType=3)
hf0<-XSAMForecastStats(FC=XSAMForecastF0,fitobj=XSAMfit,dataobj=XSAMdataobj,what="ssb",p=0.95,StatsType=3)
hblim<-XSAMForecastStats(FC=XSAMForecastF25,fitobj=XSAMfit,dataobj=XSAMdataobj,what="ssb",p=0.95,StatsType=3)


#Compare effect of options on SSB
xlim<-c(1988,2020)
PlotXSAMForecastStats(hmp,xlab="Year",ylab="SSB",xlim=xlim,ylim=c(2000,7500),lwd=2)
lines(c(0,1e4),c(5,5)*1e3,col="gray",lwd=2,lty=3)
PlotXSAMForecastStats(hmsy,add=T,col=2,lwd=2)
PlotXSAMForecastStats(hfpa,add=T,col=3,lwd=2)
PlotXSAMForecastStats(hf0,add=T,col=4,lwd=2)
PlotXSAMForecastStats(hblim,add=T,col=5,lwd=2)
PlotXSAMForecastStats(hmp,add=T,col=1,lwd=2)
legend("bottomleft",lty=1,legend=c("MP","MSY",expression(F==F[pa]),"F=0",expression(B[lim])),title="Rationale:",col=1:5,bty="n",lwd=2)
lines(c(0,1e4),c(5,5)*1e3,col="gray",lwd=2,lty=3)
#lines(c(0,1e4),c(4360,4360))

#############
#Extract some stats on F
hmp<-XSAMForecastStats(FC=XSAMForecastFMP,fitobj=XSAMfit,dataobj=XSAMdataobj,what="FbarW",p=0.95,StatsType=3)
hmsy<-XSAMForecastStats(FC=XSAMForecastFMSY,fitobj=XSAMfit,dataobj=XSAMdataobj,what="FbarW",p=0.95,StatsType=3)
hfpa<-XSAMForecastStats(FC=XSAMForecastFpa,fitobj=XSAMfit,dataobj=XSAMdataobj,what="FbarW",p=0.95,StatsType=3)
hf0<-XSAMForecastStats(FC=XSAMForecastF0,fitobj=XSAMfit,dataobj=XSAMdataobj,what="FbarW",p=0.95,StatsType=3)
hblim<-XSAMForecastStats(FC=XSAMForecastF25,fitobj=XSAMfit,dataobj=XSAMdataobj,what="FbarW",p=0.95,StatsType=3)


#Compare effect of options on F
PlotXSAMForecastStats(hmp,xlab="Year",ylab="F",xlim=xlim,ylim=c(0,0.5),lwd=2)
lines(c(0,1e4),c(.125,.125),col="gray",lwd=2,lty=3)
PlotXSAMForecastStats(hmsy,add=T,col=2,lwd=2)
PlotXSAMForecastStats(hfpa,add=T,col=3,lwd=2)
PlotXSAMForecastStats(hf0,add=T,col=4,lwd=2)
PlotXSAMForecastStats(hblim,add=T,col=5,lwd=2)
PlotXSAMForecastStats(hmp,add=T,col=1,lwd=2)
legend("bottomleft",lty=1,legend=c("MP","MSY",expression(F==F[pa]),"F=0",expression(B[lim])),title="Rationale:",col=1:5,bty="n",lwd=2)

#################
printManagementOptionsTable(ManagementOptionsTable(XSAMForecastFMP,XSAMdataobj,weightedF=T),printBasis=T)

objMP<-ManagementOptionsTable(XSAMForecastFMP,XSAMdataobj,weightedF=T)
printManagementOptionsTable(objMP,printBasis=T)

#AND THE OTHER OPTIONS
objMSY<-ManagementOptionsTable(XSAMForecastFMSY,XSAMdataobj,weightedF=T)
printManagementOptionsTable(objMSY,printBasis=T)

objFpa<-ManagementOptionsTable(XSAMForecastFpa,XSAMdataobj,weightedF=T)
printManagementOptionsTable(objFpa,printBasis=T)

objF0<-ManagementOptionsTable(XSAMForecastF0,XSAMdataobj,weightedF=T)
printManagementOptionsTable(objF0,printBasis=T)

objF25<-ManagementOptionsTable(XSAMForecastF25,XSAMdataobj,weightedF=T)
printManagementOptionsTable(objF25,printBasis=T)

printManagementOptionsTable(objMP,printBasis=T)
printManagementOptionsTable(objMSY,printBasis=T)
printManagementOptionsTable(objFpa,printBasis=T)
printManagementOptionsTable(objF0,printBasis=T)
printManagementOptionsTable(objFsq,printBasis=T)
printManagementOptionsTable(objF25,printBasis=T)	
#####################
#Selection patterns in the forecast
dim(XSAMForecastFMP$FC0$FF)
length(XSAMForecastFMP$years)
y<-XSAMForecastFMP$years
y<-y[1:(length(y)-1)]
x<-XSAMForecastFMP$FC0$FF[,y==2014]

x<-XSAMForecastFMP$FC0$FF[,y==2014]
plot(2:12,x/max(x),type="l",xlab="Age",ylab="Selection")
title("Exploitation pattern")
x<-XSAMForecastFMP$FC0$FF[,y==2015]
lines(2:12,x/max(x),type="l",col=2)
x<-XSAMForecastFMP$FC0$FF[,y==2016]
lines(2:12,x/max(x),type="l",col=3)
x<-XSAMForecastFMP$FC0$FF[,y==2017]
lines(2:12,x/max(x),type="l",col=4)
x<-XSAMForecastFMP$FC0$FF[,y==2018]
lines(2:12,x/max(x),type="l",col=5)

legend("topleft",lty=1,legend=paste(2014:2018),title="Years:",col=1:5,bty="n")

#Herring example

setwd("C:/Users/aanes/Documents/XSAMcourse/")
mypath<-getwd()

codepath<-paste(mypath,"/XSAMcode/",sep="")

#------------------------------------#
#Read the data
#------------------------------------#
#Read data the herring data

source(paste(mypath,"/Examples/ReadHerringData.R",sep=""))
#Data now in datalist
datalist
#age range in original data is 0-15, no plusgroup

#------------------------------------#
#create plusgroup 12
#------------------------------------#
datalist12<-CreateXSAMPlusGroups(datalist,12)
CheckDataList(datalist12)

#Make sure that the plusgroup sum is OK..
rowSums(datalist$caa[,paste(12:15)])
datalist12$caa[,"12"]

rowSums(datalist$SurveyIndices[[1]][,paste(12:15)])
datalist12$SurveyIndices[[1]][,"12"]

datalist$SurveyIndices[[4]]
datalist12$SurveyIndices[[4]]

rowSums(datalist$SurveyIndices[[5]][,paste(12:15)])
datalist12$SurveyIndices[[5]][,"12"]


#------------------------------------#
#Read model configuration
#------------------------------------#
setuppath<-paste(mypath,"/XSAMsetup/",sep="")
XSAMsettings<-ReadXSAMsettings(paste(setuppath,"XSAMsetupNSSherring.txt",sep=""))
XSAMsettings
#---------------------------------------------------#
# Create data object to be passed to the C template
#---------------------------------------------------#

#Create the XSAM data object
XSAMdataobj<-CreateDataObject(settings=XSAMsettings,datalist=datalist12)
names(XSAMdataobj)

#Default error structures
XSAMdataobj$sd_C
XSAMdataobj$R_C

dim(XSAMdataobj$sd_I)
dim(XSAMdataobj$R_I)
XSAMdataobj$sd_I
XSAMdataobj$R_I

#---------------------------------------------------#
#Replacing error structures with external values
#---------------------------------------------------#
sdpath<-paste(mypath,"/Data/herring/sds/",sep="")

#-------------#
# CATCH DATA
#-------------#
#load SEs for
tmp<-load(paste(sdpath,"sd_C1_12_88_16",sep=""))#smoothed values
tmp1<-load(paste(sdpath,"sd_Cemp1_12_88_16",sep=""))#empirical values
plot(sd_C1_12_88_16,sd_Cemp1_12_88_16)
plot(log(sd_C1_12_88_16),log(sd_Cemp1_12_88_16))

tmpsd<-sd_C1_12_88_16[2:12,]
#Setting last year (i.e. 2017) equal second last year (i.e. 2017) since no data from 2017 yet...
tmpsd[,30]<-tmpsd[,29]
dim(XSAMdataobj$sd_C)
dim(tmpsd)
XSAMdataobj$sd_C<-tmpsd
#Set last years sd in catch data equal to average
XSAMdataobj$sd_C[,ncol(XSAMdataobj$sd_C)]<-rowMeans(tmpsd[,1:(ncol(tmpsd)-1)])
#CVs
sqrt(exp(XSAMdataobj$sd_C^2)-1)

#-------------#
# Fleet 1
#-------------#

tmp<-load(paste(sdpath,"sd_I_fleet1NEW_3_12_88_17",sep=""))
tmp1<-load(paste(sdpath,"sd_Iemp_fleet1NEW_3_12_88_17",sep=""))
plot(sd_I_fleet1NEW_3_12_88_17,sd_Iemp_fleet1NEW_3_12_88_17)
plot(log(sd_I_fleet1NEW_3_12_88_17),log(sd_Iemp_fleet1NEW_3_12_88_17))

tmpsd<-sd_I_fleet1NEW_3_12_88_17
dim(XSAMdataobj$sd_I[,,1])
dim(tmpsd)
XSAMdataobj$sd_I[,,1]<-tmpsd
#CVs
sqrt(exp(XSAMdataobj$sd_I[,,1]^2)-1)

#-------------#
# Fleet 4
#-------------#

tmp<-load(paste(sdpath,"sd_I_fleet4NEW_1_2_88_17",sep=""))
tmp1<-load(paste(sdpath,"sd_Iemp_fleet4NEW_1_2_88_17",sep=""))
plot(sd_I_fleet4NEW_1_2_88_17,sd_Iemp_fleet4NEW_1_2_88_17)
plot(log(sd_I_fleet4NEW_1_2_88_17),log(sd_Iemp_fleet4NEW_1_2_88_17))

tmpsd<-sd_I_fleet4NEW_1_2_88_17
dim(XSAMdataobj$sd_I[,,2])
dim(tmpsd)
XSAMdataobj$sd_I[,,2][1,]<-tmpsd[1,]
#CVs
sqrt(exp(XSAMdataobj$sd_I[,,2]^2)-1)

#-------------#
# Fleet 5
#-------------#

tmp<-load(paste(sdpath,"sd_I_fleet5_3_12_88_17",sep=""))
tmp1<-load(paste(sdpath,"sd_Iemp_fleet5_3_12_88_17",sep=""))

tmpsd<-sd_I_fleet5_3_12_88_17
tmpsd1<-sd_Iemp_fleet5_3_12_88_17
plot(tmpsd,tmpsd1)
plot(log(tmpsd),log(tmpsd1))

dim(XSAMdataobj$sd_I[,,3])
dim(tmpsd)
XSAMdataobj$sd_I[,,3]<-tmpsd
#CVs
sqrt(exp(XSAMdataobj$sd_I[,,3]^2)-1)

#------------------------------------#
#Set up parameterlist
#------------------------------------#
XSAMpars<-CreateParameterList(XSAMdataobj,mapnames=NULL)
XSAMpars

#------------------------------------#
#Constraining parameters
#------------------------------------#
#Constraints on parameters
constraints<-list()
constraints[[1]]<-list()
constraints[[1]]$par<-"logs2_2"
constraints[[1]]$lower<--5
constraints[[1]]$upper<-Inf

constraints[[2]]<-list()
constraints[[2]]$par<-"logs2_1"
constraints[[2]]$lower<--5
constraints[[2]]$upper<-Inf

constraints[[3]]<-list()
constraints[[3]]$par<-"logs2_4"
constraints[[3]]$lower<--5
constraints[[3]]$upper<-Inf

#These once are not necessary to constrain but...
constraints[[4]]<-list()
constraints[[4]]$par<-"bY"
constraints[[4]]$lower<-0
constraints[[4]]$upper<-1

constraints[[5]]<-list()
constraints[[5]]$par<-"bU"
constraints[[5]]$lower<-0
constraints[[5]]$upper<-1

#------------------------------------#
#Fit the model...
#------------------------------------#

require(TMB)
compile(paste(codepath,"XSAM.cpp",sep=""))
dyn.load(dynlib(paste(codepath,"XSAM",sep="")))
XSAMfit<-FitXSAMFunction(XSAMdataobj,XSAMpars,returnObj=T,constraints=constraints,control=list(iter.max=500,eval.max=500))

#Convergence?
XSAMfit$optobj
XSAMfit$convergence
#------------------------------------#
#Table of parameters
#------------------------------------#
tmptab<-XSAMfit$tab
tmptab<-as.data.frame(tmptab)
tmptab
#insert CV's
tmptab$cv<-abs(tmptab[,2]/tmptab[,1])
round(tmptab,3)

#------------------------------------#
#Correlation of parameters
#------------------------------------#

mybubble2(cov2cor(XSAMfit$rep$cov.fixed),inches=F,ylab="",xlab="",axes=F)
axis(1,at=1:44,labels=dimnames(XSAMfit$tab)[[1]],las=2,cex.axis=0.7)
axis(2,at=1:44,labels=dimnames(XSAMfit$tab)[[1]],las=2,cex.axis=0.7)
box()

#------------------------------------#
#Weights
#------------------------------------#
years<-XSAMdataobj$years
ages<-XSAMdataobj$minAge:XSAMdataobj$maxAge
par(mfrow=c(2,2),mar=c(4,3,2,1))
mybubble2(XSAMfit$res$wmatC,xlabs=years,ylabs=ages,inches=F)
title("Catch@age")
#Residuals first fleet (Fleet 1)
mybubble2(XSAMfit$res$wmatI[,,1],xlabs=years,ylabs=ages,inches=F)
title("Fleet 1")
#Residuals second fleet (Fleet 4)
mybubble2(XSAMfit$res$wmatI[,,2],xlabs=years,ylabs=ages,inches=F)
title("Fleet 4")
#Residuals second fleet (Fleet 5)
mybubble2(XSAMfit$res$wmatI[,,3],xlabs=years,ylabs=ages,inches=F)
title("Fleet 5")

#Comparing input variance vs model scaled variances
plot(XSAMfit$res$varmatC,XSAMdataobj$sd_C^2,xlab="h*var(Chat)",ylab="var(Chat)")
plot(XSAMfit$res$varmatI[2:11,,1],XSAMdataobj$sd_I[,,1]^2,xlab="h*var(Ihat)",ylab="var(Ihat)")
plot(XSAMfit$res$varmatI[1,,2],XSAMdataobj$sd_I[1,,2]^2,xlab="h*var(Ihat)",ylab="var(Ihat)")
plot(XSAMfit$res$varmatI[2:11,,3],XSAMdataobj$sd_I[,,3]^2,xlab="h*var(Ihat)",ylab="var(Ihat)")

#Resultant model scaled CV's
round(sqrt(exp(XSAMfit$res$varmatC)-1),3)
round(sqrt(exp(XSAMfit$res$varmatI[,,1])-1),3)
round(sqrt(exp(XSAMfit$res$varmatI[,,2])-1),3)
round(sqrt(exp(XSAMfit$res$varmatI[,,3])-1),3)


#------------------------------------#
#Plot parameter estimates
#------------------------------------#
#---------------------#
# first catchabilities...
#---------------------#

#First extract estimates...
#First survey
logqs<-XSAMfit$tab[grep("logqpar",rownames(XSAMfit$tab)),]
XSAMsettings$keyLogqpar
idx<-XSAMsettings$keyLogqpar[1,]
idx<-idx[!is.na(idx)]
logq1<-logqs[idx,]
a1<-XSAMsettings$agerangeI[1,]
a1<-a1[1]:a1[2]

#Second survey
idx<-XSAMsettings$keyLogqpar[2,]
idx<-idx[!is.na(idx)]
logq2<-logqs[idx,]
a2<-XSAMsettings$agerangeI[2,]
a2<-a2[1]:a2[2]

#Third survey
idx<-XSAMsettings$keyLogqpar[3,]
idx<-idx[!is.na(idx)]
logq3<-logqs[idx,]
a3<-XSAMsettings$agerangeI[3,]
a3<-a3[1]:a3[2]

#Then plot
par(mfrow=c(1,3),pty='s')
m1<-logq1[,1]
l1<-m1-2*logq1[,2]
u1<-m1+2*logq1[,2]
plot(a1,m1,type="o",pch=16,ylim=range(c(l1,u1)))
lines(a1,l1,lty=3)
lines(a1,u1,lty=3)
title("Fleet 1")

m2<-logq2[1]
l2<-m2-2*logq2[2]
u2<-m2+2*logq2[2]
plot(a2,m2,type="o",pch=16,ylim=range(c(l2,u2)))
points(a2,l2)
points(a2,u2)
title("Fleet 4")

m3<-logq3[,1]
l3<-m3-2*logq3[,2]
u3<-m3+2*logq3[,2]
plot(a3,m3,type="o",pch=16,ylim=range(c(l3,u3)))
lines(a3,l3,lty=3)
lines(a3,u3,lty=3)
title("Fleet 5")

#---------------------#
# then selection patterns...
#---------------------#

#Extract the U's from the latent variable estimates
names(XSAMfit$obj)
summary(XSAMfit$rep,"random")
lpmat<-summary(XSAMfit$rep,"random")[,1]
lpmat<-matrix(lpmat,ncol=XSAMdataobj$noYears)

lpsdmat<-summary(XSAMfit$rep,"random")[,2]
lpsdmat<-matrix(lpsdmat,ncol=XSAMdataobj$noYears)

fmat<-lpmat[1:XSAMdataobj$Am,]
evec<-lpmat[XSAMdataobj$Am+XSAMdataobj$uam+1,]
sdevec<-lpsdmat[XSAMdataobj$Am+XSAMdataobj$uam+1,]

XSAMfit$stats$FF
umat<-lpmat[(XSAMdataobj$Am+1):(XSAMdataobj$Am+XSAMdataobj$uam),]
colSums(umat)
umat<-rbind(umat,-colSums(umat))
colSums(umat)

persp(y=2:11,x=1988:2017,z=t(umat),shade=0.5,theta=-50,phi=20,ticktype="detailed",xlab="Year",ylab="Age")
persp(y=2:11,x=1988:2017,z=t(umat),shade=0.5,theta=50,phi=40,ticktype="detailed",xlab="Year",ylab="Age")

persp(y=2:11,x=1988:2017,z=t(umat),shade=0.5,theta=30,phi=20,ticktype="detailed",xlab="Year",ylab="Age")

fmats<-exp(fmat)
for(i in 1:ncol(fmats))fmats[,i]<-fmats[,i]/max(fmats[,i])

#Selection pattern
persp(y=2:11,x=1988:2017,z=t(umat),shade=0.5,theta=50,phi=30,ticktype="detailed",xlab="Year",ylab="Age")
#Fishing pattern
persp(y=2:11,x=1988:2017,z=t(fmats),shade=0.5,theta=50,phi=30,ticktype="detailed",xlab="Year",ylab="Age")

#Another way of plotting selection
plot(1988:2017,umat[1,],ylim=range(umat),type="l")
for(i in 1:nrow(umat)){
	lines(1988:2017,umat[i,])
	text(x=1988:2017,y=umat[i,],i+1,cex=.5)
}


#Effort
plot(1988:2017,evec,type="l")
lines(1988:2017,evec-2*sdevec,lty=3)
lines(1988:2017,evec+2*sdevec,lty=3)

#-------------------------------------------------------------------#
# MORE DIAGNOSTICS
#-------------------------------------------------------------------#
#--------------------------#
#'RAW' residuals
#Residuals catch at age

ages<-XSAMdataobj$minAge:XSAMdataobj$maxAge
years<-XSAMdataobj$years
mybubble2(XSAMfit$res$cres,xlabs=years,ylabs=ages,inches=F)
title("Catch@age")
#Residuals first fleet (Fleet 1)
mybubble2(XSAMfit$res$ires[[1]],xlabs=years,ylabs=ages,inches=F)
#Residuals second fleet (Fleet 4)
mybubble2(XSAMfit$res$ires[[2]],xlabs=years,ylabs=ages,inches=F)
#Residuals second fleet (Fleet 5)
mybubble2(XSAMfit$res$ires[[3]],xlabs=years,ylabs=ages,inches=F)

#--------------------------#
#Standardised residuals
par(mfrow=c(2,2),mar=c(4,3,2,1))
#Residuals catch at age
mybubble2(XSAMfit$res$cstdres,xlabs=years,ylabs=ages,inches=F)
title("Catch@age r1")

##
source(paste(mypath,"/XSAMcode/XSAMProjectionProgs.R",sep=""))
r1<-OneSampleApproach(dataobj=XSAMdataobj,fitobj=XSAMfit,parobj=XSAMpars)
mybubble2(r1$cstdres,xlabs=years,ylabs=ages,inches=F)
title("Catch@age r2")

#Residuals first fleet (Fleet 1)
mybubble2(XSAMfit$res$istdres[[1]],xlabs=years,ylabs=ages,inches=F)
title("Fleet 1 r1")
mybubble2(r1$istdres[[1]],xlabs=years,ylabs=ages,inches=F)
title("Fleet 1 r2")

par(mfrow=c(2,2),mar=c(4,3,2,1))
#Residuals second fleet (Fleet 4)
mybubble2(XSAMfit$res$istdres[[2]],xlabs=years,ylabs=ages,inches=F)
title("Fleet 4 r1")
mybubble2(r1$istdres[[2]],xlabs=years,ylabs=ages,inches=F)
title("Fleet 4 r2")

#Residuals second fleet (Fleet 5)
mybubble2(XSAMfit$res$istdres[[3]],xlabs=years,ylabs=ages,inches=F)
title("Fleet 5 r1")
mybubble2(r1$istdres[[3]],xlabs=years,ylabs=ages,inches=F)
title("Fleet 5 r2")


#########################################################
#qq-plots
#------------------------------
#Relevant for independent variables, but be CAREFUL for correlated variables!!
#------------------------------
par(mfrow=c(2,3),pty='s',mar=c(4,4,1,1),oma=c(1,1,3,0))
#Catches
SummaryPlot0(XSAMfit$res$cp,XSAMfit$res$co)
qqnorm(XSAMfit$res$cstdres,main="Norm Q-Q Plot; std res1")
qqline(XSAMfit$res$cstdres)
title("CAA",outer=T,line=0)
qqnorm(r1$cstdres,main="Norm Q-Q Plot; std res2")
qqline(r1$cstdres)

#Fleet 1
SummaryPlot0(XSAMfit$res$ip[[1]],XSAMfit$res$io[[1]])
qqnorm(XSAMfit$res$istdres[[1]],main="Norm Q-Q Plot; std res1")
qqline(XSAMfit$res$istdres[[1]])
title("Fleet 1",outer=T,line=-25)
qqnorm(r1$istdres[[1]],main="Norm Q-Q Plot; std res2")
qqline(r1$istdres[[1]])

par(mfrow=c(2,3),pty='s',mar=c(4,4,1,1),oma=c(1,1,3,0))
#Fleet 4
SummaryPlot0(XSAMfit$res$ip[[2]],XSAMfit$res$io[[2]])
qqnorm(XSAMfit$res$istdres[[2]],main="Norm Q-Q Plot; std res1")
qqline(XSAMfit$res$istdres[[2]])
title("Fleet 4",outer=T,line=0)
qqnorm(r1$istdres[[2]],main="Norm Q-Q Plot; std res2")
qqline(r1$istdres[[2]])


#Fleet 5
SummaryPlot0(XSAMfit$res$ip[[3]],XSAMfit$res$io[[3]])
qqnorm(XSAMfit$res$istdres[[3]],main="Norm Q-Q Plot; std res1")
qqline(XSAMfit$res$istdres[[3]])
title("Fleet 5",outer=T,line=-25)
qqnorm(r1$istdres[[3]],main="Norm Q-Q Plot; std res2")
qqline(r1$istdres[[3]])

#-------------------------------------------------------#
#Compare survey estimates of biomass with assessment
#survey estimates are scaled with assessment catchabilities...
#-------------------------------------------------------#
years<-XSAMdataobj$years
#Calculate survey SSB
ssb1<-ssb5<-c()
for(i in 1:length(years)){
	#ssb fleet 1
	ssb1[i]<-sum(exp(XSAMfit$res$logniaa[[1]][,i])*XSAMdataobj$stockMeanWeight[i,]*XSAMdataobj$propMat[i,1:11],na.rm=T)
	#ssb fleet 5
	ssb5[i]<-sum(exp(XSAMfit$res$logniaa[[3]][,i])*XSAMdataobj$stockMeanWeight[i,]*XSAMdataobj$propMat[i,1:11],na.rm=T)
}
ssb1[ssb1==0]<-NA
ssb5[ssb5==0]<-NA
#Assessment SSB
plot(XSAMdataobj$years,XSAMfit$stats$ssb,type="l",ylim=c(0,10000),xlab="Year",ylab="SSB")
#lines(c(-100,3000),c(5000,5000))
#ADD survey SSBs
lines(years,ssb1,col=2)
lines(years,ssb5,col=3)
legend("topleft",lty=1,legend=c("Model","Fleet 1","Fleet 5"),col=1:3,bty="n")

##
#Adding approximate error measures 
#Recall: for a log normal distributed variable: SD on original scale=mu*cv on log scale. Therefore approximately
sd1<-exp(XSAMfit$res$logniaa[[1]])*XSAMfit$res$cvI[[1]]
sd2<-exp(XSAMfit$res$logniaa[[3]])*XSAMfit$res$cvI[[3]]

#and assuming no correlation...
v1<-v2<-c()
for(i in 1:length(years)){
	v1[i]<-sum((XSAMdataobj$stockMeanWeight[i,]*XSAMdataobj$propMat[i,1:11])^2*sd1[,i]^2,na.rm=T)
	v2[i]<-sum((XSAMdataobj$stockMeanWeight[i,]*XSAMdataobj$propMat[i,1:11])^2*sd2[,i]^2,na.rm=T)
}
v1[v1==0]<-NA

plot(years,XSAMfit$stats$ssb,type="l",ylim=c(0,13000),xlab="Year",ylab="SSB",axes=F)
axis(1)
box()
lines(years,XSAMfit$stats$ssb+qnorm(0.025)*XSAMfit$stats$ssbse,lty=3)
lines(years,XSAMfit$stats$ssb+qnorm(0.975)*XSAMfit$stats$ssbse,lty=3)
#Add SSB from Fleet 1
lines(years,ssb1,col=2)
lines(years,ssb1+qnorm(0.025)*sqrt(v1),col=2,lty=3)
lines(years,ssb1+qnorm(0.975)*sqrt(v1),col=2,lty=3)
lines(years,ssb5,col=3)
lines(years,ssb5+qnorm(0.025)*sqrt(v2),col=3,lty=3)
lines(years,ssb5+qnorm(0.975)*sqrt(v2),col=3,lty=3)
legend("topleft",lty=1,legend=c("Model","Fleet 1","Fleet 5"),col=1:3,bty="n")

#---------------------------------------------------------------#
#      RETRO RUN
#---------------------------------------------------------------#

#Prepare input to retro run
sdlist<-list(XSAMdataobj$sd_C,XSAMdataobj$sd_I[,,1],XSAMdataobj$sd_I[,,2],XSAMdataobj$sd_I[,,3])

#Number of years in retro
n<-5
retrorun<-RunXSAMRetro(datalist12,XSAMsettings,n=n,constraints=constraints,control=list(iter.max=500,eval.max=500),sdlist=sdlist,StartValFixed=XSAMfit$optobj$par)
warnings()
par(mfrow=c(1,2),mar=c(3,3,2,1),pty='s')
yrs<-XSAMdataobj$years
plot(yrs,retrorun[[1]]$stats$ssb,type="l",col=1,ylim=c(0,10000),xlab="",ylab="")
lines(XSAMdataobj$years,XSAMfit$stats$ssb,col=1)
for(i in 1:n)lines(yrs[1:(length(yrs)-i)],retrorun[[i+1]]$stats$ssb)
mtext(side=1,"Year",line=2)
mtext(side=2,"SSB",line=2)

if(0){
#include prediction in F
plot(yrs,retrorun[[1]]$stats$Fbar,type="l",col=1,ylim=c(0,.5),xlab="",ylab="")
for(i in 1:n)lines(yrs[1:(length(yrs)-i)],retrorun[[i+1]]$stats$Fbar)
mtext(side=1,"Year",line=2)
mtext(side=2,expression(bar(F)[5-11]),line=2)
}

#Remove prediction in F...
plot(yrs[1:(length(yrs)-1)],retrorun[[1]]$stats$Fbar[1:(length(yrs)-1)],type="l",col=1,ylim=c(0,.5),xlab="",ylab="")
for(i in 1:n)lines(yrs[1:(length(yrs)-i-1)],retrorun[[i+1]]$stats$Fbar[1:(length(yrs)-i-1)])
mtext(side=1,"Year",line=2)
mtext(side=2,expression(bar(F)[5-11]),line=2)

#---------------------------------------------------------------#
#      Profiling
#---------------------------------------------------------------#
#Profiling over loghvec
hsest<-XSAMfit$tab["loghvec",1]
hs<-seq(0,3,len=15)
hslike<-hs[abs(hs-hsest)==min(abs(hs-hsest))]

ProfileList<-list()

##OOBS: this takes some time.. since the model is fitted 15 times!!	
#Turn on key for returning likelihood components
XSAMdataobj$ReturnLL<-1
tmppar<-CreateParameterList(XSAMdataobj,mapnames="loghvec")
for(i in 1:length(hs)){
	tmppar$ParameterList$FixedPars[tmppar$ParameterNames=="loghvec"]<-hs[i]
	ProfileList[[i]]<-FitXSAMFunction(XSAMdataobj,tmppar,returnObj=T,constraints=constraints,control=list(iter.max=500,eval.max=500))
}

t1<-c()
for(i in 1:length(hs))t1[i]<-ProfileList[[i]]$optobj$objective
plot(hs,t1,type="l",xlab="",ylab="")
points(XSAMfit$tab["loghvec",1],XSAMfit$optobj$objective,col=2,pch=16)
title(expression(l[M]))
mtext("log(h)",side=1,line=2.5,cex=0.7)
mtext("-log likelihood",side=2,line=2,cex=0.7)

#Profile version 2
#Do not optimize likelihood for each value of the fixed parameter, but keeping the other parameters fixed and just calculate likelihood value

#Construct the setup for the parameters... (There must be a simpler way around this...)
randomvector<-XSAMfit$rep$par.random
dimLP<-dim(XSAMpars$ParameterList$LP)
LPmat<-matrix(randomvector,dimLP[1],dimLP[2])

fixedvector<-XSAMfit$rep$par.fixed
names(XSAMpars)
XSAMparsP<-XSAMpars$ParameterList
names(XSAMparsP)
XSAMparsP$FixedPars[XSAMpars$ParameterNames[!is.na(XSAMpars$MapVector)]]<-fixedvector
XSAMparsP$LP<-LPmat

#extract likelihood values
tmp1<-MakeADFun(XSAMdataobj,XSAMparsP,random=c("LP"),DLL="XSAM")
p2<-c(tmp1$fn())

p2<-c()
for(i in 1:length(hs)){
	XSAMparsP$FixedPars[names(XSAMparsP$FixedPars)=="loghvec"]<-hs[i]
	#extract likelihood values
	tmp<-MakeADFun(XSAMdataobj,XSAMparsP,random=c("LP"),DLL="XSAM")
	p2[i]<-c(tmp$fn())
}

lines(hs,p2,col=2)

#More on profiling
par(mfrow=c(1,3),mar=c(3,3,2,1),oma=c(1,1,2,0),pty='s')
t1<-c()
for(i in 1:length(hs))t1[i]<-ProfileList[[i]]$optobj$objective
plot(hs,t1,type="l",xlab="",ylab="")
points(XSAMfit$tab["loghvec",1],XSAMfit$optobj$objective,col=2,pch=16)
title(expression(l[M]))
mtext("log(h)",side=1,line=2.5,cex=0.7)
mtext("-log likelihood",side=2,line=2,cex=0.7)
for(i in 1:length(hs)){
	if(!is.null(ProfileList[[i]]$stats))t1[i]<-ProfileList[[i]]$stats$ssb[length(ProfileList[[i]]$stats$ssb)]
	else t1[i]<-NA

}
plot(hs[!is.na(t1)],t1[!is.na(t1)],type="l",ylim=c(0,10000),xlab="",ylab="")
points(hsest,XSAMfit$stats$ssb[length(XSAMfit$stats$ssb)],pch=16,col=2)
title("SSB")
mtext("log(h)",side=1,line=2.5,cex=0.7)
mtext("Thousand tonnes",side=2,line=2,cex=0.7)

for(i in 1:length(hs)){
	if(!is.null(ProfileList[[i]]$stats))t1[i]<-ProfileList[[i]]$stats$Fbar[length(ProfileList[[i]]$stats$Fbar)]
	else t1[i]<-NA

}
plot(hs[!is.na(t1)],t1[!is.na(t1)],type="l",ylim=c(0,0.3),xlab="",ylab="")
points(hsest,XSAMfit$stats$Fbar[length(XSAMfit$stats$Fbar)],pch=16,col=2)
title(expression(bar(F)[5-11]))
mtext("log(h)",side=1,line=2.5,cex=0.7)

#----------------------------------------------------------------------------------------------------------#
# SAVING RESULTS
#----------------------------------------------------------------------------------------------------------#
#Save fits?

XSAMfitlist<-list(
XSAMfit=XSAMfit,
XSAMdataobj=XSAMdataobj,
XSAMpars=XSAMpars,
retrorun=retrorun,
ProfileList=ProfileList)

savefile<-paste(mypath,"/Results/XSAMfitlist",sep="")
savefile<-"C:/Users/aanes/Documents/XSAMcourseW/"
#save(XSAMfitlist,file=savefile)
if(0){
	load(savefile)
	XSAMfit<-XSAMfitlist$XSAMfit
	XSAMdataobj<-XSAMfitlist$XSAMdataobj
	XSAMpars<-XSAMfitlist$XSAMpars
	retrorun<-XSAMfitlist$retrorun
	ProfileList<-XSAMfitlist$ProfileList
}


#----------------------------------------------------------------------------------------------------------#
# RESULTS
#----------------------------------------------------------------------------------------------------------#
par(mfrow=c(1,2),pty='s')
plot(XSAMdataobj$years,XSAMfit$stats$ssbse/XSAMfit$stats$ssb,type="l",xlab="Years",ylab="RSE SSB",ylim=c(0,.3))
plot(XSAMdataobj$years,XSAMfit$stats$FbarWse/XSAMfit$stats$Fbar,type="l",xlab="Years",ylab="RSE Fbar",ylim=c(0,.3))

#source(paste(path,"XSAMcode/myplot.R",sep=""))

#Key for adding model predicted catch in weight
add<-F
#add<-T
par(mfrow=c(2,2),mar=c(3,3,2,1),oma=c(0,1,0,0))
t1<-barplot(XSAMdataobj$caton,axes=F,ylab="",xlab="",names.arg=F,ylim=c(0,1.2*max(XSAMdataobj$caton)))
tmp<-summary(XSAMfit$rep)
tmp<-tmp[dimnames(tmp)[[1]]%in%"predW",]
n<-length(t1)
t1<-c(t1,t1[n]+t1[n]-t1[n-1])
if(add){
	points(t1,tmp[,1],pch=16,col=2)
	#Adding ~95% confidence ints
	for(i in 1:length(t1))lines(c(t1[i],t1[i]),c(tmp[i,1]-2*tmp[i,2],tmp[i,1]+2*tmp[i,2]),col=2)
}
axis(1,at=t1,labels=1988:2017,tick=T)
axis(2)
box()
mtext(side=2,"Million tonnes",line=2)
title("Landings")

PlotXSAMfitobj(years=XSAMdataobj$years,fitobj=XSAMfit,what="R",p=0.95,add=F,col=1,ylab="",xlab="",lty2=3,axes=F)
axis(1,at=1988:2017)
axis(2)
box()
mtext(side=2,"Millions",line=2)
title("Recruitment (age 2)")

PlotXSAMfitobj(years=XSAMdataobj$years,fitobj=XSAMfit,what="FbarW",p=0.95,add=F,ylim=c(0,.4),col=1,ylab="",xlab="",lty2=3,axes=F)
axis(1,at=1988:2017)
axis(2)
box()
mtext(side=2,expression(bar(F)[w5-11]),line=2)
title("Fishing mortality")
lines(c(-1000,3000),c(0.125,c(0.125)),col="green")
lines(c(-1000,3000),c(0.15,c(0.15)),col="gray")
legend("topright",lty=1,legend=c(expression(F[pa]),expression(F[mp])),col=c("gray","green"),bty="n")

if(0){#weighted F ages 5-12+
FF<-XSAMfit$stats$FF
N<-XSAMfit$stats$N
agid<-4:11
FFB<-colSums(FF[agid,]*N[agid,])/colSums(N[agid,])
lines(1988:2017,FFB,col=2)
}

PlotXSAMfitobj(years=XSAMdataobj$years,fitobj=XSAMfit,what="ssb",p=0.95,add=F,ylim=c(0,10000),col=1,ylab="",xlab="",lty2=3,axes=F)
axis(1,at=1988:2017)
axis(2)
box()
mtext(side=2,"Million tonnes",line=2)
title("SSB")
lines(c(-1000,3000),c(5000,5000),col="green")
lines(c(-1000,3000),c(2500,2500),col="red")
legend("topright",lty=1,legend=c(expression(B[mp]),expression(B[lim])),col=c("green","red"),bty="n")

#----------------------#
#Write reports
#----------------------#
#Gets code
source(paste(path,"XSAMcode/XSAMWriteReports.R",sep=""))
#Extracts N and F point estimates
tables<-Create.Reports(XSAMfit,XSAMdataobj)
View(tables$N)
View(tables$FF)

#Write reports
reportpath<-paste(path,"Results/Reports/",sep="")
write.table(tables$N,paste(reportpath,"N.csv",sep=""),sep=";")
write.table(tables$FF,paste(reportpath,"F.csv",sep=""),sep=";")

#Create summary table
summarytab<-CreateSummaryTable(fitobj=XSAMfit,dataobj=XSAMdataobj,p=0.95)
summarytab<-round(summarytab,3)
View(summarytab)
#and write to file...
write.table(summarytab,paste(reportpath,"summarytab.csv",sep=""),sep=";",dec=",")



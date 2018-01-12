#Exploratory Data Analysis of herring data

#Read data the herring data
mypath<-"C:/Users/aanes/Documents/XSAMcourse/Examples/"

source(paste(mypath,"ReadHerringData.R",sep=""))
#Data now in datalist
names(datalist)
str(datalist)

tmpobj<-datalist$natMor
what<-"M"
plot(rownames(tmpobj),tmpobj[,1],ylim=c(range(tmpobj)),type="l",xlab="Year",ylab=what)
for(i in 1:ncol(tmpobj))lines(rownames(tmpobj),tmpobj[,i])

#Weight at age in catch
tmpobj<-datalist$weca
tmpobj[tmpobj==0]<-NA
what<-"weca"
cols<-rainbow(ncol(tmpobj))
plot(rownames(tmpobj),tmpobj[,1],ylim=c(range(tmpobj,na.rm=T)),xlim=c(1950,2030),type="l",xlab="Year",ylab=what)
for(i in 1:ncol(tmpobj))points(rownames(tmpobj),tmpobj[,i],type="o",pch=16,cex=.3,col=cols[i])
legend("right",lty=1,legend=colnames(tmpobj),title="Age:",col=cols,bty="n")

#Weight at age in stock
tmpobj<-datalist$west
tmpobj[tmpobj==0]<-NA
what<-"west"
cols<-rainbow(ncol(tmpobj))
plot(rownames(tmpobj),tmpobj[,1],ylim=c(range(tmpobj,na.rm=T)),xlim=c(1950,2030),type="l",xlab="Year",ylab=what)
for(i in 1:ncol(tmpobj))points(rownames(tmpobj),tmpobj[,i],type="o",pch=16,cex=.3,col=cols[i])
legend("right",lty=1,legend=colnames(tmpobj),title="Age:",col=cols,bty="n")

#Proportion mature v 1; each age by year
tmpobj<-datalist$matprop
what<-"P(mature)"
cols<-rainbow(ncol(tmpobj))
plot(rownames(tmpobj),tmpobj[,1],ylim=c(range(tmpobj,na.rm=T)),xlim=c(1950,2030),type="l",xlab="Year",ylab=what)
for(i in 1:ncol(tmpobj))points(rownames(tmpobj),tmpobj[,i],type="o",pch=16,cex=.3,col=cols[i])
legend("right",lty=1,legend=colnames(tmpobj),title="Age:",col=cols,bty="n")

#Proportion mature v 2; each year by age for a subset
tmpobj<-datalist$matprop
tmpobj<-tmpobj[as.numeric(rownames(tmpobj))>2005,]
cols<-rainbow(nrow(tmpobj))
plot(colnames(tmpobj),tmpobj[1,],ylim=c(range(tmpobj,na.rm=T)),xlim=c(0,16),type="l",xlab="Age",ylab=what)
for(i in 1:nrow(tmpobj))points(colnames(tmpobj),tmpobj[i,],type="o",pch=16,cex=.3,col=cols[i])
legend("right",lty=1,legend=rownames(tmpobj),title="Year:",col=cols,bty="n")

#Proportion mature v 3; average over years by age
tmpobj<-datalist$matprop
tmpobj<-colMeans(tmpobj)
plot(0:15,tmpobj,ylim=c(range(tmpobj,na.rm=T)),xlim=c(0,16),type="l",xlab="Age",ylab=what)
axis(1,at=0:15)

#Total catches
h<-barplot(t(datalist$caton))
axis(1,at=h,labels=attributes(datalist$caton)$years)
box()

#
n<-datalist$caa
n<-t(n)*1e6
n<-n[paste(1:15),paste(1988:2016)]
yrs<-col(n)+1987
ags<-row(n)
chs<-yrs-ags

uchs<-sort(unique(c(chs)))
nchs<-length(uchs)
par(mfrow=c(5,5),mar=c(2,2,1,1))
for(i in 17:length(uchs)){
	#i<-17
	tmpchs<-uchs[i]
	tmpags<-ags[chs==tmpchs]
	tmpy<-n[chs==tmpchs]
	tmpx<-yrs[chs==tmpchs]
	if(length(tmpx)>2){
		#plot(tmpx,log(tmpy),type="l",pch=16,ylim=c(4,16))
		plot(tmpx,log(tmpy),type="l",pch=16,ylim=c(4,16),xlim=c(min(tmpx),min(tmpx)+15))
		#Add lines corresponding to total mortality Z=-.3
		addZlines(z=-.3)
		text(tmpx,log(tmpy),tmpags)
		mtext(side=1,tmpchs,line=-2)
	}
}
title("log catch by cohort",outer=T)

####
datalist$SurveyIndices[1]
n<-datalist$SurveyIndices[[1]]
n<-t(n)
n<-n[paste(2:15),paste(1988:2017)]
yrs<-col(n)+1987
ags<-row(n)
chs<-yrs-ags


uchs<-sort(unique(c(chs)))
nchs<-length(uchs)
par(mfrow=c(5,5),mar=c(2,2,1,1))
for(i in 17:length(uchs)){
	tmpchs<-uchs[i]
	tmpags<-ags[chs==tmpchs]
	tmpy<-n[chs==tmpchs]
	tmpx<-yrs[chs==tmpchs]
	if(length(tmpx)>2){
		plot(tmpx,log(tmpy),type="l",pch=16,ylim=c(2,10),xlim=c(min(tmpx),min(tmpx)+15))
		text(tmpx,log(tmpy),tmpags)
		mtext(side=1,tmpchs,line=-2)
	}
}

title("log abundance indices by cohort, Fleet 1",outer=T)

################3#
##################
################3#
datalist$SurveyIndices[5]
n<-datalist$SurveyIndices[[5]]
n<-t(n)
n<-n[paste(1:15),paste(1997:2017)]
yrs<-col(n)+1996
ags<-row(n)
chs<-yrs-ags


uchs<-sort(unique(c(chs)))
nchs<-length(uchs)
par(mfrow=c(5,5),mar=c(2,2,1,1))
for(i in 11:length(uchs)){
	tmpchs<-uchs[i]
	tmpags<-ags[chs==tmpchs]
	tmpy<-n[chs==tmpchs]
	tmpx<-yrs[chs==tmpchs]
	if(length(tmpx)>=2){
		plot(tmpx,log(tmpy),type="l",pch=16,ylim=c(3,10),xlim=c(min(tmpx),min(tmpx)+15))
		text(tmpx,log(tmpy),tmpags)
		mtext(side=1,tmpchs,line=-2)
	}
	else{
		plot(tmpx,log(tmpy),pch=16,ylim=c(3,10),xlim=c(min(tmpx),min(tmpx)+15))
		text(tmpx,log(tmpy),tmpags)
		mtext(side=1,tmpchs,line=-2)
	}
}

title("log abundance indices by cohort, Fleet 5",outer=T)


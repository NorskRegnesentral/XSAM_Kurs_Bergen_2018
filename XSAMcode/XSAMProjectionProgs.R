ReadLatentParameters<-function(dataobj,repobj,parobj){
	randomvector<-repobj$par.random
	dimLP<-dim(parobj$ParameterList$LP)
	randommat<-matrix(randomvector,dimLP[1],dimLP[2])
	idmat<-matrix(randomvector,dimLP[1],dimLP[2])
	Am<-dataobj$Am
	noYears<-dataobj$noYears
	#F
	logF<-matrix(NA,nrow=Am,ncol=noYears)
	start<-1
	stopp<-Am
	idmat[start:stopp,]<-"logF"
	logF<-randommat[start:stopp,]
	start<-stopp+1
	#U
	TSsel<-dataobj$TSsel
	uam<-dataobj$uam
	if(TSsel==1){
		stopp<-stopp+uam
		idmat[start:stopp,]<-"U"
		U<-randommat[start:stopp,]
		start<-stopp+1
	}
	else{
		U<-NULL
		
	}
	#V
	LatentEffort<-dataobj$LatentEffort
	if(LatentEffort==1){
		stopp<-stopp+1
		idmat[start:stopp,]<-"V"
		V<-randommat[start:stopp,]
		start<-stopp+1
	}
	else V<-NULL
	#Y
	stopp<-stopp+1
	idmat[start:stopp,]<-"Y"
	Y<-randommat[start:stopp,]
	start<-stopp+1
	#if(LatentEffort==0){V<-Y;Y<-NULL}
	
	#R
	RecruitmentProcess<-dataobj$RecruitmentProcess
	if(RecruitmentProcess==1){
		stopp<-stopp+1
		idmat[start:stopp,]<-"R"
		R<-randommat[start:stopp,]
		start<-stopp+1
	}
	else R<-NULL
	#varM<-dataobj$varM
	#if(varM==1){
	#	stopp<-stopp+(A-1)
	#	idmat[start:stopp,]<-"M"
	#	M<-exp(randommat[start:stopp,])
	#	start<-stopp+1
	#}
	#else M<-NULL
	#Vector of identifiers corresponding to 
	idvec<-c(t(idmat))
	#list(logF=logF,U=U,V=V,Y=Y,R=R,M=M,idmat=idmat)
	list(logF=logF,U=U,V=V,Y=Y,R=R,idmat=idmat)
}

ReadFixedParameters<-function(dataobj,repobj,parobj){
	fixedvector<-repobj$par.fixed
	#print(fixedvector)
	#
	#Fill in mapped parameters in fixedvector
	#Create empty vector with correct length
	tmp<-rep(NA,length(parobj$ParameterNames))
	#Fill in estimated parameters
	tmp[!is.na(parobj$MapVector)]<-fixedvector
	#Fill in fixed parameters from parobj
	tmp[is.na(parobj$MapVector)]<-parobj$ParameterList$FixedPars[is.na(parobj$MapVector)]
	fixedvector<-tmp
	#
	ParCounter<-0
	A<-dataobj$A
	noYears<-dataobj$noYears
	######################
	# START FILLING      #
	######################
	#
	start<-1
	stopp<-A-1
	initlogN<-fixedvector[start:stopp]
	start<-stopp+1
	#
	RecruitmentProcess<-dataobj$RecruitmentProcess
	if(RecruitmentProcess==0){
		stopp<-stopp+noYears
		R<-fixedvector[start:stopp]
		start<-stopp+1
	}
	else R<-NULL
	#
	nqs<-dataobj$nqs
	stopp<-stopp+nqs
	logqpar<-fixedvector[start:stopp]
	start<-stopp+1
	#
	nvf<-dataobj$nvf
	stopp<-stopp+nvf
	s2_1<-exp(fixedvector[start:stopp])
	start<-stopp+1
	#
	TSsel<-dataobj$TSsel
	if(TSsel==1){
		stopp<-stopp+1
		s2_2<-exp(fixedvector[start:stopp])
		start<-stopp+1
	}
	else s2_2<-NULL
	#
	LatentEffort<-dataobj$LatentEffort
	if(LatentEffort==1){
		stopp<-stopp+1
		s2_3<-exp(fixedvector[start:stopp])
		start<-stopp+1
	}
	else s2_3<-NULL
	#
	stopp<-stopp+1
	s2_4<-exp(fixedvector[start:stopp])
	start<-stopp+1
	#
	RecruitmentProcess<-dataobj$RecruitmentProcess
	if(RecruitmentProcess==1){
		stopp<-stopp+1
		sR2<-exp(fixedvector[start:stopp])
		start<-stopp+1
	}
	else sR2<-NULL
	#
	nObsVar<-dataobj$nObsVar
	stopp<-stopp+nObsVar
	hvec<-exp(fixedvector[start:stopp])
	start<-stopp+1
	#
	stopp<-stopp+1
	rec_loga<-fixedvector[start:stopp]
	start<-stopp+1
	#
	stopp<-stopp+1
	rec_logb<-fixedvector[start:stopp]
	start<-stopp+1
	#
	corFlagU<-dataobj$corFlagU
	if(TSsel==1 & corFlagU>0){
		stopp<-stopp+1
		rhoU<-fixedvector[start:stopp]
		start<-stopp+1
	}
	else rhoU<-NULL
	#
	corFlaglogF<-dataobj$corFlaglogF
	if(corFlaglogF>0){
		stopp<-stopp+1
		rhologF<-fixedvector[start:stopp]
		start<-stopp+1
	}
	else rhologF<-NULL
	#
	EstimateARInterceptEffort<-dataobj$EstimateARInterceptEffort
  	if(EstimateARInterceptEffort==1){
		stopp<-stopp+1
		aYpar<-fixedvector[start:stopp]
		start<-stopp+1
	}
	else aYpar<-NULL
	#
	stopp<-stopp+1
	bY<-fixedvector[start:stopp]
	start<-stopp+1
	#
	am<-dataobj$am
	if(TSsel==0){
		stopp<-stopp+am-1
		Upar<-fixedvector[start:stopp]
		start<-stopp+1
	}
	else Upar<-NULL
	#
	EstimateUInterceptEffort<-dataobj$EstimateUInterceptEffort
	uam<-dataobj$uam
	if(TSsel==1 & EstimateUInterceptEffort==1){
		stopp<-stopp+uam
		aUpar<-fixedvector[start:stopp]
		start<-stopp+1
	}
	else aUpar<-NULL
	#
	if(TSsel==1){
		stopp<-stopp+1
		bU<-fixedvector[start:stopp]
		start<-stopp+1
	}
	else bU<-NULL
	#
	stopp<-stopp+2
	#PredictCatchPars<-fixedvector[start:stopp]
	#start<-stopp+1
	#
	#varM<-dataobj$varM
	#if(varM==1){
	#	stopp<-stopp+1
	#	sM2<-exp(fixedvector[start:stopp])
	#	start<-stopp+1
	#}
	#else sM2<-NULL

	###########################################
	list(
	initlogN=initlogN,
	R=R,
	logqpar=logqpar,
	s2_1=s2_1,
	s2_2=s2_2,
	s2_3=s2_3,
	s2_4=s2_4,
	sR2=sR2,
	hvec=hvec,
	rec_loga=rec_loga,
	rec_logb=rec_logb,
	rhoU=rhoU,
	rhologF=rhologF,
	aYpar=aYpar,
	bY=bY,
	Upar=Upar,
	aUpar=aUpar,
	bU=bU)#,
	#PredictCatchPars=PredictCatchPars)#,
	#sM2=sM2)

}

MakeCovar<-function(s2,dim,r,type,key){
	if(length(key)!=dim)stop("Provided wrong dim in MakeCovar")
	cormat<-matrix(NA,dim,dim)
	s2<-s2[key]
	if(!is.matrix(r)){#i.e. r is a singel number
		if(type==0){cormat[]<-0;diag(cormat)<-1}
		if(type==1){cormat[]<-r;diag(cormat)<-1}
		if(type==3){
			cormat<-r^abs(col(cormat)-row(cormat))
		}
	}
	else{
		#Test if the dimensions are correct
		if(sum(r>1)>0)stop("Provided correlation matrix has values greater than 1!")
		if(dim(r)[1]!=dim(r)[2])stop("Provided correlation matrix is not a square matrix!")
		if(dim(r)[1]!=length(s2))stop("Provided correlation matrix has mismatching dimension with s2!")
		cormat<-r
	}
	covmat<-diag(sqrt(s2))%*%cormat%*%diag(sqrt(s2))
	if(sum(covmat)!=0){
		if(sum(round(cov2cor(covmat),10)==round(cormat,10))!=(dim*dim))stop("Something is wrong!!")
	}
	covmat
}

SetFCovariances<-function(sigma2_1,sigma2_2){
	fcor_U<-matrix(NA,nrow=uam,ncol=uam)
	if(corFlagU==0){fcor_U[]<-0;diag(fcor_U)<-1}
	if(corFlagU==1){fcor_U[]<-rhoU;diag(fcor_U)<-1}
	if(corFlagU==2){
		fcor_U<-rhoU^abs(col(fcor_U)-row(fcor_U))
	}
	fvar_U<-sigma2_2*fcor_U
	
	fcor_logF<-matrix(NA,nrow=Am,ncol=Am)
	if(corFlaglogF==0){fcor_logF[]<-0;diag(fcor_logF)<-1}
	if(corFlaglogF==1){fcor_logF[]<-rhologF;diag(fcor_logF)<-1}
	if(corFlaglogF==2){
		fcor_logF<-rhologF^abs(col(fcor_logF)-row(fcor_logF))
	}
	fvar_logF<-diag(sqrt(sigma2_1))%*%fcor_logF%*%diag(sqrt(sigma2_1))
	if(sum(round(cov2cor(fvar_logF),10)==round(fcor_logF,10))!=(Am*Am))stop("Something is wrong!!")
	
	list(fvar_U=fvar_U,fvar_logF=fvar_logF)
}


SimF<-function(inits,Fpars,n,biasCorrect=F){
	require(MASS)
	Upar<-Fpars$Upar
	aU<-Fpars$aUpar#Only for TSsel==1
	bU<-Fpars$bU#Only for TSsel==1
	aY<-Fpars$aY
	bY<-Fpars$bY
	fvar_logF<-Fpars$fvar_logF
	fvar_U<-Fpars$fvar_U
	sigma2_3<-Fpars$sigma2_3
	sigma2_4<-Fpars$sigma2_4
	A<-Fpars$A
	Am<-Fpars$Am
	am<-Fpars$am
	uam<-Fpars$uam
	TSsel<-Fpars$TSsel
	Y<-inits$Y
	U<-matrix(NA,am,n+1)
	muU<-matrix(NA,uam,n+1)
	mulogF<-matrix(NA,Am,n+1)
	V<-c()
	logF<-matrix(NA,A,n+1)
	logF[1:Am,1]<-inits$logF[1:Am]
	if(TSsel==1){
		U[1:uam,1]<-inits$U
	}
	else U[1:uam,1]<-Upar
	#print(U)
	for(y in 1:n){
		#Latent Effort
		Y[y+1]<-aY+bY*Y[y]+rnorm(1,0,sqrt(sigma2_4))
		#Realized effort
		if(is.null(sigma2_3))sigma2_3<-0
		V[y+1]<-Y[y+1]+rnorm(1,0,sqrt(sigma2_3))
		#print(Y);print(V)		
		for(a in 1:uam){#uam<-am-1
			if(TSsel==1){
				muU[a,y+1]<-aU[a]+bU*U[a,y]
			}
			else muU[a,y+1]<-Upar[a]
		}
		#print(muU)
		U[1:uam,y+1]<-mvrnorm(1,muU[,y+1],fvar_U)
		U[am,y+1]<--sum(U[1:uam,y+1])	
		#print(U)
	
		mulogF[1:am,y+1]<-U[1:am,y+1]+V[y+1]
		
		#Set constant above am
		if(am<Am){
			for(a in (am+1):Am){
				mulogF[a,y+1]<-mulogF[am,y+1]
			}
		}
		#print("mulogF")
		#print(mulogF)
		#if(biasCorrect)mulogF[,y+1]<-mulogF[,y+1]-0.5*diag(fvar_logF)
		if(biasCorrect)mulogF[,y+1]<-mulogF[,y+1]-0.5*(diag(fvar_logF)+sigma2_3+sigma2_4)
		logF[1:Am,y+1]<-mvrnorm(1,mulogF[,y+1],fvar_logF)
		#if(biasCorrect)logF[1:Am,y+1]<-logF[1:Am,y+1]-0.5*(diag(fvar_logF)[1:Am]+sigma2_3+sigma2_4)
		#if(biasCorrect)logF[1:Am,y+1]<-logF[1:Am,y+1]-0.5*(diag(fvar_logF)[1:Am])
		for(a in (Am+1):A){
			logF[a,y+1]<-logF[Am,y+1]
		}
	}
	list(Y=Y,V=V,U=U,logF=logF)
}

SetFpars<-function(fixedpars,dataobj){
	Upar<-fixedpars$Upar
	aUpar<-fixedpars$aUpar
	bU<-fixedpars$bU
	aY<-fixedpars$aYpar
	bY<-fixedpars$bY

	#varlogF<-MakeCovar(fixedpars$s2_1,dim=dataobj$Am,fixedpars$rhologF,dataobj$corFlaglogF,dataobj$keyVarF[1,1:dataobj$Am]+1)
	varlogF<-MakeCovar(fixedpars$s2_1,dim=dataobj$Am,fixedpars$rhologF,dataobj$corFlaglogF,dataobj$keyVarF[1:dataobj$Am]+1)
	varU<-MakeCovar(0,dim=dataobj$uam,matrix(0,dataobj$uam,dataobj$uam),dataobj$corFlagU,rep(1,dataobj$uam))

	fvar_logF<-varlogF
	fvar_U<-varU
	sigma2_3<-fixedpars$s2_3
	if(is.null(sigma2_3))sigma2_3<-0
	sigma2_4<-fixedpars$s2_4
	if(is.null(sigma2_4))sigma2_4<-0

	A<-dataobj$A
	Am<-dataobj$Am
	am<-dataobj$am
	uam<-dataobj$uam
	TSsel<-dataobj$TSsel
#	
	list(Upar=Upar,
	aUpar=aUpar,
	bU=bU,
	aY=aY,
	bY=bY,
	fvar_logF=varlogF,
	fvar_U=varU,
	sigma2_3=sigma2_3,
	sigma2_4=sigma2_4,
	A=A,
	Am=Am,
	am=am,
	uam=uam,
	TSsel=TSsel)

	
}

SetFtoEOLD<-function(Ftarget,Fbarrange,Fvector){
	n<-length(Fbarrange)
	Fmult<-(Ftarget/sum(Fvector[Fbarrange]))*n
	Fmult
}

SetFtoE<-function(Ftarget,Fbarrange,Fvector,w){
	n<-length(Fbarrange)
	w<-w/sum(w)
	#Fmult<-(Ftarget/sum(Fvector[Fbarrange]))*n
	Fmult<-(Ftarget/sum(Fvector[Fbarrange]*w))
	Fmult
}


SetFtoC<-function(CC,N,M,cw,Fvector,tol=1e-6,iter.max=50){
	Fpattern<-Fvector/max(Fvector)
	Fmult<-1
	lower_Fmult<-0
	upper_Fmult<-5
	FF<-Fmult*Fpattern
	Z<-FF+M
	newC<-FF/Z*(1-exp(-Z))*N
	newtotC<-sum(newC*cw)
	done<-F
	iter<-0
	while(!done){
		#print(iter)
		#print(Fmult)
		#print(newtotC)
		iter<-iter+1
		if(newtotC>CC){#decrease Fmult
			upper_Fmult<-Fmult
		}
		else{#increase Fmult
			lower_Fmult<-Fmult
		}
		Fmult<-0.5*(lower_Fmult+upper_Fmult)
		FF<-Fmult*Fpattern
		Z<-FF+M
		newC<-FF/Z*(1-exp(-Z))*N
		newtotC<-sum(newC*cw)
		if(iter<iter.max){
			if(abs((newtotC-CC)/CC)<tol)done<-T
		}
		else{
			done<-T
			warning("Maximum number of iteration is reached before convergence!!")
		}
	}
	Fmult
}

SetRpars<-function(fixedpars,dataobj){
	stockRecruitmentModelCode<-dataobj$stockRecruitmentModelCode
	rec_loga<-fixedpars$rec_loga
	rec_logb<-fixedpars$rec_logb
	sR2<-fixedpars$sR2
	list(stockRecruitmentModelCode=stockRecruitmentModelCode,rec_loga=rec_loga,rec_logb=rec_logb,sR2=sR2)
}
GenR<-function(Rpars,biasCorrect=T){
	if(Rpars$stockRecruitmentModelCode==5){
		if(biasCorrect)R<-exp(rnorm(1,Rpars$rec_loga-0.5*Rpars$sR2,sqrt(Rpars$sR2)))
		else R<-exp(rnorm(1,Rpars$rec_loga,sqrt(Rpars$sR2)))
	}
	else stop("No other recruitment models implemented yet!!")
	R
}

ProjectPoPOLD<-function(inits,Rpars,Fpars,Bpars,Fbarrange,n,TACs,Ftargets,DeterministicProjection=F,maxAgePlusGroup=1,ConstraintType=1){
	#ConstraintType=1: Catch constraint, uses TACs
	#ConstraintType=2: F constraint, uses Ftargets and Fbarrange
	#ConstraintType=3: Fs projected according to model
	A<-inits$A
	if(DeterministicProjection){
		Fpars$fvar_logF[]<-0
		Fpars$fvar_U[]<-0
		Fpars$sigma2_3<-0
		Fpars$sigma2_4<-0
		Rpars$sR2<-0
	}
	
	propMat<-Bpars$propMat
	stockMeanWeight<-Bpars$stockMeanWeight
	catchMeanWeight<-Bpars$catchMeanWeight
	M<-Bpars$M
	#NOTE: first element contains initial values
	Ctot<-c()
	Fbar<-c()
	SSB<-c()
	N<-matrix(NA,A,n+2)
	FF<-matrix(NA,A,n+1)
	caa<-matrix(NA,A,n+1)
	N[,1]<-inits$N0
	FF[,1]<-inits$F0
	print(N[,1])
	print(M[,1])
	print(FF[,1])
	SSB[1]<-sum(N[,1]*propMat[,1]*stockMeanWeight[,1])
	caa[,1]<-FF[,1]/(FF[,1]+M[,1])*(1-exp(-(FF[,1]+M[,1])))*N[,1]
	Ctot[1]<-sum(caa[,1]*catchMeanWeight[,1])
	#initsF<-list(Y=inits$Y0,logF=inits$logF0[1:Am])
	initsF<-list(Y=inits$Y0,logF=inits$logF0)
	print("-2")
	for(y in 1:n){
		N[1,y+1]<-GenR(Rpars)
		N[2:A,y+1]<-N[1:(A-1),y]*exp(-FF[1:(A-1),y]-M[1:(A-1),y])
		if(maxAgePlusGroup==1)N[A,y+1]=N[A,y+1]+N[A,y]*exp(-FF[A,y]-M[A,y])
		SSB[y+1]<-sum(N[,y+1]*propMat[,y+1]*stockMeanWeight[,y+1])
		#Simulate F one year ahead to account for changes in fishing pattern...
		print("-1")
		kk0<<-initsF
		kk1<<-Fpars
		firstF<-SimF(initsF,Fpars,n=1)
		Fpattern<-exp(firstF$logF[,2])
		Fpattern<-Fpattern/max(Fpattern)
		#But rescale to catch- or F-constraint
		print("0")
		if(ConstraintType==1){#CatchConstraint
			#Adjust to prescribed total catch in y+1
			Fmult<-SetFtoC(TACs[y+1],N[,y+1],M[,y+1],catchMeanWeight[,y+1],Fpattern)
			k1<<-TACs[y+1]
			k2<<-N[,y+1]
			k3<<-M[,y+1]
			k4<<-catchMeanWeight[,y+1]
			k5<<-Fpattern
			print(Fmult)
		}
		if(ConstraintType==2){#F-constraint
			#Adjust to prescribed Fbar in y+1
			Fmult<-SetFtoE(Ftarget=Ftargets[y+1],Fbarrange=Fbarrange,Fvector=Fpattern)
		}
		if(ConstraintType==3){#Projected by model
			Fpattern<-exp(firstF$logF[,2])
			Fmult<-1
		}
		FF[,y+1]<-Fpattern*Fmult
		print("1")

		#Update inits to scaling F
		#initsF<-list(Y=firstF$Y[2],logF=log(FF[1:Am,y+1]))
		initsF<-list(Y=firstF$Y[2],logF=log(FF[,y+1]))

		caa[,y+1]<-FF[,y+1]/(FF[,y+1]+M[,y+1])*(1-exp(-(FF[,y+1]+M[,y+1])))*N[,y+1]
		Ctot[y+1]<-sum(caa[,y+1]*catchMeanWeight[,y+1])

	}
	#Project stock one time step further
	N[1,n+2]<-GenR(Rpars)
	N[2:A,n+2]<-N[1:(A-1),n+1]*exp(-FF[1:(A-1),n+1]-M[1:(A-1),n+1])
	if(maxAgePlusGroup==1)N[A,n+2]<-N[A,n+2]+N[A,n+1]*exp(-FF[A,n+1]-M[A,n+1])
	SSB[n+2]<-sum(N[,n+2]*propMat[,n+2]*stockMeanWeight[,n+2])
	
	list(N=N,FF=FF,caa=caa,Ctot=Ctot,SSB=SSB)

}

SetBparsOLD<-function(fixedpars,dataobj,n){
	A<-dataobj$A
	minAge<-dataobj$minAge
	maxAge<-dataobj$maxAge
	tmppropMat<-dataobj$propMat[nrow(dataobj$propMat),paste(minAge:maxAge)]
	propMat<-NULL
	for(y in 1:(n+2))propMat<-cbind(propMat,tmppropMat)

	tmpstockMeanWeight<-dataobj$stockMeanWeight[nrow(dataobj$stockMeanWeight),paste(minAge:maxAge)]
	stockMeanWeight<-NULL
	for(y in 1:(n+2))stockMeanWeight<-cbind(stockMeanWeight,tmpstockMeanWeight)

	tmpcatchMeanWeight<-dataobj$catchMeanWeight[nrow(dataobj$catchMeanWeight),paste(minAge:maxAge)]
	catchMeanWeight<-NULL
	for(y in 1:(n+2))catchMeanWeight<-cbind(catchMeanWeight,tmpcatchMeanWeight)

	tmpM<-dataobj$natMor[nrow(dataobj$natMor),paste(minAge:maxAge)]
	M<-NULL
	for(y in 1:(n+2))M<-cbind(M,tmpM)

	list(propMat=propMat,
	stockMeanWeight=stockMeanWeight,
	catchMeanWeight=catchMeanWeight,
	M=M)

}

#Fs=NULL,Bmp=5000,Fmp=0.125,Fmin=0.05,Btrigger=5000,Blim=2500,Bpa=NULL,Fmsy=0.15,Flim=NULL,Fpa=0.15
ProjectPoP<-function(dataobj,repobj,parobj,Bpars,Fbarrange=NULL,n,TACs,Ftargets,DeterministicProjection=F,ConstraintType=1,biasCorrect=T,UseManagementPlan=F,ManagementPlanPars=NULL,weightedF=F,forceC=F){
	#ConstraintType=1: Catch constraint, uses TACs
	#ConstraintType=2: F constraint, uses Ftargets and Fbarrange
	#ConstraintType=3: Fs projected according to model
	A<-dataobj$A
	if(is.null(Fbarrange))Fbarrange<-c(dataobj$Fbarminage,dataobj$Fbarmaxage)
	fixedpars<-ReadFixedParameters(dataobj,repobj=repobj,parobj=parobj)
	randompars<-ReadLatentParameters(dataobj,repobj=repobj,parobj=parobj)
	Fpars<-SetFpars(fixedpars,dataobj)
	Rpars<-SetRpars(fixedpars,dataobj)
	if(DeterministicProjection){
		Fpars$fvar_logF[]<-0
		Fpars$fvar_U[]<-0
		Fpars$sigma2_3<-0
		Fpars$sigma2_4<-0
		Rpars$sR2<-0
	}
	propMat<-Bpars$propMat
	stockMeanWeight<-Bpars$stockMeanWeight
	catchMeanWeight<-Bpars$catchMeanWeight
	M<-Bpars$M
	if(UseManagementPlan){
		Fs<-ManagementPlanPars$Fs
		Bmp<-ManagementPlanPars$Bmp
		Fmp<-ManagementPlanPars$Fmp
		Fmin<-ManagementPlanPars$Fmin
		Btrigger<-ManagementPlanPars$Btrigger
		Blim<-ManagementPlanPars$Blim
		Bpa<-ManagementPlanPars$Bpa
		Fmsy<-ManagementPlanPars$Fmsy
		Flim<-ManagementPlanPars$Flim
		Fpa<-ManagementPlanPars$Fpa
		type<-ManagementPlanPars$type
		Fset<-ManagementPlanPars$Fset
		print(type)
	}
	#############
	N<-matrix(NA,nrow=dataobj$A,ncol=(dataobj$noYears+n+1))
	FF<-matrix(NA,nrow=dataobj$A,ncol=(dataobj$noYears+n))
	caa<-matrix(NA,nrow=dataobj$A,ncol=(dataobj$noYears+n))
	SSB<-Ctot<-c()
	#initialize
	N[1,1:dataobj$noYears]<-exp(randompars$R)
	N[2:dataobj$A,1]<-exp(fixedpars$initlogN)
	FF[1:dataobj$Am,1:dataobj$noYears]<-exp(randompars$logF)
	for(a in (dataobj$Am+1):dataobj$A)FF[a,1:dataobj$noYears]<-FF[dataobj$Am,1:dataobj$noYears]
	caa[,1]<-FF[,1]/(FF[,1]+M[,1])*(1-exp(-(FF[,1]+M[,1])))*N[,1]
	Ctot[1]<-sum(caa[,1]*catchMeanWeight[,1])
	for(y in 1:(dataobj$noYears-1)){
		N[2:dataobj$A,y+1]<-N[1:(dataobj$A-1),y]*exp(-FF[1:(dataobj$A-1),y]-M[1:(dataobj$A-1),y])
		if(dataobj$maxAgePlusGroup==1)N[dataobj$A,y+1]=N[dataobj$A,y+1]+N[dataobj$A,y]*exp(-FF[dataobj$A,y]-M[dataobj$A,y])
		SSB[y]<-sum(N[,y]*propMat[,y]*stockMeanWeight[,y])
		caa[,y+1]<-FF[,y+1]/(FF[,y+1]+M[,y+1])*(1-exp(-(FF[,y+1]+M[,y+1])))*N[,y+1]
		Ctot[y+1]<-sum(caa[,y+1]*catchMeanWeight[,y+1])
	}
	if(forceC){#Force catch in assessmentyear to be equal to catch-prediction and not be given by the models estimate
		Fpattern<-FF[,y+1]
		Fpattern<-Fpattern/max(Fpattern)
		Fmult<-SetFtoC(dataobj$CatchPrediction[1],N[,y+1],M[,y+1],catchMeanWeight[,y+1],Fpattern)
		FF[,y+1]<-Fpattern*Fmult
		caa[,y+1]<-FF[,y+1]/(FF[,y+1]+M[,y+1])*(1-exp(-(FF[,y+1]+M[,y+1])))*N[,y+1]
		Ctot[y+1]<-sum(caa[,y+1]*catchMeanWeight[,y+1])
	}
	SSB[y+1]<-sum(N[,y+1]*propMat[,y+1]*stockMeanWeight[,y+1])
	#############
	if(is.null(Ftargets))Ftargets<-c()
	#NOTE: first element contains initial values
	if(dataobj$TSsel==0)initsF<-list(Y=randompars$Y[dataobj$noYears],logF=randompars$logF[,dataobj$noYears])
	else initsF<-list(Y=randompars$Y[dataobj$noYears],logF=randompars$logF[,dataobj$noYears],U=randompars$U[,dataobj$noYears])
	#If forceF: must rescale such that initial value correspond to catch in assessmentyear
	for(y in dataobj$noYears:(dataobj$noYears+n-1)){
		N[1,y+1]<-GenR(Rpars,biasCorrect)
		N[2:A,y+1]<-N[1:(A-1),y]*exp(-FF[1:(A-1),y]-M[1:(A-1),y])
		if(dataobj$maxAgePlusGroup==1)N[A,y+1]=N[A,y+1]+N[A,y]*exp(-FF[A,y]-M[A,y])
		SSB[y+1]<-sum(N[,y+1]*propMat[,y+1]*stockMeanWeight[,y+1])
		#Simulate F one year ahead to account for changes in fishing pattern...
		firstF<-SimF(initsF,Fpars,n=1,biasCorrect)
		Fpattern<-exp(firstF$logF[,2])
		Fpattern<-Fpattern/max(Fpattern)
		if(UseManagementPlan){
			ConstraintType<-2
			#Ftargets[y+1]<-ManagementPlan(B=SSB[y+1],Fs=Fs,Bmp=Bmp,Fmp=Fmp,Fmin=Fmin,Btrigger=Btrigger,Blim=Blim,Bpa=Bpa,Fmsy=Fmsy,Flim=Flim,Fpa=Fpa)
			Ftargets[y+1]<-ManagementOptions(B=SSB[y+1],Fs=Fs,Bmp=Bmp,Fmp=Fmp,Fmin=Fmin,Btrigger=Btrigger,Blim=Blim,Bpa=Bpa,Fmsy=Fmsy,Flim=Flim,Fpa=Fpa,type=type,Fset=Fset)
			print(Ftargets)
		}
		#But rescale to catch- or F-constraint
		if(ConstraintType==1){#CatchConstraint
			#Adjust to prescribed total catch in y+1
			Fmult<-SetFtoC(TACs[y+1],N[,y+1],M[,y+1],catchMeanWeight[,y+1],Fpattern)
		}
		if(ConstraintType==2){#F-constraint
			#Adjust to prescribed Fbar in y+1
			aindeks<-(dataobj$Fbarminage:dataobj$Fbarmaxage)-dataobj$minAge+1
			print(aindeks)
			#Fmult<-SetFtoE(Ftarget=Ftargets[y+1],Fbarrange=aindeks,Fvector=Fpattern)#old version
			if(weightedF)w<-N[aindeks,y+1]
			else w<-rep(1,length(aindeks))
			Fmult<-SetFtoE(Ftarget=Ftargets[y+1],Fbarrange=aindeks,Fvector=Fpattern,w)
		}
		if(ConstraintType==3){#Projected by model
			Fpattern<-exp(firstF$logF[,2])
			Fmult<-1
		}
		#Fmult<-Fmult*exp(0.5*diag(Fpars$fvar_logF)[3])
		FF[,y+1]<-Fpattern*Fmult

		#Update inits to scaling F
		#initsF<-list(Y=firstF$Y[2],logF=log(FF[1:Am,y+1]))
		if(dataobj$TSsel==0)initsF<-list(Y=firstF$Y[2],logF=log(FF[,y+1]))
		else initsF<-list(Y=firstF$Y[2],logF=log(FF[,y+1]),U=firstF$U[1:dataobj$uam,2])

		caa[,y+1]<-FF[,y+1]/(FF[,y+1]+M[,y+1])*(1-exp(-(FF[,y+1]+M[,y+1])))*N[,y+1]
		Ctot[y+1]<-sum(caa[,y+1]*catchMeanWeight[,y+1])

	}
	#Project stock one time step further
	N[1,dataobj$noYears+n+1]<-GenR(Rpars)
	N[2:A,dataobj$noYears+n+1]<-N[1:(A-1),dataobj$noYears+n]*exp(-FF[1:(A-1),dataobj$noYears+n]-M[1:(A-1),dataobj$noYears+n])
	if(dataobj$maxAgePlusGroup==1)N[A,dataobj$noYears+n+1]<-N[A,dataobj$noYears+n+1]+N[A,dataobj$noYears+n]*exp(-FF[A,dataobj$noYears+n]-M[A,dataobj$noYears+n])
	SSB[dataobj$noYears+n+1]<-sum(N[,dataobj$noYears+n+1]*propMat[,dataobj$noYears+n+1]*stockMeanWeight[,dataobj$noYears+n+1])
	
	aindeks<-(dataobj$Fbarminage:dataobj$Fbarmaxage)-dataobj$minAge+1
	#yindeks<-dataobj$noYears:(dataobj$noYears+n)
	yindeks<-1:(dataobj$noYears+n)
	Fbar<-colMeans(FF[aindeks,])
	FbarW<-colSums(FF[aindeks,yindeks]*N[aindeks,yindeks])/colSums(N[aindeks,yindeks])
	years<-(1:(dataobj$noYears+n+1))+min(dataobj$years)-1
	list(N=N,FF=FF,caa=caa,Ctot=Ctot,SSB=SSB,Fbar=Fbar,FbarW=FbarW,years=years)

}

SetBpars<-function(dataobj,n,type=1){
	A<-dataobj$A
	minAge<-dataobj$minAge
	maxAge<-dataobj$maxAge
	tmppropMat<-dataobj$propMat[nrow(dataobj$propMat),paste(minAge:maxAge)]
	propMat<-t(dataobj$propMat[,paste(minAge:maxAge)])
	for(y in 1:(n+2))propMat<-cbind(propMat,tmppropMat)
	t1<-as.numeric(dimnames(propMat)[[2]])[1]
	t2<-(as.numeric(dimnames(propMat)[[2]])[1]+ncol(propMat)-1)
	nams<-t1:t2
	dimnames(propMat)[[2]]<-nams

	if(type==1)tmpstockMeanWeight<-dataobj$stockMeanWeight[nrow(dataobj$stockMeanWeight),paste(minAge:maxAge)]
	else tmpstockMeanWeight<-colMeans(dataobj$stockMeanWeight[(nrow(dataobj$stockMeanWeight)-2):nrow(dataobj$stockMeanWeight),paste(minAge:maxAge)])
	stockMeanWeight<-t(dataobj$stockMeanWeight[,paste(minAge:maxAge)])
	for(y in 1:(n+2))stockMeanWeight<-cbind(stockMeanWeight,tmpstockMeanWeight)
	t1<-as.numeric(dimnames(stockMeanWeight)[[2]])[1]
	t2<-(as.numeric(dimnames(stockMeanWeight)[[2]])[1]+ncol(stockMeanWeight)-1)
	nams<-t1:t2
	dimnames(stockMeanWeight)[[2]]<-nams


	if(type==1)tmpcatchMeanWeight<-dataobj$catchMeanWeight[nrow(dataobj$catchMeanWeight),paste(minAge:maxAge)]
	else tmpcatchMeanWeight<-colMeans(dataobj$catchMeanWeight[(nrow(dataobj$catchMeanWeight)-2):nrow(dataobj$catchMeanWeight),paste(minAge:maxAge)])
	catchMeanWeight<-t(dataobj$catchMeanWeight[,paste(minAge:maxAge)])
	for(y in 1:(n+2))catchMeanWeight<-cbind(catchMeanWeight,tmpcatchMeanWeight)
	t1<-as.numeric(dimnames(catchMeanWeight)[[2]])[1]
	t2<-(as.numeric(dimnames(catchMeanWeight)[[2]])[1]+ncol(catchMeanWeight)-1)
	nams<-t1:t2
	dimnames(catchMeanWeight)[[2]]<-nams

	tmpM<-dataobj$natMor[nrow(dataobj$natMor),paste(minAge:maxAge)]
	M<-t(dataobj$natMor[,paste(minAge:maxAge)])
	for(y in 1:(n+2))M<-cbind(M,tmpM)
	t1<-as.numeric(dimnames(M)[[2]])[1]
	t2<-(as.numeric(dimnames(M)[[2]])[1]+ncol(M)-1)
	nams<-t1:t2
	dimnames(M)[[2]]<-nams


	list(propMat=propMat,
	stockMeanWeight=stockMeanWeight,
	catchMeanWeight=catchMeanWeight,
	M=M)

}

setInits<-function(randompars,statsobj){
	A<-nrow(statsobj$N)
	Y0<-randompars$Y[length(randompars$Y)]
	logF0<-log(statsobj$FF[,ncol(statsobj$FF)])
	F0<-exp(logF0)
	N0<-statsobj$N[,ncol(statsobj$U)]
	U0<-randompars$Y[,ncol(randompars$U)]
	list(Y0=Y0,F0=F0,logF0=logF0,N0=N0,U0=U0,A=A)
}


#GenerateSimDist<-function(repobj,n=1){
#GenerateSimDist<-function(repobj,dataobj,n=1,biasCorrect=F,idvec){
GenerateSimDist<-function(repobj,dataobj,parobj,n=1,biasCorrect=T){
	#Since parobj only holds the fixed parameters, the latent parameters is extracted from rep 
	tmp<-ReadLatentParameters(dataobj=dataobj,repobj=repobj,parobj=parobj)
	idvec<-c(parobj$ParameterNames[!is.na(parobj$MapVector)],c(t(tmp$idmat)))
	if(is.null(repobj$jointPrecision))stop("Must provide report object using 'getJointPrecision=T'")
	nfixedpar<-length(repobj$par.fixed)
	nrandompar<-length(repobj$par.random)
	covmat<-solve(repobj$jointPrecision)
	#pars<-c(repobj$par.fixed,repobj$par.random)
	#newpars<-mvrnorm(n=n,mu=c(repobj$par.fixed,repobj$par.random),covmat)
	Am<-dataobj$Am
	noYears<-dataobj$noYears
	mu<-c(repobj$par.fixed,repobj$par.random)
	#print(length(idvec))
	#print(length(mu))
	#idvec<-c(rep("fixed",length(repobj$par.fixed)),idvec)
	if(biasCorrect){
		#Biascorrect F and R
		#mu[(nfixedpar+1):(nfixedpar+Am*noYears)]<-mu[(nfixedpar+1):(nfixedpar+Am*noYears)]-0.5*diag(covmat)[(nfixedpar+1):(nfixedpar+Am*noYears)]
		mu[idvec=="logF"]<-mu[idvec=="logF"]-0.5*diag(covmat)[idvec=="logF"]
		mu[idvec=="R"]<-mu[idvec=="R"]-0.5*diag(covmat)[idvec=="R"]
		mu[grep("initlogN",idvec)]<-mu[grep("initlogN",idvec)]-0.5*diag(covmat)[grep("initlogN",idvec)]
		#Biascorect R
		
	}
	newpars<-mvrnorm(n=n,mu=mu,covmat)
	if(!is.matrix(newpars))newpars<-matrix(newpars,nrow=1)
	par.fixed<-newpars[,1:nfixedpar,drop=F]
	par.random<-newpars[,(nfixedpar+1):(nfixedpar+nrandompar),drop=F]
	list(par.fixed=par.fixed,par.random=par.random)
}



StochasticProjection<-function(dataobj,repobj,parobj){
	Fbarrange<-dataobj$Fbarminage:dataobj$Fbarmaxage
	fixedpars<-ReadFixedParameters(dataobj,repobj=repobj,parobj=parobj)
	randompars<-ReadLatentParameters(dataobj,repobj=repobj,parobj=parobj)
	Fpars<-SetFpars(fixedpars,dataobj)
	Rpars<-SetRpars(fixedpars,dataobj)
	inits<-setInits(randompars,statsobj)
	Bpars<-SetBpars(fixedpars,dataobj,n=n)
	GenR(Rpars)


}
sp2<-function(dataobj,fixedpars,randompars){
	N<-matrix(NA,nrow=dataobj$A,ncol=dataobj$noYears)
	N[1,]<-exp(randompars$R)
	N[2:dataobj$A,1]<-exp(fixedpars$initlogN)
	FF<-exp(randompars$logF)
	FF<-rbind(FF,FF[nrow(FF),])
	M<-t(dataobj$natMor[,1:12])
	SSB<-c()
	for(y in 1:(dataobj$noYears-1)){
		N[2:dataobj$A,y+1]<-N[1:(dataobj$A-1),y]*exp(-FF[1:(dataobj$A-1),y]-M[1:(dataobj$A-1),y])
		if(dataobj$maxAgePlusGroup==1)N[dataobj$A,y+1]=N[dataobj$A,y+1]+N[dataobj$A,y]*exp(-FF[dataobj$A,y]-M[dataobj$A,y])
		SSB[y]<-sum(N[,y]*dataobj$propMat[y,1:12]*dataobj$stockMeanWeight[y,1:12])
	}
	SSB[y+1]<-sum(N[,y+1]*dataobj$propMat[y+1,1:12]*dataobj$stockMeanWeight[y+1,1:12])
	list(FF=FF,N=N,SSB=SSB)
}

if(0){
sp2(dataobj,fixedpars,randompars)

n<-1000
tst<-GenerateSimDist(objF51Nrep,n=n)
tmpobj<-objF51Nrep
tstlist<-list()
parobj<-pars
for(i in 1:nrow(tst[[1]])){
	tmpobj$par.fixed<-tst$par.fixed[i,]
	tmpobj$par.random<-tst$par.random[i,]
	fixedpars<-ReadFixedParameters(dataobj,repobj=tmpobj,parobj=parobj)
	randompars<-ReadLatentParameters(dataobj,repobj=tmpobj,parobj=parobj)
	tstlist[[i]]<-sp2(dataobj,fixedpars,randompars)
}

tmp<-lapply(tstlist,FUN=function(x)x$SSB)
tmp<-matrix(unlist(tmp),nrow=n,byrow=T)

colMeans(tmp)
sqrt(colVars(tmp))

cbind(statsobj$ssb,colMeans(tmp))
cbind(statsobj$ssbse,sqrt(colVars(tmp)))
}


ManagementPlan<-function(B,Fs=NULL,Bmp=5000,Fmp=0.125,Fmin=0.05,Btrigger=5000,Blim=2500,Bpa=NULL,Fmsy=0.15,Flim=NULL,Fpa=0.15){
	#B is the spawning stock biomass in the beginning of the quota year
	#Management plan
	if(B>Blim){
		if(B<Bmp){
			b<-(Fmp-Fmin)/(Bmp-Blim)
			a<-Fmp-b*Bmp
			Freal<-a+b*B
		}
		else Freal<-Fmp
	}
	else Freal<-Fmin
	#MSY
	#if(B<Btrigger)Freal<-Fmsy*(B/Btrigger)
	#else Freal<-Fmsy
	Freal
}

ManagementOptions<-function(B,Fs=NULL,Bmp=5000,Fmp=0.125,Fmin=0.05,Btrigger=5000,Blim=2500,Bpa=5000,Fmsy=0.15,Flim=NULL,Fpa=0.15,type="MP",Fset=0){
	if(!type%in%c("MP","MSY","setF"))stop("Type must be one of MP, MSY or setF!!")
	#B is the spawning stock biomass in the beginning of the quota year
	#Management plan
	if(type=="MP"){
		if(B>Blim){
			if(B<Bmp){
				b<-(Fmp-Fmin)/(Bmp-Blim)
				a<-Fmp-b*Bmp
				Freal<-a+b*B
			}
			else Freal<-Fmp
		}
		else Freal<-Fmin
	}
	if(type=="MSY"){
		if(B<Btrigger)Freal<-Fmsy*(B/Btrigger)
		else Freal<-Fmsy
	}
	if(type=="setF"){
		Freal<-Fset
	}
	#MSY
	#if(B<Btrigger)Freal<-Fmsy*(B/Btrigger)
	#else Freal<-Fmsy
	Freal
}


DoXSAMForecast<-function(dataobj,fitobj,parobj,StochasticForcast=T,Bpars=NULL,Fbarrange=NULL,np=2,n=1000,ManagementPlanPars,BparsType=1,weightedF=F,forceC=F){
	if(is.null(Bpars))Bpars<-SetBpars(dataobj,n=np,type=BparsType)
	#1) Deterministic projection using management plan
	FC0<-ProjectPoP(dataobj,fitobj$rep,parobj,Bpars,Fbarrange=Fbarrange,n=np,TACs=NULL,Ftargets=NULL,DeterministicProjection=T,ConstraintType=1,UseManagementPlan=T,ManagementPlanPars=ManagementPlanPars,weightedF=weightedF,forceC=forceC)
	#2) Deterministic projection using TAC derived from management plan
	FC1<-ProjectPoP(dataobj,fitobj$rep,parobj,Bpars,Fbarrange=Fbarrange,n=np,TACs=FC0$Ctot,Ftargets=NULL,DeterministicProjection=T,ConstraintType=1,UseManagementPlan=F,ManagementPlanPars=ManagementPlanPars,weightedF=weightedF,forceC=forceC)

	#3) Stochastic projection using TAC derived from MP
	TACs<-FC1$Ctot
	print(TACs)
	if(StochasticForcast){
		#Generate samples from simultaneous distribution of all fitted parameters including latent parameters
		parrepl<-GenerateSimDist(fitobj$rep,dataobj,parobj,n=n)
		tmprepobj<-list()
		FCrepl<-list()
		for(i in 1:nrow(parrepl[[1]])){
			tmprepobj$par.fixed<-parrepl$par.fixed[i,]
			tmprepobj$par.random<-parrepl$par.random[i,]
			#FCrepl[[i]]<-ProjectPoP(dataobj,tmprepobj,parobj=parobj,Bpars,Fbarrange=NULL,n=np,TACs=TACs,Ftargets=NULL,DeterministicProjection=F,ConstraintType=1,weightedF=weightedF,forceC=forceC)
			FCrepl[[i]]<-ProjectPoP(dataobj,tmprepobj,parobj=parobj,Bpars,Fbarrange=NULL,n=np,TACs=TACs,Ftargets=NULL,DeterministicProjection=F,ConstraintType=1,UseManagementPlan=F,ManagementPlanPars=ManagementPlanPars,weightedF=weightedF,forceC=forceC)
		}
	}
	else{
		FCrepl<-NULL
	}
	years<-FC0$years[1:length(FC0$years)]
	yearsF<-years[1:(length(years)-1)]
	list(FC0=FC0,FC1=FC1,FCrepl=FCrepl,TACs=TACs,years=years,yearsF=yearsF,type=ManagementPlanPars$type,Bpars=Bpars)
}

if(0){
PlotXSAMfitobj<-function(years,fitobj,what="ssb",p=0.95,add=F,ylim=NULL,xlim=NULL,col=1,lwd=1,plot=T){
	if(!what%in%c("ssb","R","Fbar","FbarW"))stop("what must be one of ssb,R,Fbar or FbarW!")
	if(what=="ssb"){
		y<-fitobj$stats$ssb
		yse<-fitobj$stats$ssbse
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
		if(!add)plot(years,y,ylim=ylim,xlim=xlim,type="l",col=col,lwd=lwd)
		else lines(years,y,col=col,lwd=lwd)
		lines(years,yl,col=col,lwd=lwd,lty=2)
		lines(years,yu,col=col,lwd=lwd,lty=2)
	}
	invisible(list(y=y,yl=yl,yu=yu))
}
}
XSAMForecastStats<-function(FC,fitobj,dataobj,what="ssb",p=0.95,StatsType=1){
	#if StatsType==1: Uses estimates and uncertainties directly from TMB till last year in assessment such that derived 
	#stats are based on replicates only for prediction.
	#if StatsType==2:derived stats are based on replicates only for the entire time series. MAY therefore deviate from TMB 
	#estimates due to Monte Carlo error
	#if StatsType==3: Uses estimates and uncertainties directly from TMB till last year in assessment-1 such that derived 
	#stats are based on replicates only for prediction.
	years<-dataobj$years
	y<-PlotXSAMfitobj(years,fitobj,what=what,p=0.95,plot=F)
	ym<-y$y
	yl<-y$yl
	yu<-y$yu
	if(what=="ssb"){
		tmp<-lapply(FC$FCrepl,FUN=function(x)x$SSB)
		ydet0<-FC$FC0$SSB
		ydet1<-FC$FC1$SSB
	}
	if(what=="R"){
		tmp<-lapply(FC$FCrepl,FUN=function(x)x$N[1,])
		ydet0<-FC$FC0$N[1,]
		ydet1<-FC$FC1$N[1,]
	}
	if(what=="Fbar"){
		tmp<-lapply(FC$FCrepl,FUN=function(x)x$Fbar)
		ydet0<-FC$FC0$Fbar
		ydet1<-FC$FC1$Fbar
	}
	if(what=="FbarW"){
		tmp<-lapply(FC$FCrepl,FUN=function(x)x$FbarW)
		ydet0<-FC$FC0$FbarW
		ydet1<-FC$FC1$FbarW
	}
	yearssim<-FC$years
	yearsFsim<-FC$yearsF
	tmp<-matrix(unlist(tmp),nrow=n,byrow=T)	
	ysimm<-colMeans(tmp)
	pl<-(1-p)/2
	pu<-1-pl
	ysiml<-apply(tmp,MAR=2,FUN=quantile,probs=pl)
	ysimu<-apply(tmp,MAR=2,FUN=quantile,probs=pu)
	##
	ymean<-ysimm
	ylower<-ysiml
	yupper<-ysimu
	if(StatsType==1 & what%in%c("ssb","R")){
		ymean[yearssim<=max(years)]<-ym
		ylower[yearssim<=max(years)]<-yl
		yupper[yearssim<=max(years)]<-yu
		years<-yearssim
	}
	if(StatsType==1 & what%in%c("Fbar","FbarW")){
		ymean[yearsFsim<=max(years)]<-ym
		ylower[yearsFsim<=max(years)]<-yl
		yupper[yearsFsim<=max(years)]<-yu
		years<-yearsFsim
	}
	if(StatsType==3 & what%in%c("ssb","R")){
		ymean[yearssim<=max(years)-1]<-ym[1:(length(ym)-1)]
		ylower[yearssim<=max(years)-1]<-yl[1:(length(ym)-1)]
		yupper[yearssim<=max(years)-1]<-yu[1:(length(ym)-1)]
		years<-yearssim
	}
	if(StatsType==3 & what%in%c("Fbar","FbarW")){
		ymean[yearsFsim<=max(years)-1]<-ym[1:(length(ym)-1)]
		ylower[yearsFsim<=max(years)-1]<-yl[1:(length(ym)-1)]
		yupper[yearsFsim<=max(years)-1]<-yu[1:(length(ym)-1)]
		years<-yearsFsim
	}
	#list(years=years,ym=ym,yl=yl,yu=yu,ydet0=ydet0,ydet1=ydet1,yearssim=yearssim,yearsFsim=yearsFsim,ysimm=ysimm,ysiml=ysiml,ysimu=ysimu)
	list(years=years,ymean=ymean,ylower=ylower,yupper=yupper,ydet0=ydet0,ydet1=ydet1)
}

PlotXSAMForecastStats<-function(Forecastobj,add=F,ylim=NULL,xlim=NULL,col=1,lwd=1,xlab="",ylab=""){
	years<-Forecastobj$years
	y<-Forecastobj$ymean
	yl<-Forecastobj$ylower
	yu<-Forecastobj$yupper
	if(!add & is.null(ylim))ylim<-c(0,max(yu))
	if(!add & is.null(xlim))xlim<-range(years)+c(-1,1)
	if(!add)plot(years,y,ylim=ylim,xlim=xlim,type="l",col=col,lwd=lwd,xlab=xlab,ylab=ylab)
	else lines(years,y,col=col,lwd=lwd)
	lines(years,yl,col=col,lwd=lwd,lty=2)
	lines(years,yu,col=col,lwd=lwd,lty=2)
}


ManagementOptionsTable<-function(fobj,dataobj,p=0.95,weightedF=F){
	print("basis")
	basis<-list()
	#F in assessment year
	basis$assessmentyear<-max(dataobj$years)
	if(weightedF)basis$F<-fobj$FC1$FbarW[fobj$FC1$years==max(dataobj$years)]
	else basis$F<-fobj$FC1$Fbar[fobj$FC1$years==max(dataobj$years)]
	#SSB January 1. in quota year
	basis$SSB<-fobj$FC1$SSB[fobj$FC1$years==(max(dataobj$years)+1)]
	yrs<-max(dataobj$years):(max(dataobj$years)+2)
	basis$Ryrs<-yrs
	basis$R<-fobj$FC1$N[1,fobj$FC1$years%in%yrs]
	basis$TACest<-fobj$FC1$Ctot[fobj$FC1$years==max(dataobj$years)]
	basis$TAC<-dataobj$CatchPrediction[1]
	print(basis)
	optionstab<-list()
	print("options")
	optionstab$type<-fobj$type
	optionstab$TAC<-fobj$FC1$Ctot[fobj$FC1$years==(max(dataobj$years)+1)]
	id<-(1:length(fobj$years))[fobj$FC1$years==(max(dataobj$years)+1)]
	if(weightedF) Frepl<-unlist(lapply(fobj$FCrepl,FUN=function(x,id)x$FbarW[id],id=id))
	else Frepl<-unlist(lapply(fobj$FCrepl,FUN=function(x,id)x$Fbar[id],id=id))
	koko<<-Frepl
	optionstab$F<-fobj$FC1$FbarW[fobj$FC1$years==max(dataobj$years)+1]#mean(Frepl)
	pl<-(1-p)/2
	pu<-1-pl
	optionstab$FCI<-quantile(Frepl,probs=c(pl,pu))
	optionstab$PFabF<-sum(Frepl>optionstab$F)/length(Frepl)
	optionstab$PFabFpm<-sum(Frepl>0.125)/length(Frepl)

	id<-(1:length(fobj$years))[fobj$FC1$years==(max(dataobj$years)+2)]
	Brepl<-unlist(lapply(fobj$FCrepl,FUN=function(x,id)x$SSB[id],id=id))
	koko<<-Brepl
	optionstab$SSB<-fobj$FC1$SSB[fobj$FC1$years==(max(dataobj$years)+2)]#mean(Brepl)
	pl<-(1-p)/2
	pu<-1-pl
	optionstab$SSBCI<-quantile(Brepl,probs=c(pl,pu))
	optionstab$PSSBblBpa<-sum(Brepl<5000)/length(Brepl)
	optionstab$PSSBblBlim<-sum(Brepl<2500)/length(Brepl)
	Breplyminus1<-unlist(lapply(fobj$FCrepl,FUN=function(x,id)x$SSB[id],id=id-1))
	Bdecrease<-(Brepl-Breplyminus1)/Breplyminus1
	optionstab$Bdecrease<-round(mean(Bdecrease)*100)
	optionstab$BdecreaseCI<-round(quantile(Bdecrease,probs=c(pl,pu))*100)
	optionstab$TACchange<-round((optionstab$TAC-basis$TAC)/basis$TAC*100)

	print(optionstab)
	list(basis=basis,optionstab=optionstab)
}

printManagementOptionsTable<-function(obj,ndig=3,printBasis=T,file="",append=F){
	if(file!="")sink(file=file,append=append)
	basis<-obj$basis
	optionstab<-obj$optionstab
	if(printBasis){
		cat("Basis:\n")
		cat(paste("F(",basis$assessmentyear,"): ",sep=""))
		cat(round(basis$F,ndig),"\n")
		cat(paste("SSB(",basis$assessmentyear+1,"): ",sep=""))
		cat(round(basis$SSB,ndig),"\n")
		cat(paste("Recruitment(",paste(range(basis$Ryrs),collapse="-"),"): ",sep=""))
		cat(paste(round(basis$R/1000,0),collapse=" "),"\n")
		cat(paste("Catch(",basis$assessmentyear,"): ",sep=""))
		cat(round(basis$TAC,ndig),"\n")
		cat("\n")
	}
	cat(paste("Rationale: ",optionstab$type,"\n",sep=""))
	cat(paste("Catch(",basis$assessmentyear+1,"): ",sep=""))
	cat(round(optionstab$TAC,ndig),"\n")
	cat(paste("F(",basis$assessmentyear+1,"): ",sep=""))
	cat(round(optionstab$F,ndig))
	cat(paste(" 95%CI: (",paste(round(optionstab$FCI,ndig),collapse=","),")",sep=""))
	cat("\n")
	cat(paste("SSB(",basis$assessmentyear+2,"): ",sep=""))
	cat(round(optionstab$SSB,ndig))
	cat(paste(" 95%CI: (",paste(round(optionstab$SSBCI,ndig),collapse=","),")",sep=""))
	#
	cat(",  P(SSB<Bpa)=")
	cat(round(optionstab$PSSBblBpa,ndig))
	cat(",  P(SSB<Blim)=")
	cat(round(optionstab$PSSBblBlim,ndig))
	cat("\n")
	#
	cat("%SSB Change: ")
	cat(round(optionstab$Bdecrease,ndig))
	cat(paste(" 95%CI: (",paste(round(optionstab$BdecreaseCI,ndig),collapse=","),")",sep=""))
	cat("\n")	
	cat("%TAC Change: ")
	cat(paste(optionstab$TACchange))
	cat("\n")	
	cat("\n")
	cat("\n")
	if(file!="")sink(file=NULL)
}

#############################################################################################################################################
#-The following is for creating a sample from the simultaneous distribution to obtain a sample of residuals that is not serially correlated.
# Analouge to onse step prediciton errors using Kalman filter
#############################################################################################################################################
SimPoP<-function(dataobj,repobj,parobj,Bpars,Fbarrange=NULL,n,TACs,Ftargets,DeterministicProjection=F,ConstraintType=1,biasCorrect=T,UseManagementPlan=F,ManagementPlanPars=NULL,weightedF=F){
	#ConstraintType=1: Catch constraint, uses TACs
	#ConstraintType=2: F constraint, uses Ftargets and Fbarrange
	#ConstraintType=3: Fs projected according to model
	A<-dataobj$A
	if(is.null(Fbarrange))Fbarrange<-c(dataobj$Fbarminage,dataobj$Fbarmaxage)
	fixedpars<-ReadFixedParameters(dataobj,repobj=repobj,parobj=parobj)
	k0<<-fixedpars
	randompars<-ReadLatentParameters(dataobj,repobj=repobj,parobj=parobj)
	Fpars<-SetFpars(fixedpars,dataobj)
	Rpars<-SetRpars(fixedpars,dataobj)
	if(DeterministicProjection){
		Fpars$fvar_logF[]<-0
		Fpars$fvar_U[]<-0
		Fpars$sigma2_3<-0
		Fpars$sigma2_4<-0
		Rpars$sR2<-0
	}
	propMat<-Bpars$propMat
	stockMeanWeight<-Bpars$stockMeanWeight
	catchMeanWeight<-Bpars$catchMeanWeight
	M<-Bpars$M
	if(UseManagementPlan){
		Fs<-ManagementPlanPars$Fs
		Bmp<-ManagementPlanPars$Bmp
		Fmp<-ManagementPlanPars$Fmp
		Fmin<-ManagementPlanPars$Fmin
		Btrigger<-ManagementPlanPars$Btrigger
		Blim<-ManagementPlanPars$Blim
		Bpa<-ManagementPlanPars$Bpa
		Fmsy<-ManagementPlanPars$Fmsy
		Flim<-ManagementPlanPars$Flim
		Fpa<-ManagementPlanPars$Fpa
		type<-ManagementPlanPars$type
		Fset<-ManagementPlanPars$Fset
		print(type)
	}
	#############
	N<-matrix(NA,nrow=dataobj$A,ncol=(dataobj$noYears+n+1))
	FF<-matrix(NA,nrow=dataobj$A,ncol=(dataobj$noYears+n))
	caa<-matrix(NA,nrow=dataobj$A,ncol=(dataobj$noYears+n))
	SSB<-Ctot<-c()
	#initialize
	N[1,1:dataobj$noYears]<-exp(randompars$R)
	N[2:dataobj$A,1]<-exp(fixedpars$initlogN)
	FF[1:dataobj$Am,1:dataobj$noYears]<-exp(randompars$logF)
	for(a in (dataobj$Am+1):dataobj$A)FF[a,1:dataobj$noYears]<-FF[dataobj$Am,1:dataobj$noYears]
	caa[,1]<-FF[,1]/(FF[,1]+M[,1])*(1-exp(-(FF[,1]+M[,1])))*N[,1]
	Ctot[1]<-sum(caa[,1]*catchMeanWeight[,1])
	for(y in 1:(dataobj$noYears-1)){
		N[2:dataobj$A,y+1]<-N[1:(dataobj$A-1),y]*exp(-FF[1:(dataobj$A-1),y]-M[1:(dataobj$A-1),y])
		if(dataobj$maxAgePlusGroup==1)N[dataobj$A,y+1]=N[dataobj$A,y+1]+N[dataobj$A,y]*exp(-FF[dataobj$A,y]-M[dataobj$A,y])
		SSB[y]<-sum(N[,y]*propMat[,y]*stockMeanWeight[,y])
		caa[,y+1]<-FF[,y+1]/(FF[,y+1]+M[,y+1])*(1-exp(-(FF[,y+1]+M[,y+1])))*N[,y+1]
		Ctot[y+1]<-sum(caa[,y+1]*catchMeanWeight[,y+1])
	}
	SSB[y+1]<-sum(N[,y+1]*propMat[,y+1]*stockMeanWeight[,y+1])
	#############
	if(is.null(Ftargets))Ftargets<-c()
	#NOTE: first element contains initial values
	if(dataobj$TSsel==0)initsF<-list(Y=randompars$Y[dataobj$noYears],logF=randompars$logF[,dataobj$noYears])
	else initsF<-list(Y=randompars$Y[dataobj$noYears],logF=randompars$logF[,dataobj$noYears],U=randompars$U[,dataobj$noYears])
	for(y in dataobj$noYears:(dataobj$noYears+n-1)){
		N[1,y+1]<-GenR(Rpars,biasCorrect)
		N[2:A,y+1]<-N[1:(A-1),y]*exp(-FF[1:(A-1),y]-M[1:(A-1),y])
		if(dataobj$maxAgePlusGroup==1)N[A,y+1]=N[A,y+1]+N[A,y]*exp(-FF[A,y]-M[A,y])
		SSB[y+1]<-sum(N[,y+1]*propMat[,y+1]*stockMeanWeight[,y+1])
		#Simulate F one year ahead to account for changes in fishing pattern...
		firstF<-SimF(initsF,Fpars,n=1,biasCorrect)
		Fpattern<-exp(firstF$logF[,2])
		Fpattern<-Fpattern/max(Fpattern)
		if(UseManagementPlan){
			ConstraintType<-2
			#Ftargets[y+1]<-ManagementPlan(B=SSB[y+1],Fs=Fs,Bmp=Bmp,Fmp=Fmp,Fmin=Fmin,Btrigger=Btrigger,Blim=Blim,Bpa=Bpa,Fmsy=Fmsy,Flim=Flim,Fpa=Fpa)
			Ftargets[y+1]<-ManagementOptions(B=SSB[y+1],Fs=Fs,Bmp=Bmp,Fmp=Fmp,Fmin=Fmin,Btrigger=Btrigger,Blim=Blim,Bpa=Bpa,Fmsy=Fmsy,Flim=Flim,Fpa=Fpa,type=type,Fset=Fset)
			print(Ftargets)
		}
		#But rescale to catch- or F-constraint
		if(ConstraintType==1){#CatchConstraint
			#Adjust to prescribed total catch in y+1
			Fmult<-SetFtoC(TACs[y+1],N[,y+1],M[,y+1],catchMeanWeight[,y+1],Fpattern)
		}
		if(ConstraintType==2){#F-constraint
			#Adjust to prescribed Fbar in y+1
			aindeks<-(dataobj$Fbarminage:dataobj$Fbarmaxage)-dataobj$minAge+1
			print(aindeks)
			#Fmult<-SetFtoE(Ftarget=Ftargets[y+1],Fbarrange=aindeks,Fvector=Fpattern)#old version
			if(weightedF)w<-N[aindeks,y+1]
			else w<-rep(1,length(aindeks))
			Fmult<-SetFtoE(Ftarget=Ftargets[y+1],Fbarrange=aindeks,Fvector=Fpattern,w)
		}
		if(ConstraintType==3){#Projected by model
			Fpattern<-exp(firstF$logF[,2])
			Fmult<-1
		}
		#Fmult<-Fmult*exp(0.5*diag(Fpars$fvar_logF)[3])
		FF[,y+1]<-Fpattern*Fmult

		#Update inits to scaling F
		#initsF<-list(Y=firstF$Y[2],logF=log(FF[1:Am,y+1]))
		if(dataobj$TSsel==0)initsF<-list(Y=firstF$Y[2],logF=log(FF[,y+1]))
		else initsF<-list(Y=firstF$Y[2],logF=log(FF[,y+1]),U=firstF$U[1:dataobj$uam,2])

		caa[,y+1]<-FF[,y+1]/(FF[,y+1]+M[,y+1])*(1-exp(-(FF[,y+1]+M[,y+1])))*N[,y+1]
		Ctot[y+1]<-sum(caa[,y+1]*catchMeanWeight[,y+1])

	}
	#Project stock one time step further
	N[1,dataobj$noYears+n+1]<-GenR(Rpars)
	N[2:A,dataobj$noYears+n+1]<-N[1:(A-1),dataobj$noYears+n]*exp(-FF[1:(A-1),dataobj$noYears+n]-M[1:(A-1),dataobj$noYears+n])
	if(dataobj$maxAgePlusGroup==1)N[A,dataobj$noYears+n+1]<-N[A,dataobj$noYears+n+1]+N[A,dataobj$noYears+n]*exp(-FF[A,dataobj$noYears+n]-M[A,dataobj$noYears+n])
	SSB[dataobj$noYears+n+1]<-sum(N[,dataobj$noYears+n+1]*propMat[,dataobj$noYears+n+1]*stockMeanWeight[,dataobj$noYears+n+1])
	
	aindeks<-(dataobj$Fbarminage:dataobj$Fbarmaxage)-dataobj$minAge+1
	#yindeks<-dataobj$noYears:(dataobj$noYears+n)
	yindeks<-1:(dataobj$noYears+n)
	Fbar<-colMeans(FF[aindeks,])
	FbarW<-colSums(FF[aindeks,yindeks]*N[aindeks,yindeks])/colSums(N[aindeks,yindeks])
	years<-(1:(dataobj$noYears+n+1))+min(dataobj$years)-1
	#Generate surveyobservation
	surveyobslist<-list()
	keyLogqpar<-dataobj$keyLogqpar
	print(keyLogqpar)
	nsurveys<-nrow(keyLogqpar)
	for(i in 1:nsurveys){
		logqpar<-fixedpars$logqpar[keyLogqpar[i,]+1]
		ags<-dataobj$agerangeI[i,1]:dataobj$agerangeI[i,2]
		agsidx<-(ags-dataobj$minAge+1)
		print(i)
		print(logqpar)
		print(ags)
		print(agsidx)
		tmpmat<-matrix(NA,nrow=nrow(N),ncol=ncol(N))
		for(y in 1:(ncol(N)-1))tmpmat[agsidx,y]<-(exp(logqpar)*N[agsidx,y]*exp(-dataobj$sampleTime[i+1]*(FF[agsidx,y]+M[agsidx,y])))
		surveyobslist[[i]]<-tmpmat
	}
	
	list(N=N,FF=FF,caa=caa,Ctot=Ctot,SSB=SSB,Fbar=Fbar,FbarW=FbarW,years=years,surveyobslist=surveyobslist)

}

	
DoXSAMSim<-function(dataobj,fitobj,parobj,StochasticForcast=T,Bpars=NULL,Fbarrange=NULL,np=2,n=1000,ManagementPlanPars,BparsType=1,weightedF=F){
	if(is.null(Bpars))Bpars<-SetBpars(dataobj,n=np,type=BparsType)
	#1) Deterministic projection using management plan
	FC0<-ProjectPoP(dataobj,fitobj$rep,parobj,Bpars,Fbarrange=Fbarrange,n=np,TACs=NULL,Ftargets=NULL,DeterministicProjection=T,ConstraintType=1,UseManagementPlan=T,ManagementPlanPars=ManagementPlanPars,weightedF=weightedF)
	#2) Deterministic projection using TAC derived from management plan
	FC1<-ProjectPoP(dataobj,fitobj$rep,parobj,Bpars,Fbarrange=Fbarrange,n=np,TACs=FC0$Ctot,Ftargets=NULL,DeterministicProjection=T,ConstraintType=1,UseManagementPlan=F,ManagementPlanPars=ManagementPlanPars,weightedF=weightedF)

	#3) Stochastic projection using TAC derived from MP
	TACs<-FC1$Ctot
	if(StochasticForcast){
		#Generate samples from simultaneous distribution of all fitted parameters including latent parameters
		parrepl<-GenerateSimDist(fitobj$rep,dataobj,parobj,n=n)
		tmprepobj<-list()
		FCrepl<-list()
		for(i in 1:nrow(parrepl[[1]])){
			tmprepobj$par.fixed<-parrepl$par.fixed[i,]
			tmprepobj$par.random<-parrepl$par.random[i,]
			#FCrepl[[i]]<-ProjectPoP(dataobj,tmprepobj,parobj=parobj,Bpars,Fbarrange=NULL,n=np,TACs=TACs,Ftargets=NULL,DeterministicProjection=F,ConstraintType=1,weightedF=weightedF)
			FCrepl[[i]]<-SimPoP(dataobj,tmprepobj,parobj=parobj,Bpars,Fbarrange=NULL,n=np,TACs=TACs,Ftargets=NULL,DeterministicProjection=F,UseManagementPlan=F,ConstraintType=3,weightedF=weightedF)
		}
	}
	else{
		FCrepl<-NULL
	}
	years<-FC0$years[1:length(FC0$years)]
	yearsF<-years[1:(length(years)-1)]
	list(FC0=FC0,FC1=FC1,FCrepl=FCrepl,TACs=TACs,years=years,yearsF=yearsF,type=ManagementPlanPars$type,Bpars=Bpars)
}


OneSampleApproach<-function(dataobj,fitobj,parobj){
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

	n<-2
	np<-1
	#Obtain one realisation of simultaneous distribution of the population trajectory. Uses the stochastic projection as used for predictions...
	XSAMsim<-DoXSAMSim(dataobj=dataobj,fitobj=fitobj,parobj=parobj,StochasticForcast=T,Bpars=NULL,Fbarrange=NULL,np=np,n=n,ManagementPlanPars,BparsType=2,weightedF=T)

	k<-1
	
	x<-log(t(dataobj$caa))
	x<-cbind(x,NA)
	cres<-x-log(XSAMsim$FCrepl[[k]]$caa[1:nrow(x),1:ncol(x)])
	#vmat<-XSAMf145TSU$res$varmatC
	vmat<-fitobj$res$varmatC
	cstdres<-cres/sqrt(vmat)
	
	#Combine with data
	ires<-istdres<-list()
	for(j in 1:dataobj$nIndices){
		x<-matrix(NA,nrow=length(dataobj$minAge:dataobj$maxAge),ncol=length(dataobj$years))
		ags<-dataobj$agerangeI[j,1]:dataobj$agerangeI[j,2]
		agsidx<-(ags-dataobj$minAge+1)
		ys<-as.numeric(dimnames(dataobj$SurveyList[[j]])[[1]])
		yidx<-(ys-min(dataobj$years)+1)
		x[agsidx,yidx]<-(t(dataobj$SurveyList[[j]]))
		y<-log(XSAMsim$FCrepl[[k]]$surveyobslist[[j]][1:nrow(x),1:ncol(x)])
		x[x<0]<-NA
		x<-log(x)
		ires[[j]]<-x-y
		vmat<-fitobj$res$varmatI[,,j]
		istdres[[j]]<-ires[[j]]/sqrt(vmat)
	}
	list(cres=cres,cstdres=cstdres,ires=ires,istdres=istdres)
}


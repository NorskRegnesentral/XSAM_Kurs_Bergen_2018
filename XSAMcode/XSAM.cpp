// Statistical assessment model from Aanes 2016 in ICES 2016.
// Copyright 2016 Sondre Aanes <sondre.aanes@nr.no>
// All rights reserved.
//  --------------------------------------------------------------------------
// X-version of the Statistical fish stock Assessment Model XSAM. This model 
// is a flexible version of a Time Series model inspired by Gudmundsson 1994
// and Aanes et al 2007. The model is documentd in:
// ICES. 2016. Report of the Benchmark Workshop on Pelagic stocks (WKPELA), 
// 29 February–4 March 2016, ICES Headquarters, Copenhagen, Denmark. ICES CM 2016/ACOM:34. 106pp. 

// 
// Redistribution and use in source and binary forms, with or without
// modification, are permitted provided that the following conditions are met:
//   * Redistributions of source code must retain the above copyright
//     notice, this list of conditions and the following disclaimer.
//   * Redistributions in binary form must reproduce the above copyright
//     notice, this list of conditions and the following disclaimer in the
//     documentation and/or other materials provided with the distribution.
//   * Neither the name of the assessment tool XSAM nor the
//     names of its contributors may be used to endorse or promote products
//     derived from this software without specific prior written permission.
// 
// THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS" 
// AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE 
// IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE 
// ARE DISCLAIMED. IN NO EVENT SHALL SONDRE AANES BE LIABLE FOR ANY DIRECT, 
// INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, 
// BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, 
// DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY 
// OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING 
// NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, 
// EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
//  --------------------------------------------------------------------------

// Some naming of variables, structure of data and parameters is inspired by the template from 
// Nielsen and Berg 2014, Fisheries Research which includes the SAM model. See the file sam.cpp in 
//"https://github.com/kaskr/adcomp/tree/master/tmb_examples" and copyright issues for the SAM template follows here:


// State space assessment model from Nielsen and Berg 2014, Fisheries Research.
//  --------------------------------------------------------------------------
// Copyright (c) 2014, Anders Nielsen <an@aqua.dtu.dk>, 
// Casper Berg <cbe@aqua.dtu.dk>, and Kasper Kristensen <kkr@aqua.dtu.dk>.
// All rights reserved.
// 
// Redistribution and use in source and binary forms, with or without
// modification, are permitted provided that the following conditions are met:
//   * Redistributions of source code must retain the above copyright
//     notice, this list of conditions and the following disclaimer.
//   * Redistributions in binary form must reproduce the above copyright
//     notice, this list of conditions and the following disclaimer in the
//     documentation and/or other materials provided with the distribution.
//   * Neither the name of the assessment tool SAM nor the
//     names of its contributors may be used to endorse or promote products
//     derived from this software without specific prior written permission.
// 
// THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS" 
// AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE 
// IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE 
// ARE DISCLAIMED. IN NO EVENT SHALL ANDERS NIELSEN, CASPER BERG OR KASPER 
// KRISTENSEN BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, 
// EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, 
// PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; 
// OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, 
// WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR 
// OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF 
// ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
//  --------------------------------------------------------------------------

 
#include <TMB.hpp>
#include <iostream>


/* Parameter transform */
template <class Type>
Type f(Type x){return Type(2)/(Type(1) + exp(-Type(2) * x)) - Type(1);}

template <class Type> 
Type square(Type x){return x*x;}

template <class Type> 
matrix<Type> NtoLNvar(vector<Type> mu, matrix<Type> Sigma, int dim){
	matrix<Type> vy(dim,dim);
	for(int i=0;i<dim;i++){
		for(int j=0;j<dim;j++){
			vy(i,j)=exp(mu(i)+mu(j)+0.5*(Sigma(i,i)+Sigma(j,j)))*(exp(Sigma(i,j))-1);
		}
	}
	return(vy);
}

template <class Type> 
vector<Type> NtoLNe(vector<Type> mu, matrix<Type> Sigma, int dim){
	vector<Type> ey(dim);
	for(int i=0;i<dim;i++){
		ey(i)=exp(mu(i)+0.5*Sigma(i,i));
	}
	return(ey);
}
 
template <class Type> 
matrix<Type> LNtoNSigma(vector<Type> ey, matrix<Type> vy, int dim){
	matrix<Type> Sigma(dim,dim);
	for(int i=0;i<dim;i++){
		for(int j=0;j<dim;j++){
			Sigma(i,j)=log(vy(i,j)/(ey(i)*ey(j))+1);
		}
	}
	return(vy);
}

template <class Type> 
vector<Type> LNtoNmu(vector<Type> ey, matrix<Type> vy, int dim){
	vector<Type> mu(dim);
	for(int i=0;i<dim;i++){
		mu(i)=log(ey(i))-0.5*log(vy(i,i)/(ey(i)*ey(i))+1);
	}
	return(ey);
}

template <class Type> 
Type VarSum(matrix<Type> vy, vector<Type> a,int dim){
	Type vsum=0;
	for(int i=0;i<dim;i++){
		for(int j=0;j<dim;j++){
			vsum+=a(i)*a(j)*vy(i,j);
		}
	}
	return(vsum);
}

template <class Type> 
matrix<Type> SIGMA(vector<Type> h, vector<Type> sd, matrix<Type> R, int dim){
	matrix<Type> sigma(dim,dim);
	for(int i=0;i<dim;i++){
		for(int j=0;j<dim;j++){
			sigma(i,j)=h(i)*sd(i)*sd(j)*R(i,j);
		}
	}
	return(sigma);
}

template <class Type> 
Type dimnotna(vector<Type> x, Type NAval, int dim){
	int dimuse=0;
	for(int i=0;i<dim;i++){
		if(CppAD::isnan(x(i))==0)dimuse=dimuse+1;
	}
	return(dimuse);
}

template <class Type> 
vector<Type> findnotNAvector(vector<Type> x, Type NAval, int dim){
	vector<Type> vec(dim);
	for(int i=0;i<dim;i++){
		if(CppAD::isnan(x(i))==0)vec(i)=1;
		else vec(i)=0;
	}
	return(vec);
}

template<class Type>
Type objective_function<Type>::operator() ()
{

  DATA_VECTOR(fleetIndex);
  DATA_VECTOR(sampleTimes);
  DATA_INTEGER(noYears);
  DATA_INTEGER(nobs);
  DATA_ARRAY(obs);
  DATA_ARRAY(propMat);
  DATA_ARRAY(stockMeanWeight); 
  DATA_ARRAY(catchMeanWeight);
  DATA_ARRAY(natMor);
  DATA_ARRAY(propF);
  DATA_ARRAY(propM);
  DATA_INTEGER(minAge);
  //DATA_INTEGER(maxAge);
  DATA_INTEGER(maxAgePlusGroup);
  DATA_VECTOR(keyLogFconst);

  DATA_INTEGER(corFlagU);
  DATA_INTEGER(corFlaglogF);

  DATA_INTEGER(Modelq);
  DATA_ARRAY(keyLogqpar);
  DATA_INTEGER(nqs);
  DATA_VECTOR(keyVarF);
  DATA_INTEGER(nvf);
  DATA_ARRAY(keyVarObs);
  DATA_VECTOR(agerangeC);
  DATA_ARRAY(agerangeI); 
  DATA_INTEGER(stockRecruitmentModelCode);
  DATA_INTEGER(RecruitmentProcess);
  DATA_INTEGER(LatentEffort);
  DATA_INTEGER(TSsel);
  DATA_INTEGER(EstimateARInterceptEffort);
  DATA_INTEGER(EstimateUInterceptEffort);
  DATA_INTEGER(ymeantype);
  DATA_INTEGER(am);
  DATA_INTEGER(Am);
  DATA_INTEGER(A);
  DATA_INTEGER(AT);

  DATA_INTEGER(uam);

  DATA_INTEGER(nIndices);
  DATA_INTEGER(nObsVar);
  DATA_VECTOR(SurveyIndex);
  DATA_INTEGER(dimC);
  DATA_VECTOR(dimI);
  DATA_INTEGER(Fbarminage);
  DATA_INTEGER(Fbarmaxage);
  DATA_VECTOR(caton); 
  DATA_INTEGER(CatchConstraint);
  DATA_INTEGER(UseCatchPred); 
  DATA_VECTOR(CatchPrediction);
  DATA_SCALAR(NAval);
  DATA_INTEGER(ReturnLL);

  DATA_ARRAY(sd_C);//matrix holding sd of catch by year (n,y)
  DATA_ARRAY(R_C);//3 dim array holding correlation matrices of catch by year (n,n,y)
  DATA_ARRAY(sd_I);//3 dim array holding sd of index by year and fleet (n,y,f)
  DATA_ARRAY(R_I);//4 dim array holding correlations of index by year and fleet (n,n,y,f)

  PARAMETER_VECTOR(FixedPars);
  PARAMETER_ARRAY(LP);

	
	/*-----------------------------------------------*/
	/*------------	FILL IN PARAMETERS	---------*/
	/*-----------------------------------------------*/
	int ParCounter=0;
	vector<Type> initlogN(A-1);
	for(int i=0;i<(A-1);i++){
		initlogN(i)=FixedPars(ParCounter);
		ParCounter=ParCounter+1;
	}
	vector<Type> R(noYears);
	if(RecruitmentProcess==0){
		for(int i=0;i<noYears;i++){
			R(i)=FixedPars(ParCounter);
			ParCounter=ParCounter+1;
		}
	}
	vector<Type> logqpar(nqs);
	for(int i=0;i<nqs;i++){
		logqpar(i)=FixedPars(ParCounter);
		ParCounter=ParCounter+1;
	}
  	//Type s2_1=exp(FixedPars(ParCounter));
	//ParCounter=ParCounter+1;
	vector<Type> s2_1(nvf);
	for(int i=0;i<nvf;i++){
		s2_1(i)=exp(FixedPars(ParCounter));
		ParCounter=ParCounter+1;
	}
	Type s2_2;
	if(TSsel==1){
		s2_2=exp(FixedPars(ParCounter));
		ParCounter=ParCounter+1;
	}
	Type s2_3;
	if(LatentEffort==1){
		s2_3=exp(FixedPars(ParCounter));
		ParCounter=ParCounter+1;
	}
	Type s2_4=exp(FixedPars(ParCounter));
	ParCounter=ParCounter+1;
	Type sR2;
	if(RecruitmentProcess==1){
		sR2=exp(FixedPars(ParCounter));
		ParCounter=ParCounter+1;
	}
	//vector<Type> hvec(nIndices+1);//ONLY ALLOW ONE VARCOMP FOR EACH DATA SOURCE....
	vector<Type> hvec(nObsVar);
	//for(int i=0;i<(nIndices+1);i++){
	for(int i=0;i<nObsVar;i++){
		hvec(i)=exp(FixedPars(ParCounter));
		ParCounter=ParCounter+1;
	}
	Type rec_loga=FixedPars(ParCounter);
	ParCounter=ParCounter+1;
	Type rec_logb=FixedPars(ParCounter);
	ParCounter=ParCounter+1;
  	Type rhoU;
	if(TSsel==1 & corFlagU>0){
		rhoU=FixedPars(ParCounter);
		ParCounter=ParCounter+1;
	}
	Type rhologF;
	if(corFlaglogF>0){
		rhologF=FixedPars(ParCounter);
		ParCounter=ParCounter+1;
	}
	Type aYpar;
  	if(EstimateARInterceptEffort==1){
		aYpar=FixedPars(ParCounter);
		ParCounter=ParCounter+1;
	}
  	Type bY=FixedPars(ParCounter);
	ParCounter=ParCounter+1;
	vector<Type> Upar(am-1);
	if(TSsel==0){
		for(int i=0;i<(am-1);i++){
			Upar(i)=FixedPars(ParCounter);
			ParCounter=ParCounter+1;
		}
	}
	vector<Type> aUpar(uam);
	if(TSsel==1 & EstimateUInterceptEffort==1){
		for(int i=0;i<uam;i++){
			aUpar(i)=FixedPars(ParCounter);
			ParCounter=ParCounter+1;
		}
	}
	Type bU;
	if(TSsel==1){
		bU=FixedPars(ParCounter);
		ParCounter=ParCounter+1;
	}

	/*---------------------------------------------------------*/

	/*---------------------------------------------------*/
	/*      FILL IN LATENT VARIABLES FOR U, V AND Y      */
	/*---------------------------------------------------*/
	array<Type> logN(AT,noYears);
	array<Type> N(AT,noYears);
	array<Type> logF(Am,noYears);
	array<Type> F(A,noYears);
	array<Type> logM(A-1,noYears);
	array<Type> M(A,noYears);
	array<Type> U(am,noYears);
	vector<Type> V(noYears);
	vector<Type> Y(noYears);


	int dimi=0;
	for(int y=0;y<noYears;y++){
		dimi=0;
		for(int a=0;a<Am;a++){
			logF(a,y)=LP(dimi,y);
			dimi=dimi+1;
		}
		if(TSsel==1){
			for(int a=0;a<(am-1);a++){
				U(a,y)=LP(dimi,y);
				dimi=dimi+1;
			}
		}
		if(LatentEffort==1){
			V(y)=LP(dimi,y);
			dimi=dimi+1;
		}
		Y(y)=LP(dimi,y);
		dimi=dimi+1;
		if(RecruitmentProcess==1){
			R(y)=LP(dimi,y);
			dimi=dimi+1;
		}
		for(int a=0;a<A;a++){
			M(a,y)=natMor(y,a);
		}

	}


	if(TSsel==0){
		for(int y=0;y<noYears;y++){
			for(int a=0;a<(am-1);a++){
				U(a,y)=Upar(a);
			}
		}
	}

	if(LatentEffort==0){
		for(int y=0;y<noYears;y++){
			V(y)=Y(y);
		}
	}


	/*	Take care of constraint in U	*/
	for(int y=0;y<noYears;y++){
		Type usum=0;
		for(int a=0;a<(am-1);a++){
			usum=usum+U(a,y);
		}
		U(am-1,y)=-usum;
	}

	/* Yields expected fishing mortalities...*/
	array<Type> mulogF(Am,noYears);
	for(int y=0;y<noYears;y++){
		for(int a=0;a<am;a++){
			mulogF(a,y)=U(a,y)+V(y);
		}
		if(am<Am){
			for(int a=am;a<Am;a++){
				mulogF(a,y)=U(am-1,y)+V(y);
				std::cout<<"mulogF"<<std::endl;
				std::cout<<mulogF(a,y)<<std::endl;
			}  
		}
	}
	std::cout<<"Am"<<std::endl;
	std::cout<<Am<<std::endl;
	std::cout<<"am"<<std::endl;
	std::cout<<am<<std::endl;

	array<Type> Mlogqpar(nIndices,A);
	if(Modelq==1){
		//Redefine logq; change qs in observation model
		int qindeks=0;
		for(int k=0;k<nIndices;k++){
			Type aq=exp(logqpar(qindeks));
			qindeks=qindeks+1;
			Type bq=exp(logqpar(qindeks));
			qindeks=qindeks+1;
			Type cq=exp(logqpar(qindeks));
			qindeks=qindeks+1;
			int oldkey=-1;
			int newkey=0;
			for(int a=0;a<A;a++){
				//Actual age
				Type age=minAge+a;
				//indeks of in actual age
				//int ageid=CppAD::Integer(agerangeI(k,0))-minAge+a;
				//get identificator 
				newkey=CppAD::Integer(keyLogqpar(k,a));
				std::cout<<"age, newkey"<<std::endl;
				std::cout<<age<<std::endl;
				std::cout<<newkey<<std::endl;
				if(newkey!=-1 & newkey!=oldkey){
					Mlogqpar(k,a)=-aq+bq*age-cq*pow(age,2);
					oldkey=newkey;
				}
				else{
					if(newkey==oldkey & a>0)Mlogqpar(k,a)=Mlogqpar(k,a-1);
					else Mlogqpar(k,a)=-999.0;
				}
			}
			std::cout<<"index"<<std::endl;
			std::cout<<k<<std::endl;
			std::cout<<"q"<<std::endl;
			std::cout<<Mlogqpar<<std::endl;
		}
	}

	//Likelihood components
	//if(ReturnLL==1){
		vector<Type> LLS(noYears);
		for(int y=0;y<noYears;y++)LLS(y)=0;
		vector<Type> llcaa(noYears);
		vector<Type> lliaa1(noYears);
		vector<Type> lliaa2(noYears);
		vector<Type> lliaa3(noYears);
		vector<Type> lliaa4(noYears);
	//}
	/*---------------------------------------------------*/
	/*      F-PROCESS				     */
	/*---------------------------------------------------*/
	using namespace density;
	Type ans=0;

	/*      FIRST EFFORT				     */

	//Not estimating intercept, but use the sum relationship to determine ay
	Type aY;
	if(EstimateARInterceptEffort==0){
		Type ymean=0;
		Type ymean1=0;
		Type ymean2=0;
		for(int y=0;y<(noYears-1);y++){
			ymean=ymean+Y(y);
			ymean1=ymean1+Y(y+1);
			ymean2=ymean2+Y(y);
		}
		ymean=ymean/(noYears-1);
		ymean1=ymean1/(noYears-1);
		ymean2=ymean2/(noYears-1);
		if(ymeantype==1)aY=ymean*(1-bY);
		else aY=ymean1-bY*ymean2;
	}
	else{
		aY=aYpar;
	}

	Type predY;
	for(int y=1;y<noYears;y++){
		predY=aY+bY*Y(y-1);
		ans+=-dnorm(Y(y),predY,sqrt(s2_4),true);// Y-Process likelihood 
	}
	if(ReturnLL==1)LLS(0)=ans;

	//Type temp=ans;
	//std::cout<<"temp0"<<std::endl;
	//std::cerr<<temp<<std::endl;

	if(LatentEffort==1){
		/*Process V*/
		Type predV;
		//for(int y=1;y<noYears;y++){
		for(int y=0;y<noYears;y++){
			predV=Y(y);
			ans+=-dnorm(V(y),predV,sqrt(s2_3),true);// V-Process likelihood 
		}
	}
	if(ReturnLL==1)LLS(1)=ans;

	/*      THEN SELECTIVITY				     */
	 /*TS for UU*/
	if(TSsel){
		//Covariance structure
		//int uam=am-1;//? or am?
		array<Type> fcor_U(uam,uam);
		for(int ai=0; ai<(uam); ++ai){
			for(int aj=0; aj<(uam); ++aj){
				if(corFlagU==0){
					if(ai!=aj){fcor_U(ai,aj)=0;}else{fcor_U(ai,aj)=1.0;}
				}
				if(corFlagU==1){
					if(ai!=aj){fcor_U(ai,aj)=rhoU;}else{fcor_U(ai,aj)=1.0;}
				}
				if(corFlagU==2){
					if(ai!=aj){fcor_U(ai,aj)=pow(rhoU,abs(ai-aj));}else{fcor_U(ai,aj)=1.0;}
				}
			}
		}

		matrix<Type> fvar_U(uam,uam);
		for(int ai=0; ai<uam; ++ai){
			for(int aj=0; aj<uam; ++aj){
				fvar_U(ai,aj)=s2_2*fcor_U(ai,aj);
			}
		}

		/*PREDICT U's*/
		//Not estimating intercept, but use the sum relationship to determine ay
		vector<Type> aU(uam);
		if(EstimateUInterceptEffort==0){
			vector<Type> umean(uam);
			vector<Type> umean1(uam);
			vector<Type> umean2(uam);
			for(int a=0;a<(uam);a++){
				umean(a)=0;
				umean1(a)=0;
				umean2(a)=0;
				for(int y=0;y<(noYears-1);y++){
					umean(a)=umean(a)+U(a,y);
					umean1(a)=umean1(a)+U(a,y+1);
					umean2(a)=umean2(a)+U(a,y);
				}
				umean(a)=umean(a)/(noYears-1);
				umean1(a)=umean1(a)/(noYears-1);
				umean2(a)=umean2(a)/(noYears-1);
				if(ymeantype==1)aU(a)=umean(a)*(1-bU);
				else aU(a)=umean1(a)-bU*umean2(a);
			}
		}
		else{
			aU=aUpar;
		}
		array<Type> predU(uam,noYears-1);
		for(int y=1;y<noYears;y++){
			for(int a=0;a<(uam);a++){
				predU(a,y-1)=aU(a)+bU*U(a,y-1);
			}
		}
		MVNORM_t<Type> neg_log_density_U(fvar_U);
		vector<Type> tempUcol(uam);
		for(int y=1;y<noYears;y++){    
			for(int a=0;a<uam;a++)tempUcol(a)=U(a,y);
			ans+=neg_log_density_U(tempUcol-predU.col(y-1)); // U-Process likelihood 
			//ans+=neg_log_density_U(U.col(y)-predU.col(y-1)); // U-Process likelihood 
		}
	}//END if(TSsel)
	if(ReturnLL==1)LLS(2)=ans;

	/*	FINALLY F	*/

	//Covariance structure

	array<Type> fcor_logF(Am,Am);
	for(int ai=0; ai<Am; ++ai){
		for(int aj=0; aj<Am; ++aj){
			if(corFlaglogF==0){
				if(ai!=aj){fcor_logF(ai,aj)=0;}else{fcor_logF(ai,aj)=1.0;}
			}
			if(corFlaglogF==1){
				if(ai!=aj){fcor_logF(ai,aj)=rhologF;}else{fcor_logF(ai,aj)=1.0;}
			}
			if(corFlaglogF==2){
				if(ai!=aj){fcor_logF(ai,aj)=pow(rhologF,abs(ai-aj));}else{fcor_logF(ai,aj)=1.0;}
			}
		}
	}

	matrix<Type> fvar_logF(Am,Am);
	for(int ai=0; ai<Am; ++ai){
		for(int aj=0; aj<Am; ++aj){
			fvar_logF(ai,aj)=sqrt(s2_1(CppAD::Integer(keyVarF(ai))))*sqrt(s2_1(CppAD::Integer(keyVarF(aj))))*fcor_logF(ai,aj);
		}
	}

	MVNORM_t<Type> neg_log_density_logF(fvar_logF);
	//for(int y=1;y<noYears;y++){    
	for(int y=0;y<noYears;y++){
		ans+=neg_log_density_logF(logF.col(y)-mulogF.col(y)); // F-Process likelihood 
	}
	if(ReturnLL==1)LLS(3)=ans;

	/*---------------------------------------------------*/
	/*      POPULATION-PROCESS			     */
	/*---------------------------------------------------*/
	std::cout<<"AT"<<std::endl;
	std::cout<<AT<<std::endl;
	std::cout<<"A"<<std::endl;
	std::cout<<A<<std::endl;
	logN(0,0)=R(0);

	for(int a=1;a<AT;a++){
		if(a<A)logN(a,0)=initlogN(a-1);
		else logN(a,0)=Type(log(0.00001));
		std::cout<<logN(a,0)<<std::endl;
	}


	vector<Type> ssb(noYears);
	ssb(0)=0.0;
	for(int a=0;a<A;a++){
		ssb(0)+=exp(logN(a,0))*exp(-exp(logF(CppAD::Integer(keyLogFconst(a)),0))*propF(0,a)-M(a,0)*propM(0,a))*propMat(0,a)*stockMeanWeight(0,a);
	}

	vector<Type> predN(noYears);
	for(int y=1;y<noYears;y++){
		if(RecruitmentProcess==1){
			if(stockRecruitmentModelCode==0){//Just an unknown parameter 
				predN(0)=logN(0,y);
			}
			else{
				if(stockRecruitmentModelCode==1){ // straight RW
					predN(0)=logN(0,y-1);
				}
				else{
					if(stockRecruitmentModelCode==2){//Ricker: OBS only valid if minAge==1
						predN(0)=rec_loga+log(ssb(y-1))-exp(rec_logb)*ssb(y-1);  
					}
					else{
						if(stockRecruitmentModelCode==3){//BH: OBS only valid if minAge==1
							predN(0)=rec_loga+log(ssb(y-1))-log(1.0+exp(rec_logb)*ssb(y-1));
						}
						else{
							if(stockRecruitmentModelCode==4){//AR1
								predN(0)=rec_loga+rec_logb*logN(0,y-1);
							}
							else{
								if(stockRecruitmentModelCode==5){//Common mean
									predN(0)=rec_loga;
								}
								else{
									error("SR model code not recognized");
								}
							}
						}
					}
				}
			}
			//True value is the latent variable R
			ans+=-dnorm(R(y),predN(0),sqrt(sR2),true);// Recruitment-Process likelihood
		}
		else{
			predN(0)=R(y);
		}

		logN(0,y)=R(y);
		//for(int a=1; a<A; a++){
		for(int a=1; a<AT; a++){
			if(a<A){
				predN(a)=logN(a-1,y-1)-exp(logF(CppAD::Integer(keyLogFconst(a-1)),y-1))-M(a-1,y-1);
				logN(a,y)=predN(a);
			}
			else{
				predN(a)=logN(a-1,y-1)-exp(logF(CppAD::Integer(keyLogFconst(A-1)),y-1))-M(A-1,y-1);
				logN(a,y)=predN(a);
			}
		}
		if(maxAgePlusGroup==1){
			predN(A-1)=log(exp(logN(A-2,y-1)-exp(logF(CppAD::Integer(keyLogFconst(A-2)),y-1))-M(A-2,y-1))+exp(logN(A-1,y-1)-exp(logF(CppAD::Integer(keyLogFconst(A-1)),y-1))-M(A-1,y-1))); 
			logN(A-1,y)=predN(A-1);
		}

		ssb(y)=0;
		for(int a=0; a<AT; a++){
			if(a<A)ssb(y)+=exp(logN(a,y))*exp(-exp(logF(CppAD::Integer(keyLogFconst(a)),y))*propF(y,a)-M(a,y)*propM(y,a))*propMat(y,a)*stockMeanWeight(y,a);
			else ssb(y)+=exp(logN(a,y))*exp(-exp(logF(CppAD::Integer(keyLogFconst(A-1)),y))*propF(y,A-1)-M(A-1,y)*propM(y,A-1))*propMat(y,A-1)*stockMeanWeight(y,A-1);
		}
	}
	if(ReturnLL==1)LLS(4)=ans;
	ADREPORT(ssb);

	/*------------------------------------*/
	/*		OBSERVATIONS	      */
	/*------------------------------------*/
	//Now observations; 
	int f, a, y;
	int minYear=CppAD::Integer((obs(0,0)));
	/*---------------------------------------------------*/  	
  	//first catch
	/*---------------------------------------------------*/
	/*--------------------------------------------------------------*/
	//OBS fleet given the index 0 is alwas treated as catch at age
	/*--------------------------------------------------------------*/
  	//Setup covariance matrices for observation error...
	vector<Type> h(dimC);
	for(int i=0; i<dimC; ++i)h(i)=hvec(CppAD::Integer(keyVarObs(0,i)));

	//REARRANGE OBSERVATIONS
	array<Type> logcaa(dimC,noYears);
	//Initialize
	for(int y=0;y<noYears;y++){
		for(int a=0;a<dimC;a++){
			//logcaa(a,y)=-999.0;
			logcaa(a,y)=log(NAval);
		}
	}
  
	// Insert data on matrix form AxY into logcaa

	for(int i=0;i<nobs;i++){
		y=CppAD::Integer(obs(i,0))-minYear;
		f=CppAD::Integer(obs(i,1));
		a=CppAD::Integer(obs(i,2))-minAge;
    
		if(f==0){// Catch@age
			logcaa(a,y)=log(obs(i,3));
		}

	}      
	 
	// Match catches to observations
	Type CatchPred=CatchPrediction(0);
	Type PredictCatchVar=CatchPrediction(1);
	vector<Type> predlogCAA(dimC);
  	Type zz;
	vector<Type> predW(noYears);
	Type VlogW;
	Type ElogW;
  	for(int yy=0;yy<noYears;yy++){
		vector<Type> notNAvector=findnotNAvector(tmbutils::vector<Type>(logcaa.col(yy)),NAval,dimC);
		Type dimuse=notNAvector.sum();
		//Model catch at age: 
		for(int aa=0;aa<dimC;aa++){
			zz=exp(logF(CppAD::Integer(keyLogFconst(aa)),yy))+M(aa,yy);
			predlogCAA(aa)=logN(aa,yy)-log(zz)+log(1-exp(-zz))+logF(CppAD::Integer(keyLogFconst(aa)),yy);
		}
		//IF ANY OBSERVATIONS AT ALL
		if(dimuse>0){
			vector<Type> tmppredlogCAA(CppAD::Integer(dimuse));
			vector<Type> tmplogCAA(CppAD::Integer(dimuse));
			vector<Type> sdC(CppAD::Integer(dimuse));
			matrix<Type> RC(CppAD::Integer(dimuse),CppAD::Integer(dimuse));
			vector<Type> hCu(CppAD::Integer(dimuse));
			int id=0;
     			for(int aa=0;aa<dimC;aa++){
				int ageid=CppAD::Integer(agerangeC(0))-minAge+aa;
				if(CppAD::Integer(notNAvector(aa))==1){
					tmplogCAA(id)=logcaa(aa,yy);
					tmppredlogCAA(id)=predlogCAA(ageid);
					sdC(id)=sd_C(aa,yy);
					int id2=0;
					for(int aaa=0; aaa<dimC; aaa++){
						if(CppAD::Integer(notNAvector(aaa))==1){
							RC(id,id2)=R_C(aa,aaa,yy);
							id2=id2+1;
						}
					}
					hCu(id)=h(aa);
					id=id+1;
				}	
			}
			//Set up covariance matrix
			matrix<Type> cvar=SIGMA(hCu,sdC,RC,CppAD::Integer(dimuse));

			//Define sampling distribution
			MVNORM_t<Type> neg_log_densityCAA(cvar);
			//Likelihood contribution if observations available, i.e. !=-999.0
			ans+=neg_log_densityCAA(tmplogCAA-tmppredlogCAA); // Observation likelihood 


	
		}//END if(dimuse>0)
		llcaa(yy)=ans;
	}
	if(ReturnLL==1)LLS(5)=ans;
	/*---------------------------------------------------------------*/
	/*    THE FOLLOWING CONCERNS CALCULATION OF TOTAL CATCH WEIGHTS  */
	/*---------------------------------------------------------------*/ 
	//Should estimated reported total catch weight contribute to the likelihood together with CAA? Probably not! However can be technically implemented as the following...
	vector<Type> cmw(dimC);//Temporary variable for mean catch weight
	vector<Type> ECAA(dimC); 
	matrix<Type> VCAA(dimC,dimC);
	//REUSE VARIABLES
	vector<Type> sdC(dimC);
	matrix<Type> RC(dimC,dimC);
	Type bcor;
	for(int yy=0;yy<noYears;yy++){
		predW(yy)=0;
		//Model catch at age: 
		for(int aa=0;aa<dimC;aa++){
			zz=exp(logF(CppAD::Integer(keyLogFconst(aa)),yy))+M(aa,yy);
			predlogCAA(aa)=logN(aa,yy)-log(zz)+log(1-exp(-zz))+logF(CppAD::Integer(keyLogFconst(aa)),yy);
			sdC(aa)=sd_C(aa,yy);
			for(int aaa=0; aaa<dimC; aaa++){
				RC(aa,aaa)=R_C(aa,aaa,yy);
			}
		}
		//Set up covariance matrix
		matrix<Type> cvar=SIGMA(h,sdC,RC,dimC);
		//#OOBS: bias correct inserted 30/8-2017 compared to WGWIDE 2016
		for(int aa=0;aa<dimC;aa++){
		//for(int aa=0;aa<dimuse;aa++){
			bcor=-0.5*cvar(aa,aa);
			predlogCAA(aa)=predlogCAA(aa)+bcor;
		}

		if(yy<(noYears-1))cmw=catchMeanWeight.transpose().col(yy);
		else cmw=catchMeanWeight.transpose().col(noYears-2);//SET CATCH WEIGHT IN LAST YEAR EQUAL TO THE PREVIOUS YEAR
		//First transform moments to original scale
		ECAA=NtoLNe(predlogCAA, cvar, dimC);
		VCAA=NtoLNvar(predlogCAA, cvar, dimC);
		for(int aa=0;aa<dimC;aa++){
			predW(yy)=predW(yy)+ECAA(aa)*cmw(aa);
		}
		//Find corresponding variance
		Type varW=VarSum(VCAA,cmw,dimC);
		//Back to log scale
		VlogW=log(varW/pow(predW(yy),2)+1);
		ElogW=log(predW(yy))-0.5*VlogW;
		//Require that caton is available until the second last year
		if(CatchConstraint==1 & yy<(noYears-1))ans+=-dnorm(log(caton(yy)),ElogW,sqrt(VlogW),true);
		//If use predicted total catch contribute to the likelihood for estimation of F in assessment year
	}
	if(ReturnLL==1)LLS(6)=ans;
	if(UseCatchPred==1){
		Type totVar=VlogW+PredictCatchVar;
		ans+=-dnorm(log(CatchPred),ElogW,sqrt(totVar),true);
	}
	if(ReturnLL==1)LLS(7)=ans;
	ADREPORT(predW);

	/*---------------------------------------------------*/
	// Then survey Indices
	/*---------------------------------------------------*/

	for(int k=0;k<nIndices;k++){
		vector<Type> hI(CppAD::Integer(dimI(k)));
		for(int i=0; i<dimI(k); ++i)hI(i)=hvec(CppAD::Integer(keyVarObs(k+1,i)));
		//REARRANGE OBSERVATIONS
		array<Type> logindex(CppAD::Integer(dimI(k)),noYears);
		//Initialize by setting all obs to missing
		for(int yy=0;yy<noYears;yy++){
			for(int aa=0;aa<dimI(k);aa++){
				logindex(aa,yy)=log(NAval);
			}

		}
  
		// Insert data on matrix form AxY into logcaa and logindex

		for(int i=0;i<nobs;i++){
			y=CppAD::Integer(obs(i,0))-minYear;
			f=CppAD::Integer(obs(i,1));
			a=CppAD::Integer(obs(i,2))-CppAD::Integer(agerangeI(k,0));
    
			if(f==SurveyIndex(k)){// survey
				logindex(a,y)=log(obs(i,3));
			}
		}

		Type bcor;
		Type zz;
		for(int yy=0;yy<noYears;yy++){
			std::cout<<"yy"<<std::endl;
			std::cout<<yy<<std::endl;
			//Find vector of indicators for NA's
			vector<Type> notNAvector=findnotNAvector(tmbutils::vector<Type>(logindex.col(yy)),NAval,CppAD::Integer(dimI(k)));
			//FIND DIMENSION OF NON NA'S
			Type dimuse=notNAvector.sum();
			//IF ANY OBSERVATIONS AT ALL
			if(dimuse>0){ 
				vector<Type> tmplogIndex(CppAD::Integer(dimuse));
				vector<Type> predlogIndex(CppAD::Integer(dimuse));
				vector<Type> sdI(CppAD::Integer(dimuse));
				matrix<Type> RI(CppAD::Integer(dimuse),CppAD::Integer(dimuse));
				vector<Type> hIu(CppAD::Integer(dimuse));
				int id=0;
				for(int aa=0;aa<dimI(k);aa++){
					//Identify the age index in the data relative to the age range defined by minAge:maxAge
					int ageid=CppAD::Integer(agerangeI(k,0))-minAge+aa;
					std::cout<<"ageid"<<std::endl;
					std::cout<<ageid<<std::endl;
					//Get corresponding z...
					zz=exp(logF(CppAD::Integer(keyLogFconst(ageid)),yy))+M(ageid,yy);
					//...and prediction based on the corresponding agerange for the population
					if(CppAD::Integer(notNAvector(aa))==1){
						tmplogIndex(id)=logindex(aa,yy);
						//keep only data that not is missing
						if(Modelq==0) predlogIndex(id)=logN(ageid,yy)-zz*sampleTimes(k+1)+logqpar(CppAD::Integer(keyLogqpar(k,ageid)));
						else predlogIndex(id)=logN(ageid,yy)-zz*sampleTimes(k+1)+Mlogqpar(k,ageid);
						sdI(id)=sd_I(aa,yy,k);
						int id2=0;
						for(int aaa=0; aaa<dimI(k); aaa++){
							if(CppAD::Integer(notNAvector(aaa))==1){
								RI(id,id2)=R_I(aa,aaa,yy,k);
								id2=id2+1;
							}
						}
						hIu(id)=hI(aa);
						id=id+1;
					}
				}
				std::cout<<"sdI"<<std::endl;
				std::cout<<sdI<<std::endl;
				std::cout<<"RI"<<std::endl;
				std::cout<<RI<<std::endl;
				std::cout<<"tmplogIndex"<<std::endl;
				std::cout<<tmplogIndex<<std::endl;
				matrix<Type> ivar=SIGMA(hIu,sdI,RI,CppAD::Integer(dimuse));
				MVNORM_t<Type> neg_log_densityIndex(ivar);
				ans+=neg_log_densityIndex(tmplogIndex-predlogIndex); // Observation likelihood
				//if(tmplogIndex(0)!=-999.0)ans+=neg_log_densityIndex(tmplogIndex-predlogIndex); // Observation likelihood
			}//END if(dimuse>0)
			if(ReturnLL==1){
				if(k==0)lliaa1(yy)=ans;
				if(k==1)lliaa2(yy)=ans;
				if(k==2)lliaa3(yy)=ans;
				if(k==3)lliaa4(yy)=ans;
			}
		}
		if(ReturnLL==1)LLS(8+k)=ans;
	}//end for k



	/*---------------------------------------------------*/
	//	For report
	/*---------------------------------------------------*/
	
	for(int y=0;y<noYears;y++){
		for(int a=0;a<A;a++){
			if(a<Am)F(a,y)=exp(logF(a,y));
			else F(a,y)=exp(logF(Am-1,y));
		}
	}
	for(int y=0;y<noYears;y++){
		//for(int a=0;a<A;a++){
		for(int a=0;a<AT;a++){
			N(a,y)=exp(logN(a,y));
		}
	}
	//Total biomass
	vector<Type> totb(noYears);
	for(int y=0;y<noYears;y++){
		totb(y)=0;
		for(int a=0; a<AT; a++){
			if(a<A)totb(y)+=exp(logN(a,y))*stockMeanWeight(y,a);
			else totb(y)+=exp(logN(a,y))*stockMeanWeight(y,A-1);
		}
	}

	//Calculate average F
	int nfbar=Fbarmaxage-Fbarminage+1;
	vector<Type> Fbar(noYears);
	vector<Type> FbarW(noYears);
	Type Wsum;
	for(int i=0;i<noYears;i++){
		Fbar(i)=0;
		FbarW(i)=0;
		Wsum=0;
		//for(int j=(Fbarminage-1);j<Fbarmaxage;j++){
		for(int j=(Fbarminage-minAge);j<(Fbarmaxage-minAge+1);j++){
			Fbar(i)=Fbar(i)+F(j,i);
			FbarW(i)=FbarW(i)+N(j,i)*F(j,i);
			Wsum=Wsum+N(j,i);
		}
		Fbar(i)=Fbar(i)/nfbar;
		FbarW(i)=FbarW(i)/Wsum;
	}
	ADREPORT(Fbar); 
	ADREPORT(FbarW);
	ADREPORT(logN);
	ADREPORT(logF);
	ADREPORT(N);
	ADREPORT(F);
	ADREPORT(totb);
	if(ReturnLL==1){
		ADREPORT(LLS);
		ADREPORT(llcaa);
		if(nIndices>0)ADREPORT(lliaa1);
		if(nIndices>1)ADREPORT(lliaa2);
		if(nIndices>2)ADREPORT(lliaa3);
		if(nIndices>3)ADREPORT(lliaa4);
	}
	return ans;
}


rm(list = ls())
library(lubridate)
library(mvtnorm)
library(MASS)
library(gtools)
library(MCMCpack)


# Set error level, regression error sd and conditional association between XA and XB
## er: error level. 'P'=no error, 'L'=20% error, 'H'=40% error
## er_sd: error sd in regression. 'H' = 0.5, 'L'=0.1
## a_par: conditional association between XA and XB. Can be in{0.5,1,2,5,10}, corresponding to beta_M={0.05,0.1,0.2,0.5,1}
er = "P"
er_sd = "H"
a_par = 10

# Import file that sets above parameter sizes
source("Sizes.R")

# Define metrics to export
Sensitivity=rep(0,nsim); Sensitivity_sd=rep(0,nsim)
PPV=rep(0,nsim); PPV_sd=rep(0,nsim)
F1=rep(0,nsim); F1_sd=rep(0,nsim)
B1_bias=rep(0,nsim); B1_bias_sd=rep(0,nsim); B1c_bias=rep(0,nsim); B1c_bias_sd=rep(0,nsim)
B1_RMSE=rep(0,nsim);B1c_RMSE=rep(0,nsim);cor_true=rep(0,nsim)

# Start simulation
for(s in 1:nsim){
  
  set.seed(s)
  print(s)
  #########################################################################################################
  #Simulate data for linked individuals
  
  #V1: Date of Birth, Age~N(50 yrs, 5 yrs ^2), then convert age to DOB
  V1 = rnorm(n_m,50,5)
  V1 = 2019-V1
  DOB = format(date_decimal(V1), "%d-%m-%Y")
  DOB_year = format(date_decimal(V1), "%Y")
  DOB_month = format(date_decimal(V1), "%m")
  DOB_day = format(date_decimal(V1), "%d")
  
  #V2: Zip Code, Stepwise Distribution based on digits of zip code
  Zip1 = sample(c(1,2,3), size=n_m, replace=T)
  Zip2 = sample(c(3,4,5,6), size=n_m, replace=T)
  Zip3 = sample(c(7,8,9,0,1), size=n_m, replace=T)
  
  #V3: Gender, U(1,2)
  Gender = round(runif(n_m,1,2))
  
  #Covariate X and outcome Y
  X_m=mvrnorm(n_m,mu=rep(x_mn,4),Sigma=x_var*diag(4))
  Y_1_m=rnorm(n_m,10+cond_assoc*(a_m1*X_m[,1]+a_m2*X_m[,2]+a_m3*X_m[,3]+a_m4*X_m[,4]), error_sd)
  
  #Unique match IDs
  MatchID_m = round(runif(n_m,1000000,10000000))
  
  #########################################################################################################
  #Simulate data for unlinked individuals in file A
  
  #V1: Date of Birth, Age~N(50 yrs, 5 yrs ^2), then convert age to DOB
  V1_A = rnorm(N_A-n_m,50,5)
  V1_A = 2019-V1_A
  DOB_A = format(date_decimal(V1_A), "%d-%m-%Y")
  DOB_year_A = format(date_decimal(V1_A), "%Y")
  DOB_month_A = format(date_decimal(V1_A), "%m")
  DOB_day_A = format(date_decimal(V1_A), "%d")
  
  #V2: Zip Code, Stepwise Distribution based on digits of zip code
  Zip1_A = sample(c(1,2,3), size=N_A-n_m, replace=T)
  Zip2_A = sample(c(3,4,5,6), size=N_A-n_m, replace=T)
  Zip3_A = sample(c(7,8,9,0,1), size=N_A-n_m, replace=T)
  
  #V3: Gender, U(1,2)
  Gender_A = round(runif(N_A-n_m,1,2))
  
  #Covariate X and outcome Y
  X_A=mvrnorm(N_A-n_m,mu=rep(x_mn,4),Sigma=x_var*diag(4))
  Y_1_A=rnorm(N_A-n_m,5+cond_assoc*(a_u1*X_A[,1]+a_u2*X_A[,2]+a_u3*X_A[,3]+a_u4*X_A[,4]),error_sd)
  
  
  #Unique match IDs
  MatchID_A = sample(c(1000000:10000000),N_A-n_m,replace=F)
  
  
  #########################################################################################################
  #Simulate data for unlinked individuals in file B
  
  #V1: Date of Birth, Age~N(50 yrs, 5 yrs ^2), then convert age to DOB
  V1_B = rnorm(N_B-n_m,50,5)
  V1_B = 2017-V1_B
  DOB_B = format(date_decimal(V1_B), "%d-%m-%Y")
  DOB_year_B = format(date_decimal(V1_B), "%Y")
  DOB_month_B = format(date_decimal(V1_B), "%m")
  DOB_day_B = format(date_decimal(V1_B), "%d")
  
  #V2: Zip Code, Stepwise Distribution based on digits of zip code
  Zip1_B = sample(c(1,2,3), size=N_B-n_m, replace=T)
  Zip2_B = sample(c(3,4,5,6), size=N_B-n_m, replace=T)
  Zip3_B = sample(c(7,8,9,0,1), size=N_B-n_m, replace=T)
  
  #V3: Gender, U(1,2)
  Gender_B = round(runif(N_B-n_m,1,2))
  
  #Covariate X and outcome Y
  X_B=mvrnorm(N_B-n_m,mu=rep(x_mn,4),Sigma=x_var*diag(4))
  Y_1_B=rnorm(N_B-n_m,5+cond_assoc*(a_u1*X_B[,1]+a_u2*X_B[,2]+a_u3*X_B[,3]+a_u4*X_B[,4]),error_sd)
  
  #Unique match IDs
  MatchID_B = sample(c(1000000:10000000),N_B-n_m,replace=F)
  
  
  
  
  #Merge Linked Data together
  LinkData.A = data.frame(MatchID_m, DOB_year, DOB_month, DOB_day, Zip1, Zip2, Zip3, Gender,Y_1_m)
  
  LinkData.B = data.frame(MatchID_m, DOB_year, DOB_month, DOB_day, Zip1, Zip2, Zip3, Gender,
                          B_1_m=X_m[,1],B_2_m=X_m[,2],B_3_m=X_m[,3],B_4_m=X_m[,4])
  
  #Merge Unlinked data in file A together
  N_A_unlinked = data.frame(MatchID_A, DOB_year_A, DOB_month_A, DOB_day_A, Zip1_A, Zip2_A, Zip3_A, Gender_A, Y_1_A)
  names(N_A_unlinked)=names(LinkData.A)
  
  #Merge unlinked data in file B together
  N_B_unlinked = data.frame(MatchID_B, DOB_year_B, DOB_month_B, DOB_day_B, Zip1_B, Zip2_B, Zip3_B, Gender_B, 
                            X_B[,1], X_B[,2], X_B[,3], X_B[,4])
  names(N_B_unlinked)=names(LinkData.B)
  
  #Combine data to form datasets A and B
  Data_A=rbind(LinkData.A,N_A_unlinked)
  Data_B=rbind(LinkData.B,N_B_unlinked)
  levels(Data_A$DOB_year)=unique(c(levels(Data_A$DOB_year),levels(Data_B$DOB_year)))
  levels(Data_B$DOB_year)=unique(c(levels(Data_A$DOB_year),levels(Data_B$DOB_year)))
  
  ############################################################
  #Introduce error into entries in dataset A
  #Error in DOB
  DOB_day_error=rep(0,nrow(Data_A))
  DOB_month_error=rep(0,nrow(Data_A))
  DOB_year_error=rep(0,nrow(Data_A))
  ZIP1_error=rep(0,nrow(Data_A))
  ZIP2_error=rep(0,nrow(Data_A))
  ZIP3_error=rep(0,nrow(Data_A))
  Gender_error=rep(0,nrow(Data_A))
  
  for (i in 1:N_A){
    DOBerror=runif(1)
    ZIPerror=runif(1)
    Gendererror=runif(1)
    
    if(DOBerror>=0 & DOBerror<=(prob*(1)/(3*10))){
      Data_A$DOB_day[i]=sample(Data_B$DOB_day[Data_B$DOB_day!=Data_A$DOB_day[i]],1)
      DOB_day_error[i]=1
    } else if (DOBerror>(prob*(1)/(3*10)) & DOBerror<=(prob*(2)/(3*10))){
      Data_A$DOB_month[i]=sample(Data_B$DOB_month[Data_B$DOB_month!=Data_A$DOB_month[i]],1)
      DOB_month_error[i]=1
    } else if (DOBerror>(prob*(2)/(3*10)) & DOBerror<=(prob*(3)/(3*10))){
      Data_A$DOB_year[i]=sample(Data_B$DOB_year[Data_B$DOB_year!=Data_A$DOB_year[i]],1)
      DOB_year_error[i]=1}
    
    if(ZIPerror>=0 & ZIPerror<=(prob*(1)/(3*10))){
      Data_A$Zip1[i]=sample(Data_B$Zip1[Data_B$Zip1!=Data_A$Zip1[i]],1)
      ZIP1_error[i]=1
    } else if (ZIPerror>(prob*(1)/(3*10)) & ZIPerror<=(prob*(2)/(3*10))){
      Data_A$Zip2[i]=sample(Data_B$Zip2[Data_B$Zip2!=Data_A$Zip2[i]],1)
      ZIP2_error[i]=1
    } else if (ZIPerror>(prob*(2)/(3*10)) & ZIPerror<=(prob*(3)/(3*10))){
      Data_A$Zip3[i]=sample(Data_B$Zip3[Data_B$Zip3!=Data_A$Zip3[i]],1)
      ZIP3_error[i]=1
    }
  } 
  
  
  
  ############################################################
  #Create Comparison matrix Gamma for each of the fields in N_A and N_B
  DOB_year_comparison=t(matrix(Data_B$DOB_year,nrow=length(Data_B$DOB_year),ncol=nrow(Data_A)))
  DOB_month_comparison=t(matrix(Data_B$DOB_month,nrow=length(Data_B$DOB_month),ncol=nrow(Data_A)))
  DOB_day_comparison=t(matrix(Data_B$DOB_day,nrow=length(Data_B$DOB_day),ncol=nrow(Data_A)))
  Zip1_comparison=t(matrix(Data_B$Zip1,nrow=length(Data_B$Zip1),ncol=nrow(Data_A)))
  Zip2_comparison=t(matrix(Data_B$Zip2,nrow=length(Data_B$Zip2),ncol=nrow(Data_A)))
  Zip3_comparison=t(matrix(Data_B$Zip3,nrow=length(Data_B$Zip3),ncol=nrow(Data_A)))
  Gender_comparison=t(matrix(Data_B$Gender,nrow=length(Data_B$Gender),ncol=nrow(Data_A)))
  
  #Compare each variable in dataset A with each element in dataset B element-wise
  DOB_year_gamma=apply(DOB_year_comparison,2,'==',Data_A$DOB_year)
  DOB_month_gamma=apply(DOB_month_comparison,2,'==',Data_A$DOB_month)
  DOB_day_gamma=apply(DOB_day_comparison,2,'==',Data_A$DOB_day)
  Zip1_gamma=apply(Zip1_comparison,2,'==',Data_A$Zip1)
  Zip2_gamma=apply(Zip2_comparison,2,'==',Data_A$Zip2)
  Zip3_gamma=apply(Zip3_comparison,2,'==',Data_A$Zip3)
  Gender_gamma=apply(Gender_comparison,2,'==',Data_A$Gender)
  
  ###Create Gamma matrices from element-wise comparisons
  #Date of Birth
  Gamma2=DOB_year_gamma-DOB_month_gamma
  Gamma2[Gamma2!=1]=0
  Gamma3=DOB_year_gamma+DOB_month_gamma-DOB_day_gamma-1
  Gamma3[Gamma3!=1]=0
  Gamma4=DOB_year_gamma+DOB_month_gamma+DOB_day_gamma-2
  Gamma4[Gamma4!=1]=0
  Gamma1=-(Gamma2+Gamma3+Gamma4)+1
  Gamma1[Gamma1<=0]=0
  #Zip Code
  Gamma6=Zip1_gamma-Zip2_gamma
  Gamma6[Gamma6!=1]=0
  Gamma7=Zip1_gamma+Zip2_gamma-Zip3_gamma-1
  Gamma7[Gamma7!=1]=0
  Gamma8=Zip1_gamma+Zip2_gamma+Zip3_gamma-2
  Gamma8[Gamma8!=1]=0
  Gamma5=-(Gamma6+Gamma7+Gamma8)+1
  Gamma5[Gamma5<=0]=0
  #Gender
  Gamma10=Gender_gamma+0
  Gamma9=-(Gender_gamma)+1
  
  #Create matrix for covariate and outcome
  Y_mat=matrix(Data_A$Y_1_m,nrow=nrow(Data_A),ncol=nrow(Data_B))
  B1_mat=t(matrix(Data_B$B_1_m,nrow=nrow(Data_B),ncol=nrow(Data_A)))
  B2_mat=t(matrix(Data_B$B_2_m,nrow=nrow(Data_B),ncol=nrow(Data_A)))
  B3_mat=t(matrix(Data_B$B_3_m,nrow=nrow(Data_B),ncol=nrow(Data_A)))
  B4_mat=t(matrix(Data_B$B_4_m,nrow=nrow(Data_B),ncol=nrow(Data_A)))
  
  #Create matrix for IDs
  ID1=matrix(Data_A$MatchID_m,nrow=N_A,ncol=N_B)
  ID2=t(matrix(Data_B$MatchID_m,nrow=N_B,ncol=N_A))
  
  ######################################################################
  ######################################################################
  #Collapse matrix data to data frame
  Gamma=data.frame(ID1=c(ID1),ID2=c(ID2),Gamma1=c(Gamma1),Gamma2=c(Gamma2),Gamma3=c(Gamma3),Gamma4=c(Gamma4),Gamma5=c(Gamma5),
                   Gamma6=c(Gamma6),Gamma7=c(Gamma7),Gamma8=c(Gamma8),Gamma9=c(Gamma9),Gamma10=c(Gamma10),Y=c(Y_mat),
                   B1=c(B1_mat),B2=c(B2_mat),B3=c(B3_mat),B4=c(B4_mat))
  X=matrix(cbind(rep(1,dim(Gamma)[1]),Gamma$B1,Gamma$B2,Gamma$B3,Gamma$B4),ncol=5,nrow=dim(Gamma)[1])
  
  
  #Create and initialize linkage structure
  C=rep(0,dim(Gamma)[1])
  C[c(1,502,1003,1504,2005,2506,3007,3508,4009,4510,5011,5512,6013,6514,7015)]=1
  
  
  #Specify hyperparameters for prior distribution 
  prior_DOB_M=c(1,1,1,1)
  prior_DOB_U=c(1,1,1,1)
  prior_ZIP_M=c(1,1,1,1)
  prior_ZIP_U=c(1,1,1,1)
  prior_Gender_M=c(1,1)
  prior_Gender_U=c(1,1)
  prior_pi=c(1,1)
  
  #Initialize parameters of linkage model
  theta_M_DOB=matrix(0,nrow=length(prior_DOB_M),ncol=ngibbs)
  theta_U_DOB=matrix(0,nrow=length(prior_DOB_U),ncol=ngibbs)
  theta_M_ZIP=matrix(0,nrow=length(prior_ZIP_M),ncol=ngibbs)
  theta_U_ZIP=matrix(0,nrow=length(prior_ZIP_U),ncol=ngibbs)
  theta_M_Gender=matrix(0,nrow=length(prior_Gender_M),ncol=ngibbs)
  theta_U_Gender=matrix(0,nrow=length(prior_Gender_U),ncol=ngibbs)
  pi_M=matrix(0,nrow=1,ncol=ngibbs)
  
  #Initialize regression paramters
  beta_M=matrix(0,nrow=5,ncol=ngibbs)
  beta_U=matrix(0,nrow=5,ncol=ngibbs)
  sigma_M=matrix(0,nrow=1,ncol=ngibbs)
  sigma_U=matrix(0,nrow=1,ncol=ngibbs)
  
  # Intermediary quantities
  LinkDesignation=matrix(0,nrow=N_A,ncol=ngibbs)
  LinkProbability=matrix(0,nrow=N_A,ncol=ngibbs)
  beta_1=rep(0,ngibbs);cor_1 = rep(0,ngibbs)
  mean_1=rep(0,ngibbs)
  A_link=rep(0,dim(Data_A)[1])
  A_link[Data_A$MatchID_m%in%Gamma[C==1,]$ID1]=Gamma[C==1,]$ID2
  
  #Specify starting values for parameters
  theta_M_DOB[,1]=c(.01,.04,.05,.9)
  theta_U_DOB[,1]=c(.9099,.05,.04,.0001)
  theta_M_ZIP[,1]=c(.01,.04,.05,.9)
  theta_U_ZIP[,1]=c(.9099,.05,.04,.0001)
  theta_M_Gender[,1]=c(.25,.75)
  theta_U_Gender[,1]=c(.75,.25)
  start_coef=summary(lm(Y~B1+B2+B3+B4-1,data=Gamma))$coef[,1]
  beta_M[,1]=c(rnorm(1,10,2),rnorm(1,a_m1,1),rnorm(1,a_m2,1),rnorm(1,a_m3,1),rnorm(1,a_m4,1))
  beta_U[,1]=c(rnorm(1,5,2),rnorm(1,start_coef[1],.1),rnorm(1,start_coef[2],.1),rnorm(1,start_coef[3],.1),rnorm(1,start_coef[4],.1))
  sigma_M[1]=a_m1^2+a_m2^2+a_m3^2+a_m4^2
  sigma_U[1]=summary(lm(Y~B1+B2+B3+B4,data=Gamma))$sigma
  
  
  # Start posterior sampling
  for(t in 1:(ngibbs-1)){
    
    #########################################################################################################################
    ##Iterate through the rows of file A
    for(i in 1:N_A){
      #Reset the results of C and A_link for row i
      C[Gamma$ID1==Data_A$MatchID_m[i]]=0
      A_link[i]=0
      
      #Extract the portions of Gamma that correspond to row i
      Gamma_i=Gamma[Gamma$ID1==Data_A$MatchID_m[i],]
      
      #Identify and remove any records from file B that already have a link
      Data_B_unlinked=data.frame(ID2=setdiff(Gamma_i$ID2,A_link))
      Gamma_unlinked=merge(Gamma_i,Data_B_unlinked, by="ID2")
      
      #Calculate the ratio of likelihoods for the Gamma components
      num=apply(t(theta_M_DOB[,t]^t(Gamma_unlinked[,c(3:6)])),1,prod)*apply(theta_M_ZIP[,t]^t(Gamma_unlinked[,c(7:10)]),2,prod)*
        apply(theta_M_Gender[,t]^t(Gamma_unlinked[,c(11:12)]),2,prod)
      den=apply(theta_U_DOB[,t]^t(Gamma_unlinked[,c(3:6)]),2,prod)*apply(theta_U_ZIP[,t]^t(Gamma_unlinked[,c(7:10)]),2,prod)*
        apply(theta_U_Gender[,t]^t(Gamma_unlinked[,c(11:12)]),2,prod)
      Likelihood=log(num/den)
      
      #Calculate the ratio of likelihoods for the XA and XB components
      num_Y=dnorm(Gamma_unlinked$Y, mean=matrix(cbind(rep(1,dim(Gamma_unlinked)[1]),Gamma_unlinked$B1,Gamma_unlinked$B2,Gamma_unlinked$B3,Gamma_unlinked$B4),ncol=5)%*%beta_M[,t], sd=sqrt(sigma_M[t]),log=T)
      den_Y=dnorm(Gamma_unlinked$Y, mean=matrix(cbind(rep(1,dim(Gamma_unlinked)[1]),Gamma_unlinked$B1,Gamma_unlinked$B2,Gamma_unlinked$B3,Gamma_unlinked$B4),ncol=5)%*%beta_U[,t], sd=sqrt(sigma_U[t]),log=T)
      Likelihood_Y=num_Y-den_Y
      
      #Calculate the probability of individual i not linking
      p_nolink=(N_B-sum(C))*(N_A-sum(C)+prior_pi[2]-1)/(sum(C)+prior_pi[1])
      
      #Parse togeter possible moves and move probability
      B_unlinked=c(Gamma_unlinked$ID2,0)
      Likelihood_All = Likelihood+Likelihood_Y; hh = c(Likelihood_All, log(p_nolink)) - max(c(Likelihood_All, log(p_nolink)));B_prob=exp(hh)/sum(exp(hh))
      
      #Sample new bipartite link designation for individual i
      link_designation=sample(B_unlinked,size=1,prob=B_prob)
      
      if(link_designation!=0){
        A_link[i]=link_designation
        C[Gamma$ID1==Data_A$MatchID_m[i] & Gamma$ID2==link_designation]=1
        LinkDesignation[i,t]=link_designation
      }
    }
    
    ############################################################################################################################
    #Sample the Posterior Distribution of the parameters
    theta_M_DOB[,t+1]=rdirichlet(1,c(prior_DOB_M[1]+sum(C*Gamma$Gamma1),prior_DOB_M[2]+sum(C*Gamma$Gamma2),prior_DOB_M[3]+sum(C*Gamma$Gamma3),
                                     prior_DOB_M[4]+sum(C*Gamma$Gamma4)))
    theta_U_DOB[,t+1]=rdirichlet(1,c(prior_DOB_U[1]+sum((1-C)*Gamma$Gamma1),prior_DOB_U[2]+sum((1-C)*Gamma$Gamma2),prior_DOB_U[3]+sum((1-C)*Gamma$Gamma3),
                                     prior_DOB_U[4]+sum((1-C)*Gamma$Gamma4)))
    theta_M_ZIP[,t+1]=rdirichlet(1,c(prior_ZIP_M[1]+sum(C*Gamma$Gamma5),prior_ZIP_M[2]+sum(C*Gamma$Gamma6),prior_ZIP_M[3]+sum(C*Gamma$Gamma7),
                                     prior_ZIP_M[4]+sum(C*Gamma$Gamma8)))
    theta_U_ZIP[,t+1]=rdirichlet(1,c(prior_ZIP_U[1]+sum((1-C)*Gamma$Gamma5),prior_ZIP_U[2]+sum((1-C)*Gamma$Gamma6),prior_ZIP_U[3]+sum((1-C)*Gamma$Gamma7),
                                     prior_ZIP_U[4]+sum((1-C)*Gamma$Gamma8)))
    theta_M_Gender[,t+1]=rdirichlet(1,c(prior_Gender_M[1]+sum(C*Gamma9),prior_Gender_M[2]+sum(C*Gamma10)))
    theta_U_Gender[,t+1]=rdirichlet(1,c(prior_Gender_U[1]+sum((1-C)*Gamma10),prior_Gender_U[2]+sum((1-C)*Gamma10)))
    
    if(sum(C)<1){
      beta_M[,t+1]=c(rnorm(1,a_m1,1),rnorm(1,a_m1,1),rnorm(1,a_m1,1),rnorm(1,a_m1,1))
      sigma_M[t+1]=1
    } else{
      beta_M[,t+1]=rmvnorm(1,mean=solve(t(X[C==1,])%*%X[C==1,])%*%t(X[C==1,])%*%Gamma$Y[C==1],sigma=sigma_M[t]*solve(t(X[C==1,])%*%X[C==1,]))
      sigma_M[t+1]=rinvgamma(1, sum(C)/2, (1/2)*t(Gamma$Y[C==1]-X[C==1,]%*%beta_M[,t+1])%*%(Gamma$Y[C==1]-X[C==1,]%*%beta_M[,t+1]))
    }
    beta_U[,t+1]=rmvnorm(1,mean=solve(t(X[C==0,])%*%X[C==0,])%*%t(X[C==0,])%*%Gamma$Y[C==0],sigma=sigma_U[t]*solve(t(X[C==0,])%*%X[C==0,]))
    sigma_U[t+1]=rinvgamma(1, sum(1-C)/2, (1/2)*t(Gamma$Y[C==0]-X[C==0,]%*%beta_U[,t+1])%*%(Gamma$Y[C==0]-X[C==0,]%*%beta_U[,t+1]))
    
    
    ############################################################################################################################
    #Calculate Association for linked Sample
    Data=data.frame(Data_A,A_link)
    Merge_data=merge(Data,Data_B,by.x="A_link",by.y="MatchID_m")
    beta_1[t]=summary(lm(Y_1_m~B_1_m,data=Merge_data))$coef[2,1]
    cor_1[t]=cor(Merge_data$Y_1_m,Merge_data$B_1_m)
    mean_1[t]=mean(Merge_data$Y_1_m)
    
    Sys.sleep(0.01)
    flush.console()
  }
  
  #Store linkage quality metrics
  n=apply(LinkDesignation[,postburnin:endval],2,function(x) sum(x!=0))
  TP=apply(LinkDesignation[,postburnin:endval],2,function(x) sum(x==Data_A$MatchID_m))
  FP=n-TP
  Sensitivity[s]=mean(TP/n_m)
  Sensitivity_sd[s]=sd(TP/n_m)
  PPV[s]=mean(TP/n)
  PPV_sd[s]=sd(TP/n)
  F1[s]=mean(2*((TP/n_m)*(TP/n))/(TP/n_m+TP/n))
  F1_sd[s]=sd(2*((TP/n_m)*(TP/n))/(TP/n_m+TP/n))
  B1_bias[s]=mean(abs(beta_1[postburnin:endval]-cond_assoc*a_par))
  B1_bias_sd[s]=sd(abs(beta_1[postburnin:endval]-cond_assoc*a_par))
  B1_RMSE[s]=sqrt(sum((beta_1[postburnin:endval]-cond_assoc*a_par)^2)/100)
}

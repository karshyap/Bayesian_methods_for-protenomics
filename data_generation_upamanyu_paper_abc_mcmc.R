###Om ganeshaaya namaha##

#####Parameters used for generating the model###
##Control- class_0
##Treatement- class_1

rm(list=ls())

library(class)
library(boot)
library(caret)
library(ISLR)

clrscr<-function(){
  for(i in 1:100) {cat("\n")}
}

library(MASS) 
library(plyr)
library(ellipse)
library(ggplot2)
library(MASS)
Npro_trt=100
Npro_cntrl=100

Npro=Npro_trt

#Npro_trt_test=60
#Npro_cntrl_test=60
Npep_trt=ceiling(1.5*Npro_trt)
Npep_cntrl=ceiling(1.5*Npro_cntrl)
Npep=Npep_trt

#Npep_trt_test=90
#Npep_cntrl_test=90
t=2
theta_a=1000
number_of_samples=400 #number of samples in each class
number_of_samples_test_cntrl=400 #number of test samples fr control
number_of_samples_test_trt=400 #number of test samples fr treatement

amin=1.5
amax=1.6
al=1.55

fold_change=function(flag,amin,amax){   ##Please give amin and amax in just two decimal points.
  x1=runif(1,amin,amax)
  if(flag){
    return(x1)
  }
  else{
    return(1/x1)
  }
}

euc_norm <- function(x) sqrt(sum(x^2))



eta_treatement=c()
eta_control=c()
#mean_vec=mean()


M_cal=1000 #Number of calibarations to be made in abc rejection algo.

#Generating the synthetic sample data S0####
###Only control sample used#####

###initial parameters as used in table 2
k0=2
theta0=100
phi0=0.4
al0=1.55
gamma_dist0=rgamma(1000,shape = k0,scale = theta0)

gamma_l0=numeric(Npro_cntrl)


for(i in 1:Npro_cntrl){
  gamma_l0[i]=sample(gamma_dist0,1)
}

mean_vec_cont=gamma_l0


mean_vec_cont0=gamma_l0
sig_matrix0=matrix(numeric(Npro*Npro),nrow=Npro)
sig_sq_vec0=phi0*mean_vec_cont*mean_vec_cont

for(i in 1:Npro){
  sig_matrix0[i,i]=sig_sq_vec0[i]
}


c_pro_control0=mvrnorm(n=number_of_samples,mu = mean_vec_cont0,Sigma = sig_matrix0)
c_pro_treatement0=mvrnorm(n=number_of_samples,mu = mean_vec_cont,Sigma = sig_matrix0)
	
	



ratio_vec0=colMeans(c_pro_control0)/colMeans(c_pro_treatement0)



#####With initial value of al
fold_change_vec0=numeric(Npro)

for(i in 1:Npro){
  if(ratio_vec0[i]>1)
  {
    fold_change_vec0[i]=fold_change(1,al0,al0)
  }
  else{
    fold_change_vec0[i]=fold_change(0,al0,al0)
    
  }
}

mean_vec_treatement0=mean_vec_cont0*fold_change_vec0


###Defining c_pro_treatement_0 again for convinience

c_pro_treatement0=mvrnorm(n=number_of_samples,mu = mean_vec_treatement0,Sigma = sig_matrix0)

	




k_list=list()
theta_list=list()
phi_list=list()
count_1=0
###Synthetic sample data S0 done#####

###ABC-Rejection Sampling##########

norm_array=numeric(M_cal)

for(j in 1:M_cal){
  k=sample(160:240,1)/100
  theta=sample(800:1200,1)
  gamma_dist=rgamma(1000,shape=k,scale=theta_a)
  gamma_l=numeric(Npro_cntrl)
  
  
  
  for(i in 1:Npro_cntrl){
    gamma_l[i]=sample(gamma_dist,1)
  }
  
  mean_vec_cont=gamma_l
  
  phi=runif(Npro,0.3,0.5) ###Defining the coefficient of variation
  
  
  sig_matrix=matrix(numeric(Npro*Npro),nrow=Npro)
  sig_sq_vec=phi*mean_vec_cont*mean_vec_cont
  
  for(i in 1:Npro){
    sig_matrix[i,i]=sig_sq_vec[i]
  }
  
  
  c_pro_control=mvrnorm(n=number_of_samples,mu = mean_vec_cont,Sigma = sig_matrix)
  
  diff_vec=colMeans(c_pro_control)-colMeans(c_pro_control0)
  norm_array[j]=euc_norm(diff_vec)
 # print(euc_norm(diff_vec))
  thresh_key=20000
  #print(j)
  if(norm_array[j]<thresh_key){
    #print(diff_vec)
    #print("gnasher")
    count_1=count_1+1
    k_list[count_1]=k
    theta_list[count_1]=theta
    phi_list[count_1]=phi
  }
  #cat("\n")
    
  
  
  
  }
#print(mean(norm_array))
#cat("\n")

k_vec=theta_vec=phi_vec=numeric(length(k_list))
for(i in 1:length(k_vec)){
  k_vec[i]=k_list[[i]]
  theta_vec[i]=theta_list[[i]]
  phi_vec[i]=phi_list[[i]]
}


k_opt=mean(k_vec)
theta_opt=mean(theta_vec)
phi_opt=mean(phi_vec) 

#######Upamanyu's algorithm-3...The ABC-MCMC-algorithm#####
	
#####First three steps of ABC-MCMC#####

###gamma_0 is the gamma related to S^(0)_(0) in the paper. This is generated with k_opt and theta_opt
###gamma_dist_0 is the proper gamma distribution pertaining to S_(0) in the paper. Its generated with k_0 and theta_0.
	
#Step1: Sampling gamma and generating mean vectors##
gamma_dist_0=rgamma(1000,shape=k_opt,scale=theta_opt)
ratio_vec=colMeans(c_pro_treatement0)/colMeans(c_pro_control0)
fold_change_vec=ratio_vec
c_pro_control_0=c_pro_control0
print("hoi")
c_pro_treatment_0=c_pro_treatement0
gamma_0=numeric(Npro)

for(i in 1:Npro_cntrl){
  gamma_0[i]=sample(gamma_dist_0,1)
}

mean_vec_cntrl_0_0=gamma_0
mean_vec_trt_0_0=gamma_0*fold_change_vec

sig_matrix_0_0=matrix(numeric(Npro*Npro),nrow=Npro)
sig_sq_vec_0_0=phi_opt*mean_vec_cntrl_0_0*mean_vec_cntrl_0_0
for(i in 1:Npro){
    sig_matrix_0_0[i,i]=sig_sq_vec_0_0[i]
  }
  

###Step2: Generating protein data
c_pro_control_0_0=mvrnorm(n=number_of_samples,mu = mean_vec_cntrl_0_0,Sigma = sig_matrix_0_0)
c_pro_treatment_0_0=mvrnorm(n=number_of_samples,mu = mean_vec_trt_0_0,Sigma = sig_matrix_0_0)

diff_vec_cntrl=colMeans(c_pro_control_0_0)-colMeans(c_pro_control_0)
diff_vec_trt=colMeans(c_pro_treatment_0_0)-colMeans(c_pro_treatment_0)

norm_cntrl=euc_norm(diff_vec_cntrl)

norm_trt=euc_norm(diff_vec_trt)

while(norm_cntrl>20000 | norm_trt>25000){
  print("i am here")
  print(norm_cntrl)
  print(norm_trt)
  for(i in 1:Npro_cntrl){
    gamma_0[i]=sample(gamma_dist_0,1)
  }
  
  mean_vec_cntrl_0_0=gamma_0
  mean_vec_trt_0_0=gamma_0*fold_change_vec
  
  sig_matrix_0_0=matrix(numeric(Npro*Npro),nrow=Npro)
  sig_sq_vec_0_0=phi_opt*mean_vec_cntrl_0_0*mean_vec_cntrl_0_0
  for(i in 1:Npro){
    sig_matrix_0_0[i,i]=sig_sq_vec_0_0[i]
  }
  
  
  ###Step2: Generating protein data
  c_pro_control_0_0=mvrnorm(n=number_of_samples,mu = mean_vec_cntrl_0_0,Sigma = sig_matrix_0_0)
  c_pro_treatment_0_0=mvrnorm(n=number_of_samples,mu = mean_vec_trt_0_0,Sigma = sig_matrix_0_0)
  
  diff_vec_cntrl=colMeans(c_pro_control_0_0)-colMeans(c_pro_control_0)
  diff_vec_trt=colMeans(c_pro_treatment_0_0)-colMeans(c_pro_treatment_0)
  
  norm_cntrl=euc_norm(diff_vec_cntrl)
  
  norm_trt=euc_norm(diff_vec_trt)
  
}

#####Steps 5, 6 and 7 of the markov chain####
c_pro_control_markov=c_pro_control_0_0
c_pro_trt_markov=c_pro_treatment_0_0
gamma_vec=colMeans(c_pro_control_markov)
#####Main MCMC chain#####
count_markov=0
markov_iter=100
gamma_markov=matrix(0L, nrow = markov_iter,ncol=Npro)
for(j in 1:markov_iter){
  mean_vec_cntrl_markov=gamma_vec
  mean_vec_trt_markov=gamma_vec*fold_change_vec
  #print(i)
  sig_matrix_markov=matrix(numeric(Npro*Npro),nrow=Npro)
  sig_sq_vec_markov=phi_opt*mean_vec_cntrl_markov*mean_vec_cntrl_markov
  for(i in 1:Npro){
    sig_matrix_markov[i,i]=sig_sq_vec_markov[i]
  }
  
  c_pro_control_markov=mvrnorm(n=number_of_samples,mu = mean_vec_cntrl_markov,Sigma = sig_matrix_markov)
  c_pro_trt_markov=mvrnorm(n=number_of_samples,mu = mean_vec_trt_markov,Sigma = sig_matrix_markov)
  
  diff_vec_cntrl_markov=colMeans(c_pro_control_markov)-colMeans(c_pro_control_0)
  diff_vec_trt_markov=colMeans(c_pro_trt_markov)-colMeans(c_pro_treatment_0)
  
  norm_cntrl_markov=euc_norm(diff_vec_cntrl_markov)
  
  norm_trt_markov=euc_norm(diff_vec_trt_markov)
  print("-----------")
  print(norm_cntrl_markov)
  print(norm_trt_markov)
  print("-----------")
  
  if(norm_cntrl_markov<20000 & norm_trt_markov<25000){
  
  gamma_vec=colMeans(c_pro_control_markov)
  #print(j)
  count_markov=count_markov+1
  
  }
  gamma_markov[j,]=gamma_vec


  
}
######Algorithm-3 in upamanyu paper done#####



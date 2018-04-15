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
Npro_trt=35
Npro_cntrl=35

#Npro_trt_test=60
#Npro_cntrl_test=60
Npep_trt=ceiling(1.5*Npro_trt)
Npep_cntrl=ceiling(1.5*Npro_cntrl)
Npep=Npep_trt

#Npep_trt_test=90
#Npep_cntrl_test=90
t=2
theta_a=1000
number_of_samples=50 #number of samples in each class
number_of_samples_test_cntrl=400 #number of test samples fr control
number_of_samples_test_trt=400 #number of test samples fr treatement

amin=1.5
amax=1.6

k0=2
theta0=100
phi0=0.4
M_cal=1000 #Number of calibarations to be made in abc rejection algo.


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

main_function<-function(Npro_cntrl,Npro_trt,Npep_cntrl,Npep_trt,t,theta_a,number_of_samples,number_of_samples_test_cntrl,
                        number_of_samples_test_trt,amin,amax,k0,theta0,phi0,M_cal){
  Npro=Npro_trt
  Npep=Npep_trt
  
  eta_treatement=c()
  eta_control=c()
  #mean_vec=mean()
  
  
  
  #Generating the synthetic sample data S0####
  ###Only control sample used#####
  
  ###initial parameters as used in table 2
  
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
  
  
  c_pro_control0=mvrnorm(n=number_of_samples,mu = mean_vec_cont,Sigma = sig_matrix0)
  
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
    
    phi=runif(Npro,phi0-0.1,phi0+1) ###Defining the coefficient of variation
    
    
    sig_matrix=matrix(numeric(Npro*Npro),nrow=Npro)
    sig_sq_vec=phi*mean_vec_cont*mean_vec_cont
    
    for(i in 1:Npro){
      sig_matrix[i,i]=sig_sq_vec[i]
    }
    
    
    c_pro_control=mvrnorm(n=number_of_samples,mu = mean_vec_cont,Sigma = sig_matrix)
    
    diff_vec=colMeans(c_pro_control)-colMeans(c_pro_control0)
    norm_array[j]=euc_norm(diff_vec)
    # #print(euc_norm(diff_vec))
    print(diff_vec)
    thresh_key=2*(median(norm_array)+min(norm_array))/3
    
    if(norm_array[j]<thresh_key){
      ##print(diff_vec)
      
      count_1=count_1+1
      k_list[count_1]=k
      theta_list[count_1]=theta
      phi_list[count_1]=phi
    }
    #cat("\n")
    
    
    
    
  }
  ##print(mean(norm_array))
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
  
  
  #############ABC Rejection sampling done#######
  
  ##### Data from optimal parameters being generated#######
  
  
  
  gamma_dist_opt=rgamma(1000,shape=k_opt,scale=theta_opt)
  gamma_l_opt=numeric(Npro_cntrl)
  
  for(i in 1:Npro_cntrl){
    gamma_l_opt[i]=sample(gamma_dist,1)
  }
  
  mean_vec_cont_opt=gamma_l_opt
  
  sig_matrix_opt=matrix(numeric(Npro*Npro),nrow=Npro)
  sig_sq_vec_opt=phi_opt*mean_vec_cont_opt*mean_vec_cont_opt
  
  for(i in 1:Npro){
    sig_matrix_opt[i,i]=sig_sq_vec_opt[i]
  }
  
  
  #####Generating dummy data till i understand the protein file###
  flag_vec=numeric(Npro)
  
  for(i in 1:Npro)
  {
    flag_vec[i]=sample(0:1,1)
  }
  
  ####flag_vec is used to say if the protein is overexpressed or not #####
  
  mean_vec_treatement1=gamma_l_opt
  
  fold_change_vec=numeric(Npro)
  for(i in 1:Npro){
    fold_change_vec[i]=fold_change(flag_vec[i],amin,amax)
  }
  
  
  
  mean_vec_treatement=mean_vec_treatement1*fold_change_vec
  
  
  
  c_pro_control=mvrnorm(n=number_of_samples,mu = mean_vec_cont_opt,Sigma = sig_matrix_opt)
  c_pro_treatement=mvrnorm(n=number_of_samples,mu = mean_vec_treatement,Sigma = sig_matrix_opt)
  
  c_pro_control_test=mvrnorm(n=number_of_samples_test_cntrl,mu = mean_vec_cont_opt,Sigma = sig_matrix_opt)
  c_pro_treatement_test=mvrnorm(n=number_of_samples_test_trt,mu = mean_vec_treatement,Sigma = sig_matrix_opt)
  
  
  #c_pep_
  #c_pep_control=
  #c_pep_treatement=
  
  c_pep_control=matrix(numeric(Npep_cntrl*number_of_samples),nrow = number_of_samples)
  c_pep_treatement=matrix(numeric(Npep_trt*number_of_samples),nrow = number_of_samples)
  
  c_pep_control_test=matrix(numeric(Npep_cntrl*number_of_samples_test_cntrl),nrow = number_of_samples_test_cntrl)
  c_pep_treatement_test=matrix(numeric(Npep_trt*number_of_samples_test_trt),nrow = number_of_samples_test_trt)
  
  prot_pept_list_cntrl=list()
  prot_pept_list_trt=list()
  
  prot_pept_list_cntrl_test=list()
  prot_pept_list_trt_test=list()
  
  ####Generating random peptide control data###
  for(j in 1:length(c_pep_control[1,])){
    p3=sample(1:length(c_pro_control[1,]),sample(1:3,1))
    prot_pept_list_cntrl[[j]]=p3
    for(i in 1:length(c_pep_control[,1])){
      
      #cat(i,"-",j,"-",p3,"\n")
      for(k in 1:length(p3)){
        c_pep_control[i,j]=c_pep_control[i,j]+c_pro_control[i,p3[k]]
      }
    }
    if(1){
      for(i in 1:length(c_pep_control_test[,1])){
        for(k in 1:length(p3)){
          c_pep_control_test[i,j]=c_pep_control_test[i,j]+c_pro_control_test[i,p3[k]]
        }
      }
    }
    
  }
  
  if(1){
    for(j in 1:length(c_pep_control_test[1,])){
      p3=sample(1:length(c_pro_control_test[1,]),sample(1:3,1))
      prot_pept_list_cntrl_test[[j]]=p3
      for(i in 1:length(c_pep_control_test[,1])){
        
        #cat(i,"-",j,"-",p3,"\n")
        for(k in 1:length(p3)){
          c_pep_control_test[i,j]=c_pep_control_test[i,j]+c_pro_control_test[i,p3[k]]
        }
      }
      
      
      
    }
  }
  
  ##print("grain")
  ####Generating random peptide treatement data###
  
  for(j in 1:length(c_pep_treatement[1,])){
    p3=sample(1:length(c_pro_treatement[1,]),sample(1:3,1))
    prot_pept_list_trt[[j]]=p3
    for(i in 1:length(c_pep_treatement[,1])){
      
      #cat(i,"-",j,"-",p3,"\n")
      for(k in 1:length(p3)){
        c_pep_treatement[i,j]=c_pep_treatement[i,j]+c_pro_treatement[i,p3[k]]
      }
    }
    if(1){
      for(i in 1:length(c_pep_treatement_test[,1])){
        
        #cat(i,"-",j,"-",p3,"\n")
        for(k in 1:length(p3)){
          c_pep_treatement_test[i,j]=c_pep_treatement_test[i,j]+c_pro_treatement_test[i,p3[k]]
        }
      }
    }
  }
  
  #print("hello")
  if(1){
    for(j in 1:length(c_pep_treatement_test[1,])){
      p3=sample(1:length(c_pro_treatement_test[1,]),sample(1:3,1))
      prot_pept_list_trt_test[[j]]=p3
      for(i in 1:length(c_pep_treatement_test[,1])){
        
        #cat(i,"-",j,"-",p3,"\n")
        for(k in 1:length(p3)){
          c_pep_treatement_test[i,j]=c_pep_treatement_test[i,j]+c_pro_treatement_test[i,p3[k]]
        }
      }
      
      
      
    }
    
  }
  
  prot_pept_list_cntrl_test=prot_pept_list_cntrl
  prot_pept_list_trt_test=prot_pept_list_trt
  
  ####Parameters to find mu_ij
  
  mu_matrix_cntrl=matrix(numeric(Npep_cntrl*number_of_samples),nrow = number_of_samples)
  mu_matrix_trt=matrix(numeric(Npep_trt*number_of_samples),nrow = number_of_samples)
  mu_matrix_cntrl_test=matrix(numeric(Npep_cntrl*number_of_samples_test_cntrl),nrow = number_of_samples_test_cntrl)
  mu_matrix_trt_test=matrix(numeric(Npep_trt*number_of_samples_test_trt),nrow = number_of_samples_test_trt)
  
  kappa=5
  efficiency_vector=runif(Npep,0.1,1)
  print("maxhappan")
  for(i in 1:length(mu_matrix_cntrl[,1])){
    for(j in 1:length(mu_matrix_cntrl[1,])){
      mu_matrix_cntrl[i,j]=c_pep_control[i,j]*kappa*efficiency_vector[j]
    }
  }
  
  for(i in 1:length(mu_matrix_cntrl_test[,1])){
    for(j in 1:length(mu_matrix_cntrl_test[1,])){
      mu_matrix_cntrl_test[i,j]=c_pep_control_test[i,j]*kappa*efficiency_vector[j]
    }
  }
  
  
  
  for(i in 1:length(mu_matrix_trt[,1])){
    for(j in 1:length(mu_matrix_trt[1,])){
      mu_matrix_trt[i,j]=c_pep_treatement[i,j]*kappa*efficiency_vector[j]
    }
  }
  
  
  for(i in 1:length(mu_matrix_trt_test[,1])){
    for(j in 1:length(mu_matrix_trt_test[1,])){
      mu_matrix_trt_test[i,j]=c_pep_treatement_test[i,j]*kappa*efficiency_vector[j]
    }
  }
  
  ###Generate noisy gaussian###
  alpha=0.03
  beta=3.6
  
  var_vector_noisy_gaussian_cntrl=alpha*(mu_matrix_cntrl^2)+beta*mu_matrix_cntrl
  var_vector_noisy_gaussian_trt=alpha*(mu_matrix_trt^2)+beta*mu_matrix_trt
  
  var_vector_noisy_gaussian_cntrl_test=alpha*(mu_matrix_cntrl_test^2)+beta*mu_matrix_cntrl_test
  var_vector_noisy_gaussian_trt_test=alpha*(mu_matrix_trt_test^2)+beta*mu_matrix_trt_test
  
  
  total_vector_cntrl=matrix(numeric(Npep_cntrl*number_of_samples),nrow = number_of_samples)
  total_vector_trt=matrix(numeric(Npep_cntrl*number_of_samples),nrow = number_of_samples)
  
  
  total_vector_cntrl_test=matrix(numeric(Npep_cntrl*number_of_samples_test_cntrl),nrow = number_of_samples_test_cntrl)
  total_vector_trt_test=matrix(numeric(Npep_cntrl*number_of_samples_test_trt),nrow = number_of_samples_test_trt)
  
  
  for(i in 1:length(total_vector_cntrl[,1])){
    for(j in 1:length(total_vector_cntrl[1,])){
      #cat(i,j,"\n")
      total_vector_cntrl[i,j]=mu_matrix_cntrl[i,j]+20*mvrnorm(n=1,mu = 0, Sigma = abs(var_vector_noisy_gaussian_cntrl[i,j]))
    }
  }
  
  for(i in 1:length(total_vector_trt[,1])){
    for(j in 1:length(total_vector_trt[1,])){
      #cat(i,j,"\n")
      total_vector_trt[i,j]=mu_matrix_trt[i,j]+20*mvrnorm(n=1,mu = 0, Sigma = abs(var_vector_noisy_gaussian_trt[i,j]))
    }
  }
  
  
  for(i in 1:length(total_vector_cntrl_test[,1])){
    for(j in 1:length(total_vector_cntrl_test[1,])){
      #cat(i,j,"\n")
      total_vector_cntrl_test[i,j]=mu_matrix_cntrl_test[i,j]+20*mvrnorm(n=1,mu = 0, Sigma = abs(var_vector_noisy_gaussian_cntrl_test[i,j]))
    }
  }
  
  for(i in 1:length(total_vector_trt_test[,1])){
    for(j in 1:length(total_vector_trt_test[1,])){
      #cat(i,j,"\n")
      total_vector_trt_test[i,j]=mu_matrix_trt_test[i,j]+20*mvrnorm(n=1,mu = 0, Sigma = abs(var_vector_noisy_gaussian_trt_test[i,j]))
    }
  }
  #total_vector_cntrl=data.frame(cbind(total_vector_cntrl,numeric(50)))
  #total_vector_trt=data.frame(cbind(total_vector_trt,numeric(50)+1))
  
  #total_vector_cntrl$X31=factor(total_vector_cntrl$X31)
  #total_vector_trt$X31=factor(total_vector_trt$X31)
  
  
  ###To calculate rolled up abundances###
  #print("greek4")
  
  x_pro_control=matrix(numeric(Npro_cntrl*number_of_samples),nrow = number_of_samples)
  x_pro_trt=matrix(numeric(Npro_trt*number_of_samples),nrow = number_of_samples)
  
  x_pro_control_test=matrix(numeric(Npro_cntrl*number_of_samples_test_cntrl),nrow = number_of_samples_test_cntrl)
  x_pro_trt_test=matrix(numeric(Npro_trt*number_of_samples_test_trt),nrow = number_of_samples_test_trt)
  
  
  pept_to_prot_cntrl=list()
  pept_to_prot_trt=list()
  
  pept_to_prot_cntrl_test=list()
  pept_to_prot_trt_test=list()
  
  
  for(i in 1:length(prot_pept_list_cntrl)){
    for(j in 1:length(prot_pept_list_cntrl[i][[1]])){
      ##print(prot_pept_list_cntrl[i][[1]][j])
      k1=prot_pept_list_cntrl[i][[1]][j]
      ##print(k1)
      pept_to_prot_cntrl[k1][[1]][length(pept_to_prot_cntrl[k1][[1]])+1]=i
      
    }
  }
  
  for(i in 1:length(prot_pept_list_trt)){
    for(j in 1:length(prot_pept_list_trt[i][[1]])){
      ##print(prot_pept_list_cntrl[i][[1]][j])
      k1=prot_pept_list_trt[i][[1]][j]
      ##print(k1)
      pept_to_prot_trt[k1][[1]][length(pept_to_prot_trt[k1][[1]])+1]=i
      
    }
  }
  
  #print("greek1")
  
  for(i in 1:length(prot_pept_list_cntrl_test)){
    for(j in 1:length(prot_pept_list_cntrl_test[i][[1]])){
      ##print(prot_pept_list_cntrl[i][[1]][j])
      k1=prot_pept_list_cntrl_test[i][[1]][j]
      ##print(k1)
      pept_to_prot_cntrl_test[k1][[1]][length(pept_to_prot_cntrl_test[k1][[1]])+1]=i
      
    }
  }
  
  #print("greek2")
  
  
  for(i in 1:length(prot_pept_list_trt_test)){
    for(j in 1:length(prot_pept_list_trt_test[i][[1]])){
      ##print(prot_pept_list_cntrl[i][[1]][j])
      k1=prot_pept_list_trt_test[i][[1]][j]
      ##print(k1)
      pept_to_prot_trt_test[k1][[1]][length(pept_to_prot_trt_test[k1][[1]])+1]=i
      ##print("gopo")
    }
  }
  
  ##print("greek3")
  
  
  xlj_control=matrix(numeric(Npro_cntrl*number_of_samples),nrow = number_of_samples)
  xlj_trt=matrix(numeric(Npro_trt*number_of_samples),nrow = number_of_samples)
  
  xlj_control_test=matrix(numeric(Npro_cntrl*number_of_samples_test_cntrl),nrow = number_of_samples_test_cntrl)
  xlj_trt_test=matrix(numeric(Npro_trt*number_of_samples_test_trt),nrow = number_of_samples_test_trt)
  
  ###############################
  
  for(i in 1:length(pept_to_prot_cntrl)){
    if(pept_to_prot_cntrl[2]=="NULL"){
      pept_to_prot_cntrl[2]=pept_to_prot_cntrl[3]
    }
    
    if(pept_to_prot_cntrl[1]=="NULL"){
      pept_to_prot_cntrl[1]=pept_to_prot_cntrl[2]
    }
    
    
    if(pept_to_prot_cntrl[i]=="NULL"){
      pept_to_prot_cntrl[i]=pept_to_prot_cntrl[1]
    }
    
    
    
    
  }
  
  ############################
  
  for(i in 1:length(pept_to_prot_trt)){
    if(pept_to_prot_trt[2]=="NULL"){
      pept_to_prot_trt[2]=pept_to_prot_trt[3]
    }
    
    if(pept_to_prot_trt[1]=="NULL"){
      pept_to_prot_trt[1]=pept_to_prot_trt[2]
    }
    
    
    if(pept_to_prot_trt[i]=="NULL"){
      pept_to_prot_trt[i]=pept_to_prot_trt[1]
    }
    
    
    
    
  }
  
  
  ##print("greek1")
  
  
  for(i in 1:length(pept_to_prot_cntrl_test)){
    if(pept_to_prot_cntrl_test[2]=="NULL"){
      pept_to_prot_cntrl_test[2]=pept_to_prot_cntrl_test[3]
    }
    
    if(pept_to_prot_cntrl_test[1]=="NULL"){
      pept_to_prot_cntrl_test[1]=pept_to_prot_cntrl_test[2]
    }
    
    
    if(pept_to_prot_cntrl_test[i]=="NULL"){
      pept_to_prot_cntrl_test[i]=pept_to_prot_cntrl_test[1]
    }
    
    
    
    
  }
  
  
  for(i in 1:length(pept_to_prot_trt_test)){
    if(pept_to_prot_trt_test[2]=="NULL"){
      pept_to_prot_trt_test[2]=pept_to_prot_trt_test[3]
    }
    
    if(pept_to_prot_trt_test[1]=="NULL"){
      pept_to_prot_trt_test[1]=pept_to_prot_trt_test[2]
    }
    
    
    if(pept_to_prot_trt_test[i]=="NULL"){
      pept_to_prot_trt_test[i]=pept_to_prot_trt_test[1]
    }
    
    
    
    
  }
  
  
  for(j in 1:length(xlj_control[1,])){ #for all proteins
    for(i in 1:length(xlj_control[,1])) { #for all samples
      for(k in 1:length(pept_to_prot_cntrl[j][[1]])){
        p4=pept_to_prot_cntrl[j][[1]][k]
        if(is.null(p4)){
          p4=pept_to_prot_cntrl[1][[1]][1]
        }
        xlj_control[i,j]=xlj_control[i,j]+total_vector_cntrl[i,p4]
        
      }
    }
    ##print(j)
    
  }
  
  #print("hello123")
  
  for(j in 1:length(xlj_trt[1,])){ #for all proteins
    for(i in 1:length(xlj_trt[,1])) { #for all samples
      for(k in 1:length(pept_to_prot_trt[j][[1]])){
        p4=pept_to_prot_trt[j][[1]][k]
        if(is.null(p4)){
          #print("gonu")
          p4=pept_to_prot_trt[1][[1]][1]
          
        }
        #if(is.null(p4)){
        #}
        xlj_trt[i,j]=xlj_trt[i,j]+total_vector_trt[i,p4]
        
      }
    }
    ##print(j)
    
  }
  ##print("greekza")
  ########
  for(j in 1:length(xlj_control_test[1,])){ #for all proteins
    for(i in 1:length(xlj_control_test[,1])) { #for all samples
      for(k in 1:length(pept_to_prot_cntrl_test[j][[1]])){
        p4=pept_to_prot_cntrl_test[j][[1]][k]
        if(is.null(p4)){
          p4=pept_to_prot_cntrl_test[1][[1]][1]
        }
        
        xlj_control_test[i,j]=xlj_control_test[i,j]+total_vector_cntrl_test[i,p4]
        
      }
    }
    ##print(j)
    
  }
  
  
  for(j in 1:length(xlj_trt_test[1,])){ #for all proteins
    for(i in 1:length(xlj_trt_test[,1])) { #for all samples
      for(k in 1:length(pept_to_prot_trt_test[j][[1]])){
        p4=pept_to_prot_trt_test[j][[1]][k]
        if(is.null(p4)){
          p4=pept_to_prot_trt_test[1][[1]][1]
        }
        
        xlj_trt_test[i,j]=xlj_trt_test[i,j]+total_vector_trt_test[i,p4]
        
      }
    }
    ##print(j)
    
  }
  
  
  
  #print("hi man")
  xlj_control=data.frame(xlj_control)
  xlj_control=data.frame(cbind(xlj_control,numeric(number_of_samples)))
  colnames(xlj_control)[Npro+1]="type"
  #xlj_control$type=as.factor(xlj_control$type)
  
  xlj_control_test=data.frame(xlj_control_test)
  xlj_control_test=data.frame(cbind(xlj_control_test,numeric(number_of_samples_test_cntrl)))
  colnames(xlj_control_test)[Npro+1]="type"
  
  xlj_trt=data.frame(xlj_trt)
  xlj_trt=data.frame(cbind(xlj_trt,rep(1,number_of_samples)))
  colnames(xlj_trt)[Npro+1]="type"
  
  
  xlj_trt_test=data.frame(xlj_trt_test)
  xlj_trt_test=data.frame(cbind(xlj_trt_test,rep(1,number_of_samples_test_trt)))
  colnames(xlj_trt_test)[Npro+1]="type"
  #xlj_trt$type=as.factor(xlj_trt$type)
  
  
  
  
  
  #xlj_control$"type"=numeric(50)
  #xlj_trt$"type"=(numeric(number_of_samples)+1)
  
  if(1){
    
    xlj_control=data.frame(xlj_control)
    xlj_trt=data.frame(xlj_trt)
    
    xlj_control_test=data.frame(xlj_control_test)
    xlj_trt_test=data.frame(xlj_trt_test)
    
    
    
    
    xlj_all=rbind(xlj_control,xlj_trt)
    xlj_all$type=as.factor(xlj_all$type)
    
    
    xlj_all_test=rbind(xlj_control_test,xlj_trt_test)
    xlj_all_test$type=as.factor(xlj_all_test$type)
    
    
  }
  #print("greek-pal")
  ####xlj_all is the required data as per equation 12###
  ####xlj_all_test is the required test data####
  
  ####Designing classifiers####
  
  #1. LDA classifier
  
  lda_classifier_50<-lda(type ~ .,data=xlj_all)
  predictions_lda_50=predict(lda_classifier_50,xlj_all_test[,1:Npro])$class
  table_data_lda_50=table(predictions_lda_50,xlj_all_test[,Npro+1])
  
  predictions_lda_50_app=predict(lda_classifier_50,xlj_all[,1:Npro])$class
  table_data_lda_50_app=table(predictions_lda_50_app,xlj_all[,Npro+1])
  
  
  training_labels=xlj_all$type
  knn_trained<-knn(train = xlj_all[,1:Npro]  , test =xlj_all_test[,1:Npro] , cl = training_labels, k=3)
  table_data_knn_50=table(knn_trained,xlj_all_test[,Npro+1])
  knn_trained_app<-knn(train = xlj_all[,1:Npro]  , test =xlj_all[,1:Npro] , cl = training_labels, k=3)
  table_data_knn_50_app=table(knn_trained_app,xlj_all[,Npro+1])
  
  lda_error=(table_data_lda_50[2]+table_data_lda_50[3])/sum(table_data_lda_50)
  knn_error=(table_data_knn_50[2]+table_data_knn_50[3])/sum(table_data_knn_50)
  
  #return(knn_error)
  #cat("lda-",lda_error)
  #cat("knn-",knn_error)
  #print(list(lda_error,knn_error))
  return(c(lda_error,knn_error))
  
  
}



n_vector=c(10,20,30,40,50)
Npro_vec=c(3,5,8,10)
phi_vec=c(0.2,0.4,0.6,0.8,1)
alpha_vec=c(0,0.2,0.4,0.6,0.8,1)

n_err_vec=matrix(rep(1,2*length(n_vector)),ncol = 2)
Npro_err_vec=matrix(rep(1,2*length(Npro_vec)),ncol = 2)
phi_err_vec=matrix(rep(1,2*length(phi_vec)),ncol = 2)
alpha_err_vec=matrix(rep(1,2*length(alpha_vec)),ncol = 2)




###Measuring the effect of error for lda and knn as a function of sample size
if(1){
  for(i in 1:length(n_err_vec[,1])) 
  {
    print(i)
    n_err_vec[i,]=main_function(Npro_cntrl,Npro_trt,Npep_cntrl,Npep_trt,t,theta_a,n_vector[i],number_of_samples_test_cntrl,
                                number_of_samples_test_trt,amin,amax,k0,theta0,phi0,M_cal)
  }
  #print("hi man")
  n_err_vec=data.frame(n_err_vec)
  n_err_vec=cbind(n_err_vec,n_vector)
  colnames(n_err_vec)=c("lda","knn","n")
  plot_n=ggplot(data = n_err_vec, aes(x=n, y=lda))
  plot_n=plot_n+geom_line(size=2.5,color="pink")+geom_point(size=2.5,color="blue")
}



####Measuring the effect of error for lda as a fuction of dimension
if(0){
  for(i in 1:length(Npro_vec)) 
  {
    Npro_err_vec[i,]=main_function(Npro_cntrl = Npro_vec[i],Npro_trt = Npro_vec[i],Npep_cntrl,Npep_trt,t,theta_a,30,number_of_samples_test_cntrl,
                                   number_of_samples_test_trt,amin,amax,k0,theta0,phi0,M_cal)
  }
  Npro_err_vec=data.frame(Npro_err_vec)
  Npro_err_vec=cbind(Npro_err_vec,Npro_vec)
  colnames(Npro_err_vec)=c("lda","knn","Npro")
  plot_Npro=ggplot(data = Npro_err_vec, aes(x=Npro, y=lda))
  plot_Npro=plot_Npro+geom_line(size=2.5,color="pink")+geom_point(size=2.5,color="blue")
}

####Measuring the effect of error for lda as a fuction of phi
if(0){
  for(i in 1:length(phi_err_vec[,1])) 
  {
    phi_err_vec[i,]=main_function(Npro_cntrl,Npro_trt,Npep_cntrl,Npep_trt,t,theta_a,30,number_of_samples_test_cntrl,
                                  number_of_samples_test_trt,amin,amax,k0,theta0,phi_err_vec[i],M_cal)
  }
  phi_err_vec=data.frame(phi_err_vec)
  phi_err_vec=cbind(phi_err_vec,phi_vec)
  colnames(phi_err_vec)=c("lda","knn","phi")
  plot_phi=ggplot(data = phi_err_vec, aes(x=phi, y=knn))
  plot_phi=plot_phi+geom_line(size=2.5,color="pink")+geom_point(size=2.5,color="blue")
}


#print("its-done")

####Measuring the effect of error for lda as a fuction of peptide efficiency alpha


#total_vector_cntrl=data.frame(cbind(total_vector_cntrl,numeric(50)))
#total_vector_trt=data.frame(cbind(total_vector_trt,numeric(50)+1))

#total_vector_cntrl$X31=factor(total_vector_cntrl$X31)
#total_vector_trt$X31=factor(total_vector_trt$X31)


####Add the mean and the gaussian vector####



#SNR=1/(alpha+(beta/alpha))











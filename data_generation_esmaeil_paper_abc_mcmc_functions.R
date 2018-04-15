###Om ganeshaaya namaha##

#####Parameters used for generating the model###
##Control- class_0
##Treatement- class_1
#rm(list=ls())
#rm(.Random.seed, envir=globalenv())
#time_init=proc.time()
print("New-Code")
load("variables.RData")
set.seed( as.integer((as.double(Sys.time())*1000+Sys.getpid()) %% 2^31) )

#print(runif(n=100,min=0,max=100))
#####Please check with random seed###
#rm(list=ls())
#set.seed(101)
library(class)
library(boot)
library(MASS) 

#library(caret)
#library(ISLR)
library("Biostrings")
setwd("/home/kashyap/Desktop/Masters_thesis_related/Codes")
time2=proc.time()
normalize=function(x){return(x/max(x))}



t_test_vals=function(c_pro_all)
{
  t_test_value_pro=numeric((length(c_pro_all[1,])-1))
  len_t_test_pro=length(t_test_value_pro)
  
  for(i in 1:length(t_test_value_pro)){
    k1=t.test(c_pro_all[,i]~c_pro_all[,len_t_test_pro+1])
    t_test_value_pro[i]=abs(k1$statistic)
    
  }
  
  return(t_test_value_pro)
}

clrscr<-function(){
  for(i in 1:100) {cat("\n")}
}


diff_elem_removal=function(x,y){
  ##X is an array from which the elements which are there in y should be removed
  remov_array=rep(0,(length(x)+100))
  count=0
  for(i in 1:length(x)){
    match_flag=0
    for(j in 1:length(y)){
      if(x[i]==y[j]){
        match_flag=1
      }
    }
    if(match_flag){
      count=count+1
      remov_array[count]=i
    }
  }
  remov_array=remov_array[remov_array!=0]
  
  return(x[-remov_array])
}





diff_elem_removal_total=function(p1,q1)
{
  #Here p1 is the protein to peptide list vector
  #q1 is the proteins selected by t-test
  r1=list()
  count=0
  for(i in 1:length(p1)){
    k1=p1[[i]]
    v1=setdiff(k1,q1)
    #print(k1)
    #print(v1)
    if(length(v1)){
      l1=diff_elem_removal(k1,v1)
      
      #print(l1)
      #print("----------------")
      if(length(l1)>0){
        #if(1){  
        count=count+1
        
        r1[[count]]=l1
      }
    }
    else{
      count=count+1
      r1[[count]]=k1
    }
    
  }
  return(r1)
}

#####This function is to reduce the prot_pept_list from random proteins to the ascending stuff ######
prot_pept_list_reduction=function(r1){
  max_elem=0
  count_2=0
  prot_uniq_array=rep(0,100000)
  for(i in 1:length(r1)){
    for(j in 1:length(r1[[i]])){
      count_2=count_2+1
      prot_uniq_array[count_2]=r1[[i]][j]
    }
  }
  r2=r1
  prot_uniq_array1=prot_uniq_array[prot_uniq_array>0]
  prot_uniq_array1=unique(prot_uniq_array1)
  uniq_prots=1:length(prot_uniq_array1)
  
  for(i in 1:length(r1)){
    for(j in 1:length(r1[[i]])){
      for(k in 1:length(uniq_prots)){
        if(r1[[i]][j]==prot_uniq_array1[k]){
          r2[[i]][j]=uniq_prots[k]
        }
      }
    }
  }
  return(r2)
}

#####Reduction of prot_pept_list ends######

######Function for the protein to peptide list starts####


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

if(1){
  pro_trt_file="prot_file1.fasta"
  pep_trt_file="pept_file1.fasta"
}
pro_cntrl_file="prot_file1.fasta"
pep_cntrl_file="pept_file1.fasta"

if(1){
  
  Npro_factor=0.02
  #Npro_trt=fasta_file_to_array_conv(pro_trt_file,pep_trt_file,"pro")
  
  #Npro_cntrl=fasta_file_to_array_conv(pro_cntrl_file,pep_cntrl_file,"pro")
  
  Npro_for_analysis=ceiling(Npro_factor*Npro_trt)
  
  #Npep_trt=fasta_file_to_array_conv(pro_trt_file,pep_trt_file,"pep")
  #Npep_cntrl=fasta_file_to_array_conv(pro_cntrl_file,pep_cntrl_file,"pep")
  
  
}
	
noise_factor_gauss=1
noise_factor_exp=1
Npro=Npro_trt

#Npro_trt_test=60
#Npro_cntrl_test=60
#Npep_trt=ceiling(15*Npro_trt)
#Npep_cntrl=ceiling(15*Npro_cntrl)
Npep=Npep_trt	
	
####Function for the generation of the data in the paper. This is like main function #####

peptide_generation=function(gamma_dist_c_opt,gamma_dist_a_opt,fold_change_vec,phi_opt,Npro_cntrl,perc_cntrl,Npep_cntrl,Npro_trt,perc_trt,Npep_trt,prot_pept_list_cntrl,Npro_for_analysis,noise_factor_gauss,
                            noise_factor_exp,prot_pept_list_cntrl_test,prot_pept_list_trt,prot_pept_list_trt_test,train_test_flag)
{

#gamma_dist_a_opt=rgamma(1000,shape=ka_opt,scale=thetaa_opt)
#gamma_dist_c_opt=rgamma(1000,shape=kc_opt,scale=thetac_opt)
  Npro_cntrl_a=ceiling(Npro_cntrl*frac_npro)
  Npro_cntrl_c=Npro_cntrl-Npro_cntrl_a
  
  Npro_trt_a=ceiling(Npro_trt*frac_npro)
  Npro_trt_c=Npro_trt-Npro_trt_a
  Npro=Npro_cntrl

gamma_l_opt=numeric(Npro_cntrl)

for(i in 1:(Npro_cntrl_c+Npro_cntrl_a)){
  if(i<=Npro_cntrl_c){
    gamma_l_opt[i]=sample(gamma_dist_c_opt,1)
  }
  else{
    gamma_l_opt[i]=sample(gamma_dist_a_opt,1)
  }
}

mean_vec_cont_opt=gamma_l_opt
	
sig_matrix_opt=matrix(numeric(Npro*Npro),nrow=Npro)
sig_sq_vec_opt=phi_opt*mean_vec_cont_opt*mean_vec_cont_opt

for(i in 1:Npro_cntrl){
  sig_matrix_opt[i,i]=sig_sq_vec_opt[i]
}

mean_vec_treatement_opt=mean_vec_cont_opt*fold_change_vec

c_pro_control=mvrnorm(n=number_of_samples_train_cntrl,mu = mean_vec_cont_opt,Sigma = sig_matrix_opt)
#print("gundanna")
c_pro_treatement=mvrnorm(n=number_of_samples_train_trt,mu = mean_vec_treatement_opt,Sigma = sig_matrix_opt)
#print("hindanna")
	
c_pro_control_test=mvrnorm(n=number_of_samples_test_cntrl,mu = mean_vec_cont_opt,Sigma = sig_matrix_opt)
c_pro_treatement_test=mvrnorm(n=number_of_samples_test_trt,mu = mean_vec_treatement_opt,Sigma = sig_matrix_opt)

	
c_pro_control[,Npro_cntrl_c:length(c_pro_control[1,])]=runif(1,min=0,max=10e-6)*c_pro_control[,Npro_cntrl_c:length(c_pro_control[1,])]
c_pro_treatement[,Npro_trt_c:length(c_pro_treatement[1,])]=runif(1,min=0,max=10e-6)*c_pro_treatement[,Npro_trt_c:length(c_pro_treatement[1,])]


c_pro_control_test[,Npro_cntrl_c:length(c_pro_control_test[1,])]=runif(1,min=0,max=10e-6)*c_pro_control_test[,Npro_cntrl_c:length(c_pro_control_test[1,])]
c_pro_treatement_test[,Npro_trt_c:length(c_pro_treatement_test[1,])]=runif(1,min=0,max=10e-6)*c_pro_treatement_test[,Npro_trt_c:length(c_pro_treatement_test[1,])]

c_pep_control=matrix(numeric(Npep_cntrl*number_of_samples_train_cntrl),nrow = number_of_samples_train_cntrl)
c_pep_treatement=matrix(numeric(Npep_trt*number_of_samples_train_trt),nrow = number_of_samples_train_trt)

c_pep_control_test=matrix(numeric(Npep_cntrl*number_of_samples_test_cntrl),nrow = number_of_samples_test_cntrl)
c_pep_treatement_test=matrix(numeric(Npep_trt*number_of_samples_test_trt),nrow = number_of_samples_test_trt)

for(j in 1:length(c_pep_control[1,])){
    #p3=sample(1:length(c_pro_control[1,]),sample(1:3,1))
    #prot_pept_list_cntrl[[j]]=p3
    p3=prot_pept_list_cntrl[[j]]
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
  
  if(0){
    for(j in 1:length(c_pep_control_test[1,])){
      #p3=sample(1:length(c_pro_control_test[1,]),sample(1:3,1))
      #prot_pept_list_cntrl_test[[j]]=p3
      p3=prot_pep_list_cntrl_test[[j]]
      
      for(i in 1:length(c_pep_control_test[,1])){
        
        #cat(i,"-",j,"-",p3,"\n")
        for(k in 1:length(p3)){
          c_pep_control_test[i,j]=c_pep_control_test[i,j]+c_pro_control_test[i,p3[k]]
        }
      }
      
      
      
    }
  }
  for(j in 1:length(c_pep_treatement[1,])){
    #p3=sample(1:length(c_pro_treatement[1,]),sample(1:3,1))
    #prot_pept_list_trt[[j]]=p3
    p3=prot_pept_list_trt[[j]]
    
    for(i in 1:length(c_pep_treatement[,1])){
      
      #cat(i,"-",j,"-",p3,"\n")
      for(k in 1:length(p3)){
        ###print("mada")
        c_pep_treatement[i,j]=c_pep_treatement[i,j]+c_pro_treatement[i,p3[k]]
      }
    }
    
    if(1){
      for(i in 1:length(c_pep_treatement_test[,1])){
        
        #cat(i,"-",j,"-",p3,"\n")
        for(k in 1:length(p3)){
          ###print("sheldon")
          c_pep_treatement_test[i,j]=c_pep_treatement_test[i,j]+c_pro_treatement_test[i,p3[k]]
        }
      }
    }
  }
  
  ###print("am out of this loop")
  if(0){ 
    for(j in 1:length(c_pep_treatement_test[1,])){
      #p3=sample(1:length(c_pro_treatement_test[1,]),sample(1:3,1))
      #prot_pept_list_trt_test[[j]]=p3
      p3=prot_pept_list_trt_test[[j]]
      
      for(i in 1:length(c_pep_treatement_test[,1])){
        
        #cat(i,"-",j,"-",p3,"\n")
        for(k in 1:length(p3)){
          c_pep_treatement_test[i,j]=c_pep_treatement_test[i,j]+c_pro_treatement_test[i,p3[k]]
        }
      }
      
      
      
    }
    
  }
  #print("USA")

prot_pept_list_cntrl_test=prot_pept_list_cntrl
prot_pept_list_trt_test=prot_pept_list_trt
####Parameters to find mu_ij

mu_matrix_cntrl=matrix(numeric(Npep_cntrl*number_of_samples_train_cntrl),nrow = number_of_samples_train_cntrl)
mu_matrix_trt=matrix(numeric(Npep_trt*number_of_samples_train_trt),nrow = number_of_samples_train_trt)
mu_matrix_cntrl_test=matrix(numeric(Npep_cntrl*number_of_samples_test_cntrl),nrow = number_of_samples_test_cntrl)
mu_matrix_trt_test=matrix(numeric(Npep_trt*number_of_samples_test_trt),nrow = number_of_samples_test_trt)

kappa=5
#efficiency_vector=runif(Npep,0.1,1)
efficiency_vector_cntrl=runif(Npep_cntrl,0.1,1)
efficiency_vector_trt=runif(Npep_trt,0.1,1)
#print("ginger")

for(i in 1:length(mu_matrix_cntrl[,1])){
  for(j in 1:length(mu_matrix_cntrl[1,])){
    mu_matrix_cntrl[i,j]=c_pep_control[i,j]*kappa*efficiency_vector_cntrl[j]
  }
}



for(i in 1:length(mu_matrix_cntrl_test[,1])){
  for(j in 1:length(mu_matrix_cntrl_test[1,])){
    mu_matrix_cntrl_test[i,j]=c_pep_control_test[i,j]*kappa*efficiency_vector_cntrl[j]
  }
}



for(i in 1:length(mu_matrix_trt[,1])){
  for(j in 1:length(mu_matrix_trt[1,])){
    mu_matrix_trt[i,j]=c_pep_treatement[i,j]*kappa*efficiency_vector_trt[j]
  }
}

#print("ginger")

for(i in 1:length(mu_matrix_trt_test[,1])){
  for(j in 1:length(mu_matrix_trt_test[1,])){
    #print(i)
    mu_matrix_trt_test[i,j]=c_pep_treatement_test[i,j]*kappa*efficiency_vector_trt[j]
  }
}
#print("ginger")

###Generate noisy gaussian###
alpha=0.03
beta=3.6

var_vector_noisy_gaussian_cntrl=alpha*(mu_matrix_cntrl^2)+beta*mu_matrix_cntrl
var_vector_noisy_gaussian_trt=alpha*(mu_matrix_trt^2)+beta*mu_matrix_trt

var_vector_noisy_gaussian_cntrl_test=alpha*(mu_matrix_cntrl_test^2)+beta*mu_matrix_cntrl_test
var_vector_noisy_gaussian_trt_test=alpha*(mu_matrix_trt_test^2)+beta*mu_matrix_trt_test


total_vector_cntrl=matrix(numeric(Npep_cntrl*number_of_samples_train_cntrl),nrow = number_of_samples_train_cntrl)
total_vector_trt=matrix(numeric(Npep_trt*number_of_samples_train_trt),nrow = number_of_samples_train_trt)


total_vector_cntrl_test=matrix(numeric(Npep_cntrl*number_of_samples_test_cntrl),nrow = number_of_samples_test_cntrl)
total_vector_trt_test=matrix(numeric(Npep_trt*number_of_samples_test_trt),nrow = number_of_samples_test_trt)

####Data for transition modification### ####Both equation 13 and 15 in esmaeil's paper are condensed into one.
for(i in 1:length(total_vector_cntrl[,1])){
  for(j in 1:length(total_vector_cntrl[1,])){
    #cat(i,j,"\n")
  #  if(rexp(1,rate = mu_matrix_cntrl[i,j])=="NaN"){
   #   print("NAN spotted!!!")
    #  print(mu_matrix_cntrl[i,j])
    #}
    total_vector_cntrl[i,j]=mu_matrix_cntrl[i,j]+noise_factor_gauss*mvrnorm(n=1,mu = 0, Sigma = abs(var_vector_noisy_gaussian_cntrl[i,j]))+noise_factor_exp*rexp(1,rate = abs(mu_matrix_cntrl[i,j]))
  }
}
#print("lina")
#print(length(mu_matrix_trt[1,]))
#print(length(mu_matrix_trt[,1]))
#print(length(total_vector_trt[1,]))

#print(length(total_vector_trt[,1]))

for(i in 1:length(total_vector_trt[,1])){
  for(j in 1:length(total_vector_trt[1,])){
    #cat(i,j,"\n")
    total_vector_trt[i,j]=mu_matrix_trt[i,j]+noise_factor_gauss*mvrnorm(n=1,mu = 0, Sigma = abs(var_vector_noisy_gaussian_trt[i,j]))+noise_factor_exp*rexp(1,rate = abs(mu_matrix_trt[i,j]))#+
  }
}



for(i in 1:length(total_vector_cntrl_test[,1])){
  for(j in 1:length(total_vector_cntrl_test[1,])){
    #cat(i,j,"\n")
    total_vector_cntrl_test[i,j]=mu_matrix_cntrl_test[i,j]+noise_factor_gauss*mvrnorm(n=1,mu = 0, Sigma = abs(var_vector_noisy_gaussian_cntrl_test[i,j]))+noise_factor_exp*rexp(1,rate = abs(mu_matrix_cntrl_test[i,j]))#+
  }
}

for(i in 1:length(total_vector_trt_test[,1])){
  for(j in 1:length(total_vector_trt_test[1,])){
    #cat(i,j,"\n")
    total_vector_trt_test[i,j]=mu_matrix_trt_test[i,j]+noise_factor_gauss*mvrnorm(n=1,mu = 0, Sigma = abs(var_vector_noisy_gaussian_trt_test[i,j]))+noise_factor_exp*rexp(1,rate = abs(mu_matrix_trt_test[i,j]))#
  }
}
#total_vector_cntrl=data.frame(cbind(total_vector_cntrl,numeric(50)))
#total_vector_trt=data.frame(cbind(total_vector_trt,numeric(50)+1))

#total_vector_cntrl$X31=factor(total_vector_cntrl$X31)
#total_vector_trt$X31=factor(total_vector_trt$X31)


###To calculate rolled up abundances###


x_pro_control=matrix(numeric(Npro_cntrl*number_of_samples_train_cntrl),nrow = number_of_samples_train_cntrl)
x_pro_trt=matrix(numeric(Npro_trt*number_of_samples_train_trt),nrow = number_of_samples_train_trt)

x_pro_control_test=matrix(numeric(Npro_cntrl*number_of_samples_test_cntrl),nrow = number_of_samples_test_cntrl)
x_pro_trt_test=matrix(numeric(Npro_trt*number_of_samples_test_trt),nrow = number_of_samples_test_trt)


pept_to_prot_cntrl=list()
pept_to_prot_trt=list()

pept_to_prot_cntrl_test=list()
pept_to_prot_trt_test=list()


for(i in 1:length(prot_pept_list_cntrl)){
  for(j in 1:length(prot_pept_list_cntrl[i][[1]])){
    #print(prot_pept_list_cntrl[i][[1]][j])
    k1=prot_pept_list_cntrl[i][[1]][j]
    #print(k1)
    pept_to_prot_cntrl[k1][[1]][length(pept_to_prot_cntrl[k1][[1]])+1]=i
    
  }
}

for(i in 1:length(prot_pept_list_trt)){
  for(j in 1:length(prot_pept_list_trt[i][[1]])){
    #print(prot_pept_list_cntrl[i][[1]][j])
    k1=prot_pept_list_trt[i][[1]][j]
    #print(k1)
    pept_to_prot_trt[k1][[1]][length(pept_to_prot_trt[k1][[1]])+1]=i
    
  }
}

#print("greek2")

for(i in 1:length(prot_pept_list_cntrl_test)){
  for(j in 1:length(prot_pept_list_cntrl_test[i][[1]])){
    #print(prot_pept_list_cntrl[i][[1]][j])
    k1=prot_pept_list_cntrl_test[i][[1]][j]
    #print(k1)
    pept_to_prot_cntrl_test[k1][[1]][length(pept_to_prot_cntrl_test[k1][[1]])+1]=i
    
  }
}

#print("greek1")


for(i in 1:length(prot_pept_list_trt_test)){
  for(j in 1:length(prot_pept_list_trt_test[i][[1]])){
    #print(prot_pept_list_cntrl[i][[1]][j])
    k1=prot_pept_list_trt_test[i][[1]][j]
    #print(k1)
    pept_to_prot_trt_test[k1][[1]][length(pept_to_prot_trt_test[k1][[1]])+1]=i
    #print("gopo")
  }
}

#print("greek3")


xlj_control=matrix(numeric(Npro_cntrl*number_of_samples_train_cntrl),nrow = number_of_samples_train_cntrl)
xlj_trt=matrix(numeric(Npro_trt*number_of_samples_train_trt),nrow = number_of_samples_train_trt)

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


#print("greek1")


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
  #print(j)
  
}


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
  #print(j)
  
}
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
  #print(j)
  
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
  #print(j)
  
}




xlj_control=data.frame(xlj_control)
xlj_control=data.frame(cbind(xlj_control,numeric(number_of_samples_train_cntrl)))
colnames(xlj_control)[Npro+1]="type"
#xlj_control$type=as.factor(xlj_control$type)

xlj_control_test=data.frame(xlj_control_test)
xlj_control_test=data.frame(cbind(xlj_control_test,numeric(number_of_samples_test_cntrl)))
colnames(xlj_control_test)[Npro+1]="type"

xlj_trt=data.frame(xlj_trt)
xlj_trt=data.frame(cbind(xlj_trt,rep(1,number_of_samples_train_trt)))
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
  t_test_value=numeric((length(xlj_all[1,])-1))
  len_t_test=length(t_test_value)
  for(i in 1:length(t_test_value)){
    k1=t.test(xlj_all[,i]~xlj_all[,len_t_test+1])
    t_test_value[i]=abs(k1$statistic)
  }
  ##print(t_test_value)
  ##print(order(t_test_value))
  
  p11=xlj_all$type
  p22=xlj_all_test$type
  
  if(0){
    xlj_all=xlj_all[,order(t_test_value,decreasing = TRUE)]
    xlj_all=xlj_all[,1:Npro_for_analysis]
    xlj_all$"type"=p11
    xlj_all_test=xlj_all_test[,order(t_test_value,decreasing = TRUE)]
    xlj_all_test=xlj_all_test[,1:Npro_for_analysis]
    xlj_all_test$"type"=p22
    fold_change_vec=fold_change_vec[order(t_test_value,decreasing = TRUE)]
    fold_change_vec=fold_change_vec[1:Npro_for_analysis]
  }
  ##print("Reached the end :-) ")
  #print("hina-rabbani")
  
  if(train_test_flag){
    return(xlj_all)
  }
  else{
    return(xlj_all_test)
  }
####xlj_all is the required data as per equation 12###
####xlj_all_test is the required test data####
	
##################:-)###########################	
  
	  
}

main_function=function(pro_trt_file,pep_trt_file,pro_cntrl_file,pep_cntrl_file,noise_factor_gauss,noise_factor_exp,theta_a,number_of_samples_train_cntrl,number_of_samples_train_trt,
                     number_of_samples_test_cntrl,number_of_samples_test_trt,amin,amax,ka0,kc0,thetaa0,thetac0,phi0,al0,M_cal,markov_iter,frac_npro)
{
print("I am inside")
Npro_trt=fasta_file_to_array_conv(pro_trt_file,pep_trt_file,"pro")


Npro_cntrl=fasta_file_to_array_conv(pro_cntrl_file,pep_cntrl_file,"pro")

Npro_cntrl_a=ceiling(Npro_cntrl*frac_npro)
Npro_cntrl_c=Npro_cntrl-Npro_cntrl_a

Npro_trt_a=ceiling(Npro_trt*frac_npro)
Npro_trt_c=Npro_trt-Npro_trt_a
Npro=Npro_cntrl
Npro_for_analysis=0.6*Npro_trt

Npep_trt=fasta_file_to_array_conv(pro_trt_file,pep_trt_file,"pep")
Npep_cntrl=fasta_file_to_array_conv(pro_cntrl_file,pep_cntrl_file,"pep")


eta_treatement=c()
eta_control=c()
	
gamma_dist_c0=rgamma(1000,shape = kc0,scale = thetac0)
gamma_dist_a0=rgamma(1000,shape = ka0,scale = thetaa0)
	
gamma_l0=numeric(Npro_cntrl)

print("henene")
for(i in 1:(Npro_cntrl_c+Npro_cntrl_a)){
  if(i<Npro_cntrl_c){
  gamma_l0[i]=sample(gamma_dist_c0,1)
  }
  else{
    gamma_l0[i]=sample(gamma_dist_a0,1)
  }
}

mean_vec_cont=gamma_l0
mean_vec_cont0=gamma_l0

sig_matrix0=matrix(numeric(Npro_cntrl*Npro_cntrl),nrow=Npro_cntrl)
sig_sq_vec0=phi0*mean_vec_cont*mean_vec_cont


for(i in 1:Npro_cntrl){
  sig_matrix0[i,i]=sig_sq_vec0[i]
}

	
c_pro_control0=mvrnorm(n=number_of_samples_train_cntrl,mu = mean_vec_cont,Sigma = sig_matrix0)
c_pro_treatement0=mvrnorm(n=number_of_samples_train_cntrl,mu = mean_vec_cont,Sigma = sig_matrix0)
ratio_vec0=colMeans(c_pro_control0)/colMeans(c_pro_treatement0)
	
	
k_rand=sample(160:240,1)/100
kc_rand=sample(160:240,1)/100
ka_rand=sample(400:600,1)/100

thetac_rand=sample(80:120,1)
thetaa_rand=sample(9e6:11e6,1)
	
phi_rand=sample(0.3:0.5,1)
gamma_distc_rand=rgamma(1000,shape=kc_rand,scale=thetac_rand)
gamma_dista_rand=rgamma(1000,shape=ka_rand,scale=thetaa_rand)
	
	
gamma_rand=numeric(Npro_cntrl)


for(i in 1:(Npro_cntrl_c+Npro_cntrl_a)){
  if(i<Npro_cntrl_c){
  gamma_rand[i]=sample(gamma_distc_rand,1)
  }
  else{
    gamma_rand[i]=sample(gamma_dista_rand,1)
  }
}

mean_vec_cont_rand=gamma_rand
mean_vec_cont_rand=gamma_rand
	
sig_matrix_rand=matrix(numeric(Npro_cntrl*Npro_cntrl),nrow=Npro_cntrl)
sig_sq_vec_rand=phi_rand*mean_vec_cont_rand*mean_vec_cont_rand

for(i in 1:Npro_cntrl){
  sig_matrix_rand[i,i]=sig_sq_vec_rand[i]
}


c_pro_control_rand=mvrnorm(n=number_of_samples_train_cntrl,mu = mean_vec_cont_rand,Sigma = sig_matrix_rand)
c_pro_treatement_ran=mvrnorm(n=number_of_samples_train_cntrl,mu = mean_vec_cont_rand,Sigma = sig_matrix_rand)
	
diff_rand_vec_cntrl=colMeans(c_pro_control_rand)-colMeans(c_pro_control0)
thresh_key_cntrl=1*euc_norm(diff_rand_vec_cntrl)

fold_change_vec0=numeric(Npro)
print("goofy")
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
mean_vec_treatement_rand=mean_vec_cont_rand*fold_change_vec0
	
###Defining c_pro_treatement_0 again for convinience#####

c_pro_treatement0=mvrnorm(n=number_of_samples_train_trt,mu = mean_vec_treatement0,Sigma = sig_matrix0)
c_pro_treatement_rand=mvrnorm(n=number_of_samples_train_trt,mu = mean_vec_treatement_rand,Sigma = sig_matrix_rand)

diff_rand_vec_trt=colMeans(c_pro_treatement_rand)-colMeans(c_pro_treatement0)
thresh_key_trt=1.3*euc_norm(diff_rand_vec_trt)
	
	
print("India")
k_list=list()
theta_list=list()
phi_list=list()
thetaa_list=list()
thetac_list=list()
ka_list=list()
kc_list=list()
count_1=0
	
###ABC-Rejection Sampling##########

norm_array=numeric(M_cal)

####Esmaeil parameters###
kc=sample(160:240,1)/100
ka=sample(450:550,1)/100

thetac=sample(80:120,1)
thetaa=sample(9e6:11e6,1)

print("hinj")
for(j in 1:M_cal){
  k=sample(160:240,1)/100
  theta=sample(800:1200,1)
  kc=sample(160:240,1)/100
  ka=sample(450:550,1)/100
  thetac=sample(80:120,1)
  thetaa=sample(9e6:11e6,1)
  #gamma_dist=rgamma(1000,shape=k,scale=theta_a)
  gamma_dist_c=rgamma(1000,shape = kc,scale = thetac)
  gamma_dist_a=rgamma(1000,shape = ka,scale = thetaa)
  
  gamma_l=numeric(Npro_cntrl)
  
  
  
  #for(i in 1:Npro_cntrl){
   # gamma_l[i]=sample(gamma_dist,1)
  #}
  
  for(i in 1:(Npro_cntrl_c+Npro_cntrl_a)){
    if(i<Npro_cntrl_c){
      gamma_l[i]=sample(gamma_dist_c,1)
    }
    else{
      gamma_l[i]=sample(gamma_dist_a,1)
    }
  }
  mean_vec_cont=gamma_l
  
  phi=runif(Npro,0.3,0.5) ###Defining the coefficient of variation
  
  
  sig_matrix=matrix(numeric(Npro_cntrl*Npro_cntrl),nrow=Npro_cntrl)
  sig_sq_vec=phi*mean_vec_cont*mean_vec_cont
  
  for(i in 1:Npro_cntrl){
    sig_matrix[i,i]=sig_sq_vec[i]
  }
  
  
  c_pro_control=mvrnorm(n=number_of_samples_train_cntrl,mu = mean_vec_cont,Sigma = sig_matrix)
  
  diff_vec=colMeans(c_pro_control)-colMeans(c_pro_control0)
  norm_array[j]=euc_norm(diff_vec)
  # print(euc_norm(diff_vec))
  #thresh_key=2*(median(norm_array)+min(norm_array))/3
  if(norm_array[j]<thresh_key_cntrl){
    #print(diff_vec)
    
    count_1=count_1+1
    k_list[count_1]=k
    theta_list[count_1]=theta
    phi_list[count_1]=phi
    
    ka_list[count_1]=ka
    kc_list[count_1]=kc
    
    thetaa_list[count_1]=thetaa
    thetac_list[count_1]=thetac
  }
  #cat("\n")
  
  
  
  
}
#print(mean(norm_array))
#cat("\n")

k_vec=theta_vec=phi_vec=numeric(length(k_list))

ka_vec=kc_vec=thetaa_vec=thetac_vec=phi_vec=numeric(length(ka_list))

for(i in 1:length(k_vec)){
  k_vec[i]=k_list[[i]]
  theta_vec[i]=theta_list[[i]]
  phi_vec[i]=phi_list[[i]]
}
print("french")
for(i in 1:length(ka_vec)){
  ka_vec[i]=ka_list[[i]]
  kc_vec[i]=kc_list[[i]]
  
  thetaa_vec[i]=thetaa_list[[i]]
  thetac_vec[i]=thetac_list[[i]]
  
}

k_opt=mean(k_vec)
theta_opt=mean(theta_vec)
phi_opt=mean(phi_vec)

ka_opt=mean(ka_vec)
kc_opt=mean(kc_vec)

thetaa_opt=mean(thetaa_vec)
thetac_opt=mean(thetac_vec)

#############ABC Rejection sampling done#######

##### Data from optimal parameters being generated#######
	


print("glucose")
gamma_dist_opt=rgamma(1000,shape=k_opt,scale=theta_opt)
gamma_dist_a_opt=rgamma(1000,shape=ka_opt,scale=thetaa_opt)
gamma_dist_c_opt=rgamma(1000,shape=kc_opt,scale=thetac_opt)


gamma_l_opt=numeric(Npro_cntrl)

for(i in 1:(Npro_cntrl_c+Npro_cntrl_a)){
  if(i<=Npro_cntrl_c){
    gamma_l_opt[i]=sample(gamma_dist_c_opt,1)
  }
  else{
    gamma_l_opt[i]=sample(gamma_dist_a_opt,1)
  }
}

mean_vec_cont_opt=gamma_l_opt

sig_matrix_opt=matrix(numeric(Npro*Npro),nrow=Npro)
sig_sq_vec_opt=phi_opt*mean_vec_cont_opt*mean_vec_cont_opt

for(i in 1:Npro_cntrl){
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

ratio_vec=colMeans(c_pro_treatement0)/colMeans(c_pro_control0)

for(i in 1:Npro){
  if(ratio_vec[i]>1)
  {
    fold_change_vec[i]=fold_change(1,amin,amax)
  }
  else{
    fold_change_vec[i]=fold_change(0,amin,amax)
    
  }
  
}

fold_change_vec=ratio_vec

prot_pept_list_cntrl=fasta_file_to_array_conv(pro_cntrl_file,pep_cntrl_file,1)
prot_pept_list_trt=fasta_file_to_array_conv(pro_trt_file,pep_trt_file,1)



####Generating random peptide control data###



prot_pept_list_cntrl_test=prot_pept_list_cntrl
prot_pept_list_trt_test=prot_pept_list_trt

	

##print("MCMC algorithm starting")
#######Upamanyu's algorithm-3...The ABC-MCMC-algorithm for esmaeil's paper#####

#####First three steps of ABC-MCMC#####

###gamma_0 is the gamma related to S^(0)_(0) in the paper. This is generated with k_opt and theta_opt
###gamma_dist_0 is the proper gamma distribution pertaining to S_(0) in the paper. Its generated with k_0 and theta_0.
	 
#Step1: Sampling gamma and generating mean vectors##
	  
print("goosehammer")

gamma_dist_0=rgamma(1000,shape=k_opt,scale=theta_opt)
gamma_dist_a_0=rgamma(1000,shape=ka_opt,scale=thetaa_opt)
gamma_dist_c_0=rgamma(1000,shape=kc_opt,scale=thetac_opt)
	  
c_pro_control_0=c_pro_control0
c_pro_treatment_0=c_pro_treatement0
gamma_0=numeric(Npro)

for(i in 1:(Npro_cntrl_c+Npro_cntrl_a)){
  if(i<=Npro_cntrl_c){
    gamma_0[i]=sample(gamma_dist_c_0,1)
  }
  else{
    gamma_0[i]=sample(gamma_dist_a_0,1)
  }
}

mean_vec_cntrl_0_0=gamma_0
mean_vec_trt_0_0=gamma_0*fold_change_vec

sig_matrix_0_0=matrix(numeric(Npro*Npro),nrow=Npro)
sig_sq_vec_0_0=phi_opt*mean_vec_cntrl_0_0*mean_vec_cntrl_0_0
for(i in 1:Npro){
  sig_matrix_0_0[i,i]=sig_sq_vec_0_0[i]
}
	  
###Step2: Generating protein data
c_pro_control_0_0=mvrnorm(n=number_of_samples_train_cntrl,mu = mean_vec_cntrl_0_0,Sigma = sig_matrix_0_0)
c_pro_treatment_0_0=mvrnorm(n=number_of_samples_train_trt,mu = mean_vec_trt_0_0,Sigma = sig_matrix_0_0)

#print("goi")
	  
	  
diff_vec_cntrl=colMeans(c_pro_control_0_0)-colMeans(c_pro_control_0)
diff_vec_trt=colMeans(c_pro_treatment_0_0)-colMeans(c_pro_treatment_0)

norm_cntrl=euc_norm(diff_vec_cntrl)

norm_trt=euc_norm(diff_vec_trt)

	  
while(norm_cntrl>1.8*thresh_key_cntrl | norm_trt>thresh_key_trt){
for(i in 1:(Npro_cntrl_c+Npro_cntrl_a)){
  if(i<=Npro_cntrl_c){
    gamma_0[i]=sample(gamma_dist_c_0,1)
  }
  else{
    gamma_0[i]=sample(gamma_dist_a_0,1)
  }
}
  mean_vec_cntrl_0_0=gamma_0
  mean_vec_trt_0_0=gamma_0*fold_change_vec
	
  sig_matrix_0_0=matrix(numeric(Npro*Npro),nrow=Npro)
  sig_sq_vec_0_0=phi_opt*mean_vec_cntrl_0_0*mean_vec_cntrl_0_0
  for(i in 1:Npro){
    sig_matrix_0_0[i,i]=sig_sq_vec_0_0[i]
  }
	
###Step2: Generating protein data
  c_pro_control_0_0=mvrnorm(n=number_of_samples_train_cntrl,mu = mean_vec_cntrl_0_0,Sigma = sig_matrix_0_0)
  c_pro_treatment_0_0=mvrnorm(n=number_of_samples_train_trt,mu = mean_vec_trt_0_0,Sigma = sig_matrix_0_0)
  
  diff_vec_cntrl=colMeans(c_pro_control_0_0)-colMeans(c_pro_control_0)
  diff_vec_trt=colMeans(c_pro_treatment_0_0)-colMeans(c_pro_treatment_0)
  
  norm_cntrl=euc_norm(diff_vec_cntrl)
  
  norm_trt=euc_norm(diff_vec_trt)
  
}
#####Step 3 is just an if condition which is taken care at the begining

#####Steps 5, 6 and 7 of the markov chain####
c_pro_control_markov=c_pro_control_0_0
c_pro_trt_markov=c_pro_treatment_0_0
gamma_vec=colMeans(c_pro_control_markov)
	  
	  
	  
#####Main MCMC chain#####
count_markov=0
gamma_markov=matrix(0L, nrow = markov_iter,ncol=Npro)
for(j in 1:markov_iter){
  if(j==1){gamma_vec=gamma_l_opt}
  else{
    gamma_vec=colMeans(c_pro_control_markov)
  }
  mean_vec_cntrl_markov=gamma_vec
  mean_vec_trt_markov=gamma_vec*fold_change_vec
  ###print(i)
  sig_matrix_markov=matrix(numeric(Npro*Npro),nrow=Npro)
  sig_sq_vec_markov=phi_opt*mean_vec_cntrl_markov*mean_vec_cntrl_markov
  sig_sq_vec_markov=rep(0.1,length(mean_vec_cont))
  for(i in 1:Npro){
    sig_matrix_markov[i,i]=sig_sq_vec_markov[i]
  }
  
  c_pro_control_markov=mvrnorm(n=number_of_samples_train_cntrl,mu = mean_vec_cntrl_markov,Sigma = sig_matrix_markov)
  c_pro_trt_markov=mvrnorm(n=number_of_samples_train_trt,mu = mean_vec_trt_markov,Sigma = sig_matrix_markov)
  
  diff_vec_cntrl_markov=colMeans(c_pro_control_markov)-colMeans(c_pro_control_0)
  diff_vec_trt_markov=colMeans(c_pro_trt_markov)-colMeans(c_pro_treatment_0)
  
  norm_cntrl_markov=euc_norm(diff_vec_cntrl_markov)
  
  norm_trt_markov=euc_norm(diff_vec_trt_markov)
  ###print("-----------")
  ###print(norm_cntrl_markov)
  ###print(norm_trt_markov)
  ###print("-----------")

if(norm_cntrl_markov<(0.3*thresh_key_cntrl) & norm_trt_markov<(0.3*thresh_key_trt)){
    #mean_ratio=
    gamma_vec=colMeans(c_pro_control_markov)
    ###print(j)
    count_markov=count_markov+1
    ##print("zhezh")
  }
  gamma_markov[j,]=gamma_vec
  
}		
print("MCMC algorithm done")
######Algorithm-3 in upamanyu paper done#####




#####Generation of kernel data and then classification###

print("running the function")
##print(gamma_vec)
#fold_change_vec=sample(5:10,Npro,replace = TRUE)




train_test_flag=1

print("Came till here")
#print(fold_change_vec)
xlj_all=peptide_generation(gamma_dist_c_opt,gamma_dist_a_opt,fold_change_vec,phi_opt,Npro_cntrl,perc_cntrl,Npep_cntrl,Npro_trt,perc_trt,Npep_trt,prot_pept_list_cntrl,Npro_for_analysis,noise_factor_gauss,
                            noise_factor_exp,prot_pept_list_cntrl_test,prot_pept_list_trt,prot_pept_list_trt_test,train_test_flag)
print("passed-step-one")

train_test_flag=0
xlj_all_test=peptide_generation(gamma_dist_c_opt,gamma_dist_a_opt,fold_change_vec,phi_opt,Npro_cntrl,perc_cntrl,Npep_cntrl,Npro_trt,perc_trt,Npep_trt,prot_pept_list_cntrl,Npro_for_analysis,noise_factor_gauss,
                                noise_factor_exp,prot_pept_list_cntrl_test,prot_pept_list_trt,prot_pept_list_trt_test,train_test_flag)

print("passed-step-two")


xlj_all_test_new=xlj_all_test
xlj_all_new=xlj_all
	
#1. LDA classifier#
if(1){
  print("running lda and knn on new data")
  xlj_all_test_new=xlj_all_test[sample(nrow(xlj_all_test_new)),]
  
  lda_classifier_50<-lda(type ~ .,data=xlj_all_new)
  predictions_lda_50=predict(lda_classifier_50,xlj_all_test_new[,1:Npro])$class
  table_data_lda_50_2=table(predictions_lda_50,xlj_all_test_new[,Npro+1])
  
  predictions_lda_50_app=predict(lda_classifier_50,xlj_all_new[,1:Npro])$class
  table_data_lda_50_app=table(predictions_lda_50_app,xlj_all_new[,Npro+1])
  print(table_data_lda_50_2)
  ##2. KNN classifier####
  
  training_labels=xlj_all_new$type
  knn_trained<-knn(train = xlj_all_new[,1:Npro]  , test =xlj_all_test_new[,1:Npro] , cl = training_labels, k=3)
  table_data_knn_50_2=table(knn_trained,xlj_all_test_new[,Npro+1])
  knn_trained_app<-knn(train = xlj_all_new[,1:Npro]  , test =xlj_all_new[,1:Npro] , cl = training_labels, k=3)
  table_data_knn_50_app=table(knn_trained_app,xlj_all_new[,Npro+1])
  
}


train_test_flag=1

v2=gamma_vec
f2=fold_change_vec
  
markov_Data=peptide_generation(gamma_dist_c_opt,gamma_dist_a_opt,fold_change_vec,phi_opt,Npro_cntrl,perc_cntrl,Npep_cntrl,Npro_trt,perc_trt,Npep_trt,prot_pept_list_cntrl,Npro_for_analysis,noise_factor_gauss,
                            noise_factor_exp,prot_pept_list_cntrl_test,prot_pept_list_trt,prot_pept_list_trt_test,train_test_flag)

gaussian_normal_kernel=mvrnorm(n=1000,mu=numeric(Npro),Sigma = diag(Npro))

gaussian_normal_kernel_density=density(gaussian_normal_kernel)
  
abc_mcmc_result_vector=rep(2,length(xlj_all_test[,1]))
#print("1234")
total_outside_sum_one_arr=numeric(nrow(xlj_all_test))
total_outside_sum_zero_arr=numeric(nrow(xlj_all_test))
outside_full_sum_one_arr=numeric(nrow(xlj_all_test))
outside_full_sum_zero_arr=numeric(nrow(xlj_all_test))

print("Final abc-mcmc-model... yay")
#######Prediction by the ABC MCMC model #####
if(1){
for(i in 1:nrow(xlj_all_test)){
    test_data=xlj_all_test[i,1:Npro]
    total_outside_sum_zero=0 ####outside Sum for all the markov iterations
    total_outside_sum_one=0
    outside_full_sum_zero=0
    outside_full_sum_one=1
    print(paste("Running the row",i))
    
    ##Only last 70% are considered for the calculations as first 30% are considered as the burn-in stage
    count_123=0
    ##print("here man")
    for(j in (ceiling(0.4*markov_iter)):(markov_iter-1)){
      # ##print(paste("Markov-chain_number",j))
      #total_outside_sum_zero=0
      #total_outside_sum_one=0
     
    markov_Data=peptide_generation(gamma_dist_c_opt,gamma_dist_a_opt,fold_change_vec,phi_opt,Npro_cntrl,perc_cntrl,Npep_cntrl,Npro_trt,perc_trt,Npep_trt,prot_pept_list_cntrl,Npro_for_analysis,noise_factor_gauss,
                            noise_factor_exp,prot_pept_list_cntrl_test,prot_pept_list_trt,prot_pept_list_trt_test,train_test_flag)

    markov_row_length=length(markov_Data[,1]) 
    markov_length_class=markov_row_length/2
    total_inside_sum_zero=0 ###inside sum is for sum of all the n-values of markov_data
    total_inside_sum_one=0
    inside_full_sum_zero=0
    inside_full_sum_one=0
    for(k in 1:(markov_length_class)){
      #inside_sum_zero=0
      diff_data_zero=markov_Data[k,1:Npro]-test_data
      diff_data_zero=markov_Data[k,1:Npro]-test_data
      t0=diff_data_zero
      ##print(paste("Zero-",sum(diff_data_zero)))
      #diff_data_zero=abs(normalize(markov_Data[k,1:Npro])-normalize(test_data))
        
      ##print("###########zero#########")
      ##print(diff_data_zero)
      mean_diff_data_zero=sum(diff_data_zero)/length(diff_data_zero)
      # #print(paste("Zero-",sum(diff_data_zero)))
        
      ##print(paste("zero-",mean_diff_data_zero))
      #inside_full_sum_zero=inside_full_sum_zero+mean_diff_data_zero
      inside_full_sum_zero=inside_full_sum_zero+sum(diff_data_zero)
if(1){
          #diff_data_zero=diff_data_zero/abs(max(diff_data_zero))
          diff_data_zero=diff_data_zero/5000
          
          kernel_data=approx(gaussian_normal_kernel_density$x,gaussian_normal_kernel_density$y,xout = diff_data_zero)
          kernel_data$y[is.na(kernel_data$y)]=0
          total_inside_sum_zero=total_inside_sum_zero+mean(kernel_data$y)
        }
      }
      
      for(l in (markov_length_class+1):(markov_row_length)){
        #inside_sum_one=0
        ###print("gunda")
        ###print(l)
        diff_data_one=markov_Data[l,1:Npro]-test_data
        # #print(paste("one-",sum(diff_data_one)))
        
        #diff_data_one=abs(normalize(markov_Data[l,1:Npro])-normalize(test_data))
        ###print(diff_data_one)
        t1=diff_data_one
        ##print("###########one#########")
        mean_diff_data_one=sum(diff_data_one)/length(diff_data_one)
        
        ##print(paste("one-",mean_diff_data_one))
        # #print(paste("one-",sum(diff_data_one)))
        
        #inside_full_sum_one=inside_full_sum_one+mean_diff_data_one
        inside_full_sum_one=inside_full_sum_one+sum(diff_data_one)
        #plot(diff_data_one)
        #k22=cbind(diff_data_one,diff_data_zero)
        #plot(colMeans(k22))
        ##print("#########################")
        if(1){
          # diff_data_one=diff_data_one/abs(max(diff_data_one))
          diff_data_one=diff_data_one/5000
          
          kernel_data=approx(gaussian_normal_kernel_density$x,gaussian_normal_kernel_density$y,xout = diff_data_one)
          kernel_data$y[is.na(kernel_data$y)]=0
          total_inside_sum_one=total_inside_sum_one+mean(kernel_data$y)
        }
        k22=cbind(diff_data_one,diff_data_zero)
        # plot(colMeans(k22))
        
        if(is.na(total_inside_sum_one)& count_123==0){
          ###print(kernel_data$y)
          count_123=1
          p11=t1
        }
        ###print(total_inside_sum_one)
        
      }
      #total_outside_sum_zero=total_outside_sum_zero+total_inside_sum_zero
      #total_outside_sum_one=total_outside_sum_one+total_inside_sum_one
    }
    total_outside_sum_zero=total_outside_sum_zero+total_inside_sum_zero
    total_outside_sum_one=total_outside_sum_one+total_inside_sum_one
    outside_full_sum_zero=outside_full_sum_zero+inside_full_sum_zero
    outside_full_sum_one=outside_full_sum_one+inside_full_sum_one
    
    ###print("glagla")
    ###print(total_outside_sum_one)
    ###print(total_outside_sum_zero)
    total_outside_sum_one_arr[i]=total_outside_sum_one
    total_outside_sum_zero_arr[i]=total_outside_sum_zero
    outside_full_sum_one_arr[i]=outside_full_sum_one
    
    outside_full_sum_zero_arr[i]=outside_full_sum_zero
    
    if(total_outside_sum_one<total_outside_sum_zero){
      #if(abs(outside_full_sum_one)>abs(outside_full_sum_zero)){
      abc_mcmc_result_vector[i]=0
    }
    else{
      abc_mcmc_result_vector[i]=1
    }
    
  }
  #print("------------------------")
  print("Fullll ABC-MCMC model done")
  
  #print("------------------------")
  
  table_data_abc_mcmc_50=table(abc_mcmc_result_vector,xlj_all_test[,Npro+1])

}
lda_error=(table_data_lda_50_2[2]+table_data_lda_50_2[3])/sum(table_data_lda_50_2)
knn_error=(table_data_knn_50_2[2]+table_data_knn_50_2[3])/sum(table_data_knn_50_2)

abc_mcmc_error=(table_data_abc_mcmc_50[2]+table_data_abc_mcmc_50[3])/sum(table_data_abc_mcmc_50)
print("i izzzz here")
return(c(lda_error,knn_error,abc_mcmc_error))
}
    
######Main function ends here######
final_time=proc.time()

####Paremeters of the function###

#Npro_trt=fasta_file_to_array_conv("pro_trt_file.fasta","pep_trt_il.fasta","pro")
#Npro_cntrl=fasta_file_to_array_conv("pro_cntrl.fasta","pep_cntrl.fasta","pro")
#Npro_for_analysis=0.6*Npro_trt

pro_trt_file="pro_trt.fasta"
pep_trt_file="pep_trt.fasta"

if(0){
  pro_trt_file="pro_cntrl.fasta"
  pep_trt_file="pep_cntrl.fasta"
}
pro_cntrl_file="pro_cntrl.fasta"
pep_cntrl_file="pep_cntrl.fasta"
noise_factor=50
noise_factor_gauss=1
noise_factor_exp=1
theta_a=100
number_of_samples_train_cntrl=20 #number of samples in each class
number_of_samples_train_trt=20

number_of_samples_test_cntrl=20 #number of test samples fr control
number_of_samples_test_trt=20 
amin=1.5
amax=1.6
k0=2
ka0=5
kc0=2
thetac0=100
thetaa0=10000000
theta0=100
phi0=0.4
al0=1.55
M_cal=1000 
markov_iter=10
frac_npro=0.85
kashyap="kikakika"


table_data_abc_mcmc=main_function(pro_trt_file,pep_trt_file,pro_cntrl_file,pep_cntrl_file,noise_factor_gauss,noise_factor_exp,theta_a,number_of_samples_train_cntrl,number_of_samples_train_trt,
number_of_samples_test_cntrl,number_of_samples_test_trt,amin,amax,ka0,kc0,thetaa0,thetac0,phi0,al0,M_cal,markov_iter,frac_npro)
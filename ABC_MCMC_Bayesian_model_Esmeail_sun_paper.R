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
norm_flag=0

#print(runif(n=100,min=0,max=100))
#####Please check with random seed###
#rm(list=ls())
#set.seed(101)
library(class)
library(boot)
library(MASS) 

#library(caret)
#library(ISLR)
#library("Biostrings")
#setwd("/home/kashyap/Desktop/Masters_thesis_related/Codes")
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




euc_norm <- function(x) sqrt(sum(x^2))

if(1){
  pro_trt_file="prot_file1.fasta"
  pep_trt_file="pept_file1.fasta"
}
pro_cntrl_file="prot_file1.fasta"
pep_cntrl_file="pept_file1.fasta"

if(1){
  
  Npro_factor=0.03
  #Npro_trt=fasta_file_to_array_conv(pro_trt_file,pep_trt_file,"pro")
  
  #Npro_cntrl=fasta_file_to_array_conv(pro_cntrl_file,pep_cntrl_file,"pro")
  
  Npro_for_analysis=ceiling(Npro_factor*Npro_trt)
  
  #Npep_trt=fasta_file_to_array_conv(pro_trt_file,pep_trt_file,"pep")
  #Npep_cntrl=fasta_file_to_array_conv(pro_cntrl_file,pep_cntrl_file,"pep")
  
  
}
	
noise_factor_gauss=1
noise_factor_exp=0
Npro=Npro_trt

#Npro_trt_test=60
#Npro_cntrl_test=60
#Npep_trt=ceiling(15*Npro_trt)
#Npep_cntrl=ceiling(15*Npro_cntrl)
Npep=Npep_trt	
theta_a=100
thetac_a=100
thetaa_a=10000000	
number_of_samples_train_cntrl=10 #number of samples in each class
number_of_samples_train_trt=10

number_of_samples_test_cntrl=10 #number of test samples fr control
number_of_samples_test_trt=10 #number of test samples fr treatement
	
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

	
M_cal=20 #Number of calibarations to be made in abc rejection algo.

#Generating the synthetic sample data S0####
###Only control sample used#####

###initial parameters as used in table 2
ka0=5
kc0=2
thetac0=100
thetaa0=10000000
theta0=100
phi0=0.4
al0=1.55	
gamma_dist_a0=rgamma(1000,shape=ka0,scale=thetaa0)
gamma_dist_c0=rgamma(1000,shape=kc0,scale=thetac0)

gamma_l0=numeric(Npro_cntrl)
Npro_cntrl_a=floor(0.85*Npro)
Npro_cntrl_c=Npro-Npro_cntrl_a

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

sig_matrix0=matrix(numeric(Npro*Npro),nrow=Npro)
sig_sq_vec0=phi0*mean_vec_cont*mean_vec_cont
	
for(i in 1:Npro){
  sig_matrix0[i,i]=sig_sq_vec0[i]
}
	

c_pro_control0=mvrnorm(n=number_of_samples_train_cntrl,mu = mean_vec_cont,Sigma = sig_matrix0)
c_pro_treatement0=mvrnorm(n=number_of_samples_train_cntrl,mu = mean_vec_cont,Sigma = sig_matrix0)
if(norm_flag){
c_pro_control0[,Npro_cntrl_c:length(c_pro_control0[1,])]=runif(1,min=0,max=10e-6)*c_pro_control0[,Npro_cntrl_c:length(c_pro_control0[1,])]
c_pro_treatement0[,Npro_cntrl_c:length(c_pro_treatement0[1,])]=runif(1,min=0,max=10e-6)*c_pro_treatement0[,Npro_cntrl_c:length(c_pro_treatement0[1,])]
}
ratio_vec0=colMeans(c_pro_control0)/colMeans(c_pro_treatement0)	

######For calculating the threshold for rejection sampling

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

sig_matrix_rand=matrix(numeric(Npro_cntrl*Npro_cntrl),nrow=Npro_cntrl)
sig_sq_vec_rand=phi_rand*mean_vec_cont_rand*mean_vec_cont_rand

for(i in 1:Npro_cntrl){
  sig_matrix_rand[i,i]=sig_sq_vec_rand[i]
}


c_pro_control_rand=mvrnorm(n=number_of_samples_train_cntrl,mu = mean_vec_cont_rand,Sigma = sig_matrix_rand)
c_pro_treatement_rand=mvrnorm(n=number_of_samples_train_cntrl,mu = mean_vec_cont_rand,Sigma = sig_matrix_rand)
	
if(norm_flag){
  
c_pro_control_rand[,Npro_cntrl_c:length(c_pro_control_rand[1,])]=runif(1,min=0,max=10e-6)*c_pro_control_rand[,Npro_cntrl_c:length(c_pro_control_rand[1,])]
c_pro_treatement_rand[,Npro_cntrl_c:length(c_pro_treatement_rand[1,])]=runif(1,min=0,max=10e-6)*c_pro_treatement_rand[,Npro_cntrl_c:length(c_pro_treatement_rand[1,])]
}
diff_rand_vec_cntrl=colMeans(c_pro_control_rand)-colMeans(c_pro_control0)
thresh_key_cntrl=1*Npro_factor*euc_norm(diff_rand_vec_cntrl)

#print(rgamma(100,shape=k_rand,scale=theta_rand))
#print(thresh_key_cntrl)
#print(Npro_factor)
while(0){
  x=1
}
print("ganganna")

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
	
c_pro_treatement0=mvrnorm(n=number_of_samples_train_trt,mu = mean_vec_treatement0,Sigma = sig_matrix0)
c_pro_treatement_rand=mvrnorm(n=number_of_samples_train_trt,mu = mean_vec_treatement_rand,Sigma = sig_matrix_rand)
if(norm_flag){
  
c_pro_treatement0[,Npro_cntrl_c:length(c_pro_treatement0[1,])]=runif(1,min=0,max=10e-6)*c_pro_treatement0[,Npro_cntrl_c:length(c_pro_treatement0[1,])]
c_pro_treatement_rand[,Npro_cntrl_c:length(c_pro_treatement_rand[1,])]=runif(1,min=0,max=10e-6)*c_pro_treatement_rand[,Npro_cntrl_c:length(c_pro_treatement_rand[1,])]
}  
  
  
diff_rand_vec_trt=colMeans(c_pro_treatement_rand)-colMeans(c_pro_treatement0)
thresh_key_trt=euc_norm(diff_rand_vec_trt)
thresh_key_cntrl=euc_norm(diff_rand_vec_cntrl)
	
print("India")
k_list=list()
theta_list=list()
phi_list=list()
thetaa_list=list()
thetac_list=list()
ka_list=list()
kc_list=list()
count_1=0

###Synthetic sample data S0 done#####

###ABC-Rejection Sampling#########

norm_array=numeric(M_cal)

	
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
  if(norm_flag){
    
  c_pro_control[,Npro_cntrl_c:length(c_pro_control[1,])]=runif(1,min=0,max=10e-6)*c_pro_control[,Npro_cntrl_c:length(c_pro_control[1,])]
}
  
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
	
}
	
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

mean_vec_treatement=mean_vec_treatement1*fold_change_vec

####Generating random peptide control data###
 c_pro_control=mvrnorm(n=number_of_samples_train_cntrl,mu = mean_vec_cont_opt,Sigma = sig_matrix_opt)
c_pro_treatement=mvrnorm(n=number_of_samples_train_trt,mu = mean_vec_treatement,Sigma = sig_matrix_opt)
if(norm_flag){
  
c_pro_control[,Npro_cntrl_c:length(c_pro_control[1,])]=runif(1,min=0,max=10e-6)*c_pro_control[,Npro_cntrl_c:length(c_pro_control[1,])]
c_pro_treatement[,Npro_cntrl_c:length(c_pro_treatement[1,])]=runif(1,min=0,max=10e-6)*c_pro_treatement[,Npro_cntrl_c:length(c_pro_treatement[1,])]
  
}
  
  
l1=c_pro_control
l1=data.frame(cbind(l1,rep(0,number_of_samples_train_cntrl)))
colnames(l1)[Npro+1]="type"

l2=c_pro_treatement
l2=data.frame(cbind(l2,rep(1,number_of_samples_train_trt)))
colnames(l2)[Npro+1]="type"
#colnames(xlj_control)[Npro+1]="type"

c_pro_all=rbind(l1,l2)
c_pro_all$type=as.factor(c_pro_all$type)

t_test_value_pro=numeric((length(c_pro_all[1,])-1))
len_t_test_pro=length(t_test_value_pro)

for(i in 1:length(t_test_value_pro)){
  k1=t.test(c_pro_all[,i]~c_pro_all[,len_t_test_pro+1])
  t_test_value_pro[i]=abs(k1$statistic)
}



mean_vec_cont_opt=mean_vec_cont_opt[order(t_test_value_pro,decreasing = TRUE)]
mean_vec_cont_opt=mean_vec_cont_opt[1:Npro_for_analysis]

fold_change_vec=fold_change_vec[order(t_test_value_pro,decreasing = TRUE)]
fold_change_vec=fold_change_vec[1:Npro_for_analysis]

pro_pept_listed=1:Npro_cntrl
pro_list=pro_pept_listed[order(t_test_value_pro,decreasing=TRUE)][1:Npro_for_analysis]

Npro=Npro_for_analysis
Npro_trt=Npro
Npro_cntrl=Npro
Npro_cntrl_a=floor(0.85*Npro)
Npro_cntrl_c=Npro-Npro_cntrl_a


Npep_trt=ceiling(1.5*Npro_trt)
Npep_cntrl=ceiling(1.5*Npro_cntrl)
Npep=Npep_trt



sig_matrix_opt=matrix(numeric(Npro*Npro),nrow=Npro)
sig_sq_vec_opt=phi_opt*mean_vec_cont_opt*mean_vec_cont_opt


###########prot_pept_file start#######

#prot_pept_list_cntrl=fasta_file_to_array_conv(pro_cntrl_file,pep_cntrl_file,1)
#prot_pept_list_trt=fasta_file_to_array_conv(pro_trt_file,pep_trt_file,1)


prot_pept_list_cntrl_new1=diff_elem_removal_total(prot_pept_list_cntrl,pro_list)
prot_pept_list_cntrl_new=prot_pept_list_reduction(prot_pept_list_cntrl_new1)

prot_pept_list_cntrl=prot_pept_list_cntrl_new
prot_pept_list_trt=prot_pept_list_cntrl

#print(prot_pept_list_cntrl)
Npep_cntrl=length(prot_pept_list_cntrl)
Npep_trt=Npep_cntrl
Npep=Npep_cntrl


#######Prot_pept_file_end############

for(i in 1:Npro){
  sig_matrix_opt[i,i]=sig_sq_vec_opt[i]
}

#xlj_trt=data.frame(cbind(xlj_trt,rep(1,number_of_samples_train_trt)))
#print("gulshan")
mean_vec_treatement=mean_vec_cont_opt*fold_change_vec
  
#####For threshold generation####
if(1){
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

sig_matrix_rand=matrix(numeric(Npro_cntrl*Npro_cntrl),nrow=Npro_cntrl)
sig_sq_vec_rand=phi_rand*mean_vec_cont_rand*mean_vec_cont_rand
  
 for(i in 1:Npro){ 
    sig_matrix_rand[i,i]=sig_sq_vec_rand[i]
  }
  mean_vec_treatement_rand=mean_vec_cont_rand*fold_change_vec
  #print(mean_vec_treatement_rand)
  c_pro_control_rand=mvrnorm(n=number_of_samples_train_cntrl,mu = mean_vec_cont_rand,Sigma = sig_matrix_rand)
  c_pro_treatement_rand=mvrnorm(n=number_of_samples_test_trt,mu = mean_vec_treatement_rand,Sigma = sig_matrix_rand)
  if(norm_flag){
    
  c_pro_control_rand[,Npro_cntrl_c:length(c_pro_control_rand[1,])]=runif(1,min=0,max=10e-6)*c_pro_control_rand[,Npro_cntrl_c:length(c_pro_control_rand[1,])]
  c_pro_treatement_rand[,Npro_cntrl_c:length(c_pro_treatement_rand[1,])]=runif(1,min=0,max=10e-6)*c_pro_treatement_rand[,Npro_cntrl_c:length(c_pro_treatement_rand[1,])]
}
  
  diff_rand_vec_cntrl=colMeans(c_pro_control_rand)-colMeans(c_pro_control0)
  print(mean_vec_cont_rand)
  thresh_key_cntrl=euc_norm(diff_rand_vec_cntrl)
  
  print("lellina")
  #print(thresh_key_cntrl)
}

#####Peptide generating function#####
  
  

  
peptide_generation1=function(mean_vec_cont_opt,fold_change_vec,phi_opt,Npro_cntrl,Npro_trt,prot_pept_list_cntrl,Npro_for_analysis,
                             prot_pept_list_cntrl_test,prot_pept_list_trt,prot_pept_list_trt_test,train_test_flag,number_of_samples_train_cntrl,number_of_samples_train_trt)
{
  time1=proc.time()
  #print(fold_change_vec)
  #print(number_of_samples_train_cntrl)
  #print("lellina23")
  #print(thresh_key_cntrl)
  mean_vec_treatement=mean_vec_cont_opt*fold_change_vec
  #print("missile")
  #print(mean_vec_treatement)
  c_pro_control=mvrnorm(n=number_of_samples_train_cntrl,mu = mean_vec_cont_opt,Sigma = sig_matrix_opt)
  c_pro_treatement=mvrnorm(n=number_of_samples_train_trt,mu = mean_vec_treatement,Sigma = sig_matrix_opt)
  #print("gokarna")
  #print(length(c_pro_treatement[1,]))
  #print(c_pro_treatement[,Npro_cntrl_c:length(c_pro_treatement[1,])])
  #print("MOOCS")
  if(norm_flag){
    
  c_pro_control[,Npro_cntrl_c:length(c_pro_control[1,])]=runif(1,min=0,max=10e-6)*c_pro_control[,Npro_cntrl_c:length(c_pro_control[1,])]
  c_pro_treatement[,Npro_cntrl_c:length(c_pro_treatement[1,])]=runif(1,min=0,max=10e-6)*c_pro_treatement[,Npro_cntrl_c:length(c_pro_treatement[1,])]
}
    
  c_pep_control=matrix(numeric(Npep_cntrl*number_of_samples_train_cntrl),nrow = number_of_samples_train_cntrl)
  c_pep_treatement=matrix(numeric(Npep_trt*number_of_samples_train_trt),nrow = number_of_samples_train_trt)
  
  max_length_train=max(length(c_pep_control[1,]),length(c_pep_treatement[1,]))
  
  for(j in 1:max_length_train){
    #p3=sample(1:length(c_pro_control[1,]),sample(1:3,1))
    #prot_pept_list_cntrl[[j]]=
    p3=prot_pept_list_cntrl[[j]]
    
    for(i in 1:length(c_pep_control[,1])){
      
      #cat(i,"-",j,"-",p3,"\n")
      for(k in 1:length(p3)){
        c_pep_control[i,j]=c_pep_control[i,j]+c_pro_control[i,p3[k]]
        c_pep_treatement[i,j]=c_pep_treatement[i,j]+c_pro_treatement[i,p3[k]]
      }
    }
    
    
  }
  #print(c_pep_control)
  
  
  
  ###print("grain")
  ####Generating random peptide treatement data###
  checpoint_flag=0
  ###print("am out of this loop")
  if(checpoint_flag){
    print("checkpoint1")
    print(proc.time()-time1)
  }
  
mu_matrix_cntrl=matrix(numeric(Npep_cntrl*number_of_samples_train_cntrl),nrow = number_of_samples_train_cntrl)
mu_matrix_trt=matrix(numeric(Npep_trt*number_of_samples_train_trt),nrow = number_of_samples_train_trt)
  
  
kappa=5
efficiency_vector=runif(Npep,0.1,1)
  
  
 for(i in 1:length(mu_matrix_cntrl[,1])){
    for(j in 1:length(mu_matrix_cntrl[1,])){
      mu_matrix_cntrl[i,j]=c_pep_control[i,j]*kappa*efficiency_vector[j]
      mu_matrix_trt[i,j]=c_pep_treatement[i,j]*kappa*efficiency_vector[j]
      
    }
  } 
  
  
if(checpoint_flag){
  print("checkpoint2")
  print(proc.time()-time1)
 }
###print("ginger")

alpha=0.03
beta=3.6
  
var_vector_noisy_gaussian_cntrl=alpha*(mu_matrix_cntrl^2)+beta*mu_matrix_cntrl
var_vector_noisy_gaussian_trt=alpha*(mu_matrix_trt^2)+beta*mu_matrix_trt
  
total_vector_cntrl=matrix(numeric(Npep_cntrl*number_of_samples_train_cntrl),nrow = number_of_samples_train_cntrl)
total_vector_trt=matrix(numeric(Npep_cntrl*number_of_samples_train_trt),nrow = number_of_samples_train_trt)
  
if(checpoint_flag){
    
  print("checkpoint2-a")
  print(proc.time()-time1)
} 

  p1=matrix(numeric(Npep_cntrl*number_of_samples_train_cntrl),nrow = number_of_samples_train_cntrl)
  p2=p1
  if(0){	
  for(i in 1:length(total_vector_cntrl[,1])){
    for(j in 1:length(total_vector_cntrl[1,])){
      p1[i,j]=mvrnorm(n=1,mu = 0, Sigma = abs(var_vector_noisy_gaussian_cntrl[i,j]))
      p2[i,j]=mvrnorm(n=1,mu = 0, Sigma = abs(var_vector_noisy_gaussian_cntrl[i,j]))
      
      
    }
    
  }
  total_vector_cntrl1=mu_matrix_cntrl+p1
  total_vector_trt1=mu_matrix_trt+p2
  }
  if(checpoint_flag){
    
    print("checkpoint2-a1")
    print(proc.time()-time1)
  }

  
for(i in 1:length(total_vector_cntrl[,1])){
    #print("hello-1")
    #print(length(total_vector_cntrl[1,]))
    #print(length(total_vector_cntrl[,1]))
    
    #print("hello-2")
    
    for(j in 1:length(total_vector_cntrl[1,])){
      #cat(i,j,"\n")
      total_vector_cntrl[i,j]=mu_matrix_cntrl[i,j]+noise_factor_gauss*mvrnorm(n=1,mu = 0, Sigma = abs(var_vector_noisy_gaussian_cntrl[i,j]))+noise_factor_exp*rexp(1,rate = abs(mu_matrix_cntrl[i,j]))
      total_vector_trt[i,j]=mu_matrix_trt[i,j]+noise_factor_gauss*mvrnorm(n=1,mu = 0, Sigma = abs(var_vector_noisy_gaussian_trt[i,j]))+noise_factor_exp*rexp(1,rate = abs(mu_matrix_trt[i,j]))#+

		#total_vector_cntrl[i,j]=mu_matrix_cntrl[i,j]+noise_factor*mvrnorm(n=1,mu = 0, Sigma = abs(var_vector_noisy_gaussian_cntrl[i,j]))
      #total_vector_trt[i,j]=mu_matrix_trt[i,j]+noise_factor*mvrnorm(n=1,mu = 0, Sigma = abs(var_vector_noisy_gaussian_trt[i,j]))
      #total_vector_cntrl[,j]=mu_matrix_cntrl[,j]+noise_factor*mvrnorm(n=1,mu = 0, Sigma = abs(var_vector_noisy_gaussian_cntrl[1,j]))
      #total_vector_trt[,j]=mu_matrix_trt[,j]+noise_factor*mvrnorm(n=1,mu = 0, Sigma = abs(var_vector_noisy_gaussian_trt[1,j]))
      #total_vector_cntrl[i,]=mu_matrix_cntrl[i,]+noise_factor*mvrnorm(n=1,mu = 0, Sigma = abs(var_vector_noisy_gaussian_cntrl[i,1]))
      #total_vector_trt[i,]=mu_matrix_trt[i,]+noise_factor*mvrnorm(n=1,mu = 0, Sigma = abs(var_vector_noisy_gaussian_trt[i,1]))
    }
  }
  if(checpoint_flag){
    
    print("checkpoint2-b")
    print(proc.time()-time1)  
  } 
    #total_vector_cntrl=data.frame(cbind(total_vector_cntrl,numeric(50)))
  #total_vector_trt=data.frame(cbind(total_vector_trt,numeric(50)+1))
  
  #total_vector_cntrl$X31=factor(total_vector_cntrl$X31)
  #total_vector_trt$X31=factor(total_vector_trt$X31)
  
  
  ###To calculate rolled up abundances###
  #lolol
  
  x_pro_control=matrix(numeric(Npro_cntrl*number_of_samples_train_cntrl),nrow = number_of_samples_train_cntrl)
  x_pro_trt=matrix(numeric(Npro_trt*number_of_samples_train_trt),nrow = number_of_samples_train_trt)
  
  
  
  pept_to_prot_cntrl=list()
  pept_to_prot_trt=list()
  
  if(checpoint_flag){
    
    print("checkpoint3")
    print(proc.time()-time1)
  }
    
for(i in 1:length(prot_pept_list_cntrl)){
    for(j in 1:length(prot_pept_list_cntrl[i][[1]])){
      ###print(prot_pept_list_cntrl[i][[1]][j])
      k1=prot_pept_list_cntrl[i][[1]][j]
      ###print(k1)
      pept_to_prot_cntrl[k1][[1]][length(pept_to_prot_cntrl[k1][[1]])+1]=i
      pept_to_prot_trt[k1][[1]][length(pept_to_prot_trt[k1][[1]])+1]=i
      
      
    }
  }
  
  #print("||||||||")
  #print(prot_pept_list_trt)
  #print("||||||||")
  
  
  
  ###print("greek2")
  
  
  
  ###print("greek1")
  
  
  
  ###print("greek3")
  if(checpoint_flag){
    
    print("checkpoint4")
    print(proc.time()-time1)
  }
  xlj_control=matrix(numeric(Npro_cntrl*number_of_samples_train_cntrl),nrow = number_of_samples_train_cntrl)
  xlj_trt=matrix(numeric(Npro_trt*number_of_samples_train_trt),nrow = number_of_samples_train_trt)
  
  
  #print("bozer")
  #print(total_vector_cntrl)
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
  
  pept_to_prot_trt=pept_to_prot_cntrl
  
  ############################
 #print("fine till here")
  
  for(j in 1:length(xlj_control[1,])){ #for all proteins
    for(i in 1:length(xlj_control[,1])) { #for all samples
      for(k in 1:length(pept_to_prot_cntrl[j][[1]])){
        p4=pept_to_prot_cntrl[j][[1]][k]
        #print(p4)
        if(is.null(p4)){
          p4=pept_to_prot_cntrl[1][[1]][1]
        }
        xlj_control[i,j]=xlj_control[i,j]+total_vector_cntrl[i,p4]
        xlj_control[i,j]=xlj_control[i,j]/length(p4)
        
        xlj_trt[i,j]=xlj_trt[i,j]+total_vector_trt[i,p4]
        xlj_trt[i,j]=xlj_trt[i,j]/length(p4)
        
        
        
      }
    }
  }
  # print("bozer")
  #print(xlj_trt)
  
  ###print("greekza")
  ########
  if(checpoint_flag){
    
    print("checkpoint5")
    print(proc.time()-time1)
    
  }
  
  
  
  
  xlj_control=data.frame(xlj_control)
  xlj_control=data.frame(cbind(xlj_control,numeric(number_of_samples_train_cntrl)))
  colnames(xlj_control)[Npro+1]="type"
  #xlj_control$type=as.factor(xlj_control$type)
  
  
  
  xlj_trt=data.frame(xlj_trt)
  xlj_trt=data.frame(cbind(xlj_trt,rep(1,number_of_samples_train_trt)))
  colnames(xlj_trt)[Npro+1]="type"
  
  
  #print("briggy")
  #print(xlj_control_test)
  
  
  
  #xlj_control$"type"=numeric(50)
  #xlj_trt$"type"=(numeric(number_of_samples_train_trt)+1)
  
  if(1){
    
    xlj_control=data.frame(xlj_control)
    xlj_trt=data.frame(xlj_trt)
    
    xlj_all=rbind(xlj_control,xlj_trt)
    xlj_all$type=as.factor(xlj_all$type)
    
    
    
    
    
  }
    
 t_test_value=numeric((length(xlj_all[1,])-1))
  len_t_test=length(t_test_value)
  
  
  xlj_all_orig=xlj_all
  
  
  p11=xlj_all$type
  #p22=xlj_all_test$type
  
  
  
  Npro=Npro_for_analysis
  Npro_trt=Npro
  Npro_cntrl=Npro
  
  #print("kikiki")
  #print(xlj_all_test)
  if(checpoint_flag){
    
    print("checkpoint6")
    print(proc.time()-time1)
  }  
  if(train_test_flag){
    return(xlj_all)}
  
  
}

####xlj_all is the required data as per equation 12###
####xlj_all_test is the required test data####

if(1){
  xlj_all=peptide_generation1(mean_vec_cont_opt,fold_change_vec,phi_opt,Npro_cntrl,Npro_trt,prot_pept_list_cntrl,Npro_for_analysis,
                              prot_pept_list_cntrl_test,prot_pept_list_trt,prot_pept_list_trt_test,train_test_flag=1,number_of_samples_train_cntrl,number_of_samples_train_trt)
  
  xlj_all_test=peptide_generation1(mean_vec_cont_opt,fold_change_vec,phi_opt,Npro_cntrl,Npro_trt,prot_pept_list_cntrl,Npro_for_analysis,
                                   prot_pept_list_cntrl_test,prot_pept_list_trt,prot_pept_list_trt_test,train_test_flag=1,number_of_samples_test_cntrl,number_of_samples_test_trt)
  
  #xlj_all_test= xlj_all_test[sample(nrow(xlj_all_test)),]
}  
	
#1. LDA classifier#



####Add the mean and the gaussian vector####



#SNR=1/(alpha+(beta/alpha))

##print("MCMC algorithm starting")
#######Upamanyu's algorithm-3...The ABC-MCMC-algorithm#####

#####First three steps of ABC-MCMC#####

###gamma_0 is the gamma related to S^(0)_(0) in the paper. This is generated with k_opt and theta_opt
###gamma_dist_0 is the proper gamma distribution pertaining to S_(0) in the paper. Its generated with k_0 and theta_0.

#Step1: Sampling gamma and generating mean vectors##
	
gamma_dist_0=rgamma(100,shape=k_opt,scale=theta_opt)
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
print("here maxhappan")
print(gamma_0)
	
while(norm_cntrl>2*thresh_key_cntrl | norm_trt>2*thresh_key_cntrl){
  #print("i am here")
  # print(norm_cntrl)
  #print(norm_trt)
  #print(0.8*thresh_key_cntrl)
  
  #print(0.8*thresh_key_trt)
  #print("-----------")
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
 mean_vec_cntrl_0_0=mean_vec_cont_opt
  #print("DISHA")
  mean_vec_trt_0_0=mean_vec_cont_opt*fold_change_vec
	
###Step2: Generating protein data
  c_pro_control_0_0=mvrnorm(n=number_of_samples_train_cntrl,mu = mean_vec_cntrl_0_0,Sigma = sig_matrix_0_0)
  c_pro_treatment_0_0=mvrnorm(n=number_of_samples_train_trt,mu = mean_vec_trt_0_0,Sigma = sig_matrix_0_0)
  
  diff_vec_cntrl=colMeans(c_pro_control_0_0)-colMeans(c_pro_control_0)
  diff_vec_trt=colMeans(c_pro_treatment_0_0)-colMeans(c_pro_treatment_0)
  # print("hig")
  norm_cntrl=euc_norm(diff_vec_cntrl)
  #print(norm_cntrl)
  norm_trt=euc_norm(diff_vec_trt)
  #print(norm_trt)
  #print("kig")
  
}
	
	
#####Step 3 is just an if condition which is taken care at the begining
#####Steps 5, 6 and 7 of the markov chain####
c_pro_control_markov=c_pro_control_0_0
c_pro_trt_markov=c_pro_treatment_0_0
gamma_vec=colMeans(c_pro_control_markov)
	
#####Main MCMC chain#####

count_markov=0
markov_iter=10000
gamma_markov=matrix(0L, nrow = markov_iter,ncol=Npro)	

for(j in 1:markov_iter){
  gamma_vec=colMeans(c_pro_control_markov)
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
	
  if(norm_cntrl_markov<0.6*thresh_key_cntrl & norm_trt_markov<0.6*thresh_key_cntrl){
    #mean_ratio=
    gamma_vec=colMeans(c_pro_control_markov)
    ###print(j)
    count_markov=count_markov+1
    
  }
  gamma_markov[j,]=gamma_vec
  
  
  
}
print("MCMC algorithm done")
######Algorithm-3 in upamanyu paper done#####

#####Generation of kernel data and then classification###



#print("running the function")
#print(gamma_vec)	
	
markov_Data=peptide_generation1(mean_vec_cont_opt,fold_change_vec,phi_opt,Npro_cntrl,Npro_trt,prot_pept_list_cntrl,Npro_for_analysis,
                              prot_pept_list_cntrl_test,prot_pept_list_trt,prot_pept_list_trt_test,train_test_flag=1,number_of_samples_train_cntrl,number_of_samples_train_trt)

print("mrinal")
markov_Data_test=peptide_generation1(mean_vec_cont_opt,fold_change_vec,phi_opt,Npro_cntrl,Npro_trt,prot_pept_list_cntrl,Npro_for_analysis,
                                     prot_pept_list_cntrl_test,prot_pept_list_trt,prot_pept_list_trt_test,train_test_flag = 0,number_of_samples_test_cntrl,number_of_samples_test_trt)

print("dingchaka")
gaussian_normal_kernel=mvrnorm(n=1000,mu=numeric(Npro),Sigma = diag(Npro))
gaussian_normal_kernel_density=density(gaussian_normal_kernel,bw=0.5)
xlj_all1=xlj_all
xlj_all_test1=xlj_all_test
#xlj_all_test1=xlj_all_test


	
	
####Designing classifiers####
 xlj_all1=xlj_all
xlj_all_test1=xlj_all_test
if(1){
  lda_classifier_51<-lda(type ~ .,data=xlj_all1)
  predictions_lda_51=predict(lda_classifier_51,xlj_all_test1[,1:Npro])$class
  table_data_lda_51=table(predictions_lda_51,xlj_all_test1[,Npro+1])
  predictions_lda_51_app=predict(lda_classifier_51,xlj_all1[,1:Npro])$class
  table_data_lda_51_app=table(predictions_lda_51_app,xlj_all1[,Npro+1])
  
  print(table_data_lda_51)
  
  training_labels1=xlj_all1$type
  knn_trained1<-knn(train = xlj_all1[,1:Npro]  , test =xlj_all_test1[,1:Npro] , cl = training_labels1, k=3)
  table_data_knn_51=table(knn_trained1,xlj_all_test1[,Npro+1])
  knn_trained_app1<-knn(train = xlj_all1[,1:Npro]  , test =xlj_all1[,1:Npro] , cl = training_labels1, k=3)
  table_data_knn_51_app1=table(knn_trained_app1,xlj_all1[,Npro+1])
  
  print(table_data_knn_51)
  
  
}

abc_mcmc_result_vector=rep(2,length(xlj_all_test[,1]))
abc_mcmc_result_vector1=rep(2,length(xlj_all_test[,1]))
abc_mcmc_result_vector2=rep(2,length(xlj_all_test[,1]))
abc_mcmc_result_vector3=rep(2,length(xlj_all_test[,1]))
abc_mcmc_result_vector4=rep(2,length(xlj_all_test[,1]))
abc_mcmc_result_vector5=rep(2,length(xlj_all_test[,1]))
abc_mcmc_result_vector6=rep(2,length(xlj_all_test[,1]))
abc_mcmc_result_vector7=rep(2,length(xlj_all_test[,1]))
abc_mcmc_result_vector8=rep(2,length(xlj_all_test[,1]))
abc_mcmc_result_vector9=rep(2,length(xlj_all_test[,1]))
abc_mcmc_result_vector10=rep(2,length(xlj_all_test[,1]))
abc_mcmc_result_vector11=rep(2,length(xlj_all_test[,1]))
abc_mcmc_result_vector12=rep(2,length(xlj_all_test[,1]))
abc_mcmc_result_vector13=rep(2,length(xlj_all_test[,1]))
abc_mcmc_result_vector14=rep(2,length(xlj_all_test[,1]))

print("1234")
#print(table_data_lda_50)
#print(table_data_knn_50)
	
total_outside_sum_one_arr=numeric(nrow(xlj_all_test))
total_outside_sum_zero_arr=numeric(nrow(xlj_all_test))

total_outside_sum_one_arr1=numeric(nrow(xlj_all_test))
total_outside_sum_zero_arr1=numeric(nrow(xlj_all_test))

total_outside_sum_one_arr2=numeric(nrow(xlj_all_test))
total_outside_sum_zero_arr2=numeric(nrow(xlj_all_test))

total_outside_sum_one_arr3=numeric(nrow(xlj_all_test))
total_outside_sum_zero_arr3=numeric(nrow(xlj_all_test))

total_outside_sum_one_arr4=numeric(nrow(xlj_all_test))
total_outside_sum_zero_arr4=numeric(nrow(xlj_all_test))

total_outside_sum_one_arr5=numeric(nrow(xlj_all_test))
total_outside_sum_zero_arr5=numeric(nrow(xlj_all_test))

total_outside_sum_one_arr6=numeric(nrow(xlj_all_test))
total_outside_sum_zero_arr6=numeric(nrow(xlj_all_test))

total_outside_sum_one_arr7=numeric(nrow(xlj_all_test))
total_outside_sum_zero_arr7=numeric(nrow(xlj_all_test))

total_outside_sum_one_arr8=numeric(nrow(xlj_all_test))
total_outside_sum_zero_arr8=numeric(nrow(xlj_all_test))

total_outside_sum_one_arr9=numeric(nrow(xlj_all_test))
total_outside_sum_zero_arr9=numeric(nrow(xlj_all_test))

total_outside_sum_one_arr10=numeric(nrow(xlj_all_test))
total_outside_sum_zero_arr10=numeric(nrow(xlj_all_test))

total_outside_sum_one_arr11=numeric(nrow(xlj_all_test))
total_outside_sum_zero_arr11=numeric(nrow(xlj_all_test))

total_outside_sum_one_arr12=numeric(nrow(xlj_all_test))
total_outside_sum_zero_arr12=numeric(nrow(xlj_all_test))

total_outside_sum_one_arr13=numeric(nrow(xlj_all_test))
total_outside_sum_zero_arr13=numeric(nrow(xlj_all_test))

total_outside_sum_one_arr14=numeric(nrow(xlj_all_test))
total_outside_sum_zero_arr14=numeric(nrow(xlj_all_test))



outside_full_sum_one_arr=numeric(nrow(xlj_all_test))
outside_full_sum_zero_arr=numeric(nrow(xlj_all_test))	

	
#######Prediction by the ABC MCMC model #####
if(1){
for(i in 1:nrow(xlj_all_test)){
	test_data=xlj_all_test[i,1:Npro]
    total_outside_sum_zero=0 ####outside Sum for all the markov iterations
    total_outside_sum_one=0
    total_outside_sum_zero=0 ####outside Sum for all the markov iterations
    total_outside_sum_one=0
    outside_full_sum_zero=0
    outside_full_sum_one=0
    
    total_outside_sum_zero1=0 ####outside Sum for all the markov iterations
    total_outside_sum_one1=0
    total_outside_sum_zero2=0 ####outside Sum for all the markov iterations
    total_outside_sum_one2=0
    
    total_outside_sum_zero1=0 ####outside Sum for all the markov iterations
    total_outside_sum_one1=0
    total_outside_sum_zero2=0 ####outside Sum for all the markov iterations
    total_outside_sum_one2=0
    
    total_outside_sum_zero3=0 ####outside Sum for all the markov iterations
    total_outside_sum_one3=0
    
    
    total_outside_sum_zero4=0 ####outside Sum for all the markov iterations
    total_outside_sum_one4=0
    
    total_outside_sum_zero5=0 ####outside Sum for all the markov iterations
    total_outside_sum_one5=0
    
    total_outside_sum_zero6=0 ####outside Sum for all the markov iterations
    total_outside_sum_one6=0
    
    total_outside_sum_zero7=0 ####outside Sum for all the markov iterations
    total_outside_sum_one7=0
    
    total_outside_sum_zero8=0 ####outside Sum for all the markov iterations
    total_outside_sum_one8=0
    
    total_outside_sum_zero9=0 ####outside Sum for all the markov iterations
    total_outside_sum_one9=0
    
    total_outside_sum_zero10=0 ####outside Sum for all the markov iterations
    total_outside_sum_one10=0
    
    total_outside_sum_zero11=0 ####outside Sum for all the markov iterations
    total_outside_sum_one11=0
    
    
    total_outside_sum_zero12=0 ####outside Sum for all the markov iterations
    total_outside_sum_one12=0
    
    total_outside_sum_zero13=0 ####outside Sum for all the markov iterations
    total_outside_sum_one13=0
    
    total_outside_sum_zero14=0 ####outside Sum for all the markov iterations
    total_outside_sum_one14=0
    
    
    
    print(paste("Running the row",i))
    
    ##Only last 70% are considered for the calculations as first 30% are considered as the burn-in stage
    count_123=0
    print("here man")
 ##Only last 70% are considered for the calculations as first 30% are considered as the burn-in stage
    count_123=0
    print("here man")
    
    for(j in (ceiling(0.9985*markov_iter)):(markov_iter-1)){
      # ##print(paste("Markov-chain_number",j))
      #total_outside_sum_zero=0
      #total_outside_sum_one=0
      #ptm <- proc.time()
      
      markov_Data=peptide_generation1(gamma_markov[j,],fold_change_vec,phi_opt,Npro_cntrl,Npro_trt,prot_pept_list_cntrl,Npro_for_analysis,
                                      prot_pept_list_cntrl_test,prot_pept_list_trt,prot_pept_list_trt_test,train_test_flag = 1,number_of_samples_train_cntrl,number_of_samples_train_trt)
      
      ###print(proc.time() - ptm)
      #print("time-1")
      markov_row_length=length(markov_Data[,1]) 
      markov_length_class=markov_row_length/2
      total_inside_sum_zero=0 ###inside sum is for sum of all the n-values of markov_data
      total_inside_sum_one=0
      inside_full_sum_zero=0
      inside_full_sum_one=0
      
      total_inside_sum_zero1=0 ###inside sum is for sum of all the n-values of markov_data
      total_inside_sum_one1=0
      inside_full_sum_zero1=0
      inside_full_sum_one1=0
      
      total_inside_sum_zero2=0 ###inside sum is for sum of all the n-values of markov_data
      total_inside_sum_one2=0
      inside_full_sum_zero2=0
      inside_full_sum_one2=0
      
      total_inside_sum_zero3=0 ###inside sum is for sum of all the n-values of markov_data
      total_inside_sum_one3=0
      inside_full_sum_zero3=0
      inside_full_sum_one3=0
      
      total_inside_sum_zero4=0 ###inside sum is for sum of all the n-values of markov_data
      total_inside_sum_one4=0
      inside_full_sum_zero4=0
      inside_full_sum_one4=0
      
      total_inside_sum_zero5=0 ###inside sum is for sum of all the n-values of markov_data
      total_inside_sum_one5=0
      inside_full_sum_zero5=0
      inside_full_sum_one5=0
      
      total_inside_sum_zero6=0 ###inside sum is for sum of all the n-values of markov_data
      total_inside_sum_one6=0
      inside_full_sum_zero6=0
      inside_full_sum_one6=0
      
      total_inside_sum_zero7=0 ###inside sum is for sum of all the n-values of markov_data
      total_inside_sum_one7=0
      inside_full_sum_zero7=0
      inside_full_sum_one7=0
      
      total_inside_sum_zero8=0 ###inside sum is for sum of all the n-values of markov_data
      total_inside_sum_one8=0
      inside_full_sum_zero8=0
      inside_full_sum_one8=0
      
      total_inside_sum_zero9=0 ###inside sum is for sum of all the n-values of markov_data
      total_inside_sum_one9=0
      inside_full_sum_zero9=0
      inside_full_sum_one9=0
      
      total_inside_sum_zero10=0 ###inside sum is for sum of all the n-values of markov_data
      total_inside_sum_one10=0
      inside_full_sum_zero10=0
      inside_full_sum_one10=0
      
      total_inside_sum_zero11=0 ###inside sum is for sum of all the n-values of markov_data
      total_inside_sum_one11=0
      inside_full_sum_zero11=0
      inside_full_sum_one11=0
      
      total_inside_sum_zero12=0 ###inside sum is for sum of all the n-values of markov_data
      total_inside_sum_one12=0
      inside_full_sum_zero12=0
      inside_full_sum_one12=0
      
      total_inside_sum_zero13=0 ###inside sum is for sum of all the n-values of markov_data
      total_inside_sum_one13=0
      inside_full_sum_zero13=0
      inside_full_sum_one13=0
      
      total_inside_sum_zero14=0 ###inside sum is for sum of all the n-values of markov_data
      total_inside_sum_one14=0
      inside_full_sum_zero14=0
      inside_full_sum_one14=0
	
  
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
        
        
        ##print("#########################")
        if(1){
          #diff_data_zero=diff_data_zero/abs(max(diff_data_zero))
          #print(diff_data_zero)
          k0=diff_data_zero
          k00=k0
          k0=k00/4000
          diff_data_zero=k0/12000
          diff_data_zero1=k0/14000
          diff_data_zero2=k0/16000
          diff_data_zero3=k0/18000
          diff_data_zero4=k0/20000
          diff_data_zero5=k0/22000
          diff_data_zero6=k0/24000
          diff_data_zero7=k0/26000
          diff_data_zero8=k0/28000
          diff_data_zero9=k0/30000
          diff_data_zero10=k0/32000
          diff_data_zero11=k0/34000
          diff_data_zero12=k0/36000
          diff_data_zero13=k0/38000
          diff_data_zero14=k0/40000
          
          
          kernel_data=approx(gaussian_normal_kernel_density$x,gaussian_normal_kernel_density$y,xout = diff_data_zero)
          kernel_data$y[is.na(kernel_data$y)]=0
          total_inside_sum_zero=total_inside_sum_zero+mean(kernel_data$y)
          
          kernel_data1=approx(gaussian_normal_kernel_density$x,gaussian_normal_kernel_density$y,xout = diff_data_zero1)
          kernel_data1$y[is.na(kernel_data1$y)]=0
          total_inside_sum_zero1=total_inside_sum_zero1+mean(kernel_data1$y)
          
          kernel_data2=approx(gaussian_normal_kernel_density$x,gaussian_normal_kernel_density$y,xout = diff_data_zero2)
          kernel_data2$y[is.na(kernel_data2$y)]=0
          total_inside_sum_zero2=total_inside_sum_zero2+mean(kernel_data2$y)
          
          kernel_data3=approx(gaussian_normal_kernel_density$x,gaussian_normal_kernel_density$y,xout = diff_data_zero3)
          kernel_data3$y[is.na(kernel_data3$y)]=0
          total_inside_sum_zero3=total_inside_sum_zero3+mean(kernel_data3$y)
          
          kernel_data4=approx(gaussian_normal_kernel_density$x,gaussian_normal_kernel_density$y,xout = diff_data_zero4)
          kernel_data4$y[is.na(kernel_data4$y)]=0
          total_inside_sum_zero4=total_inside_sum_zero4+mean(kernel_data4$y)
          
          kernel_data5=approx(gaussian_normal_kernel_density$x,gaussian_normal_kernel_density$y,xout = diff_data_zero5)
          kernel_data5$y[is.na(kernel_data5$y)]=0
          total_inside_sum_zero5=total_inside_sum_zero5+mean(kernel_data5$y)
          
          kernel_data6=approx(gaussian_normal_kernel_density$x,gaussian_normal_kernel_density$y,xout = diff_data_zero6)
          kernel_data6$y[is.na(kernel_data6$y)]=0
          total_inside_sum_zero6=total_inside_sum_zero6+mean(kernel_data6$y)
          
          kernel_data7=approx(gaussian_normal_kernel_density$x,gaussian_normal_kernel_density$y,xout = diff_data_zero7)
          kernel_data7$y[is.na(kernel_data7$y)]=0
          total_inside_sum_zero7=total_inside_sum_zero7+mean(kernel_data7$y)
          
          kernel_data8=approx(gaussian_normal_kernel_density$x,gaussian_normal_kernel_density$y,xout = diff_data_zero8)
          kernel_data8$y[is.na(kernel_data8$y)]=0
          total_inside_sum_zero8=total_inside_sum_zero8+mean(kernel_data8$y)
          
          kernel_data9=approx(gaussian_normal_kernel_density$x,gaussian_normal_kernel_density$y,xout = diff_data_zero9)
          kernel_data9$y[is.na(kernel_data9$y)]=0
          total_inside_sum_zero9=total_inside_sum_zero9+mean(kernel_data9$y)
          
          
          kernel_data10=approx(gaussian_normal_kernel_density$x,gaussian_normal_kernel_density$y,xout = diff_data_zero10)
          kernel_data10$y[is.na(kernel_data10$y)]=0
          total_inside_sum_zero10=total_inside_sum_zero10+mean(kernel_data10$y)
          
          kernel_data11=approx(gaussian_normal_kernel_density$x,gaussian_normal_kernel_density$y,xout = diff_data_zero11)
          kernel_data11$y[is.na(kernel_data11$y)]=0
          total_inside_sum_zero11=total_inside_sum_zero11+mean(kernel_data11$y)
          
          kernel_data12=approx(gaussian_normal_kernel_density$x,gaussian_normal_kernel_density$y,xout = diff_data_zero12)
          kernel_data12$y[is.na(kernel_data12$y)]=0
          total_inside_sum_zero12=total_inside_sum_zero12+mean(kernel_data12$y)
          
          kernel_data13=approx(gaussian_normal_kernel_density$x,gaussian_normal_kernel_density$y,xout = diff_data_zero13)
          kernel_data13$y[is.na(kernel_data13$y)]=0
          total_inside_sum_zero13=total_inside_sum_zero13+mean(kernel_data13$y)
          
          kernel_data14=approx(gaussian_normal_kernel_density$x,gaussian_normal_kernel_density$y,xout = diff_data_zero14)
          kernel_data14$y[is.na(kernel_data14$y)]=0
          total_inside_sum_zero14=total_inside_sum_zero14+mean(kernel_data14$y)
          
          
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
          # print(diff_data_one)
          k1=diff_data_one
          k11=k1
          k1=k11/4000
          diff_data_one=k1/12000
          diff_data_one1=k1/14000
          diff_data_one2=k1/16000
          diff_data_one3=k1/18000
          diff_data_one4=k1/20000
          diff_data_one5=k1/22000
          diff_data_one6=k1/24000
          diff_data_one7=k1/26000
          diff_data_one8=k1/28000
          diff_data_one9=k1/30000
          diff_data_one10=k1/32000
          diff_data_one11=k1/34000
          diff_data_one12=k1/36000
          diff_data_one13=k1/38000
          diff_data_one14=k1/40000
          
          #print(diff_data_one)
          
          kernel_data=approx(gaussian_normal_kernel_density$x,gaussian_normal_kernel_density$y,xout = diff_data_one)
          kernel_data$y[is.na(kernel_data$y)]=0
          total_inside_sum_one=total_inside_sum_one+mean(kernel_data$y)
          
          kernel_data1=approx(gaussian_normal_kernel_density$x,gaussian_normal_kernel_density$y,xout = diff_data_one1)
          kernel_data1$y[is.na(kernel_data1$y)]=0
          total_inside_sum_one1=total_inside_sum_one1+mean(kernel_data1$y)
          
          kernel_data2=approx(gaussian_normal_kernel_density$x,gaussian_normal_kernel_density$y,xout = diff_data_one2)
          kernel_data2$y[is.na(kernel_data2$y)]=0
          total_inside_sum_one2=total_inside_sum_one2+mean(kernel_data2$y)
          
          kernel_data3=approx(gaussian_normal_kernel_density$x,gaussian_normal_kernel_density$y,xout = diff_data_one3)
          kernel_data3$y[is.na(kernel_data3$y)]=0
          total_inside_sum_one3=total_inside_sum_one3+mean(kernel_data3$y)
          
          kernel_data4=approx(gaussian_normal_kernel_density$x,gaussian_normal_kernel_density$y,xout = diff_data_one4)
          kernel_data4$y[is.na(kernel_data4$y)]=0
          total_inside_sum_one4=total_inside_sum_one4+mean(kernel_data4$y)
          
          kernel_data5=approx(gaussian_normal_kernel_density$x,gaussian_normal_kernel_density$y,xout = diff_data_one5)
          kernel_data5$y[is.na(kernel_data5$y)]=0
          total_inside_sum_one5=total_inside_sum_one5+mean(kernel_data5$y)
          
          kernel_data6=approx(gaussian_normal_kernel_density$x,gaussian_normal_kernel_density$y,xout = diff_data_one6)
          kernel_data6$y[is.na(kernel_data6$y)]=0
          total_inside_sum_one6=total_inside_sum_one6+mean(kernel_data6$y)
          
          kernel_data7=approx(gaussian_normal_kernel_density$x,gaussian_normal_kernel_density$y,xout = diff_data_one7)
          kernel_data7$y[is.na(kernel_data7$y)]=0
          total_inside_sum_one7=total_inside_sum_one7+mean(kernel_data7$y)
          
          
          kernel_data8=approx(gaussian_normal_kernel_density$x,gaussian_normal_kernel_density$y,xout = diff_data_one8)
          kernel_data8$y[is.na(kernel_data8$y)]=0
          total_inside_sum_one8=total_inside_sum_one8+mean(kernel_data8$y)
          
          kernel_data9=approx(gaussian_normal_kernel_density$x,gaussian_normal_kernel_density$y,xout = diff_data_one9)
          kernel_data9$y[is.na(kernel_data9$y)]=0
          total_inside_sum_one9=total_inside_sum_one9+mean(kernel_data9$y)
          
          
          kernel_data10=approx(gaussian_normal_kernel_density$x,gaussian_normal_kernel_density$y,xout = diff_data_one10)
          kernel_data10$y[is.na(kernel_data10$y)]=0
          total_inside_sum_one10=total_inside_sum_one10+mean(kernel_data10$y)
          
          kernel_data11=approx(gaussian_normal_kernel_density$x,gaussian_normal_kernel_density$y,xout = diff_data_one11)
          kernel_data11$y[is.na(kernel_data11$y)]=0
          total_inside_sum_one11=total_inside_sum_one11+mean(kernel_data11$y)
          
          kernel_data12=approx(gaussian_normal_kernel_density$x,gaussian_normal_kernel_density$y,xout = diff_data_one12)
          kernel_data12$y[is.na(kernel_data12$y)]=0
          total_inside_sum_one12=total_inside_sum_one12+mean(kernel_data12$y)
          
          kernel_data13=approx(gaussian_normal_kernel_density$x,gaussian_normal_kernel_density$y,xout = diff_data_one13)
          kernel_data13$y[is.na(kernel_data13$y)]=0
          total_inside_sum_one13=total_inside_sum_one13+mean(kernel_data13$y)
          
          kernel_data14=approx(gaussian_normal_kernel_density$x,gaussian_normal_kernel_density$y,xout = diff_data_one14)
          kernel_data14$y[is.na(kernel_data14$y)]=0
          total_inside_sum_one14=total_inside_sum_one14+mean(kernel_data14$y)
          
          
        }
        k22=cbind(diff_data_one,diff_data_zero)
        # plot(colMeans(k22))
        
        if(is.na(total_inside_sum_one)& count_123==0){
          ###print(kernel_data$y)
          count_123=1
          p11=t1
        }
	 }
      #total_outside_sum_zero=total_outside_sum_zero+total_inside_sum_zero
      #total_outside_sum_one=total_outside_sum_one+total_inside_sum_one
      ##print(proc.time() - ptm)
      #print("time 1-a")
    }
	
total_outside_sum_zero=total_outside_sum_zero+total_inside_sum_zero
    total_outside_sum_one=total_outside_sum_one+total_inside_sum_one
    
    total_outside_sum_zero1=total_outside_sum_zero1+total_inside_sum_zero1
    total_outside_sum_one1=total_outside_sum_one1+total_inside_sum_one1
    
    total_outside_sum_zero2=total_outside_sum_zero2+total_inside_sum_zero2
    total_outside_sum_one2=total_outside_sum_one2+total_inside_sum_one2
    
    total_outside_sum_zero3=total_outside_sum_zero3+total_inside_sum_zero3
    total_outside_sum_one3=total_outside_sum_one3+total_inside_sum_one3
    
    total_outside_sum_zero4=total_outside_sum_zero4+total_inside_sum_zero4
    total_outside_sum_one4=total_outside_sum_one4+total_inside_sum_one4
    
    total_outside_sum_zero5=total_outside_sum_zero5+total_inside_sum_zero5
    total_outside_sum_one5=total_outside_sum_one5+total_inside_sum_one5
    
    total_outside_sum_zero6=total_outside_sum_zero6+total_inside_sum_zero6
    total_outside_sum_one6=total_outside_sum_one6+total_inside_sum_one6
    
    total_outside_sum_zero7=total_outside_sum_zero7+total_inside_sum_zero7
    total_outside_sum_one7=total_outside_sum_one7+total_inside_sum_one7
    
    total_outside_sum_zero8=total_outside_sum_zero8+total_inside_sum_zero8
    total_outside_sum_one8=total_outside_sum_one8+total_inside_sum_one8
    
    total_outside_sum_zero9=total_outside_sum_zero9+total_inside_sum_zero9
    total_outside_sum_one9=total_outside_sum_one9+total_inside_sum_one9
    
    total_outside_sum_zero10=total_outside_sum_zero10+total_inside_sum_zero10
    total_outside_sum_one10=total_outside_sum_one10+total_inside_sum_one10
    
    total_outside_sum_zero11=total_outside_sum_zero11+total_inside_sum_zero11
    total_outside_sum_one11=total_outside_sum_one11+total_inside_sum_one11
    
    total_outside_sum_zero12=total_outside_sum_zero12+total_inside_sum_zero12
    total_outside_sum_one12=total_outside_sum_one12+total_inside_sum_one12
    
    total_outside_sum_zero13=total_outside_sum_zero13+total_inside_sum_zero13
    total_outside_sum_one13=total_outside_sum_one13+total_inside_sum_one13
    
    total_outside_sum_zero14=total_outside_sum_zero14+total_inside_sum_zero14
    total_outside_sum_one14=total_outside_sum_one14+total_inside_sum_one14
    
    
    
    outside_full_sum_zero=outside_full_sum_zero+inside_full_sum_zero
    outside_full_sum_one=outside_full_sum_one+inside_full_sum_one
    
    ###print("glagla")
    ###print(total_outside_sum_one)
    ###print(total_outside_sum_zero)
    total_outside_sum_one_arr[i]=total_outside_sum_one
    total_outside_sum_zero_arr[i]=total_outside_sum_zero
    
    total_outside_sum_one_arr1[i]=total_outside_sum_one1
    total_outside_sum_zero_arr1[i]=total_outside_sum_zero1
    
    total_outside_sum_one_arr2[i]=total_outside_sum_one2
    total_outside_sum_zero_arr2[i]=total_outside_sum_zero2
    
    total_outside_sum_one_arr3[i]=total_outside_sum_one3
    total_outside_sum_zero_arr3[i]=total_outside_sum_zero3
    
    total_outside_sum_one_arr4[i]=total_outside_sum_one4
    total_outside_sum_zero_arr4[i]=total_outside_sum_zero4
    
    total_outside_sum_one_arr5[i]=total_outside_sum_one5
    total_outside_sum_zero_arr5[i]=total_outside_sum_zero5
    
    total_outside_sum_one_arr6[i]=total_outside_sum_one6
    total_outside_sum_zero_arr6[i]=total_outside_sum_zero6
    
    total_outside_sum_one_arr7[i]=total_outside_sum_one7
    total_outside_sum_zero_arr7[i]=total_outside_sum_zero7
    
    total_outside_sum_one_arr8[i]=total_outside_sum_one8
    total_outside_sum_zero_arr8[i]=total_outside_sum_zero8
    
    total_outside_sum_one_arr9[i]=total_outside_sum_one9
    total_outside_sum_zero_arr9[i]=total_outside_sum_zero9
    
    total_outside_sum_one_arr10[i]=total_outside_sum_one10
    total_outside_sum_zero_arr10[i]=total_outside_sum_zero10
    
    total_outside_sum_one_arr11[i]=total_outside_sum_one11
    total_outside_sum_zero_arr11[i]=total_outside_sum_zero11
    
    total_outside_sum_one_arr12[i]=total_outside_sum_one12
    total_outside_sum_zero_arr12[i]=total_outside_sum_zero12
    
    total_outside_sum_one_arr13[i]=total_outside_sum_one13
    total_outside_sum_zero_arr13[i]=total_outside_sum_zero13
    
    total_outside_sum_one_arr14[i]=total_outside_sum_one14
    total_outside_sum_zero_arr14[i]=total_outside_sum_zero14
    
    outside_full_sum_one_arr[i]=outside_full_sum_one
    
    outside_full_sum_zero_arr[i]=outside_full_sum_zero
    
    if(total_outside_sum_one<total_outside_sum_zero){
      #if(abs(outside_full_sum_one)>abs(outside_full_sum_zero)){
      abc_mcmc_result_vector[i]=0
    }
    else{
      abc_mcmc_result_vector[i]=1
    }
    
    if(total_outside_sum_one1<total_outside_sum_zero1){
      #if(abs(outside_full_sum_one)>abs(outside_full_sum_zero)){
      abc_mcmc_result_vector1[i]=0
    }
    else{
      abc_mcmc_result_vector1[i]=1
    }
    
    if(total_outside_sum_one2<total_outside_sum_zero2){
      #if(abs(outside_full_sum_one)>abs(outside_full_sum_zero)){
      abc_mcmc_result_vector2[i]=0
    }
    else{
      abc_mcmc_result_vector2[i]=1
    }
    
    if(total_outside_sum_one3<total_outside_sum_zero3){
      #if(abs(outside_full_sum_one)>abs(outside_full_sum_zero)){
      abc_mcmc_result_vector3[i]=0
    }
    else{
      abc_mcmc_result_vector3[i]=1
    }
    
    if(total_outside_sum_one4<total_outside_sum_zero4){
      #if(abs(outside_full_sum_one)>abs(outside_full_sum_zero)){
      abc_mcmc_result_vector4[i]=0
    }
    else{
      abc_mcmc_result_vector4[i]=1
    }
    
    if(total_outside_sum_one5<total_outside_sum_zero5){
      #if(abs(outside_full_sum_one)>abs(outside_full_sum_zero)){
      abc_mcmc_result_vector5[i]=0
    }
    else{
      abc_mcmc_result_vector5[i]=1
    }
    
    if(total_outside_sum_one6<total_outside_sum_zero6){
      #if(abs(outside_full_sum_one)>abs(outside_full_sum_zero)){
      abc_mcmc_result_vector6[i]=0
    }
    else{
      abc_mcmc_result_vector6[i]=1
    }
    
    if(total_outside_sum_one7<total_outside_sum_zero7){
      #if(abs(outside_full_sum_one)>abs(outside_full_sum_zero)){
      abc_mcmc_result_vector7[i]=0
    }
    else{
      abc_mcmc_result_vector7[i]=1
    }
    
    if(total_outside_sum_one8<total_outside_sum_zero8){
      #if(abs(outside_full_sum_one)>abs(outside_full_sum_zero)){
      abc_mcmc_result_vector8[i]=0
    }
    else{
      abc_mcmc_result_vector8[i]=1
    }
    
    if(total_outside_sum_one9<total_outside_sum_zero9){
      #if(abs(outside_full_sum_one)>abs(outside_full_sum_zero)){
      abc_mcmc_result_vector9[i]=0
    }
    else{
      abc_mcmc_result_vector9[i]=1
    }
    
    if(total_outside_sum_one10<total_outside_sum_zero10){
      #if(abs(outside_full_sum_one)>abs(outside_full_sum_zero)){
      abc_mcmc_result_vector10[i]=0
    }
    else{
      abc_mcmc_result_vector10[i]=1
    }
    
    if(total_outside_sum_one11<total_outside_sum_zero11){
      #if(abs(outside_full_sum_one)>abs(outside_full_sum_zero)){
      abc_mcmc_result_vector11[i]=0
    }
    else{
      abc_mcmc_result_vector11[i]=1
    }
    
    if(total_outside_sum_one12<total_outside_sum_zero12){
      #if(abs(outside_full_sum_one)>abs(outside_full_sum_zero)){
      abc_mcmc_result_vector12[i]=0
    }
    else{
      abc_mcmc_result_vector12[i]=1
    }
    
    if(total_outside_sum_one13<total_outside_sum_zero13){
      #if(abs(outside_full_sum_one)>abs(outside_full_sum_zero)){
      abc_mcmc_result_vector13[i]=0
    }
    else{
      abc_mcmc_result_vector13[i]=1
    }
    
    if(total_outside_sum_one14<total_outside_sum_zero14){
      #if(abs(outside_full_sum_one)>abs(outside_full_sum_zero)){
      abc_mcmc_result_vector14[i]=0
    }
    else{
      abc_mcmc_result_vector14[i]=1
    }
    
    #print("time-2")
    ###print(proc.time() - ptm)
    
    print(abc_mcmc_result_vector)
    print("-----------1---------")
    print(abc_mcmc_result_vector1)
    print("-----------2---------")
    print(abc_mcmc_result_vector2)
    print("-----------3---------")
    print(abc_mcmc_result_vector3)
    print("-----------4---------")
    print(abc_mcmc_result_vector4)
    print("-----------5---------")
    print(abc_mcmc_result_vector5)
    print("-----------6---------")
    print(abc_mcmc_result_vector6)
    print("-----------7---------")
    print(abc_mcmc_result_vector7)
    print("-----------8---------")
    print(abc_mcmc_result_vector8)
    print("-----------9---------")
    print(abc_mcmc_result_vector9)
    print("-----------10---------")
    print(abc_mcmc_result_vector10)
    print("-----------11---------")
    print(abc_mcmc_result_vector11)
    print("-----------12---------")
    print(abc_mcmc_result_vector12)
    print("-----------13---------")
    print(abc_mcmc_result_vector13)
    print("-----------14---------")
    print(abc_mcmc_result_vector14)
    
  }
  #print("------------------------")
  #print("Fullll ABC-MCMC model done")
  
  #print("------------------------")
  
  table_data_abc_mcmc_0=table(abc_mcmc_result_vector,xlj_all_test[,Npro+1])
  table_data_abc_mcmc_1=table(abc_mcmc_result_vector1,xlj_all_test[,Npro+1])
  table_data_abc_mcmc_2=table(abc_mcmc_result_vector2,xlj_all_test[,Npro+1])
  table_data_abc_mcmc_3=table(abc_mcmc_result_vector3,xlj_all_test[,Npro+1])
  table_data_abc_mcmc_4=table(abc_mcmc_result_vector4,xlj_all_test[,Npro+1])
  table_data_abc_mcmc_5=table(abc_mcmc_result_vector5,xlj_all_test[,Npro+1])
  table_data_abc_mcmc_6=table(abc_mcmc_result_vector6,xlj_all_test[,Npro+1])
  table_data_abc_mcmc_7=table(abc_mcmc_result_vector7,xlj_all_test[,Npro+1])
  table_data_abc_mcmc_8=table(abc_mcmc_result_vector8,xlj_all_test[,Npro+1])
  table_data_abc_mcmc_9=table(abc_mcmc_result_vector9,xlj_all_test[,Npro+1])
  table_data_abc_mcmc_10=table(abc_mcmc_result_vector10,xlj_all_test[,Npro+1])
  table_data_abc_mcmc_11=table(abc_mcmc_result_vector11,xlj_all_test[,Npro+1])
  table_data_abc_mcmc_12=table(abc_mcmc_result_vector12,xlj_all_test[,Npro+1])
  table_data_abc_mcmc_13=table(abc_mcmc_result_vector13,xlj_all_test[,Npro+1])
  table_data_abc_mcmc_14=table(abc_mcmc_result_vector14,xlj_all_test[,Npro+1])
  
  
  
}
final_time=proc.time()
#print("K balahander")
#print(table_data_abc_mcmc_50)
abc_mcmc_result_vector_final1=rbind(abc_mcmc_result_vector,abc_mcmc_result_vector1,abc_mcmc_result_vector2,abc_mcmc_result_vector3)
abc_mcmc_result_vector_final2=rbind(abc_mcmc_result_vector4,abc_mcmc_result_vector5,abc_mcmc_result_vector6,abc_mcmc_result_vector7)
abc_mcmc_result_vector_final3=rbind(abc_mcmc_result_vector8,abc_mcmc_result_vector9,abc_mcmc_result_vector10,abc_mcmc_result_vector11)
abc_mcmc_result_vector_final4=rbind(abc_mcmc_result_vector12,abc_mcmc_result_vector13,abc_mcmc_result_vector14)

abc_mcmc_result_vector_final=rbind(abc_mcmc_result_vector_final1,abc_mcmc_result_vector_final2,abc_mcmc_result_vector_final3,abc_mcmc_result_vector_final4)

max_vector_return=function(abc_mcmc_result_vector123){
  nrows=nrow(abc_mcmc_result_vector123)
  abc_mcmc_result_vector_total=rep(2,ncol(abc_mcmc_result_vector123))
  for(i in 1:ncol(abc_mcmc_result_vector123)){
    print(sum(abc_mcmc_result_vector123[,i]))
    if(sum(abc_mcmc_result_vector123[,i])>nrows/2)
    {
      abc_mcmc_result_vector_total[i]=1
    }
    else{
      abc_mcmc_result_vector_total[i]=0
      
    }
  }
  return(abc_mcmc_result_vector_total)
}
abc_mcmc_result_vector_total=max_vector_return(abc_mcmc_result_vector_final)
table_data_abc_mcmc_total=table(abc_mcmc_result_vector_total,xlj_all_test[,Npro+1])
print("------All tables-----")
print(table_data_lda_51)
print(table_data_knn_51)
print(table_data_abc_mcmc_total)

print("%%%%%")
print(proc.time()-time2)

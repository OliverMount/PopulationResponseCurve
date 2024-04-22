
# Combining to a single data frame of preferred direction of neurons and their p-val (for whether it is a tuned to one direction)  
rm(list = ls())

library(R.matlab)
library(pracma)
library(tidyverse)
library(gridExtra)
library(reshape2)
library(ggpubr)

prefDir_path='/media/olive/Research/oliver/prefDir'  # To save
pval_path='/media/olive/Research/oliver/pvals'
save_path='/media/olive/Research/oliver/IEMdecodingForCalciumData/neuron_counts/' 

conds<- c('V1_45','V1_90','V1_135','PPC_45','PPC_90','PPC_135')
percent<- c(10,20,40,60,100)
pval_threshold<-0.05
nos<-8 # number of directions

paradigm<-'task'
setwd(file.path(pval_path,paradigm))

conds=list.files(pattern = '.mat')

 # run this only once to get the  
for (cond in conds){   # for each condition
  cond_name<- gsub("\\.mat", "", cond)
   
  # loading the p-value file
  B<- readMat(file.path(pval_path,paradigm,cond))
  #B$homo<- Filter(Negate(is.null), B$homo)  
  #B$hetero<- Filter(Negate(is.null),B$hetero)
  
  L<- length(B$homo)  # No. of animals here
  
   
  df<-data.frame(Sub=rep(NA,1),
                 Condition=rep(NA,1),
                 Group=rep(NA,1),
                 Pvalue=rep(NA,1),
                 Preference=rep(NA,1))
  
  cat('Cond name  : ', cond_name, '\n')
  for (k in 1:L){  # for each animal   
      #ho_p<- as.numeric(unlist(B$homo[[k]][[1]][[1]]))
      ho_p<- as.numeric(unlist(B$homo[[k]]))
      Lf_homo<- length(ho_p)
      
      
      #he_p<-  as.numeric(unlist(B$hetero[[k]][[1]][[1]]))
      he_p<-  as.numeric(unlist(B$hetero[[k]]))
      Lf_hetero<- length(he_p) 
      
      
      cat(' Homo:  ' , Lf_homo , ' Hetero: ' ,Lf_hetero, ' \n' )
      
      Lf<-min(c(Lf_homo,Lf_hetero))
      ho_p<-ho_p[1:Lf]
      he_p<-he_p[1:Lf] 
      
      
      homo<-hetero<- rep(0,Lf) 
      
      # For homo
      
      temp<-data.frame(Sub=rep(paste0("Animal.",k),Lf),
                       Condition=rep(cond_name,Lf),
                       Group=rep("homo",Lf),
                       Pvalue=ho_p,
                       Preference=homo) 
      
      df<-rbind(df,temp)  
      
      # for hetero 
      temp<-data.frame(Sub=rep(paste0("Animal.",k),Lf),
                       Condition=rep(cond_name,Lf),
                       Group=rep("hetero",Lf),
                       Pvalue=he_p,
                       Preference=hetero) 
      
      df<-rbind(df,temp) 
      df<- na.omit(df) 
  } 
  row.names(df)<- 1:nrow(df)
  write.csv(df,file=file.path(prefDir_path,paradigm,paste0(cond_name,'.csv')))
  
} 


########  passive data ########
paradigm<-'passive'
setwd(file.path(pval_path,paradigm))

conds=list.files(pattern = '.mat')

# run this only once to get the  
for (cond in conds){   # for each condition
  cond_name<- gsub("\\.mat", "", cond)
  cat('Cond name  : ', cond_name, '\n')
  # loading the p-value file
  B<- readMat(cond)
  B$homo<- Filter(Negate(is.null), B$homo)
  B$hetero<- Filter(Negate(is.null),B$hetero) 
  
  L<- length(B$homo)  # No. of animals here
  
  df<-data.frame(Sub=rep(NA,1),
                 Condition=rep(NA,1),
                 Group=rep(NA,1),
                 Pvalue=rep(NA,1),
                 Preference=rep(NA,1))
  
  for (k in 1:L){  # for each animal   
    
    ho_p<- as.numeric(unlist(B$homo[[k]][[1]][[1]]))
    Lf_homo<- length(ho_p) 
    
    he_p<-  as.numeric(unlist(B$hetero[[k]][[1]][[1]]))
    Lf_hetero<- length(he_p) 
    
    cat(' Homo:  ' , Lf_homo , ' Hetero: ' ,Lf_hetero, ' \n' )
     
    Lf<-min(c(Lf_homo,Lf_hetero))
    ho_p<-ho_p[1:Lf]
    he_p<-he_p[1:Lf] 
    
    
    homo<-hetero<- rep(0,Lf) 
    
    # For homo
    
    temp<-data.frame(Sub=rep(paste0("Animal.",k),Lf),
                     Condition=rep(cond_name,Lf),
                     Group=rep("homo",Lf),
                     Pvalue=ho_p,
                     Preference=homo) 
    
    df<-rbind(df,temp)  
    
    # for hetero 
    temp<-data.frame(Sub=rep(paste0("Animal.",k),Lf),
                     Condition=rep(cond_name,Lf),
                     Group=rep("hetero",Lf),
                     Pvalue=he_p,
                     Preference=hetero) 
    
    df<-rbind(df,temp) 
    df<- na.omit(df) 
    
  }
  
  row.names(df)<- 1:nrow(df)
  write.csv(df,file=file.path(prefDir_path,paradigm,paste0(cond_name,'.csv')))
} 

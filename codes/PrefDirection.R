# Combining to a single data frame of preferred direction of neurons and their p-val (for whether it is a tuned to one direction)  
rm(list = ls())

library(R.matlab)
library(pracma)
library(tidyverse)
library(gridExtra)
library(reshape2)
library(ggpubr)


base_path='/media/olive/Research/oliver/prefDir' 
pval_path='/media/olive/Research/oliver/pvals'
save_path='/media/olive/Research/oliver/IEMdecodingForCalciumData/neuron_counts/' 

conds<- c('V1_45','V1_90','V1_135','PPC_45','PPC_90','PPC_135')
percent<- c(10,20,40,60,100)
pval_threshold<-0.05
nos<-8 # number of directions

paradigm<-'task'
setwd(file.path(base_path,paradigm))

conds=list.files(pattern = '.mat')

 # run this only once to get the  
for (cond in conds){   # for each condition
  cond_name<- gsub("\\.mat", "", cond)
  
  # loading preferred direction  file
  A<-readMat(cond)    
  A$prefDir.homo<- Filter(Negate(is.null), A$prefDir.homo)  # remove the empty list values
  L<- length(A$prefDir.homo) 
  A$prefDir.hetero<- Filter(Negate(is.null), A$prefDir.hetero)
  
  # loading the p-value file
  B<- readMat(file.path(pval_path,paradigm,cond))
  
   
  df<-data.frame(Sub=rep(NA,1),Condition=rep(NA,1),Group=rep(NA,1),Pvalue=rep(NA,1),Preference=rep(NA,1))
  
  for (k in 1:L){  # for each animal   
      homo<- as.numeric(unlist(A$prefDir.homo[[k]][[1]][[1]]))
      hetero<-  as.numeric(unlist(A$prefDir.hetero[[k]][[1]][[1]]))
      
      ho_p<-as.numeric(unlist(B$homo[[k]]))
      he_p<-as.numeric(unlist(B$hetero[[k]]))
      
      Lf<- length(homo)
      
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
  write.csv(df,file=file.path(base_path,paradigm,paste0(cond_name,'_prefer.csv')))
  
} 


########  passive data ########
paradigm<-'passive'
setwd(file.path(pval_path,paradigm))

conds=list.files(pattern = '.mat')

# run this only once to get the  
for (cond in conds){   # for each condition
  cond_name<- gsub("\\.mat", "", cond)
  
  # loading the p-value file
  B<- readMat(cond)
  B$pVal.homo<- Filter(Negate(is.null), B$pVal.homo)
  B$pVal.hetero<- Filter(Negate(is.null),B$pVal.hetero)
  
  L<- length(B$pVal.homo)  # No. of animals here
  
  df<-data.frame(Sub=rep(NA,1),
                 Condition=rep(NA,1),
                 Group=rep(NA,1),
                 Pvalue=rep(NA,1),
                 Preference=rep(NA,1))
  
  for (k in 1:L){  # for each animal   
    
    
    ho_p<-as.numeric(unlist(B$pVal.homo[[k]]))
    he_p<-as.numeric(unlist(B$pVal.hetero[[k]]))
    
    homo<- rep(0,length(ho_p))
    hetero<-  rep(0,length(he_p))
    
    Lf<- length(ho_p)  # number of units in each animal
    
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
    
    df<-rbind(df,temp)   # combine homo and hetero
    df<- na.omit(df)
    
  }
  
  row.names(df)<- 1:nrow(df)
  write.csv(df,file=file.path(base_path,paradigm,paste0(cond_name,'_prefer_passive.csv')))
} 

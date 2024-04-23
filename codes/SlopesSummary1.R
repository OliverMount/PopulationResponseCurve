############### Plotting  summary of slopes ############

# Combining all the csv files that contains average slopes across
# significant time points in to a single file. The csv files are made during plotting
# in the pop_decoding.py. Also, the slope summaries are plotted in the python pop_decoding.py file at the end

rm(list = ls())
library(R.matlab)
library(pracma)
library(tidyverse)
library(gridExtra)
library(reshape2)
library(ggpubr) 

save_path='/media/olive/Research/Gitrepo/PopulationResponseCurve/fig_data/' 
setwd('/media/olive/Research/oliver/pop_slopes/task')
flist<-list.files(getwd())

rois<-c('V1_45','V1_90','V1_135','PPC_45','PPC_90','PPC_135') 
percents<-c('0','10','20','40','60','100')

df_task<-data.frame(Condition=rep(NA,1),
               Percent=rep(NA,1),
               Homo=rep(NA,1), 
               Hetero=rep(NA,1))

for (roi in rois){ 
  for (percent in percents){ 
    fname<-paste0(roi,'_',percent,'.csv') 
    temp<-read.csv(fname,header = FALSE) 
    colnames(temp)<-c("Homo","Hetero")
    temp$Condition<- rep(roi,nrow(temp))
    temp$Percent<-rep(percent,nrow(temp))
    df_task<- rbind(df_task,temp)   
  } 
}

df_task<-na.omit(df_task)
df_task$paradigm<- rep("task",nrow(df_task))


### mean_data making for passive 
setwd('/media/olive/Research/oliver/pop_slopes/passive/')
flist<-list.files(getwd())

rois<-c('V1_45','V1_90','V1_135','PPC_45','PPC_90','PPC_135') 
percents<-c('0','10','20','40','60','100')

df_passive<-data.frame(Condition=rep(NA,1),
                       Percent=rep(NA,1),
                       Homo=rep(NA,1), 
                       Hetero=rep(NA,1))

for (roi in rois){ 
  for (percent in percents){ 
    
    fname<-paste0(roi,'_',percent,'.csv') 
    temp<-read.csv(fname,header = FALSE) 
    colnames(temp)<-c("Homo","Hetero")
    temp$Condition<- rep(roi,nrow(temp))
    temp$Percent<-rep(percent,nrow(temp))
    df_passive<- rbind(df_passive,temp)  
    
  } 
}

df_passive<-na.omit(df_passive)
df_passive$paradigm<- rep("passive",nrow(df_passive))


df<- rbind(df_task,df_passive)

df_new <-  df %>% pivot_longer(cols = c("Homo","Hetero"), 
                               values_to = "slopes") 

write.csv(df_new ,file=paste0(save_path,'significant.csv')) 

 
df_new$Condition<- factor(df_new$Condition,levels = c("V1_45","V1_90","V1_135","PPC_45","PPC_90","PPC_135"),
                      labels = c("V1 45","V1 90","V1 135","PPC 45","PPC 90","PPC 135"),
                      ordered = TRUE)
df_new$paradigm<- factor(df_new$paradigm,levels = c("task","passive"),
                     labels =  c("task","passive"),
                     ordered = TRUE)
df_new$name<- factor(df_new$name,levels = c("Homo","Hetero"),
                     labels =  c("Homo","Hetero"),
                     ordered = TRUE)
df_new$Percent<- factor(df_new$Percent,levels = c("0","10","20","40","60","100"))



# all
df_all <- df_new %>%
  group_by(Condition, Percent) %>%
  summarize(mean_value = mean(slopes),
            sd_value = sd(slopes),
            se_value = sd_value / sqrt(n()))
write.csv(df_all ,file=paste0(save_path,'all.csv')) 



#task passive
df_task_passive <- df_new %>%
  group_by(paradigm,Condition, Percent) %>%
  summarize(mean_value = mean(slopes),
            sd_value = sd(slopes),
            se_value = sd_value / sqrt(n()))
write.csv(df_task_passive ,file=paste0(save_path,'task_passive.csv')) 


#homo- hetero
df_ho_het <- df_new %>%
  group_by(name,Condition, Percent) %>%
  summarize(mean_value = mean(slopes),
            sd_value = sd(slopes),
            se_value = sd_value / sqrt(n()))


write.csv(df_ho_het ,file=paste0(save_path,'homo_hetero.csv'))  
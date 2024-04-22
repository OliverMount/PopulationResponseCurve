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



my_colors <- c( "red",  "blue")

# Plot using ggplot
p <- ggplot(mean_data_task, aes(x = Percent, y = mean_value, group=variable, color = variable)) +
  geom_line()+
  geom_point()+
  geom_errorbar(aes(ymin=mean_value-se_value, ymax=mean_value+se_value), width=.2,
                position=position_dodge(0.05))+
  theme_classic() +
  theme(legend.position = "top",
        axis.ticks.length.x = unit(3, 'mm'),
        axis.ticks.length.y = unit(3, 'mm'),
        axis.text = element_text(size = 20),
        axis.title.x = element_text(size = 24),
        axis.title.y = element_text(size = 24),
        plot.title = element_text(size = 24, hjust = 0.5),
        legend.title = element_blank(),
        legend.text = element_text(size=24),
        strip.text.x = element_text(size = 24)) +
  facet_wrap(~Condition, scales = "free_y") +
  labs(x = "Percentage of untuned cells added", y = "Average slope") +
  scale_y_continuous(breaks = seq(0, 0.28, 0.05), limits = c(-0.001, 0.28),
                     expand = c(0, 0))+
  scale_color_manual(values = my_colors)

# Print the plot
print(p)

### mean_data making for passive 
setwd('/media/olive/Research/oliver/pop_slopes/passive/')
flist<-list.files(getwd())

rois<-c('V1_45','V1_90','V1_135','PPC_45','PPC_90','PPC_135') 
percents<-c('0','10','20','40','60','100')

df<-data.frame(Condition=rep(NA,1),Percent=rep(NA,1),Homo=rep(NA,1), Hetero=rep(NA,1))

for (roi in rois){ 
  for (percent in percents){ 
    
    fname<-paste0(roi,'_',percent,'.csv') 
    temp<-read.csv(fname,header = FALSE) 
    colnames(temp)<-c("Homo","Hetero")
    temp$Condition<- rep(roi,nrow(temp))
    temp$Percent<-rep(percent,nrow(temp))
    df<- rbind(df,temp)  
    
  } 
}

df<-na.omit(df)

df<-melt(df)
df$Condition<- factor(df$Condition,levels = c("V1_45","V1_90","V1_135","PPC_45","PPC_90","PPC_135"),
                      labels = c("V1 45","V1 90","V1 135","PPC 45","PPC 90","PPC 135"),
                      ordered = TRUE)
df$Percent<- factor(df$Percent,levels = c("0","10","20","40","60","100"))

#Calculate mean and standard error
mean_data_passive <- df %>%
  group_by(Condition, Percent,variable) %>%
  summarize(mean_value = mean(value),
            sd_value = sd(value),
            se_value = sd_value / sqrt(n()))


mean_data <- rbind(mean_data_task,mean_data_passive)

Nt<- nrow(mean_data_task)
Np<- nrow(mean_data_passive)

temp<-data.frame(Paradigm=rep(c("Task","Passive"),each=c(Nt))) 
final<- cbind(temp,mean_data)
final$Paradigm<- factor(final$Paradigm,levels=c("Task","Passive"),ordered = TRUE)


# may be PPC and V1 has to be plotted separately
# Plot using ggplot
p <- ggplot(final, aes(x = Percent, y = mean_value, 
                       group=interaction(variable, Paradigm), color = variable,
                       linetype=Paradigm)) +
  geom_line()+
  geom_point()+
  geom_errorbar(aes(ymin=mean_value-se_value, ymax=mean_value+se_value), width=.2,
                position=position_dodge(0.05))+
  theme_classic() +
  theme(legend.position = "top",
        axis.ticks.length.x = unit(3, 'mm'),
        axis.ticks.length.y = unit(3, 'mm'),
        axis.text = element_text(size = 20),
        axis.title.x = element_text(size = 24),
        axis.title.y = element_text(size = 24),
        plot.title = element_text(size = 24, hjust = 0.5),
        legend.title = element_blank(),
        legend.text = element_text(size=24),
        strip.text.x = element_text(size = 24)) +
  facet_wrap(~Condition, scales = "free_y") +
  labs(x = "Percentage of untuned cells added", y = "Average slope") +
  scale_y_continuous(breaks = seq(0, 0.28, 0.05), limits = c(-0.001, 0.28),
                     expand = c(0, 0))+
  scale_color_manual(values = my_colors)+
  scale_linetype_manual(values=c("solid","dashed"))

# Print the plot
print(p) 
write.csv(final,file=paste0(save_path,'final.csv')) 

### IEM decoding

setwd('~/Desktop/res/task/')
flist<-list.files(getwd())

rois<-c('V1','PPC')
conds<-c('45','90','135')
percents<-c('10','20','40','60','100')


df<-data.frame(Condition=rep(NA,1),Percent=rep(NA,1),Homo=rep(NA,1), Hetero=rep(NA,1))

for (roi in rois){
  for(cond in conds){
    for (percent in percents){
      
      cond_name=paste0(roi,'_',cond)
      fname<-paste0(roi,'_',cond,'_',percent,'.csv') 
      temp<-read.csv(fname) 
      colnames(temp)<-c("Homo","Hetero")
      temp$Condition<- rep(cond_name,nrow(temp))
      temp$Percent<-rep(percent,nrow(temp))
      df<- rbind(df,temp)  
      
    }  
  } 
}

df<-na.omit(df)
write.csv(df,file=paste0(save_path,'slopes.csv')) 


df<-melt(df)
df$Condition<- factor(df$Condition,levels = c("V1_45","V1_90","V1_135","PPC_45","PPC_90","PPC_135"),
                      labels = c("V1 45","V1 90","V1 135","PPC 45","PPC 90","PPC 135"),
                      ordered = TRUE)
df$Percent<- factor(df$Percent,levels = c("10","20","40","60","100"))


#Calculate mean and standard error
mean_data <- df %>%
  group_by(Condition, Percent,variable) %>%
  summarize(mean_value = mean(value),
            sd_value = sd(value),
            se_value = sd_value / sqrt(n()))


p <- ggplot(mean_data, aes(x = Percent, y = mean_value, group=variable,color= variable)) +
  geom_line()+
  geom_point()+
  #geom_boxplot(position = "dodge",outlier.color = "red",
  #             outlier.shape = NA) + 
  #geom_jitter(width=0.1,alpha=1,size=1,color="black",shape=21,fill="grey")+
  #stat_summary(fun=mean, geom='point', shape=23, size=3,
  #             color="black",fill="magenta",alpha=0.7)+
  geom_errorbar(aes(ymin=mean_value-se_value,ymax=mean_value+se_value),width=.2,
                position=position_dodge(0.05))+
  theme_classic()+
  theme(legend.position = "top",
        axis.ticks.length.x = unit(3, 'mm'),
        axis.ticks.length.y = unit(3, 'mm'),
        axis.text = element_text(size = 20),
        axis.title.x = element_text(size = 24),
        axis.title.y = element_text(size = 24),
        plot.title = element_text(size = 24, hjust = 0.5),
        legend.title = element_blank(),
        legend.text = element_text(size=24),
        strip.text.x = element_text(size = 24)) +
  facet_wrap(~Condition, scales = "free_y") +
  labs(x = "Percentage of cells for decoding after sorting (tuned-->untuned)", y = "Average slope") +
  scale_y_continuous(breaks = seq(0, 0.28, 0.05), limits = c(-0.001, 0.28),
                     expand = c(0, 0))+
  scale_color_manual(values = my_colors)

# Print the plot
print(p) 


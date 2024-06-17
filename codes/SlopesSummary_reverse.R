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

save_path='/media/olive/Research/oliver/IEMdecodingForCalciumData/neuron_counts/' 
setwd('/media/olive/Research/oliver/pop_slopes/task')
flist<-list.files(getwd())

rois<-c('V1_45','V1_90','V1_135','PPC_45','PPC_90','PPC_135') 
#percents<-c('0','10','20','40','60','100')
#percents<-c('0','10','25','50','100')
percents<-c('0','0.05','0.5','0.9')

df<-data.frame(Condition=rep(NA,1),
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
    df<- rbind(df,temp)   
  } 
}

df<-na.omit(df)
write.csv(df,file=paste0(save_path,'slopes_popdecoding.csv')) 


df<-melt(df)
df$Condition<- factor(df$Condition,levels = c("V1_45","V1_90","V1_135","PPC_45","PPC_90","PPC_135"),
                      labels = c("V1 45","V1 90","V1 135","PPC 45","PPC 90","PPC 135"),
                      ordered = TRUE)
df$Percent<- factor(df$Percent,levels =percents)

#Calculate mean and standard error
mean_data_task <- df %>%
  group_by(Condition, Percent,variable) %>%
  summarize(mean_value = mean(value),
            sd_value = sd(value),
            se_value = sd_value / sqrt(n()))
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
df$Percent<- factor(df$Percent,levels = percents)

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
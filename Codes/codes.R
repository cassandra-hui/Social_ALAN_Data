#Analysis

#Set working Directory
setwd("~/UNR/Social/ALAN/Data/R")


#Install Libraries 
library(tidyverse)
library(dplyr)
library(ggplot2)
library(lme4)
library(plotrix)
library(cosinor)
library(circacompare)
library(readr)
library(lmerTest)




##load data
social_df <- read.csv('Social_ALAN.csv')


#Normalize Data
df <- social_df %>%  filter(!is.na(time)) #Remove dead birds


################# Normalized data
################
sp <- subset(df, Date_In != "9/1/2023") #spring data
fa <- subset(df, Date_In == "9/1/2023") #fal data
sp <- subset(sp, Date_In != "10/13/2023") #remove mel only birds

sp$norm_CRY <- (sp$CRY - min(sp$CRY, na.rm = TRUE)) / (max(sp$CRY, na.rm = TRUE) - min(sp$CRY, na.rm = TRUE))
fa$norm_CRY <- (fa$CRY - min(fa$CRY, na.rm = TRUE)) / (max(fa$CRY, na.rm = TRUE) - min(fa$CRY, na.rm = TRUE))

# For PER3
sp$norm_PER3 <- (sp$PER3 - min(sp$PER3, na.rm = TRUE)) / (max(sp$PER3, na.rm = TRUE) - min(sp$PER3, na.rm = TRUE))
fa$norm_PER3 <- (fa$PER3 - min(fa$PER3, na.rm = TRUE)) / (max(fa$PER3, na.rm = TRUE) - min(fa$PER3, na.rm = TRUE))

# For BMAL
sp$norm_BMAL <- (sp$BMAL - min(sp$BMAL, na.rm = TRUE)) / (max(sp$BMAL, na.rm = TRUE) - min(sp$BMAL, na.rm = TRUE))
fa$norm_BMAL <- (fa$BMAL - min(fa$BMAL, na.rm = TRUE)) / (max(fa$BMAL, na.rm = TRUE) - min(fa$BMAL, na.rm = TRUE))

# For L_CRY
sp$norm_L_CRY <- (sp$L_CRY - min(sp$L_CRY, na.rm = TRUE)) / (max(sp$L_CRY, na.rm = TRUE) - min(sp$L_CRY, na.rm = TRUE))
fa$norm_L_CRY <- (fa$L_CRY - min(fa$L_CRY, na.rm = TRUE)) / (max(fa$L_CRY, na.rm = TRUE) - min(fa$L_CRY, na.rm = TRUE))


# For L_BMAL
sp$norm_L_BMAL <- (sp$L_BMAL - min(sp$L_BMAL, na.rm = TRUE)) / (max(sp$L_BMAL, na.rm = TRUE) - min(sp$L_BMAL, na.rm = TRUE))
fa$norm_L_BMAL <- (fa$L_BMAL - min(fa$L_BMAL, na.rm = TRUE)) / (max(fa$L_BMAL, na.rm = TRUE) - min(fa$L_BMAL, na.rm = TRUE))

#For melatonin
sp$norm_mel <- (sp$melatonin - min(sp$melatonin, na.rm = TRUE)) / (max(sp$melatonin, na.rm = TRUE) - min(sp$melatonin, na.rm = TRUE))
fa$norm_mel <- (fa$melatonin - min(fa$melatonin, na.rm = TRUE)) / (max(fa$melatonin, na.rm = TRUE) - min(fa$melatonin, na.rm = TRUE))


df <- rbind(sp, fa)

df$norm_L_PER3 <- (df$L_PER3 - min(df$L_PER3, na.rm = TRUE)) / (max(df$L_PER3, na.rm = TRUE) - min(df$L_PER3, na.rm = TRUE))

###
#Fall PER2 is so different I do not want to use it
###


df$norm_PER2 <- NULL

df$s_PER2[df$Date_In != "9/1/2023"] <- df$PER2[df$Date_In != "9/1/2023"]
df$norm_PER2 <- (df$s_PER2 - min(df$s_PER2, na.rm = TRUE)) / (max(df$s_PER2, na.rm = TRUE) - min(df$s_PER2, na.rm = TRUE))
# 

df$norm_L_PER3 <- (df$L_PER3 - min(df$L_PER3, na.rm = TRUE)) / (max(df$L_PER3, na.rm = TRUE) - min(df$L_PER3, na.rm = TRUE))

df$norm_L_PER2 <- NULL

df$s_L_PER2[df$Date_In != "9/1/2023"] <- df$L_PER2[df$Date_In != "9/1/2023"]
df$norm_L_PER2 <- (df$s_L_PER2 - min(df$s_L_PER2, na.rm = TRUE)) / (max(df$s_L_PER2, na.rm = TRUE) - min(df$s_L_PER2, na.rm = TRUE))
# 

#remove excess columns
# Remove columns with specified names
df <- df %>%
  select(-BMAL, -CRY, -PER2, -PER3, -L_BMAL, -L_CRY, -L_PER2, -L_PER3, -s_PER2, -s_L_PER2, melatonin)

###################

isolated <- subset(df, treatment == "LD" | treatment == "A")

social <- subset(df, treatment == "SLD" | treatment == "SA")

alan <- subset(df, treatment == "SA" | treatment == "A")

dark <- subset(df, treatment == "SLD" | treatment == "LD")


#######################
### Gene Expression ###
#######################

iso_gene <- isolated %>% 
  group_by(time, treatment) %>% 
  summarise(cry_mean=mean(norm_CRY, na.rm=TRUE), cry_se=std.error(norm_CRY, na.rm=TRUE),  
            arnt_mean=mean(norm_BMAL, na.rm=TRUE), arnt_se=std.error(norm_BMAL, na.rm=TRUE), 
            per2_mean=mean(norm_PER2, na.rm=TRUE), per2_se=std.error(norm_PER2, na.rm=TRUE),
            per3_mean=mean(norm_PER3, na.rm=TRUE), per3_se=std.error(norm_PER3, na.rm=TRUE), n = n())


soc_gene <- social %>% 
  group_by(time, treatment) %>% 
  summarise(cry_mean=mean(norm_CRY, na.rm=TRUE), cry_se=std.error(norm_CRY, na.rm=TRUE),  
            arnt_mean=mean(norm_BMAL, na.rm=TRUE), arnt_se=std.error(norm_BMAL, na.rm=TRUE), 
            per2_mean=mean(norm_PER2, na.rm=TRUE), per2_se=std.error(norm_PER2, na.rm=TRUE), 
            per3_mean=mean(norm_PER3, na.rm=TRUE), per3_se=std.error(norm_PER3, na.rm=TRUE), n = n())


dk_gene <- dark %>% 
  group_by(time, treatment) %>% 
  summarise(cry_mean=mean(norm_CRY, na.rm=TRUE), cry_se=std.error(norm_CRY, na.rm=TRUE),  
            arnt_mean=mean(norm_BMAL, na.rm=TRUE), arnt_se=std.error(norm_BMAL, na.rm=TRUE), 
            per2_mean=mean(norm_PER2, na.rm=TRUE), per2_se=std.error(norm_PER2, na.rm=TRUE),
            per3_mean=mean(norm_PER3, na.rm=TRUE), per3_se=std.error(norm_PER3, na.rm=TRUE), n = n())



al_gene <- alan %>% 
  group_by(time, treatment) %>% 
  summarise(cry_mean=mean(norm_CRY, na.rm=TRUE), cry_se=std.error(norm_CRY, na.rm=TRUE),  
            arnt_mean=mean(norm_BMAL, na.rm=TRUE), arnt_se=std.error(norm_BMAL, na.rm=TRUE), 
            per2_mean=mean(norm_PER2, na.rm=TRUE), per2_se=std.error(norm_PER2, na.rm=TRUE),
            per3_mean=mean(norm_PER3, na.rm=TRUE), per3_se=std.error(norm_PER3, na.rm=TRUE), n = n())


#######
#Circacompare
#######


#BMAL
#Iso
res <- circacompare(x=isolated, col_time="time", col_group="treatment", col_outcome="norm_BMAL", alpha_threshold = "0.1")
res

#Soc
res <- circacompare(x=social, col_time="time", col_group="treatment", col_outcome="norm_BMAL")
res


#CRY
#Iso
res <- circacompare(x=isolated, col_time="time", col_group="treatment", col_outcome="norm_CRY")
res

#Soc
res <- circacompare(x=social, col_time="time", col_group="treatment", col_outcome="norm_CRY")
res


#PER2
#Iso
res <- circacompare(x=isolated, col_time="time", col_group="treatment", col_outcome="norm_PER2")
res

#Soc
res <- circacompare(x=social, col_time="time", col_group="treatment", col_outcome="norm_PER2")
res


#PER3
#Iso
res <- circacompare(x=isolated, col_time="time", col_group="treatment", col_outcome="norm_PER3")
res

#Soc
res <- circacompare(x=social, col_time="time", col_group="treatment", col_outcome="norm_PER3")
res

########################
### Individual Times ###
########################

#############
#Isolated
#############
######

#time 1 
time_point <- 1
# Subset the data for the specific time point
subset_data <- isolated %>% filter(time == time_point)
# Perform t-test for the "A" and "LD" treatments on the 'cry_mean' variable
t.test(norm_BMAL ~ treatment, data = subset_data)

#time 7 
time_point <- 7
# Subset the data for the specific time point
subset_data <- isolated %>% filter(time == time_point)
# Perform t-test for the "A" and "LD" treatments on the 'cry_mean' variable
t.test(norm_BMAL ~ treatment, data = subset_data)

#time 13
time_point <- 13
# Subset the data for the specific time point
subset_data <- isolated %>% filter(time == time_point)
# Perform t-test for the "A" and "LD" treatments on the 'cry_mean' variable
t.test(norm_BMAL ~ treatment, data = subset_data)

#time 19
time_point <- 19
# Subset the data for the specific time point
subset_data <- isolated %>% filter(time == time_point)
# Perform t-test for the "A" and "LD" treatments on the 'cry_mean' variable
t.test(norm_BMAL ~ treatment, data = subset_data)

#######
#time 1 
time_point <- 1
# Subset the data for the specific time point
subset_data <- isolated %>% filter(time == time_point)
# Perform t-test for the "A" and "LD" treatments on the 'cry_mean' variable
t.test(norm_CRY ~ treatment, data = subset_data)

#time 7 
time_point <- 7
# Subset the data for the specific time point
subset_data <- isolated %>% filter(time == time_point)
# Perform t-test for the "A" and "LD" treatments on the 'cry_mean' variable
t.test(norm_CRY ~ treatment, data = subset_data)

#time 13
time_point <- 13
# Subset the data for the specific time point
subset_data <- isolated %>% filter(time == time_point)
# Perform t-test for the "A" and "LD" treatments on the 'cry_mean' variable
t.test(norm_CRY ~ treatment, data = subset_data)

#time 19
time_point <- 19
# Subset the data for the specific time point
subset_data <- isolated %>% filter(time == time_point)
# Perform t-test for the "A" and "LD" treatments on the 'cry_mean' variable
t.test(norm_CRY ~ treatment, data = subset_data)

#######
#time 1 
time_point <- 1
# Subset the data for the specific time point
subset_data <- isolated %>% filter(time == time_point)
# Perform t-test for the "A" and "LD" treatments on the 'cry_mean' variable
t.test(norm_PER2 ~ treatment, data = subset_data)

#time 7 
time_point <- 7
# Subset the data for the specific time point
subset_data <- isolated %>% filter(time == time_point)
# Perform t-test for the "A" and "LD" treatments on the 'cry_mean' variable
t.test(norm_PER2 ~ treatment, data = subset_data)

#time 13
time_point <- 13
# Subset the data for the specific time point
subset_data <- isolated %>% filter(time == time_point)
# Perform t-test for the "A" and "LD" treatments on the 'cry_mean' variable
t.test(norm_PER2 ~ treatment, data = subset_data)

#time 19
time_point <- 19
# Subset the data for the specific time point
subset_data <- isolated %>% filter(time == time_point)
# Perform t-test for the "A" and "LD" treatments on the 'cry_mean' variable
t.test(norm_PER2 ~ treatment, data = subset_data)

#######
#time 1 
time_point <- 1
# Subset the data for the specific time point
subset_data <- isolated %>% filter(time == time_point)
# Perform t-test for the "A" and "LD" treatments on the 'cry_mean' variable
t.test(norm_PER3 ~ treatment, data = subset_data)

#time 7 
time_point <- 7
# Subset the data for the specific time point
subset_data <- isolated %>% filter(time == time_point)
# Perform t-test for the "A" and "LD" treatments on the 'cry_mean' variable
t.test(norm_PER3 ~ treatment, data = subset_data)

#time 13
time_point <- 13
# Subset the data for the specific time point
subset_data <- isolated %>% filter(time == time_point)
# Perform t-test for the "A" and "LD" treatments on the 'cry_mean' variable
t.test(norm_PER3 ~ treatment, data = subset_data)

#time 19
time_point <- 19
# Subset the data for the specific time point
subset_data <- isolated %>% filter(time == time_point)
# Perform t-test for the "A" and "LD" treatments on the 'cry_mean' variable
t.test(norm_PER3 ~ treatment, data = subset_data)

#####################
#Social
#################


#time 1 
time_point <- 1
# Subset the data for the specific time point
subset_data <- social %>% filter(time == time_point)
# Perform t-test for the "A" and "LD" treatments on the 'cry_mean' variable
t.test(norm_BMAL ~ treatment, data = subset_data)

#time 7 
time_point <- 7
# Subset the data for the specific time point
subset_data <- social %>% filter(time == time_point)
# Perform t-test for the "A" and "LD" treatments on the 'cry_mean' variable
t.test(norm_BMAL ~ treatment, data = subset_data)

#time 13
time_point <- 13
# Subset the data for the specific time point
subset_data <- social %>% filter(time == time_point)
# Perform t-test for the "A" and "LD" treatments on the 'cry_mean' variable
t.test(norm_BMAL ~ treatment, data = subset_data)

#time 19
time_point <- 19
# Subset the data for the specific time point
subset_data <- social %>% filter(time == time_point)
# Perform t-test for the "A" and "LD" treatments on the 'cry_mean' variable
t.test(norm_BMAL ~ treatment, data = subset_data)

#######
#time 1 
time_point <- 1
# Subset the data for the specific time point
subset_data <- social %>% filter(time == time_point)
# Perform t-test for the "A" and "LD" treatments on the 'cry_mean' variable
t.test(norm_CRY ~ treatment, data = subset_data)

#time 7 
time_point <- 7
# Subset the data for the specific time point
subset_data <- social %>% filter(time == time_point)
# Perform t-test for the "A" and "LD" treatments on the 'cry_mean' variable
t.test(norm_CRY ~ treatment, data = subset_data)

#time 13
time_point <- 13
# Subset the data for the specific time point
subset_data <- social %>% filter(time == time_point)
# Perform t-test for the "A" and "LD" treatments on the 'cry_mean' variable
t.test(norm_CRY ~ treatment, data = subset_data)

#time 19
time_point <- 19
# Subset the data for the specific time point
subset_data <- social %>% filter(time == time_point)
# Perform t-test for the "A" and "LD" treatments on the 'cry_mean' variable
t.test(norm_CRY ~ treatment, data = subset_data)

#######
#time 1 
time_point <- 1
# Subset the data for the specific time point
subset_data <- social %>% filter(time == time_point)
# Perform t-test for the "A" and "LD" treatments on the 'cry_mean' variable
t.test(norm_PER2 ~ treatment, data = subset_data)

#time 7 
time_point <- 7
# Subset the data for the specific time point
subset_data <- social %>% filter(time == time_point)
# Perform t-test for the "A" and "LD" treatments on the 'cry_mean' variable
t.test(norm_PER2 ~ treatment, data = subset_data)

#time 13
time_point <- 13
# Subset the data for the specific time point
subset_data <- social %>% filter(time == time_point)
# Perform t-test for the "A" and "LD" treatments on the 'cry_mean' variable
t.test(norm_PER2 ~ treatment, data = subset_data)

#time 19
time_point <- 19
# Subset the data for the specific time point
subset_data <- social %>% filter(time == time_point)
# Perform t-test for the "A" and "LD" treatments on the 'cry_mean' variable
t.test(norm_PER2 ~ treatment, data = subset_data)

#######
#time 1 
time_point <- 1
# Subset the data for the specific time point
subset_data <- social %>% filter(time == time_point)
# Perform t-test for the "A" and "LD" treatments on the 'cry_mean' variable
t.test(norm_PER3 ~ treatment, data = subset_data)

#time 7 
time_point <- 7
# Subset the data for the specific time point
subset_data <- social %>% filter(time == time_point)
# Perform t-test for the "A" and "LD" treatments on the 'cry_mean' variable
t.test(norm_PER3 ~ treatment, data = subset_data)

#time 13
time_point <- 13
# Subset the data for the specific time point
subset_data <- social %>% filter(time == time_point)
# Perform t-test for the "A" and "LD" treatments on the 'cry_mean' variable
t.test(norm_PER3 ~ treatment, data = subset_data)

#time 19
time_point <- 19
# Subset the data for the specific time point
subset_data <- social %>% filter(time == time_point)
# Perform t-test for the "A" and "LD" treatments on the 'cry_mean' variable
t.test(norm_PER3 ~ treatment, data = subset_data)


##############
#Dark
#################


#time 1 
time_point <- 1
# Subset the data for the specific time point
subset_data <- dark %>% filter(time == time_point)
# Perform t-test for the "A" and "LD" treatments on the 'cry_mean' variable
t.test(norm_BMAL ~ treatment, data = subset_data)

#time 7 
time_point <- 7
# Subset the data for the specific time point
subset_data <- dark %>% filter(time == time_point)
# Perform t-test for the "A" and "LD" treatments on the 'cry_mean' variable
t.test(norm_BMAL ~ treatment, data = subset_data)

#time 13
time_point <- 13
# Subset the data for the specific time point
subset_data <- dark %>% filter(time == time_point)
# Perform t-test for the "A" and "LD" treatments on the 'cry_mean' variable
t.test(norm_BMAL ~ treatment, data = subset_data)

#time 19
time_point <- 19
# Subset the data for the specific time point
subset_data <- dark %>% filter(time == time_point)
# Perform t-test for the "A" and "LD" treatments on the 'cry_mean' variable
t.test(norm_BMAL ~ treatment, data = subset_data)

#######
#time 1 
time_point <- 1
# Subset the data for the specific time point
subset_data <- dark %>% filter(time == time_point)
# Perform t-test for the "A" and "LD" treatments on the 'cry_mean' variable
t.test(norm_CRY ~ treatment, data = subset_data)

#time 7 
time_point <- 7
# Subset the data for the specific time point
subset_data <- dark %>% filter(time == time_point)
# Perform t-test for the "A" and "LD" treatments on the 'cry_mean' variable
t.test(norm_CRY ~ treatment, data = subset_data)

#time 13
time_point <- 13
# Subset the data for the specific time point
subset_data <- dark %>% filter(time == time_point)
# Perform t-test for the "A" and "LD" treatments on the 'cry_mean' variable
t.test(norm_CRY ~ treatment, data = subset_data)

#time 19
time_point <- 19
# Subset the data for the specific time point
subset_data <- dark %>% filter(time == time_point)
# Perform t-test for the "A" and "LD" treatments on the 'cry_mean' variable
t.test(norm_CRY ~ treatment, data = subset_data)

#######
#time 1 
time_point <- 1
# Subset the data for the specific time point
subset_data <- dark %>% filter(time == time_point)
# Perform t-test for the "A" and "LD" treatments on the 'cry_mean' variable
t.test(norm_PER2 ~ treatment, data = subset_data)

#time 7 
time_point <- 7
# Subset the data for the specific time point
subset_data <- dark %>% filter(time == time_point)
# Perform t-test for the "A" and "LD" treatments on the 'cry_mean' variable
t.test(norm_PER2 ~ treatment, data = subset_data)

#time 13
time_point <- 13
# Subset the data for the specific time point
subset_data <- dark %>% filter(time == time_point)
# Perform t-test for the "A" and "LD" treatments on the 'cry_mean' variable
t.test(norm_PER2 ~ treatment, data = subset_data)

#time 19
time_point <- 19
# Subset the data for the specific time point
subset_data <- dark %>% filter(time == time_point)
# Perform t-test for the "A" and "LD" treatments on the 'cry_mean' variable
t.test(norm_PER2 ~ treatment, data = subset_data)

#######
#time 1 
time_point <- 1
# Subset the data for the specific time point
subset_data <- dark %>% filter(time == time_point)
# Perform t-test for the "A" and "LD" treatments on the 'cry_mean' variable
t.test(norm_PER3 ~ treatment, data = subset_data)

#time 7 
time_point <- 7
# Subset the data for the specific time point
subset_data <- dark %>% filter(time == time_point)
# Perform t-test for the "A" and "LD" treatments on the 'cry_mean' variable
t.test(norm_PER3 ~ treatment, data = subset_data)

#time 13
time_point <- 13
# Subset the data for the specific time point
subset_data <- dark %>% filter(time == time_point)
# Perform t-test for the "A" and "LD" treatments on the 'cry_mean' variable
t.test(norm_PER3 ~ treatment, data = subset_data)

#time 19
time_point <- 19
# Subset the data for the specific time point
subset_data <- dark %>% filter(time == time_point)
# Perform t-test for the "A" and "LD" treatments on the 'cry_mean' variable
t.test(norm_PER3 ~ treatment, data = subset_data)

################
#ALAN
################


#time 1 
time_point <- 1
# Subset the data for the specific time point
subset_data <- alan %>% filter(time == time_point)
# Perform t-test for the "A" and "LD" treatments on the 'cry_mean' variable
t.test(norm_L_BMAL ~ treatment, data = subset_data)

#time 7 
time_point <- 7
# Subset the data for the specific time point
subset_data <- alan %>% filter(time == time_point)
# Perform t-test for the "A" and "LD" treatments on the 'cry_mean' variable
t.test(norm_L_BMAL ~ treatment, data = subset_data)

#time 13
time_point <- 13
# Subset the data for the specific time point
subset_data <- alan %>% filter(time == time_point)
# Perform t-test for the "A" and "LD" treatments on the 'cry_mean' variable
t.test(norm_L_BMAL ~ treatment, data = subset_data)

#time 19
time_point <- 19
# Subset the data for the specific time point
subset_data <- alan %>% filter(time == time_point)
# Perform t-test for the "A" and "LD" treatments on the 'cry_mean' variable
t.test(norm_L_BMAL ~ treatment, data = subset_data)

#######
#time 1 
time_point <- 1
# Subset the data for the specific time point
subset_data <- alan %>% filter(time == time_point)
# Perform t-test for the "A" and "LD" treatments on the 'cry_mean' variable
t.test(norm_CRY ~ treatment, data = subset_data)

#time 7 
time_point <- 7
# Subset the data for the specific time point
subset_data <- alan %>% filter(time == time_point)
# Perform t-test for the "A" and "LD" treatments on the 'cry_mean' variable
t.test(norm_CRY ~ treatment, data = subset_data)

#time 13
time_point <- 13
# Subset the data for the specific time point
subset_data <- alan %>% filter(time == time_point)
# Perform t-test for the "A" and "LD" treatments on the 'cry_mean' variable
t.test(norm_CRY ~ treatment, data = subset_data)

#time 19
time_point <- 19
# Subset the data for the specific time point
subset_data <- alan %>% filter(time == time_point)
# Perform t-test for the "A" and "LD" treatments on the 'cry_mean' variable
t.test(norm_CRY ~ treatment, data = subset_data)

#######
#time 1 
time_point <- 1
# Subset the data for the specific time point
subset_data <- alan %>% filter(time == time_point)
# Perform t-test for the "A" and "LD" treatments on the 'cry_mean' variable
t.test(norm_PER2 ~ treatment, data = subset_data)

#time 7 
time_point <- 7
# Subset the data for the specific time point
subset_data <- alan %>% filter(time == time_point)
# Perform t-test for the "A" and "LD" treatments on the 'cry_mean' variable
t.test(norm_PER2 ~ treatment, data = subset_data)

#time 13
time_point <- 13
# Subset the data for the specific time point
subset_data <- alan %>% filter(time == time_point)
# Perform t-test for the "A" and "LD" treatments on the 'cry_mean' variable
t.test(norm_PER2 ~ treatment, data = subset_data)

#time 19
time_point <- 19
# Subset the data for the specific time point
subset_data <- alan %>% filter(time == time_point)
# Perform t-test for the "A" and "LD" treatments on the 'cry_mean' variable
t.test(norm_PER2 ~ treatment, data = subset_data)

#######
#time 1 
time_point <- 1
# Subset the data for the specific time point
subset_data <- alan %>% filter(time == time_point)
# Perform t-test for the "A" and "LD" treatments on the 'cry_mean' variable
t.test(norm_L_PER3 ~ treatment, data = subset_data)

#time 7 
time_point <- 7
# Subset the data for the specific time point
subset_data <- alan %>% filter(time == time_point)
# Perform t-test for the "A" and "LD" treatments on the 'cry_mean' variable
t.test(norm_L_PER3 ~ treatment, data = subset_data)

#time 13
time_point <- 13
# Subset the data for the specific time point
subset_data <- alan %>% filter(time == time_point)
# Perform t-test for the "A" and "LD" treatments on the 'cry_mean' variable
t.test(norm_L_PER3 ~ treatment, data = subset_data)

#time 19
time_point <- 19
# Subset the data for the specific time point
subset_data <- alan %>% filter(time == time_point)
# Perform t-test for the "A" and "LD" treatments on the 'cry_mean' variable
t.test(norm_L_PER3 ~ treatment, data = subset_data)


##########################


####################
## Gene Expression Figures
###################

#bmal1
iso_bmal <- ggplot() +
  geom_rect(aes(xmin = 12, xmax = 24, ymin = -Inf, ymax = Inf),fill = "grey1", alpha = .2)+
  geom_line(data=iso_gene, size=1.5, 
            aes(x=time, y=arnt_mean, group=treatment, color=treatment, linetype=treatment)) +
  geom_errorbar(data=iso_gene, 
                aes(x=time, ymin = arnt_mean-arnt_se, ymax = arnt_mean+arnt_se, color=treatment), 
                position = position_dodge(width=0.8), size = 1.5, width = 0)+
  geom_point(data=iso_gene, aes(x=time, y=arnt_mean, fill=treatment, colour=treatment), 
             position = position_dodge(width=0.8), size = 5)+
  #geom_line(data=iso_plotdata, aes(x=time, y=fitY, color=factor(fitX), linetype = factor(fitX)))+
  xlab("Zeitgeber time")+
  ylab(expression(paste(italic("bmal1"), " (Normalized)"))) + # Italicize "cry1"
  #scale_y_continuous(breaks=c(0,75,150), limits=c(0,150))+
  scale_color_manual(values=c("#FFD700", "black"))+
  scale_x_continuous(breaks=c(1,7,13,19), limits=c(0,24))+
  scale_y_continuous(breaks=c(0.2, 0.4, 0.6, 0.8), limits=c(0,0.8))+
  scale_linetype_manual(values=c("solid", "dashed"))+
  theme_classic()+
  theme(
    axis.title.x = element_text( size=24, family= "Times", vjust=-0.2),
    axis.text.x  = element_text( size=16, family= "Times"),
    axis.title.y = element_text( size=24, family= "Times", vjust=2), #vjust is the distance from the text
    axis.text.y  = element_text( size=16, family= "Times"),
    legend.title = element_blank(),
    legend.text = element_text(""),
    legend.position = "none" 
  )


iso_bmal
# Save the legend as a separate image
ggsave("Gene Figures/Iso bmal.png", iso_bmal, width = 4.8, height = 3.8, dpi = 300)


#Cry
iso_cry <- ggplot() +
  geom_rect(aes(xmin = 12, xmax = 24, ymin = -Inf, ymax = Inf),fill = "grey1", alpha = .2)+
  geom_line(data=iso_gene, size=1.5, aes(x=time, y=cry_mean, group=treatment, color=treatment, linetype=treatment)) +
  geom_errorbar(data=iso_gene, 
                aes(x=time, ymin = cry_mean-cry_se, ymax = cry_mean+cry_se, color=treatment), 
                position = position_dodge(width=0.8), size = 1.5, width = 0)+
  geom_point(data=iso_gene, aes(x=time, y=cry_mean, fill=treatment, colour=treatment), 
             position = position_dodge(width=0.8), size = 5)+
  #geom_line(data=iso_plotdata, aes(x=time, y=fitY, color=factor(fitX), linetype = factor(fitX)))+
  xlab("Zeitgeber time")+
  ylab(expression(paste(italic("cry1"), " (Normalized)"))) + # Italicize "cry1"
  #scale_y_continuous(breaks=c(0,75,150), limits=c(0,150))+
  scale_color_manual(values=c("#FFD700", "black"))+
  scale_x_continuous(breaks=c(1,7,13,19), limits=c(0,24))+
  scale_y_continuous(breaks=c(0.2, 0.4, 0.6, 0.8), limits=c(0,0.8))+
  scale_linetype_manual(values=c("solid", "dashed"))+
  theme_classic()+
  theme(
    axis.title.x = element_text( size=24, family= "Times", vjust=-0.2),
    axis.text.x  = element_text( size=16, family= "Times"),
    axis.title.y = element_text( size=24, family= "Times", vjust=2), #vjust is the distance from the text
    axis.text.y  = element_text( size=16, family= "Times"),
    legend.title = element_blank(),
    legend.text = element_text(""),
    legend.position = "none" 
  )

iso_cry
# Save the legend as a separate image
ggsave("Gene Figures/Iso cry.png", iso_cry, width = 4.8, height = 3.8, dpi = 300)



#per2
iso_per2 <- ggplot() +
  geom_rect(aes(xmin = 12, xmax = 24, ymin = -Inf, ymax = Inf),fill = "grey1", alpha = .2)+
  geom_line(data=iso_gene, size=1.5, 
            aes(x=time, y=per2_mean, group=treatment, color=treatment, linetype=treatment)) +
  geom_errorbar(data=iso_gene, 
                aes(x=time, ymin = per2_mean-per2_se, ymax = per2_mean+per2_se, color=treatment), 
                position = position_dodge(width=0.8), size = 1.5, width = 0)+
  geom_point(data=iso_gene, aes(x=time, y=per2_mean, fill=treatment, colour=treatment), 
             position = position_dodge(width=0.8), size = 5)+
  #geom_line(data=iso_plotdata, aes(x=time, y=fitY, color=factor(fitX), linetype = factor(fitX)))+
  xlab("Zeitgeber time")+
  ylab(expression(paste(italic("per2"), " (Normalized)"))) + # Italicize "cry1"
  #scale_y_continuous(breaks=c(0,75,150), limits=c(0,150))+
  scale_color_manual(values=c("#FFD700", "black"))+
  scale_x_continuous(breaks=c(1,7,13,19), limits=c(0,24))+
  scale_y_continuous(breaks=c(0.2, 0.4, 0.6, 0.8), limits=c(0,0.8))+
  scale_linetype_manual(values=c("solid", "dashed"))+
  theme_classic()+
  theme(
    axis.title.x = element_text( size=24, family= "Times", vjust=-0.2),
    axis.text.x  = element_text( size=16, family= "Times"),
    axis.title.y = element_text( size=24, family= "Times", vjust=2), #vjust is the distance from the text
    axis.text.y  = element_text( size=16, family= "Times"),
    legend.title = element_blank(),
    legend.text = element_text(""),
    legend.position = "none" 
  )

iso_per2
# Save the legend as a separate image
ggsave("Gene Figures/Iso per2.png", iso_per2, width = 4.8, height = 3.8, dpi = 300)

#per3
iso_per3 <- ggplot() +
  geom_rect(aes(xmin = 12, xmax = 24, ymin = -Inf, ymax = Inf),fill = "grey1", alpha = .2)+
  geom_line(data=iso_gene, size=1.5, 
            aes(x=time, y=per3_mean, group=treatment, color=treatment, linetype=treatment)) +
  geom_errorbar(data=iso_gene, 
                aes(x=time, ymin = per3_mean-per3_se, ymax = per3_mean+per3_se, color=treatment), 
                position = position_dodge(width=0.8), size = 1.5, width = 0)+
  geom_point(data=iso_gene, aes(x=time, y=per3_mean, fill=treatment, colour=treatment), 
             position = position_dodge(width=0.8), size = 5)+
  #geom_line(data=iso_plotdata, aes(x=time, y=fitY, color=factor(fitX), linetype = factor(fitX)))+
  xlab("Zeitgeber time")+
  ylab(expression(paste(italic("per3"), " (Normalized)"))) + # Italicize "cry1"
  #scale_y_continuous(breaks=c(0,75,150), limits=c(0,150))+
  scale_color_manual(values=c("#FFD700", "black"))+
  scale_x_continuous(breaks=c(1,7,13,19), limits=c(0,24))+
  scale_y_continuous(breaks=c(0.2, 0.4, 0.6, 0.8), limits=c(0,0.8))+
  scale_linetype_manual(values=c("solid", "dashed"))+
  theme_classic()+
  theme(
    axis.title.x = element_text( size=24, family= "Times", vjust=-0.2),
    axis.text.x  = element_text( size=16, family= "Times"),
    axis.title.y = element_text( size=24, family= "Times", vjust=2), #vjust is the distance from the text
    axis.text.y  = element_text( size=16, family= "Times"),
    legend.title = element_blank(),
    legend.text = element_text(""),
    legend.position = "none" 
  )

iso_per3
# Save the legend as a separate image
ggsave("Gene Figures/Iso per3.png", iso_per3, width = 4.8, height = 3.8, dpi = 300)


#Social
##################

#bmal1
soc_bmal <- ggplot() +
  geom_rect(aes(xmin = 12, xmax = 24, ymin = -Inf, ymax = Inf),fill = "grey1", alpha = .2)+
  geom_line(data=soc_gene, size=1.5, 
            aes(x=time, y=arnt_mean, group=treatment, color=treatment, linetype=treatment)) +
  geom_errorbar(data=soc_gene, 
                aes(x=time, ymin = arnt_mean-arnt_se, ymax = arnt_mean+arnt_se, color=treatment), 
                position = position_dodge(width=0.8), size = 1.5, width = 0)+
  geom_point(data=soc_gene, aes(x=time, y=arnt_mean, shape = treatment, fill=treatment, colour=treatment), 
             position = position_dodge(width=0.8), size = 5)+
  scale_shape_manual(values = c("A" = 21, "LD" = 21, "SA" = 22, "SLD" = 22)) +
  #geom_line(data=iso_plotdata, aes(x=time, y=fitY, color=factor(fitX), linetype = factor(fitX)))+
  xlab("Zeitgeber time")+
  ylab(expression(paste(italic("bmal1"), " (Normalized)"))) + # Italicize "cry1"
  #scale_y_continuous(breaks=c(0,75,150), limits=c(0,150))+
  scale_color_manual(values=c("orange", "black"))+
  scale_fill_manual(values=c("orange", "black"))+
  scale_x_continuous(breaks=c(1,7,13,19), limits=c(0,24))+
  scale_y_continuous(breaks=c(0.2, 0.4, 0.6, 0.8), limits=c(0,0.8))+
  scale_linetype_manual(values=c("solid", "dashed"))+
  theme_classic()+
  theme(
    axis.title.x = element_text( size=24, family= "Times", vjust=-0.2),
    axis.text.x  = element_text( size=16, family= "Times"),
    axis.title.y = element_text( size=24, family= "Times", vjust=2), #vjust is the distance from the text
    axis.text.y  = element_text( size=16, family= "Times"),
    legend.title = element_blank(),
    legend.text = element_text(""),
    legend.position = "none" 
  )

soc_bmal
# Save the legend as a separate image
ggsave("Gene Figures/Soc bmal.png", soc_bmal, width = 4.8, height = 3.8, dpi = 300)



#Cry
soc_cry <- ggplot() +
  geom_rect(aes(xmin = 12, xmax = 24, ymin = -Inf, ymax = Inf),fill = "grey1", alpha = .2)+
  geom_line(data=soc_gene, size=1.5, aes(x=time, y=cry_mean, group=treatment, color=treatment, linetype=treatment)) +
  geom_errorbar(data=soc_gene, 
                aes(x=time, ymin = cry_mean-cry_se, ymax = cry_mean+cry_se, color=treatment), 
                position = position_dodge(width=0.8), size = 1.5, width = 0)+
  geom_point(data=soc_gene, aes(x=time, y=cry_mean, fill=treatment, colour=treatment, shape = treatment), 
             position = position_dodge(width=0.8), size = 5)+
  scale_shape_manual(values = c("A" = 21, "LD" = 21, "SA" = 22, "SLD" = 22)) +
  #geom_line(data=iso_plotdata, aes(x=time, y=fitY, color=factor(fitX), linetype = factor(fitX)))+
  xlab("Zeitgeber time")+
  ylab(expression(paste(italic("cry1"), " (Normalized)"))) + # Italicize "cry1"
  #scale_y_continuous(breaks=c(0,75,150), limits=c(0,150))+
  scale_color_manual(values=c("orange", "black"))+
  scale_fill_manual(values=c("orange", "black"))+
  scale_x_continuous(breaks=c(1,7,13,19), limits=c(0,24))+
  scale_y_continuous(breaks=c(0.2, 0.4, 0.6, 0.8), limits=c(0,0.8))+
  scale_linetype_manual(values=c("solid", "dashed"))+
  theme_classic()+
  theme(
    axis.title.x = element_text( size=24, family= "Times", vjust=-0.2),
    axis.text.x  = element_text( size=16, family= "Times"),
    axis.title.y = element_text( size=24, family= "Times", vjust=2), #vjust is the distance from the text
    axis.text.y  = element_text( size=16, family= "Times"),
    legend.title = element_blank(),
    legend.text = element_text(""),
    legend.position = "none" 
  )

soc_cry
# Save the legend as a separate image
ggsave("Gene Figures/Soc cryl.png", soc_cry, width = 4.8, height = 3.8, dpi = 300)


#per2
soc_per2 <- ggplot() +
  geom_rect(aes(xmin = 12, xmax = 24, ymin = -Inf, ymax = Inf),fill = "grey1", alpha = .2)+
  geom_line(data=soc_gene, size=1.5, 
            aes(x=time, y=per2_mean, group=treatment, color=treatment, linetype=treatment)) +
  geom_errorbar(data=soc_gene, 
                aes(x=time, ymin = per2_mean-per2_se, ymax = per2_mean+per2_se, color=treatment), 
                position = position_dodge(width=0.8), size = 1.5, width = 0)+
  geom_point(data=soc_gene, aes(x=time, y=per2_mean, fill=treatment, colour=treatment, shape=treatment), 
             position = position_dodge(width=0.8), size = 5)+
  scale_shape_manual(values = c("A" = 21, "LD" = 21, "SA" = 22, "SLD" = 22)) +
  #geom_line(data=iso_plotdata, aes(x=time, y=fitY, color=factor(fitX), linetype = factor(fitX)))+
  xlab("Zeitgeber time")+
  ylab(expression(paste(italic("per2"), " (Normalized)"))) + # Italicize "cry1"
  #scale_y_continuous(breaks=c(0,75,150), limits=c(0,150))+
  scale_color_manual(values=c("orange", "black"))+
  scale_fill_manual(values=c("orange", "black"))+
  scale_x_continuous(breaks=c(1,7,13,19), limits=c(0,24))+
  scale_y_continuous(breaks=c(0.2, 0.4, 0.6, 0.8), limits=c(0,0.8))+
  scale_linetype_manual(values=c("solid", "dashed"))+
  theme_classic()+
  theme(
    axis.title.x = element_text( size=24, family= "Times", vjust=-0.2),
    axis.text.x  = element_text( size=16, family= "Times"),
    axis.title.y = element_text( size=24, family= "Times", vjust=2), #vjust is the distance from the text
    axis.text.y  = element_text( size=16, family= "Times"),
    legend.title = element_blank(),
    legend.text = element_text(""),
    legend.position = "none" 
  )


soc_per2
# Save the legend as a separate image
ggsave("Gene Figures/Soc per2.png", soc_per2, width = 4.8, height = 3.8, dpi = 300)


#per3
soc_per3 <- ggplot() +
  geom_rect(aes(xmin = 12, xmax = 24, ymin = -Inf, ymax = Inf),fill = "grey1", alpha = .2)+
  geom_line(data=soc_gene, size=1.5, 
            aes(x=time, y=per3_mean, group=treatment, color=treatment, linetype=treatment)) +
  geom_errorbar(data=soc_gene, 
                aes(x=time, ymin = per3_mean-per3_se, ymax = per3_mean+per3_se, color=treatment), 
                position = position_dodge(width=0.8), size = 1.5, width = 0)+
  geom_point(data=soc_gene, aes(x=time, y=per3_mean, shape = treatment, fill=treatment, colour=treatment), 
             position = position_dodge(width=0.8), size = 5)+
  scale_shape_manual(values = c("A" = 21, "LD" = 21, "SA" = 22, "SLD" = 22)) +
  #geom_line(data=iso_plotdata, aes(x=time, y=fitY, color=factor(fitX), linetype = factor(fitX)))+
  xlab("Zeitgeber time")+
  ylab(expression(paste(italic("per3"), " (Normalized)"))) + # Italicize "cry1"
  #scale_y_continuous(breaks=c(0,75,150), limits=c(0,150))+
  scale_color_manual(values=c("orange", "black"))+
  scale_fill_manual(values=c("orange", "black"))+
  scale_x_continuous(breaks=c(1,7,13,19), limits=c(0,24))+
  scale_y_continuous(breaks=c(0.2, 0.4, 0.6, 0.8), limits=c(0,0.8))+
  scale_linetype_manual(values=c("solid", "dashed"))+
  theme_classic()+
  theme(
    axis.title.x = element_text( size=24, family= "Times", vjust=-0.2),
    axis.text.x  = element_text( size=16, family= "Times"),
    axis.title.y = element_text( size=24, family= "Times", vjust=2), #vjust is the distance from the text
    axis.text.y  = element_text( size=16, family= "Times"),
    legend.title = element_blank(),
    legend.text = element_text(""),
    legend.position = "none" 
  )


soc_per3
# Save the legend as a separate image
ggsave("Gene Figures/Soc per3.png", soc_per3, width = 4.8, height = 3.8, dpi = 300)


########################
#ALAN
############

#bmal
ggplot() +
  geom_rect(aes(xmin = 12, xmax = 24, ymin = -Inf, ymax = Inf),fill = "grey1", alpha = .2)+
  geom_line(data=al_gene, size=1.5, 
            aes(x=time, y=arnt_mean, group=treatment, color=treatment, linetype=treatment)) +
  geom_errorbar(data=al_gene, 
                aes(x=time, ymin = arnt_mean-arnt_se, ymax = arnt_mean+arnt_se, color=treatment), 
                position = position_dodge(width=0.8), size = 1.5, width = 0)+
  geom_point(data=al_gene, aes(x=time, y=arnt_mean, shape = treatment, fill=treatment, colour=treatment), 
             position = position_dodge(width=0.8), size = 5)+
  scale_shape_manual(values = c("A" = 21, "LD" = 21, "SA" = 22, "SLD" = 22)) +
  #geom_line(data=iso_plotdata, aes(x=time, y=fitY, color=factor(fitX), linetype = factor(fitX)))+
  xlab("Zeitgeber time")+
  ylab(expression(paste(italic("bmal1"), " (Normalized)"))) + # Italicize "cry1"
  #scale_y_continuous(breaks=c(0,75,150), limits=c(0,150))+
  scale_color_manual(values=c("A" = "yellow", "SA" = "orange"))+
  scale_fill_manual(values=c("A" = "yellow", "SA" = "orange"))+
  scale_x_continuous(breaks=c(1,7,13,19), limits=c(0,24))+
  #scale_linetype_manual(values=c("solid", "soild"))+
  theme_classic()+
  theme(
    axis.title.x = element_text( size=24, family= "Times", vjust=-0.2),
    axis.text.x  = element_text( size=16, family= "Times"),
    axis.title.y = element_text( size=24, family= "Times", vjust=2), #vjust is the distance from the text
    axis.text.y  = element_text( size=16, family= "Times"),
    legend.title = element_blank(),
    legend.text = element_text(""),
    legend.position = "none" 
  )





#per3
ggplot() +
  geom_rect(aes(xmin = 12, xmax = 24, ymin = -Inf, ymax = Inf),fill = "grey1", alpha = .2)+
  geom_line(data=al_gene, size=1.5, 
            aes(x=time, y=per3_mean, group=treatment, color=treatment, linetype=treatment)) +
  geom_errorbar(data=al_gene, 
                aes(x=time, ymin = per3_mean-per3_se, ymax = per3_mean+per3_se, color=treatment), 
                position = position_dodge(width=0.8), size = 1.5, width = 0)+
  geom_point(data=al_gene, aes(x=time, y=per3_mean, shape = treatment, fill=treatment, colour=treatment), 
             position = position_dodge(width=0.8), size = 5)+
  scale_shape_manual(values = c("A" = 21, "LD" = 21, "SA" = 22, "SLD" = 22)) +
  #geom_line(data=iso_plotdata, aes(x=time, y=fitY, color=factor(fitX), linetype = factor(fitX)))+
  xlab("Zeitgeber time")+
  ylab(expression(paste(italic("per3"), " (Normalized)"))) + # Italicize "cry1"
  #scale_y_continuous(breaks=c(0,75,150), limits=c(0,150))+
  scale_color_manual(values=c("A" = "yellow", "SA" = "orange"))+
  scale_fill_manual(values=c("A" = "yellow", "SA" = "orange"))+
  scale_x_continuous(breaks=c(1,7,13,19), limits=c(0,24))+
  #scale_linetype_manual(values=c("A" = "dashed", "SA" = "soild"))+
  theme_classic()+
  theme(
    axis.title.x = element_text( size=24, family= "Times", vjust=-0.2),
    axis.text.x  = element_text( size=16, family= "Times"),
    axis.title.y = element_text( size=24, family= "Times", vjust=2), #vjust is the distance from the text
    axis.text.y  = element_text( size=16, family= "Times"),
    legend.title = element_blank(),
    legend.text = element_text(""),
    legend.position = "none" 
  )




#cry
ggplot() +
  geom_rect(aes(xmin = 12, xmax = 24, ymin = -Inf, ymax = Inf),fill = "grey1", alpha = .2)+
  geom_line(data=al_gene, size=1.5, 
            aes(x=time, y=cry_mean, group=treatment, color=treatment, linetype=treatment)) +
  geom_errorbar(data=al_gene, 
                aes(x=time, ymin = cry_mean-cry_se, ymax = cry_mean+cry_se, color=treatment), 
                position = position_dodge(width=0.8), size = 1.5, width = 0)+
  geom_point(data=al_gene, aes(x=time, y=cry_mean, shape = treatment, fill=treatment, colour=treatment), 
             position = position_dodge(width=0.8), size = 5)+
  scale_shape_manual(values = c("A" = 21, "LD" = 21, "SA" = 22, "SLD" = 22)) +
  #geom_line(data=iso_plotdata, aes(x=time, y=fitY, color=factor(fitX), linetype = factor(fitX)))+
  xlab("Zeitgeber time")+
  ylab(expression(paste(italic("cry1"), " (Normalized)"))) + # Italicize "cry1"
  #scale_y_continuous(breaks=c(0,75,150), limits=c(0,150))+
  scale_color_manual(values=c("A" = "yellow", "SA" = "orange"))+
  scale_fill_manual(values=c("A" = "yellow", "SA" = "orange"))+
  scale_x_continuous(breaks=c(1,7,13,19), limits=c(0,24))+
  #scale_linetype_manual(values=c("A" = "dashed", "SA" = "soild"))+
  theme_classic()+
  theme(
    axis.title.x = element_text( size=24, family= "Times", vjust=-0.2),
    axis.text.x  = element_text( size=16, family= "Times"),
    axis.title.y = element_text( size=24, family= "Times", vjust=2), #vjust is the distance from the text
    axis.text.y  = element_text( size=16, family= "Times"),
    legend.title = element_blank(),
    legend.text = element_text(""),
    legend.position = "none" 
  )

#per2
ggplot() +
  geom_rect(aes(xmin = 12, xmax = 24, ymin = -Inf, ymax = Inf),fill = "grey1", alpha = .2)+
  geom_line(data=al_gene, size=1.5, 
            aes(x=time, y=per2_mean, group=treatment, color=treatment, linetype=treatment)) +
  geom_errorbar(data=al_gene, 
                aes(x=time, ymin = per2_mean-per2_se, ymax = per2_mean+per2_se, color=treatment), 
                position = position_dodge(width=0.8), size = 1.5, width = 0)+
  geom_point(data=al_gene, aes(x=time, y=per2_mean, shape = treatment, fill=treatment, colour=treatment), 
             position = position_dodge(width=0.8), size = 5)+
  scale_shape_manual(values = c("A" = 21, "LD" = 21, "SA" = 22, "SLD" = 22)) +
  #geom_line(data=iso_plotdata, aes(x=time, y=fitY, color=factor(fitX), linetype = factor(fitX)))+
  xlab("Zeitgeber time")+
  ylab(expression(paste(italic("per2"), " (Normalized)"))) + # Italicize "cry1"
  #scale_y_continuous(breaks=c(0,75,150), limits=c(0,150))+
  scale_color_manual(values=c("A" = "yellow", "SA" = "orange"))+
  scale_fill_manual(values=c("A" = "yellow", "SA" = "orange"))+
  scale_x_continuous(breaks=c(1,7,13,19), limits=c(0,24))+
  #scale_linetype_manual(values=c("A" = "dashed", "SA" = "soild"))+
  theme_classic()+
  theme(
    axis.title.x = element_text( size=24, family= "Times", vjust=-0.2),
    axis.text.x  = element_text( size=16, family= "Times"),
    axis.title.y = element_text( size=24, family= "Times", vjust=2), #vjust is the distance from the text
    axis.text.y  = element_text( size=16, family= "Times"),
    legend.title = element_blank(),
    legend.text = element_text(""),
    legend.position = "none" 
  )


########################


#Livers
###################


iso_liv <- isolated %>% 
  group_by(time, treatment) %>% 
  summarise(cry_mean=mean(norm_L_CRY, na.rm=TRUE), cry_se=std.error(norm_L_CRY, na.rm=TRUE),  
            arnt_mean=mean(norm_L_BMAL, na.rm=TRUE), arnt_se=std.error(norm_L_BMAL, na.rm=TRUE), 
            per2_mean=mean(norm_L_PER2, na.rm=TRUE), per2_se=std.error(norm_L_PER2, na.rm=TRUE),
            per3_mean=mean(norm_L_PER3, na.rm=TRUE), per3_se=std.error(norm_L_PER3, na.rm=TRUE), n = n())


soc_liv <- social %>% 
  group_by(time, treatment) %>% 
  summarise(cry_mean=mean(norm_L_CRY, na.rm=TRUE), cry_se=std.error(norm_L_CRY, na.rm=TRUE),  
            arnt_mean=mean(norm_L_BMAL, na.rm=TRUE), arnt_se=std.error(norm_L_BMAL, na.rm=TRUE), 
            per2_mean=mean(norm_L_PER2, na.rm=TRUE), per2_se=std.error(norm_L_PER2, na.rm=TRUE), 
            per3_mean=mean(norm_L_PER3, na.rm=TRUE), per3_se=std.error(norm_L_PER3, na.rm=TRUE), n = n())



al_liv <- alan %>% 
  group_by(time, treatment) %>% 
  summarise(cry_mean=mean(norm_L_CRY, na.rm=TRUE), cry_se=std.error(norm_L_CRY, na.rm=TRUE),  
            arnt_mean=mean(norm_L_BMAL, na.rm=TRUE), arnt_se=std.error(norm_L_BMAL, na.rm=TRUE), 
            per2_mean=mean(norm_L_PER2, na.rm=TRUE), per2_se=std.error(norm_L_PER2, na.rm=TRUE), 
            per3_mean=mean(norm_L_PER3, na.rm=TRUE), per3_se=std.error(norm_L_PER3, na.rm=TRUE), n = n())



#####################

#bmal1
iso_bmal_l <- ggplot() +
  geom_rect(aes(xmin = 12, xmax = 24, ymin = -Inf, ymax = Inf),fill = "grey1", alpha = .2)+
  geom_line(data=iso_liv, size=1.5, 
            aes(x=time, y=arnt_mean, group=treatment, color=treatment, linetype=treatment)) +
  geom_errorbar(data=iso_liv, 
                aes(x=time, ymin = arnt_mean-arnt_se, ymax = arnt_mean+arnt_se, color=treatment), 
                position = position_dodge(width=0.8), size = 1.5, width = 0)+
  geom_point(data=iso_liv, aes(x=time, y=arnt_mean, fill=treatment, colour=treatment), 
             position = position_dodge(width=0.8), size = 5)+
  #geom_line(data=iso_plotdata, aes(x=time, y=fitY, color=factor(fitX), linetype = factor(fitX)))+
  xlab("Zeitgeber time")+
  ylab(expression(paste(italic("bmal1"), " (Normalized)"))) + # Italicize "cry1"
  #scale_y_continuous(breaks=c(0,75,150), limits=c(0,150))+
  scale_color_manual(values=c("#FFD700", "black"))+
  scale_x_continuous(breaks=c(1,7,13,19), limits=c(0,24))+
  scale_y_continuous(breaks=c(0.2, 0.4, 0.6, 0.8), limits=c(0,0.8))+
  scale_linetype_manual(values=c("solid", "dashed"))+
  theme_classic()+
  theme(
    axis.title.x = element_text( size=24, family= "Times", vjust=-0.2),
    axis.text.x  = element_text( size=16, family= "Times"),
    axis.title.y = element_text( size=24, family= "Times", vjust=2), #vjust is the distance from the text
    axis.text.y  = element_text( size=16, family= "Times"),
    legend.title = element_blank(),
    legend.text = element_text(""),
    legend.position = "none" 
  )

iso_bmal_l
# Save the legend as a separate image
ggsave("Gene Figures/Iso bmal liver.png", iso_bmal_l, width = 4.8, height = 3.8, dpi = 300)


##############################


#Cry
iso_cry_l <- ggplot() +
  geom_rect(aes(xmin = 12, xmax = 24, ymin = -Inf, ymax = Inf),fill = "grey1", alpha = .2)+
  geom_line(data=iso_liv, size=1.5, aes(x=time, y=cry_mean, group=treatment, color=treatment, linetype=treatment)) +
  geom_errorbar(data=iso_liv, 
                aes(x=time, ymin = cry_mean-cry_se, ymax = cry_mean+cry_se, color=treatment), 
                position = position_dodge(width=0.8), size = 1.5, width = 0)+
  geom_point(data=iso_liv, aes(x=time, y=cry_mean, fill=treatment, colour=treatment), 
             position = position_dodge(width=0.8), size = 5)+
  #geom_line(data=iso_plotdata, aes(x=time, y=fitY, color=factor(fitX), linetype = factor(fitX)))+
  xlab("Zeitgeber time")+
  ylab(expression(paste(italic("cry1"), " (Normalized)"))) + # Italicize "cry1"
  #scale_y_continuous(breaks=c(0,75,150), limits=c(0,150))+
  scale_color_manual(values=c("#FFD700", "black"))+
  scale_x_continuous(breaks=c(1,7,13,19), limits=c(0,24))+
  scale_y_continuous(breaks=c(0.2, 0.4, 0.6, 0.8), limits=c(0,0.8))+
  scale_linetype_manual(values=c("solid", "dashed"))+
  theme_classic()+
  theme(
    axis.title.x = element_text( size=24, family= "Times", vjust=-0.2),
    axis.text.x  = element_text( size=16, family= "Times"),
    axis.title.y = element_text( size=24, family= "Times", vjust=2), #vjust is the distance from the text
    axis.text.y  = element_text( size=16, family= "Times"),
    legend.title = element_blank(),
    legend.text = element_text(""),
    legend.position = "none" 
  )

iso_cry_l
# Save the legend as a separate image
ggsave("Gene Figures/Iso cry liver.png", iso_cry_l, width = 4.8, height = 3.8, dpi = 300)



#per2
iso_per2_l <- ggplot() +
  geom_rect(aes(xmin = 12, xmax = 24, ymin = -Inf, ymax = Inf),fill = "grey1", alpha = .2)+
  geom_line(data=iso_liv, size=1.5, 
            aes(x=time, y=per2_mean, group=treatment, color=treatment, linetype=treatment)) +
  geom_errorbar(data=iso_liv, 
                aes(x=time, ymin = per2_mean-per2_se, ymax = per2_mean+per2_se, color=treatment), 
                position = position_dodge(width=0.8), size = 1.5, width = 0)+
  geom_point(data=iso_liv, aes(x=time, y=per2_mean, fill=treatment, colour=treatment), 
             position = position_dodge(width=0.8), size = 5)+
  #geom_line(data=iso_plotdata, aes(x=time, y=fitY, color=factor(fitX), linetype = factor(fitX)))+
  xlab("Zeitgeber time")+
  ylab(expression(paste(italic("per2"), " (Normalized)"))) + # Italicize "cry1"
  #scale_y_continuous(breaks=c(0,75,150), limits=c(0,150))+
  scale_color_manual(values=c("#FFD700", "black"))+
  scale_x_continuous(breaks=c(1,7,13,19), limits=c(0,24))+
  scale_y_continuous(breaks=c(0.2, 0.4, 0.6, 0.8), limits=c(0,0.8))+
  scale_linetype_manual(values=c("solid", "dashed"))+
  theme_classic()+
  theme(
    axis.title.x = element_text( size=24, family= "Times", vjust=-0.2),
    axis.text.x  = element_text( size=16, family= "Times"),
    axis.title.y = element_text( size=24, family= "Times", vjust=2), #vjust is the distance from the text
    axis.text.y  = element_text( size=16, family= "Times"),
    legend.title = element_blank(),
    legend.text = element_text(""),
    legend.position = "none" 
  )

iso_per2_l
# Save the legend as a separate image
ggsave("Gene Figures/Iso per2 liver.png", iso_per2_l, width = 4.8, height = 3.8, dpi = 300)



####################################

#per3
iso_per3_l <- ggplot() +
  geom_rect(aes(xmin = 12, xmax = 24, ymin = -Inf, ymax = Inf),fill = "grey1", alpha = .2)+
  geom_line(data=iso_liv, size=1.5, 
            aes(x=time, y=per3_mean, group=treatment, color=treatment, linetype=treatment)) +
  geom_errorbar(data=iso_liv, 
                aes(x=time, ymin = per3_mean-per3_se, ymax = per3_mean+per3_se, color=treatment), 
                position = position_dodge(width=0.8), size = 1.5, width = 0)+
  geom_point(data=iso_liv, aes(x=time, y=per3_mean, fill=treatment, colour=treatment), 
             position = position_dodge(width=0.8), size = 5)+
  #geom_line(data=iso_plotdata, aes(x=time, y=fitY, color=factor(fitX), linetype = factor(fitX)))+
  xlab("Zeitgeber time")+
  ylab(expression(paste(italic("per3"), " (Normalized)"))) + # Italicize "cry1"
  #scale_y_continuous(breaks=c(0,75,150), limits=c(0,150))+
  scale_color_manual(values=c("#FFD700", "black"))+
  scale_x_continuous(breaks=c(1,7,13,19), limits=c(0,24))+
  scale_y_continuous(breaks=c(0.2, 0.4, 0.6, 0.8), limits=c(0,0.8))+
  scale_linetype_manual(values=c("solid", "dashed"))+
  theme_classic()+
  theme(
    axis.title.x = element_text( size=24, family= "Times", vjust=-0.2),
    axis.text.x  = element_text( size=16, family= "Times"),
    axis.title.y = element_text( size=24, family= "Times", vjust=2), #vjust is the distance from the text
    axis.text.y  = element_text( size=16, family= "Times"),
    legend.title = element_blank(),
    legend.text = element_text(""),
    legend.position = "none" 
  )


iso_per3_l
# Save the legend as a separate image
ggsave("Gene Figures/Iso per3 liver.png", iso_per3_l, width = 4.8, height = 3.8, dpi = 300)



#bmal1
soc_bmal_l <- ggplot() +
  geom_rect(aes(xmin = 12, xmax = 24, ymin = -Inf, ymax = Inf),fill = "grey1", alpha = .2)+
  geom_line(data=soc_liv, size=1.5, 
            aes(x=time, y=arnt_mean, group=treatment, color=treatment, linetype=treatment)) +
  geom_errorbar(data=soc_liv, 
                aes(x=time, ymin = arnt_mean-arnt_se, ymax = arnt_mean+arnt_se, color=treatment), 
                position = position_dodge(width=0.8), size = 1.5, width = 0)+
  geom_point(data=soc_liv, aes(x=time, y=arnt_mean, shape = treatment, fill=treatment, colour=treatment), 
             position = position_dodge(width=0.8), size = 5)+
  scale_shape_manual(values = c("A" = 21, "LD" = 21, "SA" = 22, "SLD" = 22)) +
  #geom_line(data=iso_plotdata, aes(x=time, y=fitY, color=factor(fitX), linetype = factor(fitX)))+
  xlab("Zeitgeber time")+
  ylab(expression(paste(italic("bmal1"), " (Normalized)"))) + # Italicize "cry1"
  #scale_y_continuous(breaks=c(0,75,150), limits=c(0,150))+
  scale_color_manual(values=c("orange", "black"))+
  scale_fill_manual(values=c("orange", "black"))+
  scale_x_continuous(breaks=c(1,7,13,19), limits=c(0,24))+
  scale_y_continuous(breaks=c(0.2, 0.4, 0.6, 0.8), limits=c(0,0.8))+
  scale_linetype_manual(values=c("solid", "dashed"))+
  theme_classic()+
  theme(
    axis.title.x = element_text( size=24, family= "Times", vjust=-0.2),
    axis.text.x  = element_text( size=16, family= "Times"),
    axis.title.y = element_text( size=24, family= "Times", vjust=2), #vjust is the distance from the text
    axis.text.y  = element_text( size=16, family= "Times"),
    legend.title = element_blank(),
    legend.text = element_text(""),
    legend.position = "none" 
  )

soc_bmal_l
# Save the legend as a separate image
ggsave("Gene Figures/Soc bmal liver.png", soc_bmal_l, width = 4.8, height = 3.8, dpi = 300)



#Cry
soc_cry_l <- ggplot() +
  geom_rect(aes(xmin = 12, xmax = 24, ymin = -Inf, ymax = Inf),fill = "grey1", alpha = .2)+
  geom_line(data=soc_liv, size=1.5, aes(x=time, y=cry_mean, group=treatment, color=treatment, linetype=treatment)) +
  geom_errorbar(data=soc_liv, 
                aes(x=time, ymin = cry_mean-cry_se, ymax = cry_mean+cry_se, color=treatment), 
                position = position_dodge(width=0.8), size = 1.5, width = 0)+
  geom_point(data=soc_liv, aes(x=time, y=cry_mean, fill=treatment, colour=treatment, shape = treatment), 
             position = position_dodge(width=0.8), size = 5)+
  scale_shape_manual(values = c("A" = 21, "LD" = 21, "SA" = 22, "SLD" = 22)) +
  #geom_line(data=iso_plotdata, aes(x=time, y=fitY, color=factor(fitX), linetype = factor(fitX)))+
  xlab("Zeitgeber time")+
  ylab(expression(paste(italic("cry1"), " (Normalized)"))) + # Italicize "cry1"
  #scale_y_continuous(breaks=c(0,75,150), limits=c(0,150))+
  scale_color_manual(values=c("orange", "black"))+
  scale_fill_manual(values=c("orange", "black"))+
  scale_x_continuous(breaks=c(1,7,13,19), limits=c(0,24))+
  scale_y_continuous(breaks=c(0.2, 0.4, 0.6, 0.8), limits=c(0,0.8))+
  scale_linetype_manual(values=c("solid", "dashed"))+
  theme_classic()+
  theme(
    axis.title.x = element_text( size=24, family= "Times", vjust=-0.2),
    axis.text.x  = element_text( size=16, family= "Times"),
    axis.title.y = element_text( size=24, family= "Times", vjust=2), #vjust is the distance from the text
    axis.text.y  = element_text( size=16, family= "Times"),
    legend.title = element_blank(),
    legend.text = element_text(""),
    legend.position = "none" 
  )

soc_cry_l
# Save the legend as a separate image
ggsave("Gene Figures/Soc cry liver.png", soc_cry_l, width = 4.8, height = 3.8, dpi = 300)


#per2
soc_per2_l <- ggplot() +
  geom_rect(aes(xmin = 12, xmax = 24, ymin = -Inf, ymax = Inf),fill = "grey1", alpha = .2)+
  geom_line(data=soc_liv, size=1.5, 
            aes(x=time, y=per2_mean, group=treatment, color=treatment, linetype=treatment)) +
  geom_errorbar(data=soc_liv, 
                aes(x=time, ymin = per2_mean-per2_se, ymax = per2_mean+per2_se, color=treatment), 
                position = position_dodge(width=0.8), size = 1.5, width = 0)+
  geom_point(data=soc_liv, aes(x=time, y=per2_mean, fill=treatment, colour=treatment, shape = treatment), 
             position = position_dodge(width=0.8), size = 5)+
  scale_shape_manual(values = c("A" = 21, "LD" = 21, "SA" = 22, "SLD" = 22)) +
  #geom_line(data=iso_plotdata, aes(x=time, y=fitY, color=factor(fitX), linetype = factor(fitX)))+
  xlab("Zeitgeber time")+
  ylab(expression(paste(italic("per2"), " (Normalized)"))) + # Italicize "cry1"
  #scale_y_continuous(breaks=c(0,75,150), limits=c(0,150))+
  scale_color_manual(values=c("orange", "black"))+
  scale_fill_manual(values=c("orange", "black"))+
  scale_x_continuous(breaks=c(1,7,13,19), limits=c(0,24))+
  scale_y_continuous(breaks=c(0.2, 0.4, 0.6, 0.8), limits=c(0,0.8))+
  scale_linetype_manual(values=c("solid", "dashed"))+
  theme_classic()+
  theme(
    axis.title.x = element_text( size=24, family= "Times", vjust=-0.2),
    axis.text.x  = element_text( size=16, family= "Times"),
    axis.title.y = element_text( size=24, family= "Times", vjust=2), #vjust is the distance from the text
    axis.text.y  = element_text( size=16, family= "Times"),
    legend.title = element_blank(),
    legend.text = element_text(""),
    legend.position = "none" 
  )

soc_per2_l
# Save the legend as a separate image
ggsave("Gene Figures/Soc per2 liver.png", soc_per2_l, width = 4.8, height = 3.8, dpi = 300)


#per3
soc_per3_l <- ggplot() +
  geom_rect(aes(xmin = 12, xmax = 24, ymin = -Inf, ymax = Inf),fill = "grey1", alpha = .2)+
  geom_line(data=soc_liv, size=1.5, 
            aes(x=time, y=per3_mean, group=treatment, color=treatment, linetype=treatment)) +
  geom_errorbar(data=soc_liv, 
                aes(x=time, ymin = per3_mean-per3_se, ymax = per3_mean+per3_se, color=treatment), 
                position = position_dodge(width=0.8), size = 1.5, width = 0)+
  geom_point(data=soc_liv, aes(x=time, y=per3_mean, shape = treatment, fill=treatment, colour=treatment), 
             position = position_dodge(width=0.8), size = 5)+
  scale_shape_manual(values = c("A" = 21, "LD" = 21, "SA" = 22, "SLD" = 22)) +
  #geom_line(data=iso_plotdata, aes(x=time, y=fitY, color=factor(fitX), linetype = factor(fitX)))+
  xlab("Zeitgeber time")+
  ylab(expression(paste(italic("per3"), " (Normalized)"))) + # Italicize "cry1"
  #scale_y_continuous(breaks=c(0,75,150), limits=c(0,150))+
  scale_color_manual(values=c("orange", "black"))+
  scale_fill_manual(values=c("orange", "black"))+
  scale_x_continuous(breaks=c(1,7,13,19), limits=c(0,24))+
  scale_y_continuous(breaks=c(0.2, 0.4, 0.6, 0.8), limits=c(0,0.8))+
  scale_linetype_manual(values=c("solid", "dashed"))+
  theme_classic()+
  theme(
    axis.title.x = element_text( size=24, family= "Times", vjust=-0.2),
    axis.text.x  = element_text( size=16, family= "Times"),
    axis.title.y = element_text( size=24, family= "Times", vjust=2), #vjust is the distance from the text
    axis.text.y  = element_text( size=16, family= "Times"),
    legend.title = element_blank(),
    legend.text = element_text(""),
    legend.position = "none" 
  )


soc_per3_l
# Save the legend as a separate image
ggsave("Gene Figures/Soc per3 liver.png", soc_per3_l, width = 4.8, height = 3.8, dpi = 300)



#ALAN
#########

#bmal
ggplot() +
  geom_rect(aes(xmin = 12, xmax = 24, ymin = -Inf, ymax = Inf),fill = "grey1", alpha = .2)+
  geom_line(data=al_liv, size=1.5, 
            aes(x=time, y=arnt_mean, group=treatment, color=treatment, linetype=treatment)) +
  geom_errorbar(data=al_liv, 
                aes(x=time, ymin = arnt_mean-arnt_se, ymax = arnt_mean+arnt_se, color=treatment), 
                position = position_dodge(width=0.8), size = 1.5, width = 0)+
  geom_point(data=al_liv, aes(x=time, y=arnt_mean, shape = treatment, fill=treatment, colour=treatment), 
             position = position_dodge(width=0.8), size = 5)+
  scale_shape_manual(values = c("A" = 21, "LD" = 21, "SA" = 22, "SLD" = 22)) +
  #geom_line(data=iso_plotdata, aes(x=time, y=fitY, color=factor(fitX), linetype = factor(fitX)))+
  xlab("Zeitgeber time")+
  ylab(expression(paste(italic("bmal1"), " (Normalized)"))) + # Italicize "cry1"
  #scale_y_continuous(breaks=c(0,75,150), limits=c(0,150))+
  scale_color_manual(values=c("A" = "yellow", "SA" = "orange"))+
  scale_fill_manual(values=c("A" = "yellow", "SA" = "orange"))+
  scale_x_continuous(breaks=c(1,7,13,19), limits=c(0,24))+
  #scale_linetype_manual(values=c("solid", "soild"))+
  theme_classic()+
  theme(
    axis.title.x = element_text( size=24, family= "Times", vjust=-0.2),
    axis.text.x  = element_text( size=16, family= "Times"),
    axis.title.y = element_text( size=24, family= "Times", vjust=2), #vjust is the distance from the text
    axis.text.y  = element_text( size=16, family= "Times"),
    legend.title = element_blank(),
    legend.text = element_text(""),
    legend.position = "none" 
  )





#per3
ggplot() +
  geom_rect(aes(xmin = 12, xmax = 24, ymin = -Inf, ymax = Inf),fill = "grey1", alpha = .2)+
  geom_line(data=al_liv, size=1.5, 
            aes(x=time, y=per3_mean, group=treatment, color=treatment, linetype=treatment)) +
  geom_errorbar(data=al_liv, 
                aes(x=time, ymin = per3_mean-per3_se, ymax = per3_mean+per3_se, color=treatment), 
                position = position_dodge(width=0.8), size = 1.5, width = 0)+
  geom_point(data=al_liv, aes(x=time, y=per3_mean, shape = treatment, fill=treatment, colour=treatment), 
             position = position_dodge(width=0.8), size = 5)+
  scale_shape_manual(values = c("A" = 21, "LD" = 21, "SA" = 22, "SLD" = 22)) +
  #geom_line(data=iso_plotdata, aes(x=time, y=fitY, color=factor(fitX), linetype = factor(fitX)))+
  xlab("Zeitgeber time")+
  ylab(expression(paste(italic("per3"), " (Normalized)"))) + # Italicize "cry1"
  #scale_y_continuous(breaks=c(0,75,150), limits=c(0,150))+
  scale_color_manual(values=c("A" = "yellow", "SA" = "orange"))+
  scale_fill_manual(values=c("A" = "yellow", "SA" = "orange"))+
  scale_x_continuous(breaks=c(1,7,13,19), limits=c(0,24))+
  #scale_linetype_manual(values=c("A" = "dashed", "SA" = "soild"))+
  theme_classic()+
  theme(
    axis.title.x = element_text( size=24, family= "Times", vjust=-0.2),
    axis.text.x  = element_text( size=16, family= "Times"),
    axis.title.y = element_text( size=24, family= "Times", vjust=2), #vjust is the distance from the text
    axis.text.y  = element_text( size=16, family= "Times"),
    legend.title = element_blank(),
    legend.text = element_text(""),
    legend.position = "none" 
  )


#cry
ggplot() +
  geom_rect(aes(xmin = 12, xmax = 24, ymin = -Inf, ymax = Inf),fill = "grey1", alpha = .2)+
  geom_line(data=al_liv, size=1.5, 
            aes(x=time, y=cry_mean, group=treatment, color=treatment, linetype=treatment)) +
  geom_errorbar(data=al_liv, 
                aes(x=time, ymin = cry_mean-cry_se, ymax = cry_mean+cry_se, color=treatment), 
                position = position_dodge(width=0.8), size = 1.5, width = 0)+
  geom_point(data=al_liv, aes(x=time, y=cry_mean, shape = treatment, fill=treatment, colour=treatment), 
             position = position_dodge(width=0.8), size = 5)+
  scale_shape_manual(values = c("A" = 21, "LD" = 21, "SA" = 22, "SLD" = 22)) +
  #geom_line(data=iso_plotdata, aes(x=time, y=fitY, color=factor(fitX), linetype = factor(fitX)))+
  xlab("Zeitgeber time")+
  ylab(expression(paste(italic("cry1"), " (Normalized)"))) + # Italicize "cry1"
  #scale_y_continuous(breaks=c(0,75,150), limits=c(0,150))+
  scale_color_manual(values=c("A" = "yellow", "SA" = "orange"))+
  scale_fill_manual(values=c("A" = "yellow", "SA" = "orange"))+
  scale_x_continuous(breaks=c(1,7,13,19), limits=c(0,24))+
  #scale_linetype_manual(values=c("A" = "dashed", "SA" = "soild"))+
  theme_classic()+
  theme(
    axis.title.x = element_text( size=24, family= "Times", vjust=-0.2),
    axis.text.x  = element_text( size=16, family= "Times"),
    axis.title.y = element_text( size=24, family= "Times", vjust=2), #vjust is the distance from the text
    axis.text.y  = element_text( size=16, family= "Times"),
    legend.title = element_blank(),
    legend.text = element_text(""),
    legend.position = "none" 
  )


#per2
ggplot() +
  geom_rect(aes(xmin = 12, xmax = 24, ymin = -Inf, ymax = Inf),fill = "grey1", alpha = .2)+
  geom_line(data=al_liv, size=1.5, 
            aes(x=time, y=per2_mean, group=treatment, color=treatment, linetype=treatment)) +
  geom_errorbar(data=al_liv, 
                aes(x=time, ymin = per2_mean-per2_se, ymax = per2_mean+per2_se, color=treatment), 
                position = position_dodge(width=0.8), size = 1.5, width = 0)+
  geom_point(data=al_liv, aes(x=time, y=per2_mean, shape = treatment, fill=treatment, colour=treatment), 
             position = position_dodge(width=0.8), size = 5)+
  scale_shape_manual(values = c("A" = 21, "LD" = 21, "SA" = 22, "SLD" = 22)) +
  #geom_line(data=iso_plotdata, aes(x=time, y=fitY, color=factor(fitX), linetype = factor(fitX)))+
  xlab("Zeitgeber time")+
  ylab(expression(paste(italic("per2"), " (Normalized)"))) + # Italicize "cry1"
  #scale_y_continuous(breaks=c(0,75,150), limits=c(0,150))+
  scale_color_manual(values=c("A" = "yellow", "SA" = "orange"))+
  scale_fill_manual(values=c("A" = "yellow", "SA" = "orange"))+
  scale_x_continuous(breaks=c(1,7,13,19), limits=c(0,24))+
  #scale_linetype_manual(values=c("A" = "dashed", "SA" = "soild"))+
  theme_classic()+
  theme(
    axis.title.x = element_text( size=24, family= "Times", vjust=-0.2),
    axis.text.x  = element_text( size=16, family= "Times"),
    axis.title.y = element_text( size=24, family= "Times", vjust=2), #vjust is the distance from the text
    axis.text.y  = element_text( size=16, family= "Times"),
    legend.title = element_blank(),
    legend.text = element_text(""),
    legend.position = "none" 
  )


####################
### Melatonin    ###
####################


##load data
fa <- read.csv('Mel/mel_nov_27.csv')
fa <- fa %>%  filter(!is.na(treatment))
##load data
sp <- read.csv('Mel/mel_spr.csv')
sp <- sp %>%  filter(!is.na(treatment))


fa <- subset(fa, var == "high")


sp$norm_mel <- (sp$melatonin - min(sp$melatonin, na.rm = TRUE)) / (max(sp$melatonin, na.rm = TRUE) - min(sp$melatonin, na.rm = TRUE))
fa$norm_mel <- (fa$melatonin - min(fa$melatonin, na.rm = TRUE)) / (max(fa$melatonin, na.rm = TRUE) - min(fa$melatonin, na.rm = TRUE))


df2 <- rbind(sp, fa)
df2 <- df2 %>%  filter(!is.na(timepoint))


###############



iso <- subset(df2, treatment == "LD" | treatment == "A")

soc <- subset(df2, treatment == "SLD" | treatment == "SA")

al <- subset(df2, treatment == "SA" | treatment == "A")

dk <- subset(df2, treatment == "SLD" | treatment == "LD")
#############


###########


iso_mel <- iso %>% 
  group_by(timepoint, treatment) %>% 
  summarise(mel_mean=mean(norm_mel, na.rm=TRUE), mel_se=std.error(norm_mel, na.rm=TRUE),  
            n = n())


soc_mel <- soc %>% 
  group_by(timepoint, treatment) %>% 
  summarise(mel_mean=mean(norm_mel, na.rm=TRUE), mel_se=std.error(norm_mel, na.rm=TRUE),  
            n = n())


dk_mel <- dk %>% 
  group_by(timepoint, treatment) %>% 
  summarise(mel_mean=mean(norm_mel, na.rm=TRUE), mel_se=std.error(norm_mel, na.rm=TRUE),  
            n = n())


al_mel <- al %>% 
  group_by(timepoint, treatment) %>% 
  summarise(mel_mean=mean(norm_mel, na.rm=TRUE), mel_se=std.error(norm_mel, na.rm=TRUE),  
            n = n())


##########


melmod <- lmer(norm_mel ~ treatment * timepoint + (1|Bird_Id), data = iso)
summary(melmod)


##########

cosinedat<-as.data.frame(iso)

cosinedat$X[which(cosinedat$treatment=="LD")]<-1
cosinedat$X[which(cosinedat$treatment=="A")]<-0

cosinedat$X<-as.integer(as.character(cosinedat$X))

fit <- cosinor.lm(norm_mel ~ time(timepoint) + X + amp.acro(X), data = cosinedat, period = 24)
summary(fit)

test_cosinor(fit, "X", param = "amp")

# Global test: 
#   Statistic: 
#   [1] 0.12
# 
# 
# P-value: 
#   [1] 0.7308

test_cosinor(fit, "X", param = "acr")


# Global test: 
#   Statistic: 
#   [1] 0.06
# 
# 
# P-value: 
#   [1] 0.8128


#creating dataset for plot
fitplot1<-ggplot_cosinor.lm(fit, x_str = "X") 
fitplot1

plotdata1<-data.frame(
  time=fitplot1$data$time,
  fitY=fitplot1$data$Y.hat,
  fitX=fitplot1$data$X) 

#Social
##################



cosinedat<-as.data.frame(soc)

cosinedat$X[which(cosinedat$treatment=="SLD")]<-1
cosinedat$X[which(cosinedat$treatment=="SA")]<-0

cosinedat$X<-as.integer(as.character(cosinedat$X))

fit <- cosinor.lm(norm_mel ~ time(timepoint) + X + amp.acro(X), data = cosinedat, period = 24)
summary(fit)

test_cosinor(fit, "X", param = "amp")
#not sig

# Global test: 
#   Statistic: 
#   [1] 0.04
# 
# 
# P-value: 
#   [1] 0.8402

test_cosinor(fit, "X", param = "acr")
#not sig


# Global test: 
#   Statistic: 
#   [1] 0.16
# 
# 
# P-value: 
#   [1] 0.6916


#creating dataset for plot
fitplot2<-ggplot_cosinor.lm(fit, x_str = "X") 
fitplot2

plotdata2<-data.frame(
  time=fitplot2$data$time,
  fitY=fitplot2$data$Y.hat,
  fitX=fitplot2$data$X) 

##################


#Figures for Melatonin
#######################


### Plot with raw data (4x4 for MS)
iso_melFIG <- ggplot() +
  geom_rect(aes(xmin = 12, xmax = 24, ymin = -Inf, ymax = Inf),fill = "grey85")+
  geom_point(data=iso, size = 3, position = position_dodge(width = 0.4), 
             aes(x=timepoint, y=norm_mel, shape = treatment, color = treatment, fill = treatment)) +
  scale_shape_manual(values = c("A" = 21, "LD" = 21, "SA" = 22, "SLD" = 22)) +
  # geom_point(data=iso, 
  #            aes(x=timepoint, y=norm_mel, fill=treatment, colour=treatment), 
  #            position = position_dodge(width=0.3), size = 4)+
  geom_line(data=plotdata1, aes(x=time, y=fitY, color=factor(fitX), linetype = factor(fitX)), size = 1.05)+
  #geom_line(data=plotdata3, aes(x=time, y=fitY, color=factor(fitX), linetype = factor(fitX)))+
  xlab("Zeitgeber time")+
  ylab("Melatonin (Normalized)")+
  #scale_y_continuous(breaks=c(0,50,100,150), limits=c(0,2000))+
  #scale_x_continuous(breaks=c(1,5,9,13,17,21), labels=c("08", "12", "16", "20", "00", "04"))+
  scale_x_continuous(breaks=c(1,7,13,19), limits=c(0,24))+
  scale_color_manual(values=c("#FFD700", "black", "#FFD700", "black")) +
  scale_fill_manual(values = c("#FFD700", "black", "#FFD700", "black")) +
  scale_linetype_manual(values=c("solid", "dashed")) +
  theme_classic()+
  theme(
    axis.title.x = element_text( size=24, family= "Times", vjust=-0.2),
    axis.text.x  = element_text( size=16, family= "Times"),
    axis.title.y = element_text( size=24, family= "Times", vjust=2), #vjust is the distance from the text
    axis.text.y  = element_text( size=16, family= "Times"),
    legend.title = element_blank(),
    legend.text = element_text(""),
    legend.position = "none" 
  )

iso_melFIG
# Save the legend as a separate image
ggsave("Gene Figures/Iso Mel.png", iso_melFIG, width = 6.3, height = 4.5, dpi = 300)



soc_melFIG <- ggplot() +
  geom_rect(aes(xmin = 12, xmax = 24, ymin = -Inf, ymax = Inf),fill = "grey85")+
  geom_point(data=iso, size = 3, position = position_dodge(width = 0.4), 
             aes(x=timepoint, y=norm_mel, shape = treatment, color = treatment, fill = treatment)) +
  scale_shape_manual(values = c(22, 22)) +
  # geom_point(data=soc, 
  #            aes(x=timepoint, y=norm_mel, fill=treatment, colour=treatment), 
  #            position = position_dodge(width=0.3), size = 4)+
  geom_line(data=plotdata2, aes(x=time, y=fitY, color=factor(fitX), linetype = factor(fitX)), size = 1.05)+
  xlab("Zeitgeber time")+
  ylab("Melatonin (Normalized)")+
  #scale_y_continuous(breaks=c(0,50,100,150), limits=c(0,2000))+
  #scale_x_continuous(breaks=c(1,5,9,13,17,21), labels=c("08", "12", "16", "20", "00", "04"))+
  scale_x_continuous(breaks=c(1,7,13,19), limits=c(0,24))+
  scale_color_manual(values=c("orange", "black", "orange", "black")) +
  scale_fill_manual(values=c("orange", "black", "orange", "black")) +
  scale_linetype_manual(values=c("solid", "dashed")) +
  theme_classic()+
  theme(
    axis.title.x = element_text( size=24, family= "Times", vjust=-0.2),
    axis.text.x  = element_text( size=16, family= "Times"),
    axis.title.y = element_text( size=24, family= "Times", vjust=2), #vjust is the distance from the text
    axis.text.y  = element_text( size=16, family= "Times"),
    legend.title = element_blank(),
    legend.text = element_text(""),
    legend.position = "none" 
  )

soc_melFIG
ggsave("Gene Figures/Soc Mel.png", soc_melFIG, width = 6.3, height = 4.5, dpi = 300)


#########################




#########################
#################
#Activity Onset and offset 


# Gather the columns into key-value pairs
dat_long <- dat %>%
  gather(key, value, starts_with("on"), starts_with("off")) %>%
  separate(key, into = c("type", "day_num"), sep = "(?<=\\D)(?=\\d)") %>% # Separate key into type and day_num
  spread(type, value) %>%  # Spread into separate columns for onset and offset
  mutate(day = as.integer(gsub("\\D", "", day_num))) %>%  # Extract day number
  select(-day_num) %>%  # Remove day_num column
  arrange(Bird_Id, day)  # Arrange by Bird_Id and day

# Convert from wide to long format
dat_long <- subset(dat_long, select = c("Bird_Id", "Cage", "Sex", "treatment", "on", "off", "day"))
dat_long$day <- as.numeric(dat_long$day)

dat_long$treatment <- relevel(dat_long$treatment, ref = "SLD")  # Set LD as reference
slm2.1 <- lmer(on ~ treatment * day + (1 | Cage), na.action = na.exclude, data = dat_long)
summary(slm2.1)

# compared to LD
# Estimate Std. Error        df t value Pr(>|t|)    
# treatmentA        -96.4603    22.4720   83.0121  -4.292 4.76e-05 ***

# Compared to SLD
#  treatmentSA     -299.0911    16.7089  715.5162 -17.900  < 2e-16 ***


# Compared to A
#  treatmentSA      -183.414     27.874   72.233  -6.580 6.43e-09 ***


anova(slm2.1)

library(multcomp)
summary(glht(slm2.1, linfct = mcp(treatment = "Tukey")), test = adjusted("bonferroni"))

dat_long$treatment <- relevel(dat_long$treatment, ref = "LD")  # Set LD as reference
offslm2.1 <- lmer(off ~ treatment + day + (1 | Cage), na.action = na.exclude, data = dat_long)
summary(offslm2.1)

## Treatment + day in LD
# Estimate Std. Error       df t value Pr(>|t|)    
#   treatmentA     65.626     20.856   24.656   3.147  0.00428 ** 

## Compared to SLD
#   treatmentSA  100.550     15.206 692.865   6.612 7.55e-11 ***

## Compared to A
#   treatmentSA    59.103     28.718  25.168   2.058  0.05008 . 

anova(offslm2.1) # of treatment * day
# the interaction is not significant so removing and re running


# Day 9
last_day <- dat_long %>%
  filter(day %in% c("9"))


m0.1 <- lmer(on ~ treatment + (1| Cage), data = last_day)
summary(m0.1)


#Figures


sfig2C <- ggplot(onset_m, aes(x = day, y = on_mean, fill = treatment)) + 
  geom_rect(aes(xmin = -Inf, xmax = Inf, ymin = -Inf, ymax = 0), fill = "grey85", alpha = 0.2) +
  scale_x_continuous("Day", limits = c(0.5, 9.5), breaks = c(1, 3, 5, 7, 9), labels = c("1", "3", "5", "7", "9")) +
  scale_y_continuous("Activity onset", breaks = c(0, -100, -200)) +
  geom_abline(slope = 0, intercept = 0, linetype = "dashed") +
  scale_colour_manual(values = c("yellow", "black", "orange", "black")) +
  scale_fill_manual(values = c("yellow", "black", "orange", "black")) +
  geom_errorbar(data = onset_m, aes(ymin = on_mean - on_ste, ymax = on_mean + on_ste, 
                                    color = treatment), 
                size = 1, width = 0, position = position_dodge(width = 0.4)) +
  geom_point(size = 3, position = position_dodge(width = 0.4), aes(shape = treatment, 
                                                                   color = treatment)) +
  scale_shape_manual(values = c("A" = 21, "LD" = 21, "SA" = 22, "SLD" = 22)) +
  annotate("text", x = 18.5, y = 30, label = "light-on", size = 8) +
  theme(axis.line = element_line(colour = "black"),
        panel.background = element_blank(),
        legend.position = "none",
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        axis.text.x = element_text(family = "Times", size = 16, colour = "black"),
        axis.text.y = element_text(family = "Times", size = 16, colour = "black"),
        axis.title.x = element_text(family = "Times", size = 24, colour = "black"),
        axis.title.y = element_text(family = "Times", size = 24, angle = 90))

sfig2C




sfig2B <- ggplot(sactidayfull2, aes(x = day, y = off.x, fill = treatment)) + 
  geom_rect(aes(xmin = -Inf, xmax = Inf, ymin = 0, ymax = Inf), fill = "grey85", alpha = 0.2) +
  scale_x_continuous("", limits = c(0.5, 9.5), breaks = c(1, 3, 5, 7, 9), labels = c("1", "3", "5", "7", "9")) +
  scale_y_continuous("Activity offset") +
  geom_abline(slope = 0, intercept = 0, linetype = "dashed") +
  scale_fill_manual(values = c("yellow", "black", "orange", "black")) +
  scale_colour_manual(values = c("yellow", "black", "orange", "black")) +
  geom_errorbar(data = sactidayfull2, aes(ymin = off.x - Offset_sem, ymax = off.x + Offset_sem, 
                                          color = treatment), 
                size = 1, width = 0, position = position_dodge(width = 0.4)) +
  geom_point(size = 3, position = position_dodge(width = 0.4), aes(shape = treatment, 
                                                                   color = treatment)) +
  scale_shape_manual(values = c("A" = 21, "LD" = 21, "SA" = 22, "SLD" = 22)) +
  theme(axis.line = element_line(colour = "black"),
        panel.background = element_blank(),
        legend.position = "none",
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        axis.text.x = element_blank(),
        axis.text.y = element_text(family = "Times", size = 16, colour = "black"),
        axis.title.x = element_blank(),
        axis.title.y = element_text(family = "Times", size = 24, angle = 90))

sfig2B


sfig2BC<-ggarrange(sfig2B, sfig2C, labels = c("", ""),ncol = 1, nrow = 2)
sfig2BC


########################
#
# Gene Expression Predicts Activity Onset
#
#######################

# Filter out ALAN from tp1
time_point <- 1
tp1 <- df %>% filter(time == time_point)

Atp1 <- subset(tp1, Light_Condition == "ALAN")

actbmod <- lm(on9 ~ norm_BMAL *  treatment, data = Atp1)
summary(actbmod)

# Coefficients:
#   Estimate Std. Error t value Pr(>|t|)  
# (Intercept)              1.442     21.216   0.068   0.9475  
# norm_BMAL             -149.413     57.259  -2.609   0.0312 *
#   treatmentSA            -12.004     27.687  -0.434   0.6761  
# norm_BMAL:treatmentSA   16.655     84.603   0.197   0.8488 

actbmod.2 <- lm(on9 ~ norm_BMAL +  treatment, data = Atp1)
summary(actbmod.2)

# Coefficients:
#   Estimate Std. Error t value Pr(>|t|)   
# (Intercept)   -0.8678    16.7047  -0.052  0.95971   
# norm_BMAL   -141.7842    39.8377  -3.559  0.00613 **
#   treatmentSA   -7.5951    15.3870  -0.494  0.63341


new.b_act <- ggplot(Atp1, aes(x = norm_BMAL, y = on9, color = treatment)) +
  geom_point(size = 5, aes(shape = treatment)) +
  geom_smooth(method = "lm", color = "black", 
              #se = FALSE, 
              aes(linetype = treatment)) +
  scale_color_manual(values = c("SA" = "orange", "SLD" = "black", "LD" = "dark grey", "A" = "#FFD700")) +
  scale_shape_manual(values = c("A" = 16, "LD" = 21, "SA" =15, "SLD" = 22)) +
  scale_linetype_manual(values = c("SA" = "dashed", "A" = "solid")) +
  labs(
    x = expression(italic("bmal1")),
    y = "") +
  theme_classic() +
  theme(
    axis.title.x = element_text( size=28, family= "Times", vjust=-0.2),
    axis.text.x  = element_text( size=22, family= "Times"),
    axis.title.y = element_text( size=24, family= "Times", vjust=2), #vjust is the distance from the text
    axis.text.y  = element_text( size=22, family= "Times"),
    legend.title = element_blank(),
    legend.text = element_text(""),
    legend.position = "none" 
  )
new.b_act
ggsave("Gene Figures/bmal_onSE.png", new.b_act, width = 4.4, height = 3.8, dpi = 600)

actcmod <- lm(on9 ~ norm_CRY *  treatment, data = Atp1)
summary(actcmod)

# Coefficients:
#   Estimate Std. Error t value Pr(>|t|)   
# (Intercept)             47.72      23.69   2.014  0.07873 . 
# norm_CRY              -173.75      41.31  -4.206  0.00297 **
#   treatmentSA            -50.24      28.82  -1.743  0.11947   
# norm_CRY:treatmentSA    90.23      51.57   1.750  0.11827 

actcmod.2 <- lm(on9 ~ norm_CRY + treatment, data = Atp1)
summary(actcmod.2)

# Coefficients:
#   Estimate Std. Error t value Pr(>|t|)   
# (Intercept)   17.217     17.782   0.968  0.35823   
# norm_CRY    -115.843     27.412  -4.226  0.00222 **
#   treatmentSA   -4.673     13.691  -0.341  0.74072

c_act <- ggplot(Atp1, aes(x = norm_CRY, y = on9)) +
  geom_point(size = 5, aes(color = treatment, shape = treatment)) +
  geom_smooth(method = "lm", 
              #se = FALSE, 
              color = "black", aes(linetype = treatment)) +
  scale_color_manual(values = c("SA" = "orange", "SLD" = "black", "LD" = "dark grey", "A" = "#FFD700")) +
  scale_shape_manual(values = c("A" = 16, "LD" = 21, "SA" =15, "SLD" = 22)) +
  scale_linetype_manual(values = c("SA" = "dashed", "A" = "solid")) +
  labs(
    x = expression(italic("cry1")),
    y = "") +
  theme_classic() +
  theme(
    axis.title.x = element_text( size=28, family= "Times", vjust=-0.2),
    axis.text.x  = element_text( size=22, family= "Times"),
    axis.title.y = element_text( size=24, family= "Times", vjust=2), #vjust is the distance from the text
    axis.text.y  = element_text( size=22, family= "Times"),
    legend.title = element_blank(),
    legend.text = element_text(""),
    legend.position = "none" 
  )
c_act
ggsave("Gene Figures/cry_onSE.png", c_act, width = 4.4, height = 3.8, dpi = 600)


actpmod <- lm(on9 ~ norm_PER2 + treatment, data = Atp1)
summary(actpmod)

# Coefficients:
#   Estimate Std. Error t value Pr(>|t|)  
# (Intercept) -133.164     26.591  -5.008   0.0153 *
#   norm_PER2    208.997     47.418   4.408   0.0217 *
#   treatmentSA   -3.858     11.806  -0.327   0.7653 

# Not enough points for isolated A
Atp1.1 <- subset(Atp1, treatment == "SA")


actpmod <- lm(on9 ~ norm_PER2, data = Atp1.1)
summary(actpmod)

# Coefficients:
#   Estimate Std. Error t value Pr(>|t|)
# (Intercept)  -107.03      38.33  -2.793    0.108
# norm_PER2     133.72      95.06   1.407    0.295



#####################
## Correlation Matrix
#####################

#separate by time
#############
#A
A_1 <- subset(data_A, time == "1")
A_13 <- subset(data_A, time == "13")


#SA
SA_1 <- subset(data_SA, time == "1")
SA_13 <- subset(data_SA, time == "13")





process_data <- function(data) {
  data_processed <- data %>%
    select(-treatment, -time, -com_gene) %>%
    mutate(gene = sub("^norm_", "", gene)) %>%
    mutate(expression = coalesce(expression, 0)) %>%
    pivot_wider(names_from = gene, values_from = expression) %>%
    column_to_rownames(var = "Bird_Id") %>%
    as.matrix()
  
  return(data_processed)
}



# Plots
A_1_m <- process_data(A_1)
res <- cor(A_1_m)

A1_cor <- corrplot(res, method = "color")
# Open a PNG device
png("Gene Figures/A1_cor.png", width = 5, height = 4.6, units = "in", res = 300)

# Create the plot
corrplot(res, method = "color")

# Close the device
dev.off()



A_13_m <- process_data(A_13)
cor_res13 <- cor(A_13_m)

corrplot(cor_res13, method = "color")
# Open a PNG device
png("Gene Figures/A13_cor.png", width = 5, height = 4.6, units = "in", res = 300)

# Create the plot
corrplot(cor_res13, method = "color")

# Close the device
dev.off()

#################################


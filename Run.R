require(cowplot)
require(dplyr)
require(stringr)
library(read.dbc)
library(xts)
require(ggplot2)
library(scales)

# reading data
load("dataSIHRJ-respiratoria.RData")

levels(dataRJ$SEXO) <- c('1','0')
data.deaths <- as.xts(dataRJ$MORTE,order.by=dataRJ$DT_INTER)
weekly.deaths <- apply.weekly(data.deaths,sum)
data.inter <- as.xts(rep(1,dim(dataRJ)[1]),order.by=dataRJ$DT_INTER)
weekly.inter <- apply.weekly(data.inter,sum)
data.cost <- as.xts(dataRJ$US_TOT,order.by=dataRJ$DT_INTER)
weekly.cost <- apply.weekly(data.cost,mean)
data.pac_diarias <- as.xts(dataRJ$DIAS_PERM,order.by=dataRJ$DT_INTER)
weekly.pac_diarias <- apply.weekly(data.pac_diarias,mean)
data.sexo <- as.xts(as.numeric(as.character(dataRJ$SEXO)),order.by=dataRJ$DT_INTER)
weekly.sexo <- apply.weekly(data.sexo,mean)
data.idade <- as.xts(dataRJ$IDADE,order.by=dataRJ$DT_INTER)
weekly.idade <- apply.weekly(data.idade,mean)
dataRJ$MUCIP_COMUM <- as.numeric(as.numeric(as.character(dataRJ$MUNIC_RES)) != as.numeric(as.character(dataRJ$MUNIC_MOV)))
data.mun_comum <- as.xts(dataRJ$MUCIP_COMUM,order.by=dataRJ$DT_INTER)
weekly.mun_comum <- apply.weekly(data.mun_comum,mean)
data.carater_intern <- as.xts(as.numeric(as.character(dataRJ$CAR_INT)),order.by=dataRJ$DT_INTER)
weekly.carater_intern <- apply.weekly(data.carater_intern,mean)

# data frame
data.week <- as.data.frame((as.vector(weekly.deaths)))
names(data.week) <- 'deaths'
data.week$hosp.week <- as.vector(weekly.inter)
data.week$cost <- as.vector((weekly.cost))
data.week$hosp_days <- as.vector((weekly.pac_diarias))
data.week$gender <- as.vector((weekly.sexo))
data.week$age <- as.vector((weekly.idade))
data.week$mun_comum <- as.vector((weekly.mun_comum))
data.week$hosp_elective <- as.vector((weekly.carater_intern))
data.week$week <- index(weekly.deaths)

data.week <- data.week[order(as.Date(data.week$week)),]
T1 <- which(data.week$week > as.Date('2020-03-11'))[1] - 1
n <- dim(data.week)[1]

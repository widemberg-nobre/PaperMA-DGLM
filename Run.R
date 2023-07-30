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


# Predicting Mediator
omega1 <- 2*pi/52
omega2 <- 2*pi/26

Gm.ajuste <- matrix(0,7,7)

Gm.ajuste[1,1]     <- 1  
Gm.ajuste[2:3,2:3] <- diag(2)
Gm.ajuste[4:5,4:5] <- G.sazonal(omega1)
Gm.ajuste[6:7,6:7] <- G.sazonal(omega2)

res1.mediator <- dglmNORMAL_FUN(F_med = cbind(1,data.week$gender,data.week$mun_comum,1,0,1,0),
                          G_med = Gm.ajuste,
                          discount = c(.98,.995),
                          idx_discount = list(2,1,2:7))

# Predicting Y[1,M(0)]
resultado1.Y_10 <- dglmBIN_FUN_1.0(F_out = cbind(1,data.week$age,data.week$gender,data.week$mun_comum),
                                    F_out_2 = cbind(1,c(res1.mediator[1:T1,1],res1.mediator[(T1+1):n,2]),
                                                    data.week$gender,data.week$mun_comum),
                          G_out = diag(4),
                          discount = c(.95),
                          idx_discount = list(1,c(1:4)))

# Predicting Y[0,M(0)]
resultado1.Y_00 <- dglmBIN_FUN(F_out = cbind(1,age.POST_INT,data.week$gender,data.week$mun_comum),
                          G_out = diag(4),
                          discount = c(.95),
                          idx_discount = list(1,c(1:4)))

# Predicting Y[0,M(1)]
resultado1.Y_01 <- dglmBIN_FUN(F_out = cbind(1,data.week$age,data.week$gender,data.week$mun_comum),
                               G_out = diag(4),
                               discount = c(.95),
                               idx_discount = list(1,c(1:4)))

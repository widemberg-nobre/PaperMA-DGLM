dglmBIN_FUN_1.0 <- function(G_out,F_out,F_out.0,discount,idx_discount,Y=data.week$deaths,N=data.week$hosp.week){
  st_dim <- dim(G_out)[1]
  dim_disc <- idx_discount[[1]]
  
  Ypred <- matrix(NA,1000,n)
  
  # initial values
  a <- m <- matrix(0,n,st_dim)
  a[1,1] <- log(mean(data.week$deaths/data.week$hosp.week)/(1-mean(data.week$deaths/data.week$hosp.week)))
  C0 <- diag(100,st_dim)
  s <- r <- f <- q <- f.star <- q.star <- c()
  R <- C <- replicate(n,data.frame())
  R[[1]] <- G_out%*%C0%*%t(G_out)
  
  # discount
  for(dd in 1:dim_disc){
    R[[1]][idx_discount[[dd+1]],idx_discount[[dd+1]]] <- 
      G_out[idx_discount[[dd+1]],idx_discount[[dd+1]]]%*%C0[idx_discount[[dd+1]],idx_discount[[dd+1]]]%*%t(G_out[idx_discount[[dd+1]],idx_discount[[dd+1]]])/discount[dd]}
  
  R[[1]] <- (t(R[[1]]) + R[[1]])/2
  f[1] <- F_out[1,]%*%a[1,]
  q[1] <- t(F_out[1,])%*%R[[1]]%*%F_out[1,]
  r[1] <- (q[1])^{-1}*(1 + exp(f[1]))
  s[1] <- (q[1])^{-1}*(1 + exp(-f[1]))
  f.star[1] <- log((r[1]+Y[1])/(s[1]+N[1]-Y[1]))
  q.star[1] <- 1/(r[1]+Y[1]) + 1/(s[1]+N[1]-Y[1])
  m[1,] <- c(a[1,]) + c(R[[1]]%*%F_out[1,]*(f.star[1] - f[1])/q[1])
  C[[1]] <- R[[1]] - (R[[1]]%*%F_out[1,]%*%t(F_out[1,])%*%R[[1]])*((1-q.star[1]/q[1])/q[1])
  
  for(ll in 1:1000){
    Theta <- mvrnorm(1,m[1,],C[[1]])
    Ypred[ll,1] <- rbinom(1,N[1],expit(c(F_out[1,]%*%Theta)))
  }
  
  # Updating
  for(t in 2:n){
    a[t,] <- G_out%*%m[t-1,]
    R[[t]] <- G_out%*%C[[t-1]]%*%t(G_out)
    
    # discount
    for(dd in 1:dim_disc){
      R[[t]][idx_discount[[dd+1]],idx_discount[[dd+1]]] <- 
        G_out[idx_discount[[dd+1]],idx_discount[[dd+1]]]%*%C[[t-1]][idx_discount[[dd+1]],idx_discount[[dd+1]]]%*%t(G_out[idx_discount[[dd+1]],idx_discount[[dd+1]]])/discount[dd]}
    f[t] <- F_out[t,]%*%a[t,]
    q[t] <- t(F_out[t,])%*%R[[t]]%*%F_out[t,]
    r[t] <- (q[t])^{-1}*(1 + exp(f[t]))
    s[t] <- (q[t])^{-1}*(1 + exp(-f[t]))
    f.star[t] <- log((r[t]+Y[t])/(s[t]+N[t]-Y[t]))
    q.star[t] <- 1/(r[t]+Y[t]) + 1/(s[t]+N[t]-Y[t])
    m[t,] <- a[t,] + R[[t]]%*%F_out[t,]*((f.star[t] - f[t])/q[t])
    C[[t]] <- R[[t]] - (R[[t]]%*%(F_out[t,]%*%t(F_out[t,]))%*%R[[t]])*((1-q.star[t]/q[t])/q[t])
    for(ll in 1:1000){
      Theta <- mvrnorm(1,m[t,],C[[t]])
      Ypred[ll,t] <- rbinom(1,N[t],expit(c(F_out.0[t,]%*%Theta)))
    }
  }
  return(Ypred)
}

dglmBIN_FUN <- function(G_out,F_out,discount,idx_discount,Y=data.week$deaths,N=data.week$hosp.week){
  st_dim <- dim(G_out)[1]
  dim_disc <- idx_discount[[1]]
  
  Ypred <- matrix(NA,1000,n)
  
  # initial values
  a <- m <- matrix(0,n,st_dim)
  a[1,1] <- log(mean(data.week$deaths/data.week$hosp.week)/(1-mean(data.week$deaths/data.week$hosp.week)))
  C0 <- diag(100,st_dim)
  s <- r <- f <- q <- f.star <- q.star <- c()
  R <- C <- replicate(n,data.frame())
  R[[1]] <- G_out%*%C0%*%t(G_out)
  
  # discount
  for(dd in 1:dim_disc){
    R[[1]][idx_discount[[dd+1]],idx_discount[[dd+1]]] <- 
      G_out[idx_discount[[dd+1]],idx_discount[[dd+1]]]%*%C0[idx_discount[[dd+1]],idx_discount[[dd+1]]]%*%t(G_out[idx_discount[[dd+1]],idx_discount[[dd+1]]])/discount[dd]}
  
  R[[1]] <- (t(R[[1]]) + R[[1]])/2
  f[1] <- F_out[1,]%*%a[1,]
  q[1] <- t(F_out[1,])%*%R[[1]]%*%F_out[1,]
  r[1] <- (q[1])^{-1}*(1 + exp(f[1]))
  s[1] <- (q[1])^{-1}*(1 + exp(-f[1]))
  f.star[1] <- log((r[1]+Y[1])/(s[1]+N[1]-Y[1]))
  q.star[1] <- 1/(r[1]+Y[1]) + 1/(s[1]+N[1]-Y[1])
  m[1,] <- c(a[1,]) + c(R[[1]]%*%F_out[1,]*(f.star[1] - f[1])/q[1])
  C[[1]] <- R[[1]] - (R[[1]]%*%F_out[1,]%*%t(F_out[1,])%*%R[[1]])*((1-q.star[1]/q[1])/q[1])
  
  for(ll in 1:1000){
    Theta <- mvrnorm(1,m[1,],C[[1]])
    Ypred[ll,1] <- rbinom(1,N[1],expit(c(F_out[1,]%*%Theta)))
  }
  
  # Updating
  for(t in 2:T1){
    a[t,] <- G_out%*%m[t-1,]
    R[[t]] <- G_out%*%C[[t-1]]%*%t(G_out)
    
    # discount
    for(dd in 1:dim_disc){
      R[[t]][idx_discount[[dd+1]],idx_discount[[dd+1]]] <- 
        G_out[idx_discount[[dd+1]],idx_discount[[dd+1]]]%*%C[[t-1]][idx_discount[[dd+1]],idx_discount[[dd+1]]]%*%t(G_out[idx_discount[[dd+1]],idx_discount[[dd+1]]])/discount[dd]}
    f[t] <- F_out[t,]%*%a[t,]
    q[t] <- t(F_out[t,])%*%R[[t]]%*%F_out[t,]
    r[t] <- (q[t])^{-1}*(1 + exp(f[t]))
    s[t] <- (q[t])^{-1}*(1 + exp(-f[t]))
    f.star[t] <- log((r[t]+Y[t])/(s[t]+N[t]-Y[t]))
    q.star[t] <- 1/(r[t]+Y[t]) + 1/(s[t]+N[t]-Y[t])
    m[t,] <- a[t,] + R[[t]]%*%F_out[t,]*((f.star[t] - f[t])/q[t])
    C[[t]] <- R[[t]] - (R[[t]]%*%(F_out[t,]%*%t(F_out[t,]))%*%R[[t]])*((1-q.star[t]/q[t])/q[t])
    for(ll in 1:1000){
      Theta <- mvrnorm(1,m[t,],C[[t]])
      Ypred[ll,t] <- rbinom(1,N[t],expit(c(F_out[t,]%*%Theta)))
    }
  }
  # Discount factor (prediction Y[0,M(0)] or Y[0,M(1)])
  P_disc <- G_out%*%C[[t]]%*%t(G_out)
  W.y <- matrix(0,st_dim,st_dim)
  for(dd in 1:dim_disc)
    W.y[idx_discount[[dd+1]],idx_discount[[dd+1]]] <- 
    P_disc[idx_discount[[dd+1]],idx_discount[[dd+1]]]*((1-discount[dd])/discount[dd])
  
  # Prediction - Y[1,M(0)]
  r2 <- s2 <- c()
  f.k <- q.k <-  c()
  a.k <- matrix(NA,n-T1,st_dim) 
  R.k <- C.k <- replicate(n-T1,data.frame())
  a.k[1,] <- G_out%*%m[T1,]
  f.k[1] <- F_out[t,]%*%a.k[1,]
  R.k[[1]] <- G_out%*%C[[T1]]%*%t(G_out) + W.y 
  q.k[1]  <- F_out[T1,]%*%R.k[[1]]%*%F_out[T1,]
  r2[1] <- (q.k[1])^{-1}*(1 + exp(f.k[1]))
  s2[1] <- (q.k[1])^{-1}*(1 + exp(-f.k[1]))
  
  k <- 1
  for(t in (T1+2):n){
    k <- k+1
    a.k[k,] <- c(G_out%*%a.k[k-1,])
    R.k[[k]] <- G_out%*%R.k[[k-1]]%*%t(G_out) + W.y
    f.k[k] <- F_out[t,]%*%a.k[k,]
    q.k[k]  <- F_out[t,]%*%R.k[[k]]%*%F_out[t,]
    r2[k] <- (q.k[k])^{-1}*(1 + exp(f.k[k]))
    s2[k] <- (q.k[k])^{-1}*(1 + exp(-f.k[k]))
  }
  
  # Predicted values
  mu.aux <- matrix(NA,1000,n)
  r.post <- c(r[1:T1],r2)
  s.post <- c(s[1:T1],s2)
  
  for(lll in 1:1000){
    for(ppp in (T1+1):n){
      mu.aux[lll,ppp] <- rbeta(1,r.post[ppp]+Ypred[lll,ppp-1],s.post[ppp] + N[ppp-1] - Ypred[lll,ppp-1])
      Ypred[lll,ppp] <- rbinom(1,N[ppp],mu.aux[lll,ppp])
    }
  }
  return(Ypred)
}

# G specification

G.sazonal <- function(x){
  return(matrix(c(cos(x),sin(x),-sin(x),cos(x)),2,2,byrow=TRUE))
}
omega1 <- 2*pi/52
omega2 <- 2*pi/26
omega3 <- 2*pi/13

G5 <- matrix(0,12,12)
G5[1:4,1:4] <- diag(4)
G5[5:6,5:6] <- matrix(c(1,1,0,1),2,2,byrow=T)
G5[7:8,7:8] <- G.sazonal(omega1)
G5[9:10,9:10] <- G.sazonal(omega2)
G5[11:12,11:12] <- G.sazonal(omega3)

# results
age.10 <- c(data.week$age[1:T1],apply(res5.m,2,mean)[(T1+1):n])
resultado1.Y_10 <- dglmBIN_FUN_1.0(F_out.0 = cbind(1,age.10,data.week$gender,data.week$mun_comum,1,0,1,0,1,0,1,0),
                                   F_out =cbind(1,data.week$age,data.week$gender,data.week$mun_comum,1,0,1,0,1,0,1,0),
                                   G_out = G5,
                                   discount = c(.98),
                                   idx_discount = list(1,c(1:12)))


resultado1.Y_01 <- dglmBIN_FUN(F_out = cbind(1,data.week$age,data.week$gender,data.week$mun_comum,1,0,1,0,1,0,1,0),
                               G_out = G5,
                               discount = c(.98),
                               idx_discount = list(1,1:12))

resultado1.Y_00 <- dglmBIN_FUN(F_out = cbind(1,age.10,data.week$gender,data.week$mun_comum,1,0,1,0,1,0,1,0),
                               G_out = G5,
                               discount = c(.98),
                               idx_discount = list(1,1:12))

library(scales)
data.plot <- as.data.frame(data.week$deaths/N)
names(data.plot) <- 'outcome' 
data.plot$time <- rep(index(weekly.deaths))
data.plot$age_obs <- c(as.vector(data.week$age))


# saving data M(0)
data.plot$M0_mean <- apply(res5.m,2,mean)
data.plot$M0_q025 <- apply(res5.m,2,quantile,.025)
data.plot$M0_q975 <- apply(res5.m,2,quantile,.975)

# saving data Y(0,M(0))
data.plot$Y00_mean <- apply(resultado1.Y_00,2,mean)/N
data.plot$Y00_q025 <- apply(resultado1.Y_00,2,quantile,.025)/N
data.plot$Y00_q975 <- apply(resultado1.Y_00,2,quantile,.975)/N

# saving data Y(0,M(1))
data.plot$Y01_mean <- apply(resultado1.Y_01,2,mean)/N
data.plot$Y01_q025 <- apply(resultado1.Y_01,2,quantile,.025)/N
data.plot$Y01_q975 <- apply(resultado1.Y_01,2,quantile,.975)/N

# saving data Y(1,M(0))
data.plot$Y10_mean <- apply(resultado1.Y_10,2,mean)/N
data.plot$Y10_q025 <- apply(resultado1.Y_10,2,quantile,.025)/N
data.plot$Y10_q975 <- apply(resultado1.Y_10,2,quantile,.975)/N

# plot predictive
theme_set(theme_classic() ) # Change the theme to my preference
Sys.setlocale("LC_TIME", "English")

graph.M0 <- ggplot(aes(x = time, y = age_obs), data = data.plot[16:n,]) +  geom_point(size=.85,color='grey1') + 
  geom_line(aes(x = time, y = M0_mean),color='grey70',size=.8) +
  xlab('') + ylab(expression(M[t])) + 
  scale_x_date(date_breaks = "4 month", date_labels =  "%b %Y",
               limits = as.Date(c("2015-03-01", "2020-10-30")))+
  scale_color_manual(values = c('grey66','grey1'))+
  coord_cartesian(ylim = c(20, 60))+
  theme(plot.title = element_text(hjust = 0.5,size = 22),
        axis.text=element_text(size=16),
        axis.title.y=element_text(size=20),
        #axis.text.x=element_text(angle=60, hjust=1),
        axis.title.x=element_blank(),
        axis.text.x=element_blank(),
        axis.ticks.x=element_blank(),
        axis.title=element_text(size=20,face="bold"),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank())+
  geom_ribbon(aes(ymin=M0_q025, ymax=M0_q975),color=NA, alpha=0.2, fill = "grey10")+
  geom_vline(xintercept = as.Date('2020-03-11'), linetype="dashed",color = "grey77", size=1.5)

library(scales)
graph.Y <- ggplot(aes(x = time, y = outcome), data = data.plot[16:n,]) +  geom_point(size=.85,color='grey1') + 
  geom_line(aes(x = time, y = Y00_mean),linetype = "solid",color='grey80',size=.8) +
  geom_line(aes(x = time, y = Y10_mean),linetype = "longdash",color='grey55',size=.8) +
  geom_line(aes(x = time, y = Y01_mean),linetype = "twodash",color='grey30',size=.8) +
  xlab('') + ylab(expression(paste(Y[t]/N[t],"   (%)"))) + 
  scale_y_continuous(labels = function(x) paste0(x*100)) +
  scale_x_date(date_breaks = "4 month", date_labels =  "%b %Y",
               limits = as.Date(c("2015-03-01", "2020-10-30")))+
  scale_color_manual(values = c('grey66','grey1'))+
  coord_cartesian(ylim = c(0.07, 0.32))+
  theme(plot.title = element_text(hjust = 0.5,size = 22),
        axis.text=element_text(size=16),
        axis.title.y=element_text(size=18),
        axis.text.x=element_text(angle=60, hjust=1),
        axis.title.x=element_text(size=16),
        #axis.text.x=element_blank(),
        axis.ticks.x=element_blank(),
        axis.title=element_text(size=20,face="bold"),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank())+
  geom_ribbon(aes(ymin=Y00_q025, ymax=Y00_q975),color=NA, alpha=0.2, fill = "grey80")+
  geom_ribbon(aes(ymin=Y10_q025, ymax=Y10_q975),color=NA, alpha=0.2, fill = "grey55")+
  geom_ribbon(aes(ymin=Y01_q025, ymax=Y01_q975),color=NA, alpha=0.2, fill = "grey30")+
  geom_vline(xintercept = as.Date('2020-03-11'), linetype="dashed",color = "grey77", size=1.5)

require(cowplot)
figure <- plot_grid(graph.M0,graph.Y, 
                    labels = c('Mediator','Outcome') ,nrow = 2,ncol = 1,hjust = -.85,vjust = 1,
                    label_size = 18, rel_heights = c(0.42, 0.58))


# Figure PW effect
                     
Y <- data.week$deaths
####
delta.IE.1 <- matrix(rep(Y,1000),1000,301,byrow=TRUE) -  resultado1.Y_10
delta.DE.1 <- matrix(rep(Y,1000),1000,301,byrow=TRUE) -  resultado1.Y_01
delta.IE.0 <- resultado1.Y_01 -  resultado1.Y_00
delta.DE.0 <- resultado1.Y_10 -  resultado1.Y_00
delta.TE <- matrix(rep(Y,1000),1000,301,byrow=TRUE) -  resultado1.Y_00


data.plot$delta.TE.mean   <- apply(delta.TE  ,2,mean)
data.plot$delta.TE.q025   <- apply(delta.TE  ,2,quantile,.025)
data.plot$delta.TE.q975   <- apply(delta.TE  ,2,quantile,.975)

data.plot$delta.IE.1.mean <- apply(delta.IE.1,2,mean)
data.plot$delta.IE.1.q025   <- apply(delta.IE.1,2,quantile,.025)
data.plot$delta.IE.1.q975   <- apply(delta.IE.1,2,quantile,.975)

data.plot$delta.DE.1.mean <- apply(delta.DE.1,2,mean)
data.plot$delta.DE.1.q025   <- apply(delta.DE.1,2,quantile,.025)
data.plot$delta.DE.1.q975   <- apply(delta.DE.1,2,quantile,.975)

data.plot$delta.IE.0.mean <- apply(delta.IE.0,2,mean)
data.plot$delta.IE.0.q025   <- apply(delta.IE.0,2,quantile,.025)
data.plot$delta.IE.0.q975   <- apply(delta.IE.0,2,quantile,.975)

data.plot$delta.DE.0.mean <- apply(delta.DE.0,2,mean)
data.plot$delta.DE.0.q025   <- apply(delta.DE.0,2,quantile,.025)
data.plot$delta.DE.0.q975   <- apply(delta.DE.0,2,quantile,.975)


theme_set(theme_classic() ) # Change the theme to my preference
Sys.setlocale("LC_TIME", "English")

graph.TE <- ggplot(aes(x = time, y = delta.TE.mean), data = data.plot[(T1+1):n,]) +  geom_line(size=.6,color='grey1') + 
  #  geom_line(aes(x = time, y = M0_mean),color='grey70',size=.8) +
  xlab('') + ylab('PW Effect') + 
  scale_x_date(date_breaks = "3 month", date_labels =  "%b %Y")+
  scale_color_manual(values = c('grey66','grey1'))+
  coord_cartesian(ylim = c(-60, 210))+
  theme(plot.title = element_text(hjust = 0.5,size = 22),
        axis.text=element_text(size=16),
        axis.title.y=element_text(size=14),
        #axis.text.x=element_text(angle=60, hjust=1),
        axis.title.x=element_blank(),
        axis.text.x=element_blank(),
        axis.ticks.x=element_blank(),
        axis.title=element_text(size=20,face="bold"),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank())+
  geom_ribbon(aes(ymin=delta.TE.q025, ymax=delta.TE.q975),color=NA, alpha=0.4, fill = "grey70")+
  geom_hline(yintercept = 0, linetype="dashed",color = "grey77", size=1.5)


graph.IE.1 <- ggplot(aes(x = time, y = delta.IE.1.mean), data = data.plot[(T1+1):n,]) +  geom_line(size=.6,color='grey1') + 
  #  geom_line(aes(x = time, y = M0_mean),color='grey70',size=.8) +
  xlab('') + ylab('PW Effect') + 
  scale_x_date(date_breaks = "3 month", date_labels =  "%b %Y")+
  scale_color_manual(values = c('grey66','grey1'))+
  coord_cartesian(ylim = c(-60, 210))+
  theme(plot.title = element_text(hjust = 0.5,size = 22),
        axis.text=element_text(size=16),
        axis.title.y=element_text(size=14),
        #axis.text.x=element_text(angle=60, hjust=1),
        axis.title.x=element_blank(),
        axis.text.x=element_blank(),
        axis.ticks.x=element_blank(),
        axis.title=element_text(size=20,face="bold"),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank())+
  geom_ribbon(aes(ymin=delta.IE.1.q025, ymax=delta.IE.1.q975),color=NA, alpha=0.4, fill = "grey70")+
  geom_hline(yintercept = 0, linetype="dashed",color = "grey77", size=1.5)


graph.DE.1 <- ggplot(aes(x = time, y = delta.DE.1.mean), data = data.plot[(T1+1):n,]) +  geom_line(size=.6,color='grey1') + 
  #  geom_line(aes(x = time, y = M0_mean),color='grey70',size=.8) +
  xlab('') + ylab('PW Effect') + 
  scale_x_date(date_breaks = "3 month", date_labels =  "%b %Y")+
  scale_color_manual(values = c('grey66','grey1'))+
  coord_cartesian(ylim = c(-60, 210))+
  theme(plot.title = element_text(hjust = 0.5,size = 22),
        axis.text=element_text(size=16),
        axis.title.y=element_text(size=14),
        #axis.text.x=element_text(angle=60, hjust=1),
        axis.title.x=element_blank(),
        axis.text.x=element_blank(),
        axis.ticks.x=element_blank(),
        axis.title=element_text(size=20,face="bold"),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank())+
  geom_ribbon(aes(ymin=delta.DE.1.q025, ymax=delta.DE.1.q975),color=NA, alpha=0.4, fill = "grey70")+
  geom_hline(yintercept = 0, linetype="dashed",color = "grey77", size=1.5)

graph.IE.0 <- ggplot(aes(x = time, y = delta.IE.0.mean), data = data.plot[(T1+1):n,]) +  geom_line(size=.6,color='grey1') + 
  #  geom_line(aes(x = time, y = M0_mean),color='grey70',size=.8) +
  xlab('') + ylab('PW Effect') + 
  scale_x_date(date_breaks = "3 month", date_labels =  "%b %Y")+
  scale_color_manual(values = c('grey66','grey1'))+
  coord_cartesian(ylim = c(-60, 210))+
  theme(plot.title = element_text(hjust = 0.5,size = 22),
        axis.text=element_text(size=16),
        axis.title.y=element_text(size=14),
        #axis.text.x=element_text(angle=60, hjust=1),
        axis.title.x=element_blank(),
        axis.text.x=element_blank(),
        axis.ticks.x=element_blank(),
        axis.title=element_text(size=20,face="bold"),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank())+
  geom_ribbon(aes(ymin=delta.IE.0.q025, ymax=delta.IE.0.q975),color=NA, alpha=0.4, fill = "grey70")+
  geom_hline(yintercept = 0, linetype="dashed",color = "grey77", size=1.5)


graph.DE.0 <- ggplot(aes(x = time, y = delta.DE.0.mean), data = data.plot[(T1+1):n,]) +  geom_line(size=.6,color='grey1') + 
  #  geom_line(aes(x = time, y = M0_mean),color='grey70',size=.8) +
  xlab('') + ylab('PW Effect') + 
  scale_x_date(date_breaks = "1 month", date_labels =  "%b %Y")+
  scale_color_manual(values = c('grey66','grey1'))+
  coord_cartesian(ylim = c(-60, 210))+
  theme(plot.title = element_text(hjust = 0.5,size = 22),
        axis.text=element_text(size=16),
        axis.title.y=element_text(size=14),
        axis.text.x=element_text(angle=60, hjust=1),
        axis.title.x=element_blank(),
        #axis.text.x=element_blank(),
        axis.ticks.x=element_blank(),
        axis.title=element_text(size=16,face="bold"),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank())+
  geom_ribbon(aes(ymin=delta.DE.0.q025, ymax=delta.DE.0.q975),color=NA, alpha=0.4, fill = "grey70")+
  geom_hline(yintercept = 0, linetype="dashed",color = "grey77", size=1.5)


require(cowplot)

figure2 <- plot_grid(graph.TE,graph.IE.0,graph.DE.1,graph.IE.1,graph.DE.0, 
                     labels = c('   TE',' IE(0)','DE(1)',' IE(1)','DE(0)') ,
                     nrow = 5,ncol = 1,hjust = -1.5,vjust = 1,
                     label_size = 18, rel_heights = c(0.17, 0.17,0.17,0.17,0.3))

# Figure Cumulative effect

aux1 <- apply(delta.TE[,(T1+1):n],1,cumsum)
data.plot$delta.TEcum.mean   <- c(rep(0,T1),apply(aux1,1,mean))
data.plot$delta.TEcum.q025   <- c(rep(0,T1),apply(aux1,1,quantile,.025))
data.plot$delta.TEcum.q975   <- c(rep(0,T1),apply(aux1,1,quantile,.975))

aux2 <- apply(delta.IE.1[,(T1+1):n],1,cumsum)
data.plot$delta.IE.1cum.mean   <- c(rep(0,T1),apply(aux2,1,mean))
data.plot$delta.IE.1cum.q025   <- c(rep(0,T1),apply(aux2,1,quantile,.025))
data.plot$delta.IE.1cum.q975   <- c(rep(0,T1),apply(aux2,1,quantile,.975))

aux3 <- apply(delta.IE.0[,(T1+1):n],1,cumsum)
data.plot$delta.IE.0cum.mean   <- c(rep(0,T1),apply(aux3,1,mean))
data.plot$delta.IE.0cum.q025   <- c(rep(0,T1),apply(aux3,1,quantile,.025))
data.plot$delta.IE.0cum.q975   <- c(rep(0,T1),apply(aux3,1,quantile,.975))

aux4 <- apply(delta.DE.1[,(T1+1):n],1,cumsum)
data.plot$delta.DE.1cum.mean   <- c(rep(0,T1),apply(aux4,1,mean))
data.plot$delta.DE.1cum.q025   <- c(rep(0,T1),apply(aux4,1,quantile,.025))
data.plot$delta.DE.1cum.q975   <- c(rep(0,T1),apply(aux4,1,quantile,.975))

aux5 <- apply(delta.DE.0[,(T1+1):n],1,cumsum)
data.plot$delta.DE.0cum.mean   <- c(rep(0,T1),apply(aux5,1,mean))
data.plot$delta.DE.0cum.q025   <- c(rep(0,T1),apply(aux5,1,quantile,.025))
data.plot$delta.DE.0cum.q975   <- c(rep(0,T1),apply(aux5,1,quantile,.975))

graph.TEcum <- ggplot(aes(x = time, y = delta.TEcum.mean), data = data.plot[(T1+1):n,]) +  geom_line(size=.6,color='grey1') + 
  #  geom_line(aes(x = time, y = M0_mean),color='grey70',size=.8) +
  xlab('') + ylab('Cum. Effect') + scale_y_continuous(breaks = c(0, 1000, 2000)) + 
  scale_x_date(date_breaks = "3 month", date_labels =  "%b %Y")+
  scale_color_manual(values = c('grey66','grey1'))+
  coord_cartesian(ylim = c(-60, 2500))+
  theme(plot.title = element_text(hjust = 0.5,size = 22),
        axis.text=element_text(size=16),
        axis.title.y=element_text(size=14),
        #axis.text.x=element_text(angle=60, hjust=1),
        axis.title.x=element_blank(),
        axis.text.x=element_blank(),
        axis.ticks.x=element_blank(),
        axis.title=element_text(size=16,face="bold"),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank())+
  geom_ribbon(aes(ymin=delta.TEcum.q025, ymax=delta.TEcum.q975),color=NA, alpha=0.4, fill = "grey70") +
  geom_hline(yintercept = 0, linetype="dashed",color = "grey77", size=1)


graph.IE.1cum <- ggplot(aes(x = time, y = delta.IE.1cum.mean), data = data.plot[(T1+1):n,]) +  geom_line(size=.6,color='grey1') + 
  #  geom_line(aes(x = time, y = M0_mean),color='grey70',size=.8) +
  xlab('') + ylab('Cum. Effect') + scale_y_continuous(breaks = c(0, 1000, 2000)) + 
  scale_x_date(date_breaks = "3 month", date_labels =  "%b %Y")+
  scale_color_manual(values = c('grey66','grey1'))+
  coord_cartesian(ylim = c(-60, 2500))+
  theme(plot.title = element_text(hjust = 0.5,size = 22),
        axis.text=element_text(size=16),
        axis.title.y=element_text(size=14),
        #axis.text.x=element_text(angle=60, hjust=1),
        axis.title.x=element_blank(),
        axis.text.x=element_blank(),
        axis.ticks.x=element_blank(),
        axis.title=element_text(size=16,face="bold"),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank())+
  geom_ribbon(aes(ymin=delta.IE.1cum.q025, ymax=delta.IE.1cum.q975),color=NA, alpha=0.4, fill = "grey70") +
  geom_hline(yintercept = 0, linetype="dashed",color = "grey77", size=1)


graph.DE.1cum <- ggplot(aes(x = time, y = delta.DE.1cum.mean), data = data.plot[(T1+1):n,]) +  geom_line(size=.6,color='grey1') + 
  #  geom_line(aes(x = time, y = M0_mean),color='grey70',size=.8) +
  xlab('') + ylab('Cum. Effect') + scale_y_continuous(breaks = c(0, 1000, 2000)) + 
  scale_x_date(date_breaks = "3 month", date_labels =  "%b %Y")+
  scale_color_manual(values = c('grey66','grey1'))+
  coord_cartesian(ylim = c(-60, 2500))+
  theme(plot.title = element_text(hjust = 0.5,size = 22),
        axis.text=element_text(size=16),
        axis.title.y=element_text(size=14),
        #axis.text.x=element_text(angle=60, hjust=1),
        axis.title.x=element_blank(),
        axis.text.x=element_blank(),
        axis.ticks.x=element_blank(),
        axis.title=element_text(size=16,face="bold"),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank())+
  geom_ribbon(aes(ymin=delta.DE.1cum.q025, ymax=delta.DE.1cum.q975),color=NA, alpha=0.4, fill = "grey70")+
  geom_hline(yintercept = 0, linetype="dashed",color = "grey77", size=1)

graph.IE.0cum <- ggplot(aes(x = time, y = delta.IE.0cum.mean), data = data.plot[(T1+1):n,]) +  geom_line(size=.6,color='grey1') + 
  #  geom_line(aes(x = time, y = M0_mean),color='grey70',size=.8) +
  xlab('') + ylab('Cum. Effect') + scale_y_continuous(breaks = c(0, 1000, 2000)) + 
  scale_x_date(date_breaks = "3 month", date_labels =  "%b %Y")+
  scale_color_manual(values = c('grey66','grey1'))+
  coord_cartesian(ylim = c(-60, 2500))+
  theme(plot.title = element_text(hjust = 0.5,size = 22),
        axis.text=element_text(size=16),
        axis.title.y=element_text(size=14),
        #axis.text.x=element_text(angle=60, hjust=1),
        axis.title.x=element_blank(),
        axis.text.x=element_blank(),
        axis.ticks.x=element_blank(),
        axis.title=element_text(size=16,face="bold"),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank())+
  geom_ribbon(aes(ymin=delta.IE.0cum.q025, ymax=delta.IE.0cum.q975),color=NA, alpha=0.4, fill = "grey70")+
  geom_hline(yintercept = 0, linetype="dashed",color = "grey77", size=1)


graph.DE.0cum <- ggplot(aes(x = time, y = delta.DE.0cum.mean), data = data.plot[(T1+1):n,]) +  geom_line(size=.6,color='grey1') + 
  #  geom_line(aes(x = time, y = M0_mean),color='grey70',size=.8) +
  xlab('') + ylab('Cum. Effect') + scale_y_continuous(breaks = c(0, 1000, 2000)) + 
  scale_x_date(date_breaks = "1 month", date_labels =  "%b %Y")+
  scale_color_manual(values = c('grey66','grey1'))+
  coord_cartesian(ylim = c(-60, 2500))+
  theme(plot.title = element_text(hjust = 0.5,size = 22),
        axis.text=element_text(size=16),
        axis.title.y=element_text(size=14),
        axis.text.x=element_text(angle=60, hjust=1),
        axis.title.x=element_blank(),
        #axis.text.x=element_blank(),
        axis.ticks.x=element_blank(),
        axis.title=element_text(size=16,face="bold"),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank())+
  geom_ribbon(aes(ymin=delta.DE.0cum.q025, ymax=delta.DE.0cum.q975),color=NA, alpha=0.4, fill = "grey70") +
  geom_hline(yintercept = 0, linetype="dashed",color = "grey77", size=1)


figure3 <- plot_grid(graph.TEcum,graph.IE.0cum,graph.DE.1cum,graph.IE.1cum,graph.DE.0cum, 
                     labels = c('   TE',' IE(0)','DE(1)',' IE(1)','DE(0)') ,
                     nrow = 5,ncol = 1,hjust = -1.5,vjust = 1,
                     label_size = 18, rel_heights = c(0.17, 0.17,0.17,0.17,0.3))


expit <- function(x){(1+exp(-x))^(-1)}

dglmBIN_FUN <- function(G_out,F_out,discount,idx_discount,Y=data.week$deaths,N=data.week$hosp.week){
  st_dim <- dim(G_out)[1]
  dim_disc <- idx_discount[[1]]
  
  Ypred <- matrix(NA,1000,n)
  
  # initial values
  a <- m <- matrix(0,T1,st_dim)
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
  return(Ypred)
}

G.sazonal <- function(x){
  return(matrix(c(cos(x),sin(x),-sin(x),cos(x)),2,2,byrow=TRUE))
}
omega1 <- 2*pi/52
omega2 <- 2*pi/26
omega3 <- 2*pi/13

G1 <- diag(4)

G2 <- matrix(0,6,6)
G2[1:4,1:4] <- diag(4)
G2[5:6,5:6] <- matrix(c(1,1,0,1),2,2,byrow=T)

G3 <- matrix(0,8,8)
G3[1:6,1:6] <- G2
G3[7:8,7:8] <- G.sazonal(omega1)

G4 <- matrix(0,10,10)
G4[1:8,1:8] <- G3
G4[9:10,9:10] <- G.sazonal(omega2)

G5 <- matrix(0,12,12)
G5[1:10,1:10] <- G4
G5[11:12,11:12] <- G.sazonal(omega3)

G6 <- G5[-c(5:6),-c(5:6)]
dim(G9)
G7 <- G6[-c(9:10),-c(9:10)]
G8 <- G7[-c(7:8),-c(7:8)]

G9 <- G6[-c(5:6),-c(5:6)]
G10 <- G9[-c(5:6),-c(5:6)]




resultado1 <- dglmBIN_FUN(F_out = cbind(1,data.week$age,data.week$gender,data.week$mun_comum),
                          G_out = G1,
                          discount = c(.95),
                          idx_discount = list(1,1:4))

resultado2 <- dglmBIN_FUN(F_out = cbind(1,data.week$age,data.week$gender,data.week$mun_comum,1,0),
                          G_out = G2,
                          discount = c(.95),
                          idx_discount = list(1,1:6))

resultado3 <- dglmBIN_FUN(F_out = cbind(1,data.week$age,data.week$gender,data.week$mun_comum,1,0,1,0),
                          G_out = G3,
                          discount = c(.95),
                          idx_discount = list(1,1:8))

resultado4 <- dglmBIN_FUN(F_out = cbind(1,data.week$age,data.week$gender,data.week$mun_comum,1,0,1,0,1,0),
                          G_out = G4,
                          discount = c(.95),
                          idx_discount = list(1,1:10))

resultado5 <- dglmBIN_FUN(F_out = cbind(1,data.week$age,data.week$gender,data.week$mun_comum,1,0,1,0,1,0,1,0),
                          G_out = G5,
                          discount = c(.95),
                          idx_discount = list(1,1:12))

resultado6 <- dglmBIN_FUN(F_out = cbind(1,data.week$age,data.week$gender,data.week$mun_comum,1,0,1,0,1,0),
                          G_out = G6,
                          discount = c(.95),
                          idx_discount = list(1,1:10))

resultado7 <- dglmBIN_FUN(F_out = cbind(1,data.week$age,data.week$gender,data.week$mun_comum,1,0,1,0),
                          G_out = G7,
                          discount = c(.95),
                          idx_discount = list(1,1:8))

resultado8 <- dglmBIN_FUN(F_out = cbind(1,data.week$age,data.week$gender,data.week$mun_comum,1,0),
                          G_out = G8,
                          discount = c(.95),
                          idx_discount = list(1,1:6))

resultado9 <- dglmBIN_FUN(F_out = cbind(1,data.week$age,data.week$gender,data.week$mun_comum,1,0,1,0),
                          G_out = G9,
                          discount = c(.95),
                          idx_discount = list(1,1:8))

resultado10 <- dglmBIN_FUN(F_out = cbind(1,data.week$age,data.week$gender,data.week$mun_comum,1,0),
                          G_out = G10,
                          discount = c(.95),
                          idx_discount = list(1,1:6))

# root mean square predictive error

# media do rMSPE
rmspe<- c(sqrt(mean(apply((resultado1[,16:T1] - matrix(rep(data.week$deaths[16:T1],1000),ncol=T1-15,byrow = TRUE))^2,2,mean))),
          sqrt(mean(apply((resultado2[,16:T1] - matrix(rep(data.week$deaths[16:T1],1000),ncol=T1-15,byrow = TRUE))^2,2,mean))),
          sqrt(mean(apply((resultado3[,16:T1] - matrix(rep(data.week$deaths[16:T1],1000),ncol=T1-15,byrow = TRUE))^2,2,mean))),
          sqrt(mean(apply((resultado4[,16:T1] - matrix(rep(data.week$deaths[16:T1],1000),ncol=T1-15,byrow = TRUE))^2,2,mean))),
          sqrt(mean(apply((resultado5[,16:T1] - matrix(rep(data.week$deaths[16:T1],1000),ncol=T1-15,byrow = TRUE))^2,2,mean))),
          sqrt(mean(apply((resultado6[,16:T1] - matrix(rep(data.week$deaths[16:T1],1000),ncol=T1-15,byrow = TRUE))^2,2,mean))),
          sqrt(mean(apply((resultado7[,16:T1] - matrix(rep(data.week$deaths[16:T1],1000),ncol=T1-15,byrow = TRUE))^2,2,mean))),
          sqrt(mean(apply((resultado8[,16:T1] - matrix(rep(data.week$deaths[16:T1],1000),ncol=T1-15,byrow = TRUE))^2,2,mean))),
          sqrt(mean(apply((resultado9[,16:T1] - matrix(rep(data.week$deaths[16:T1],1000),ncol=T1-15,byrow = TRUE))^2,2,mean))),
          sqrt(mean(apply((resultado10[,16:T1] - matrix(rep(data.week$deaths[16:T1],1000),ncol=T1-15,byrow = TRUE))^2,2,mean))))

# media das medias a posteriori do erro relativo absoluto
rae <- c(mean((apply(abs(resultado1[,16:T1] - matrix(rep(data.week$deaths[16:T1],1000),ncol=T1-15,byrow = TRUE))/matrix(rep(data.week$deaths[16:T1],1000),ncol=T1-15,byrow = TRUE),2,mean))),
         mean((apply(abs(resultado2[,16:T1] - matrix(rep(data.week$deaths[16:T1],1000),ncol=T1-15,byrow = TRUE))/matrix(rep(data.week$deaths[16:T1],1000),ncol=T1-15,byrow = TRUE),2,mean))),
         mean((apply(abs(resultado3[,16:T1] - matrix(rep(data.week$deaths[16:T1],1000),ncol=T1-15,byrow = TRUE))/matrix(rep(data.week$deaths[16:T1],1000),ncol=T1-15,byrow = TRUE),2,mean))),
         mean((apply(abs(resultado4[,16:T1] - matrix(rep(data.week$deaths[16:T1],1000),ncol=T1-15,byrow = TRUE))/matrix(rep(data.week$deaths[16:T1],1000),ncol=T1-15,byrow = TRUE),2,mean))),
         mean((apply(abs(resultado5[,16:T1] - matrix(rep(data.week$deaths[16:T1],1000),ncol=T1-15,byrow = TRUE))/matrix(rep(data.week$deaths[16:T1],1000),ncol=T1-15,byrow = TRUE),2,mean))),
         mean((apply(abs(resultado6[,16:T1] - matrix(rep(data.week$deaths[16:T1],1000),ncol=T1-15,byrow = TRUE))/matrix(rep(data.week$deaths[16:T1],1000),ncol=T1-15,byrow = TRUE),2,mean))),
         mean((apply(abs(resultado7[,16:T1] - matrix(rep(data.week$deaths[16:T1],1000),ncol=T1-15,byrow = TRUE))/matrix(rep(data.week$deaths[16:T1],1000),ncol=T1-15,byrow = TRUE),2,mean))),
         mean((apply(abs(resultado8[,16:T1] - matrix(rep(data.week$deaths[16:T1],1000),ncol=T1-15,byrow = TRUE))/matrix(rep(data.week$deaths[16:T1],1000),ncol=T1-15,byrow = TRUE),2,mean))),
         mean((apply(abs(resultado9[,16:T1] - matrix(rep(data.week$deaths[16:T1],1000),ncol=T1-15,byrow = TRUE))/matrix(rep(data.week$deaths[16:T1],1000),ncol=T1-15,byrow = TRUE),2,mean))),
         mean((apply(abs(resultado10[,16:T1] - matrix(rep(data.week$deaths[16:T1],1000),ncol=T1-15,byrow = TRUE))/matrix(rep(data.week$deaths[16:T1],1000),ncol=T1-15,byrow = TRUE),2,mean))))

# absolute relative error
cbind(round(rmspe.m,3),round(rae.m*100,3))

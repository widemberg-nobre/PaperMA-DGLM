library(MASS)
# auxiliar functions
expit <- function(x){(1+exp(-x))^(-1)}
integrand <- function(x,f00,q00) {log(1/(1+exp(x)))*dnorm(x,f00,sqrt(q00))}
func <- function(inf,sup,leng,f0,q0){
  c1 <- seq(-inf+1.01,sup,len=leng)
  c2 <- seq(-inf+.01,sup,len=leng)
  c3 <- f0
  c4 <- integrate(integrand,f00=f0,q00=q0, lower = c3-4*sqrt(q0), upper = c3+4*sqrt(q0))$value
  aa1 <- matrix(NA,leng,leng)
  for(i in 1:leng){
    for(j in 1:leng){
      if (c1[i] >= c2[j]) {aa1[i,j] <- -Inf}
      if (c1[i] <  c2[j]) {aa1[i,j] <- -(sum(c(abs(digamma(c1[i] + 1) - digamma(c2[j] - c1[i] + 1) - c3),
                                               abs(digamma(c2[j] - c1[i] + 1) - digamma(c2[j] + 2) - c4))))}
    }
  }
  return(c(c1[which(apply(aa1,1,max) == max(aa1))],c2[which(apply(aa1,2,max) == max(aa1))]))
}

# function for outcome model adustment
dglmBIN_FUN <- function(G_out,F_out,discount,idx_discount,Y=data.week$deaths,N=data.week$hosp.week){
  st_dim <- dim(G_out)[1]
  dim_disc <- idx_discount[[1]]
  
  Ypred <- expit.lamb <- matrix(NA,5000,n)
  # initial values
  a <- m <- matrix(0,n,st_dim)
  a[1,1] <- log(mean(data.week$deaths/data.week$hosp.week)/(1-mean(data.week$deaths/data.week$hosp.week)))
  C0 <- diag(0.2,st_dim)
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
  
  tauSTAR <- func(-2,1000,400,f[1],q[1]) + c(Y[1],N[1])
  f.star[1] <- digamma(tauSTAR[1]+1) -  digamma(tauSTAR[2]-tauSTAR[1]+1)
  q.star[1] <- trigamma(tauSTAR[1]+1) +  digamma(tauSTAR[2]-tauSTAR[1]+1)
  
  m[1,] <- c(a[1,]) + c(R[[1]]%*%F_out[1,]*(f.star[1] - f[1])/q[1])
  C[[1]] <- R[[1]] - (R[[1]]%*%F_out[1,]%*%t(F_out[1,])%*%R[[1]])*((1-q.star[1]/q[1])/q[1])
  
  expit.lamb[,1] <- rbeta(1000,tauSTAR[1]+1,tauSTAR[2]-tauSTAR[1]+1)
  Ypred[,1] <- rbinom(1000,N[1],expit.lamb[,1])
  
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
    
    tauSTAR <- func(-2,1000,400,f[t],q[t]) + c(Y[t],N[t])
    f.star[t] <- digamma(tauSTAR[1]+1) -  digamma(tauSTAR[2]-tauSTAR[1]+1)
    q.star[t] <- trigamma(tauSTAR[1]+1) +  digamma(tauSTAR[2]-tauSTAR[1]+1)
    
    m[t,] <- a[t,] + R[[t]]%*%F_out[t,]*((f.star[t] - f[t])/q[t])
    C[[t]] <- R[[t]] - (R[[t]]%*%(F_out[t,]%*%t(F_out[t,]))%*%R[[t]])*((1-q.star[t]/q[t])/q[t])
    
    expit.lamb[,t] <-  rbeta(1000,tauSTAR[1]+1,tauSTAR[2]-tauSTAR[1]+1)
    Ypred[,t] <- rbinom(1000,N[t],expit.lamb[,t])
  }
  
  # Discount factor (prediction Y[0,M(0)] or Y[0,M(1)])
  P_disc <- G_out%*%C[[t]]%*%t(G_out)
  W.y <- matrix(0,st_dim,st_dim)
  for(dd in 1:dim_disc)
    W.y[idx_discount[[dd+1]],idx_discount[[dd+1]]] <- 
    P_disc[idx_discount[[dd+1]],idx_discount[[dd+1]]]*((1-discount[dd])/discount[dd])
  
  # Prediction - Y[1,M(0)]
  f.k <- q.k <-  c()
  a.k <- matrix(NA,n-T1,st_dim) 
  R.k <- C.k <- replicate(n-T1,data.frame())
  a.k[1,] <- G_out%*%m[T1,]
  f.k[1] <- F_out[t,]%*%a.k[1,]
  R.k[[1]] <- G_out%*%C[[T1]]%*%t(G_out) + W.y 
  q.k[1]  <- F_out[T1,]%*%R.k[[1]]%*%F_out[T1,]
  tauSTAR <- func(-2,125000,1250,f.k[1],q.k[1])
  
  expit.lamb[,T1+1] <-  rbeta(1000,tauSTAR[1]+1,tauSTAR[2]-tauSTAR[1]+1)
  Ypred[,T1+1] <- rbinom(1000,N[T1+1],expit.lamb[,T1+1])
  
  k <- 1
  for(t in (T1+2):n){
    k <- k+1
    a.k[k,] <- c(G_out%*%a.k[k-1,])
    R.k[[k]] <- G_out%*%R.k[[k-1]]%*%t(G_out) + W.y
    f.k[k] <- F_out[t,]%*%a.k[k,]
    q.k[k]  <- F_out[t,]%*%R.k[[k]]%*%F_out[t,]
    tauSTAR <- func(-2,125000,1250,f.k[k],q.k[k])
    
    expit.lamb[,T1+k] <-  rbeta(5000,tauSTAR[1]+1,tauSTAR[2]-tauSTAR[1]+1)
    Ypred[,T1+k] <- rbinom(5000,N[T1+k],expit.lamb[,T1+k])
  }
  return(Ypred)
}

# Specification of the transition matrix G_t
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


# Results

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

resultado1 <- dglmBIN_FUN(F_out = cbind(1,data.week$age,data.week$gender,data.week$mun_comum),
                          G_out = G1,
                          discount = c(.99),
                          idx_discount = list(1,1:4))

resultado2 <- dglmBIN_FUN(F_out = cbind(1,data.week$age,data.week$gender,data.week$mun_comum,1,0),
                          G_out = G2,
                          discount = c(.99),
                          idx_discount = list(1,1:6))

resultado3 <- dglmBIN_FUN(F_out = cbind(1,data.week$age,data.week$gender,data.week$mun_comum,1,0,1,0),
                          G_out = G3,
                          discount = c(.99),
                          idx_discount = list(1,1:8))

resultado4 <- dglmBIN_FUN(F_out = cbind(1,data.week$age,data.week$gender,data.week$mun_comum,1,0,1,0,1,0),
                          G_out = G4,
                          discount = c(.99),
                          idx_discount = list(1,1:10))

resultado5 <- dglmBIN_FUN(F_out = cbind(1,data.week$age,data.week$gender,data.week$mun_comum,1,0,1,0,1,0,1,0),
                          G_out = G5,
                          discount = c(.99),
                          idx_discount = list(1,1:12))

resultado6 <- dglmBIN_FUN(F_out = cbind(1,data.week$age,data.week$gender,data.week$mun_comum,1,0,1,0,1,0),
                          G_out = G6,
                          discount = c(.99),
                          idx_discount = list(1,1:10))

resultado7 <- dglmBIN_FUN(F_out = cbind(1,data.week$age,data.week$gender,data.week$mun_comum,1,0,1,0),
                          G_out = G7,
                          discount = c(.99),
                          idx_discount = list(1,1:8))

resultado8 <- dglmBIN_FUN(F_out = cbind(1,data.week$age,data.week$gender,data.week$mun_comum,1,0),
                          G_out = G8,
                          discount = c(.99),
                          idx_discount = list(1,1:6))

resultado9 <- dglmBIN_FUN(F_out = cbind(1,data.week$age,data.week$gender,data.week$mun_comum,1,0,1,0),
                          G_out = G9,
                          discount = c(.99),
                          idx_discount = list(1,1:8))

resultado10 <- dglmBIN_FUN(F_out = cbind(1,data.week$age,data.week$gender,data.week$mun_comum,1,0),
                          G_out = G10,
                          discount = c(.99),
                          idx_discount = list(1,1:6))

# root mean square predictive error
rmspe<- c(sqrt(mean(apply((resultado1[,16:T1] - matrix(rep(data.week$deaths[16:T1],5000),ncol=T1-15,byrow = TRUE))^2,2,mean))),
          sqrt(mean(apply((resultado2[,16:T1] - matrix(rep(data.week$deaths[16:T1],5000),ncol=T1-15,byrow = TRUE))^2,2,mean))),
          sqrt(mean(apply((resultado3[,16:T1] - matrix(rep(data.week$deaths[16:T1],5000),ncol=T1-15,byrow = TRUE))^2,2,mean))),
          sqrt(mean(apply((resultado4[,16:T1] - matrix(rep(data.week$deaths[16:T1],5000),ncol=T1-15,byrow = TRUE))^2,2,mean))),
          sqrt(mean(apply((resultado5[,16:T1] - matrix(rep(data.week$deaths[16:T1],5000),ncol=T1-15,byrow = TRUE))^2,2,mean))),
          sqrt(mean(apply((resultado6[,16:T1] - matrix(rep(data.week$deaths[16:T1],5000),ncol=T1-15,byrow = TRUE))^2,2,mean))),
          sqrt(mean(apply((resultado7[,16:T1] - matrix(rep(data.week$deaths[16:T1],5000),ncol=T1-15,byrow = TRUE))^2,2,mean))),
          sqrt(mean(apply((resultado8[,16:T1] - matrix(rep(data.week$deaths[16:T1],5000),ncol=T1-15,byrow = TRUE))^2,2,mean))),
          sqrt(mean(apply((resultado9[,16:T1] - matrix(rep(data.week$deaths[16:T1],5000),ncol=T1-15,byrow = TRUE))^2,2,mean))),
          sqrt(mean(apply((resultado10[,16:T1] - matrix(rep(data.week$deaths[16:T1],5000),ncol=T1-15,byrow = TRUE))^2,2,mean))))


# Relative absolute error

rae <- c(mean((apply(abs(resultado1[,16:T1] - matrix(rep(data.week$deaths[16:T1],5000),ncol=T1-15,byrow = TRUE))/matrix(rep(data.week$deaths[16:T1],5000),ncol=T1-15,byrow = TRUE),2,mean))),
         mean((apply(abs(resultado2[,16:T1] - matrix(rep(data.week$deaths[16:T1],5000),ncol=T1-15,byrow = TRUE))/matrix(rep(data.week$deaths[16:T1],5000),ncol=T1-15,byrow = TRUE),2,mean))),
         mean((apply(abs(resultado3[,16:T1] - matrix(rep(data.week$deaths[16:T1],5000),ncol=T1-15,byrow = TRUE))/matrix(rep(data.week$deaths[16:T1],5000),ncol=T1-15,byrow = TRUE),2,mean))),
         mean((apply(abs(resultado4[,16:T1] - matrix(rep(data.week$deaths[16:T1],5000),ncol=T1-15,byrow = TRUE))/matrix(rep(data.week$deaths[16:T1],5000),ncol=T1-15,byrow = TRUE),2,mean))),
         mean((apply(abs(resultado5[,16:T1] - matrix(rep(data.week$deaths[16:T1],5000),ncol=T1-15,byrow = TRUE))/matrix(rep(data.week$deaths[16:T1],5000),ncol=T1-15,byrow = TRUE),2,mean))),
         mean((apply(abs(resultado6[,16:T1] - matrix(rep(data.week$deaths[16:T1],5000),ncol=T1-15,byrow = TRUE))/matrix(rep(data.week$deaths[16:T1],5000),ncol=T1-15,byrow = TRUE),2,mean))),
         mean((apply(abs(resultado7[,16:T1] - matrix(rep(data.week$deaths[16:T1],5000),ncol=T1-15,byrow = TRUE))/matrix(rep(data.week$deaths[16:T1],5000),ncol=T1-15,byrow = TRUE),2,mean))),
         mean((apply(abs(resultado8[,16:T1] - matrix(rep(data.week$deaths[16:T1],5000),ncol=T1-15,byrow = TRUE))/matrix(rep(data.week$deaths[16:T1],5000),ncol=T1-15,byrow = TRUE),2,mean))),
         mean((apply(abs(resultado9[,16:T1] - matrix(rep(data.week$deaths[16:T1],5000),ncol=T1-15,byrow = TRUE))/matrix(rep(data.week$deaths[16:T1],5000),ncol=T1-15,byrow = TRUE),2,mean))),
         mean((apply(abs(resultado10[,16:T1] - matrix(rep(data.week$deaths[16:T1],5000),ncol=T1-15,byrow = TRUE))/matrix(rep(data.week$deaths[16:T1],5000),ncol=T1-15,byrow = TRUE),2,mean))))

# Comparison
cbind(round(rmspe,3),round(rae*100,3))

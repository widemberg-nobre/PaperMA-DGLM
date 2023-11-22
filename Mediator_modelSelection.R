dglmNORMAL_FUN <- function(G_med,F_med,discount,idx_discount,M=data.week$age){
  st_dim <- dim(G_med)[1]
  dim_disc <- idx_discount[[1]]
  
  Mpred <- matrix(NA,500,n)
  
  # initial values
  S <- e <- f <- Q <- c()
  a <-  m <- matrix(NA,n,st_dim)
  P <- A <- replicate(n,data.frame())
  S0 <- 0.1
  C0 <- diag(10,st_dim)
  a[1,] <- c(mean(M),rep(0,st_dim-1))
  f[1] <- c((F_med[1,])%*%a[1,])
  e[1] <- (M[1]-f[1])
  R <- C <- replicate(n,data.frame())
  R[[1]] <- G_med%*%C0%*%t(G_med)
  for(dd in 1:dim_disc){
    R[[1]][idx_discount[[dd+1]],idx_discount[[dd+1]]] <- 
      G_med[idx_discount[[dd+1]],idx_discount[[dd+1]]]%*%C0[idx_discount[[dd+1]],idx_discount[[dd+1]]]%*%t(G_med[idx_discount[[dd+1]],idx_discount[[dd+1]]])/discount[dd]}
  R[[1]] <- (t(R[[1]]) + R[[1]])/2
  Q[1] <- S0 + t(F_med[1,])%*%R[[1]]%*%F_med[1,]
  A[[1]] <- c(R[[1]]%*%F_med[1,]/Q[1])
  m[1,] <- a[1,] + A[[1]]*e[1]
  S[1] <- S0 + (S0)*(e[1]^2/Q[1] - 1) 
  C[[1]] <- (R[[1]] - (A[[1]]*Q[1])%*%t(A[[1]]))*(S[1]/S0)
  nt <- (1:n)[1]
  
  for(jj in 1:500){
    V.up <- 1/rgamma(1,nt/2,nt*S[1]/2)
    Theta <- mvrnorm(1,m[1,],V.up*(C[[1]]))
    Mpred[jj,1] <- rnorm(1,c(F_med[1,])%*%Theta,sqrt(S[1]))
  }
  
  # Kalman-filter
  for(t in 2:T1){
    a[t,] <- c(G_med%*%m[t-1,])
    f[t] <- c((F_med[t,])%*%a[t,])
    e[t] <- (M[t]-f[t])
    R[[t]] <- G_med%*%C[[t-1]]%*%t(G_med)
    for(dd in 1:dim_disc){
      R[[t]][idx_discount[[dd+1]],idx_discount[[dd+1]]] <- 
        G_med[idx_discount[[dd+1]],idx_discount[[dd+1]]]%*%C[[t-1]][idx_discount[[dd+1]],idx_discount[[dd+1]]]%*%t(G_med[idx_discount[[dd+1]],idx_discount[[dd+1]]])/discount[dd]}
    Q[t] <- S[[t-1]] + t(F_med[t,])%*%R[[t]]%*%F_med[t,]
    nt <- (1:n)[t]
    S[t] <- S[t-1] + (S[t-1]/nt)*(e[t]^2/Q[t] - 1) 
    A[[t]] <- c(R[[t]]%*%F_med[t,])/Q[t]
    m[t,] <- a[t,] + A[[t]]*e[t]
    C[[t]] <- (R[[t]] - (A[[t]]*Q[t])%*%t(A[[t]]))*(S[t]/S[t-1])
    
    for(jj in 1:500){
      V.up <- 1/rgamma(1,nt/2,nt*S[t]/2)
      Theta <- mvrnorm(1,m[t,],V.up*C[[t]])
      Mpred[jj,t] <- rnorm(1,c((F_med[t,]))%*%Theta,sqrt(V.up))
    }
  }
  
  # Discount factor (prediction)
  P_disc <- G_med%*%C[[t]]%*%t(G_med)
  W.m <- matrix(0,st_dim,st_dim)
  for(dd in 1:dim_disc)
    W.m[idx_discount[[dd+1]],idx_discount[[dd+1]]] <- 
    P_disc[idx_discount[[dd+1]],idx_discount[[dd+1]]]*((1-discount[dd])/discount[dd])
  
  # Prediction
  f.k <-  c()
  a.k <-  matrix(NA,n-T1,st_dim) 
  R.k <- C.k <- replicate(n-T1,data.frame())
  a.k[1,] <- G_med%*%m[T1,]
  f.k[1] <- F_med[t,]%*%a.k[1,]
  R.k[[1]] <- G_med%*%C[[T1]]%*%t(G_med) + W.m
  Q[T1+1] <- F_med[T1+1,]%*%R.k[[1]]%*%F_med[T1+1,]
  for(jj in 1:500){
    Mpred[jj,T1+1] <- rnorm(1,f.k[1],sqrt(Q[T1+1]))
  }
  
  k <- 1
  for(t in (T1+2):n){
    k <- k+1
    a.k[k,] <- c(G_med%*%a.k[k-1,])
    R.k[[k]] <- G_med%*%R.k[[k-1]]%*%t(G_med) + W.m
    f.k[k] <- F_med[t,]%*%a.k[k,]
    Q[T1+k] <- F_med[t,]%*%R.k[[k]]%*%F_med[t,]
    for(jj in 1:500){
      Mpred[jj,t] <- rnorm(1,f.k[k],sqrt(Q[T1+k]))
    }
  }
  return(Mpred)
}



# Seasonal function
G.sazonal <- function(x){
  return(matrix(c(cos(x),sin(x),-sin(x),cos(x)),2,2,byrow=TRUE))
}
omega1 <- 2*pi/52
omega2 <- 2*pi/26
omega3 <- 2*pi/13

G1m <- diag(3)

G2m <- matrix(0,5,5)
G2m[1:3,1:3] <- diag(3)
G2m[4:5,4:5] <- matrix(c(1,1,0,1),2,2,byrow=T)

G3m <- matrix(0,7,7)
G3m[1:5,1:5] <- G2m
G3m[6:7,6:7] <- G.sazonal(omega1)

G4m <- matrix(0,9,9)
G4m[1:7,1:7] <- G3m
G4m[8:9,8:9] <- G.sazonal(omega2)

G5m <- matrix(0,11,11)
G5m[1:9,1:9] <- G4m
G5m[10:11,10:11] <- G.sazonal(omega3)

G6m <- G5m[-c(4:5),-c(4:5)]

G7m <- G6m[-c(8:9),-c(8:9)]
G8m <- G7m[-c(6:7),-c(6:7)]

G9m <- G6m[-c(4:5),-c(4:5)]
G10m <- G9m[-c(4:5),-c(4:5)]


res1.m <- dglmNORMAL_FUN(F_med = cbind(1,data.week$gender,data.week$mun_comum),
                         G_med = G1m,
                         discount = c(.98),
                         idx_discount = list(1,1:3))

res2.m <- dglmNORMAL_FUN(F_med = cbind(1,data.week$gender,data.week$mun_comum,1,0),
                         G_med = G2m,
                         discount = c(.98),
                         idx_discount = list(1,1:5))

res3.m <- dglmNORMAL_FUN(F_med = cbind(1,data.week$gender,data.week$mun_comum,1,0,1,0),
                         G_med = G3m,
                         discount = c(.98),
                         idx_discount = list(1,1:7))

res4.m <- dglmNORMAL_FUN(F_med = cbind(1,data.week$gender,data.week$mun_comum,1,0,1,0,1,0),
                         G_med = G4m,
                         discount = c(.98),
                         idx_discount = list(1,1:9))

res5.m <- dglmNORMAL_FUN(F_med = cbind(1,data.week$gender,data.week$mun_comum,1,0,1,0,1,0,1,0),
                         G_med = G5m,
                         discount = c(.98),
                         idx_discount = list(1,1:11))

res6.m <- dglmNORMAL_FUN(F_med = cbind(1,data.week$gender,data.week$mun_comum,1,0,1,0,1,0),
                         G_med = G6m,
                         discount = c(.98),
                         idx_discount = list(1,1:9))

res7.m <- dglmNORMAL_FUN(F_med = cbind(1,data.week$gender,data.week$mun_comum,1,0,1,0),
                         G_med = G7m,
                         discount = c(.98),
                         idx_discount = list(1,1:7))

res8.m <- dglmNORMAL_FUN(F_med = cbind(1,data.week$gender,data.week$mun_comum,1,0),
                         G_med = G8m,
                         discount = c(.98),
                         idx_discount = list(1,1:5))

res9.m <- dglmNORMAL_FUN(F_med = cbind(1,data.week$gender,data.week$mun_comum,1,0,1,0),
                         G_med = G9m,
                         discount = c(.98),
                         idx_discount = list(1,1:7))

res10.m <- dglmNORMAL_FUN(F_med = cbind(1,data.week$gender,data.week$mun_comum,1,0),
                         G_med = G10m,
                         discount = c(.98),
                         idx_discount = list(1,1:5))


# media do rMSPE
rmspe.m <- c(sqrt(mean(apply((res1.m[,16:T1] - matrix(rep(data.week$age[16:T1],500),ncol=T1-15,byrow = TRUE))^2,2,mean))),
             sqrt(mean(apply((res2.m[,16:T1] - matrix(rep(data.week$age[16:T1],500),ncol=T1-15,byrow = TRUE))^2,2,mean))),
             sqrt(mean(apply((res3.m[,16:T1] - matrix(rep(data.week$age[16:T1],500),ncol=T1-15,byrow = TRUE))^2,2,mean))),
             sqrt(mean(apply((res4.m[,16:T1] - matrix(rep(data.week$age[16:T1],500),ncol=T1-15,byrow = TRUE))^2,2,mean))),
             sqrt(mean(apply((res5.m[,16:T1] - matrix(rep(data.week$age[16:T1],500),ncol=T1-15,byrow = TRUE))^2,2,mean))),
             sqrt(mean(apply((res6.m[,16:T1] - matrix(rep(data.week$age[16:T1],500),ncol=T1-15,byrow = TRUE))^2,2,mean))),
             sqrt(mean(apply((res7.m[,16:T1] - matrix(rep(data.week$age[16:T1],500),ncol=T1-15,byrow = TRUE))^2,2,mean))),
             sqrt(mean(apply((res8.m[,16:T1] - matrix(rep(data.week$age[16:T1],500),ncol=T1-15,byrow = TRUE))^2,2,mean))),
             sqrt(mean(apply((res9.m[,16:T1] - matrix(rep(data.week$age[16:T1],500),ncol=T1-15,byrow = TRUE))^2,2,mean))),
             sqrt(mean(apply((res10.m[,16:T1] - matrix(rep(data.week$age[16:T1],500),ncol=T1-15,byrow = TRUE))^2,2,mean))))

# media das medias a posteriori do erro relativo absoluto
rae.m <-c(mean((apply(abs(res1.m[,16:T1] - matrix(rep(data.week$age[16:T1],500),ncol=T1-15,byrow = TRUE))/matrix(rep(data.week$age[16:T1],500),ncol=T1-15,byrow = TRUE),2,mean))),
          mean((apply(abs(res2.m[,16:T1] - matrix(rep(data.week$age[16:T1],500),ncol=T1-15,byrow = TRUE))/matrix(rep(data.week$age[16:T1],500),ncol=T1-15,byrow = TRUE),2,mean))),
          mean((apply(abs(res3.m[,16:T1] - matrix(rep(data.week$age[16:T1],500),ncol=T1-15,byrow = TRUE))/matrix(rep(data.week$age[16:T1],500),ncol=T1-15,byrow = TRUE),2,mean))),
          mean((apply(abs(res4.m[,16:T1] - matrix(rep(data.week$age[16:T1],500),ncol=T1-15,byrow = TRUE))/matrix(rep(data.week$age[16:T1],500),ncol=T1-15,byrow = TRUE),2,mean))),
          mean((apply(abs(res5.m[,16:T1] - matrix(rep(data.week$age[16:T1],500),ncol=T1-15,byrow = TRUE))/matrix(rep(data.week$age[16:T1],500),ncol=T1-15,byrow = TRUE),2,mean))),
          mean((apply(abs(res6.m[,16:T1] - matrix(rep(data.week$age[16:T1],500),ncol=T1-15,byrow = TRUE))/matrix(rep(data.week$age[16:T1],500),ncol=T1-15,byrow = TRUE),2,mean))),
          mean((apply(abs(res7.m[,16:T1] - matrix(rep(data.week$age[16:T1],500),ncol=T1-15,byrow = TRUE))/matrix(rep(data.week$age[16:T1],500),ncol=T1-15,byrow = TRUE),2,mean))),
          mean((apply(abs(res8.m[,16:T1] - matrix(rep(data.week$age[16:T1],500),ncol=T1-15,byrow = TRUE))/matrix(rep(data.week$age[16:T1],500),ncol=T1-15,byrow = TRUE),2,mean))),
          mean((apply(abs(res9.m[,16:T1] - matrix(rep(data.week$age[16:T1],500),ncol=T1-15,byrow = TRUE))/matrix(rep(data.week$age[16:T1],500),ncol=T1-15,byrow = TRUE),2,mean))),
          mean((apply(abs(res10.m[,16:T1] - matrix(rep(data.week$age[16:T1],500),ncol=T1-15,byrow = TRUE))/matrix(rep(data.week$age[16:T1],500),ncol=T1-15,byrow = TRUE),2,mean))))

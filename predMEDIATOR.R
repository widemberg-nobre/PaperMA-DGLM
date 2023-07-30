dglmNORMAL_FUN <- function(G_med,F_med,discount,idx_discount,M=data.week$age){
  st_dim <- dim(G_med)[1]
  dim_disc <- idx_discount[[1]]

  # initial values
  S <- e <- f <- Q <- c()
  a <-  m <- matrix(NA,n,st_dim)
  P <- A <- replicate(n,data.frame())
  S0 <- 0.1
  C0 <- diag(1000,7)
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
  k <- 1
  for(t in (T1+2):n){
    k <- k+1
    a.k[k,] <- c(G_med%*%a.k[k-1,])
    R.k[[k]] <- G_med%*%R.k[[k-1]]%*%t(G_med) + W.m
    f.k[k] <- F_med[t,]%*%a.k[k,]
    Q[T1+k] <- F_med[t,]%*%R.k[[k]]%*%F_med[t,]
  }
  M_0.pred <- cbind(M,c(f[1:T1],f.k),Q)
  return(M_0.pred)
}

# Seasonal function
G.sazonal <- function(x){
  return(matrix(c(cos(x),sin(x),-sin(x),cos(x)),2,2,byrow=TRUE))
}
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

plot.ts(data.week$age,ylim=c(25,60))
lines(res1.mediator[,2],col=2,lty=2,lwd=2)
lines(sqrt(res1.mediator[,3])*1.96 + res1.mediator[,2],col=4,lty=1,lwd=1)
lines(-sqrt(res1.mediator[,3])*1.96 + res1.mediator[,2],col=4,lty=1,lwd=1)


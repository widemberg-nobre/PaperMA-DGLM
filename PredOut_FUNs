dglmBIN_FUN_1.0 <- function(G_out,F_out,F_out_2,discount,idx_discount,Y=data.week$deaths,N=data.week$hosp.week){
  st_dim <- dim(G_out)[1]
  dim_disc <- idx_discount[[1]]
  
  # initial values
  a <- m <- matrix(0,n,st_dim)
  a[1,1] <- log(mean(data.week$deaths/data.week$hosp.week)/(1-mean(data.week$deaths/data.week$hosp.week)))
  C0 <- diag(100,st_dim)
  s <- r <- f <- q <- f.star <- q.star <- rep(NA,n)
  f_1.0 <- q_1.0 <- r_1.0 <- s_1.0 <- rep(NA,n)
  R <- C <- replicate(n,data.frame())
  R[[1]] <- G_out%*%C0%*%t(G_out)
  
  # discount
  for(dd in 1:dim_disc){
    R[[1]][idx_discount[[dd+1]],idx_discount[[dd+1]]] <- 
      G_out[idx_discount[[dd+1]],idx_discount[[dd+1]]]%*%C0[idx_discount[[dd+1]],idx_discount[[dd+1]]]%*%t(G_out[idx_discount[[dd+1]],idx_discount[[dd+1]]])/discount[dd]}
  
  # correcting numerical error 
  R[[1]] <- (t(R[[1]]) + R[[1]])/2
  f[1] <- F_out[1,]%*%a[1,]
  q[1] <- t(F_out[1,])%*%R[[1]]%*%F_out[1,]
  r[1] <- (q[1])^{-1}*(1 + exp(f[1]))
  s[1] <- (q[1])^{-1}*(1 + exp(-f[1]))
  f.star[1] <- log((r[1]+Y[1])/(s[1]+N[1]-Y[1]))
  q.star[1] <- 1/(r[1]+Y[1]) + 1/(s[1]+N[1]-Y[1])
  m[1,] <- c(a[1,]) + c(R[[1]]%*%F_out[1,]*(f.star[1] - f[1])/q[1])
  C[[1]] <- R[[1]] - (R[[1]]%*%F_out[1,]%*%t(F_out[1,])%*%R[[1]])*((1-q.star[1]/q[1])/q[1])
  
  # updating
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
    
 # updating LB factors for Y[1,M(0)]
    if(t > T1){
      f_1.0[t] <- F_out_2[t,]%*%a[t,]
      q_1.0[t] <- t(F_out_2[t,])%*%R[[t]]%*%F_out_2[t,]
      r_1.0[t] <- (q_1.0[t])^{-1}*(1 + exp(f_1.0[t]))
      s_1.0[t] <- (q_1.0[t])^{-1}*(1 + exp(-f_1.0[t]))
    }
  }
  
  # Predicted values
  Ypred <- mu.aux <- matrix(NA,10000,n)
  r.post <- c(r[1:T1],r_1.0[(1+T1):n])
  s.post <- c(s[1:T1],s_1.0[(1+T1):n])
  
  for(lll in 1:10000){
    for(ppp in 1:T1){
      mu.aux[lll,ppp] <- rbeta(1,r.post[ppp]+Y[ppp],s.post[ppp] + N[ppp] - Y[ppp])
      Ypred[lll,ppp] <- rbinom(1,N[ppp],mu.aux[lll,ppp])
    }
    for(ppp in (T1+1):n){
      mu.aux[lll,ppp] <- rbeta(1,r.post[ppp]+Ypred[lll,ppp-1],s.post[ppp] + N[ppp-1] - Ypred[lll,ppp-1])
      Ypred[lll,ppp] <- rbinom(1,N[ppp],mu.aux[lll,ppp])
    }
  }
  return(Ypred)
}


dglmBIN_FUN <- function(G_out,F_out,discount,idx_discount,Y=data.week$deaths,N=data.week$hosp.week){
  st_dim <- dim(G_out)[1]
  dim_disc <- idx_discount[[1]]
  
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
  Ypred <- mu.aux <- matrix(NA,10000,n)
  r.post <- c(r[1:T1],r2)
  s.post <- c(s[1:T1],s2)

  for(lll in 1:10000){
    for(ppp in 1:T1){
      mu.aux[lll,ppp] <- rbeta(1,r.post[ppp]+Y[ppp],s.post[ppp] + N[ppp] - Y[ppp])
      Ypred[lll,ppp] <- rbinom(1,N[ppp],mu.aux[lll,ppp])
    }
    for(ppp in (T1+1):n){
      mu.aux[lll,ppp] <- rbeta(1,r.post[ppp]+Ypred[lll,ppp-1],s.post[ppp] + N[ppp-1] - Ypred[lll,ppp-1])
      Ypred[lll,ppp] <- rbinom(1,N[ppp],mu.aux[lll,ppp])
    }
  }
  return(Ypred)
}

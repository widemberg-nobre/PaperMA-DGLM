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

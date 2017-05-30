KMoEM <- function(dim.s,dim.r,dim.x,dim.y,TT,p,initpar = NULL,FSL = NULL,Expect = NULL, 
                  parfD = NULL, eta = NULL,alpha=NULL,maxit = 1000, tol = 1e-2)
{
  #Now only suitable for lag <= 2
  if(is.null(initpar)| is.null(FSL) | is.null(parfD) | is.null(Expect))
  {
    stop("Need Parameters or Data")
  }
  if(is.null(eta)|is.null(alpha))
  {
    stop("Need to set the learning rate")
  }
  #Fixed parts and free parts of each parameter
  fb <- as.matrix(parfD$B$fixed[,,1])
  Db <- as.matrix(parfD$B$free[,,1])
  fz <- as.matrix(parfD$Z$fixed[,,1])
  Dz <- as.matrix(parfD$Z$free[,,1])
  fq <- as.matrix(parfD$augQ$fixed[,,1])
  Dq <- as.matrix(parfD$augQ$free[,,1])
  fr <- as.matrix(parfD$R$fixed[,,1])
  Dr <- as.matrix(parfD$R$free[,,1])
  #Initial values of each parameter
  #B0 stands for the old estimate
  #B would vary it's the new estimate
  B0 <- B <- initpar$B
  Z0 <- Z <- initpar$Z
  Q0 <- Q <- initpar$Q
  R0 <- R <- initpar$R
  # H <- initpar$H
  G <- initpar$G
  x00 <- initpar$x00
  V00 <- initpar$V00
  #Matrices to extract free parts
  #Should be modified to fit situations where p>2
  OMGp_q <- trans(G)
  # OMGp_r <- trans(H)
  #Get the predicted, filtered, smoothed values
  #and lag on covariance matrix
  xtT <- FSL$xtT
  Vtt_1 <- FSL$Vtt_1
  Vtt_1T <- FSL$Vtt_1T
  VtT <- FSL$VtT
  #Extract some important expectation
  # list(Eytyt = Eytyt,Eytxt = Eytxt,Extxt = Extxt, Extxt_1 = Extxt_1,
  #      sumytyt = sumytyt, sumytxt = sumytxt, sumxtxt = sumxtxt,
  #      sumxt_1xt_1 = sumxt_1xt_1, sumxtxt_1 = sumxtxt_1)
  Eytyt <- Expect$Eytyt
  Eytxt <- Expect$Eytxt
  Extxt <- Expect$Extxt
  Extxt_1 <- Expect$Extxt_1
  sumytyt <- Expect$sumytyt 
  sumytxt <- Expect$sumytxt
  sumxtxt <- Expect$sumxtxt
  sumxt_1xt_1 <- Expect$sumxt_1xt_1
  sumxtxt_1 <- Expect$sumxtxt_1
  #
  
  dif <- Inf
  iter <- 0
  while((dif > tol)&(iter <= maxit))
  {
    #UPDATE R
    W <- sumytyt - tcrossprod(sumytxt ,Z) - 
      tcrossprod(Z,sumytxt) + tcrossprod(Z %*% sumxtxt, Z)
    rpre  <- r <- R[lower.tri(R,diag = TRUE)]
    r0 <- R0[lower.tri(R0,diag = TRUE)]
    DrTDr <- crossprod(Dr)
    r <- crossprod(Dr,vec(W)) + 4*eta* DrTDr%*% r0
    r <- 1/(TT+4*eta)*pcholinv(DrTDr) %*% r
    R[lower.tri(R,diag = TRUE)] <- r
    R[upper.tri(R)] <- trans(R)[upper.tri(trans(R))]
    diftemp <- norm(r-rpre,type="F")
    if(diftemp < dif)
    {
      dif <- diftemp
    }
    # Riter <- 0
    # Rdif <- Inf
    # Rpre <- Rtemp <- R
    # while((Riter <= maxit) & (Rdif > tol))
    # {
    #   dR <- -TT/2 * pcholinv(R) + 0.5 * 
    #     tcrossprod(pcholinv(R) %*%  W, pcholinv(R) )-
    #     2 * eta * (R-R0)
    #   Rtemp <- R + alpha * dR
    #   Riter <- Riter + 1
    #   Rdif <- norm(alpha * dR,type="F")
    #   R <- Rtemp
    # }
    # diftemp <- norm(R-Rpre,type="F")
    # if(diftemp < dif)
    # {
    #   dif <- diftemp
    # }
    cat("R is updated\n")
    #UPDATE Q
    S <- sumxtxt - tcrossprod(sumxtxt_1,B) - 
      tcrossprod(B,sumxtxt_1) + tcrossprod(B%*%sumxt_1xt_1,B)
    qpre  <- q <- Q[lower.tri(Q,diag = TRUE)]
    q0 <- Q0[lower.tri(Q0,diag = TRUE)]
    DqTDq <- crossprod(Dq)
    q <- crossprod(Dq,vec(S)) + 4*eta* DqTDq%*% q0
    q <- 1/(TT+4*eta)*pcholinv(DqTDq) %*% q
    Q[lower.tri(Q,diag = TRUE)] <- q
    Q[upper.tri(Q)] <- trans(Q)[upper.tri(trans(Q))]
    diftemp <- norm(q-qpre,type="F")
    if(diftemp < dif)
    {
      dif <- diftemp
    }
    # Qiter <- 0
    # Qdif <- Inf
    # Qpre <- Qtemp <- Q
    # while((Qiter <= maxit) & (Qdif > tol))
    # {
    #   dQ <- -TT/2 * pcholinv(Q) + 0.5 * tcrossprod(pcholinv(Q) %*% OMGp_q %*%S,
    #                                   pcholinv(Q) %*% OMGp_q)- 2 * eta * (Q-Q0)
    #   Qtemp <- Q + alpha * dQ
    #   Qiter <- Qiter + 1
    #   Qdif <- norm(alpha * dQ,type="F")
    #   Q <- Qtemp
    # }
    # diftemp <- norm(Q-Qpre,type="F")
    # if(diftemp < dif)
    # {
    #   dif <- diftemp
    # }
    cat("Q is updated\n")
    #UPDATE B
    betapre <- beta <- vec(OMGp_q %*% B)
    beta0 <- vec(OMGp_q %*% B0)
    Qstar <- crossprod(crossprod(pcholinv(Q),OMGp_q),OMGp_q)
    A <- matrix(0,ncol = 1,nrow = nrow(fb))
    A1 <- matrix(0, nrow = (dim.x)^2 , ncol = (dim.x)^2)
    for(t in 1:TT)
    {
      A <- A + vec(Qstar %*% Extxt_1[,,t])-kronecker(Extxt[,,t],Qstar)%*%fb
      A1 <- A1 + kronecker(Extxt[,,t], Qstar)
    }
    beta <- crossprod(Db, A+2*eta * Db %*% beta0)
    beta <- pcholinv((crossprod(Db, A1)%*% Db +2*eta*crossprod(Db))) %*% beta
    diftemp <- norm(beta-betapre,type="F")
    if(diftemp < dif)
    {
      dif <- diftemp
    }
    B[1:dim.s,] <- matrix(beta, nrow = dim.s, ncol = dim.x) 
    cat("B is updated\n")
    #UPDATE Z
    zetapre <- zeta <- vec(Z)
    zeta0 <- vec(Z0)
    Rstar <- pcholinv(R)
    L <- matrix(0,ncol = 1,nrow = nrow(fz))
    L1 <- matrix(0, nrow = dim.y *dim.x , ncol = dim.y *dim.x)
    for(t in 1:TT)
    {
      L <- L + vec(Rstar %*% Eytxt[,,t]) - kronecker(Extxt[,,t+1], Rstar) %*% fz
      L1 <- L1 + kronecker(Extxt[,,t+1] , Rstar)
    }
    zeta <- crossprod(Dz, L+2*eta * Dz %*% zeta0)
    zeta <- pcholinv((crossprod(Dz, L1)%*% Dz +2*eta*crossprod(Dz))) %*% zeta
    diftemp <- norm(zeta-zetapre,type="F")
    if(diftemp < dif)
    {
      dif <- diftemp
    }
    Z <- matrix(zeta, nrow = dim.r , ncol = dim.x)
    # Z[(dim.s+1):(dim.y),] <- matrix(zeta, nrow = dim.r , ncol = dim.x)
    cat("Z is updated\n")
    #
    iter <- iter + 1
    cat(iter,"iteration", "max difference is ", dif,"\n")
  }
  x_0 <- xtT[,1]
  V_0 <- VtT[,,1]
  return(list(B=B,Z=Z,Q=Q,R=R,x00=x_0, V00 = V_0,
              G=G,TT=TT,augQ = G %*% Q %*% trans(G)))
}
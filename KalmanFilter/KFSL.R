KFSL <- function(y,dim.s,dim.r,dim.x,dim.y,TT,p, parameter=NULL)
{
  # (s_t,s_t-1,...,s_t-p+1)^T = (phi_1,...,phi_p
  #                                I,...
  #                                  I,...)(s_t-1,s_t-2,...,s_t-p)^T+ (I,0,...,0)^Tw_t
  # Let (s_t,s_t-1,...,s_t-p+1)^T= x_t^T
  # (phi_1,...,phi_p
  # I,...
  # I,...) be B, (I,0,...,0)^T be G
  # w_t~N(0,Q)
  # (s_t,r_t)^T = (I,0,...,0 \\ psi_0,...,psi_p-1) x_t^T + (0,I)^T v_t
  # Let (s_t,r_t)^T = y_t^T, (I,0,...,0 \\ psi_0,...,psi_p-1) be Z, (0,I)^T be H
  # v_t~N(0,R). Notice that 0 in H is 32-by_420 and I is 420-by-420
  # We fix the initial value, x_0 to be ks, which means that F=0
  
  ####
  #Given s and r. Both should be centered
  #both s and r should be in long form: rows are obs
  ####
  if(is.null(parameter))
  {
    stop("No Parameters")
  }
  ####
  B <- parameter$B
  Z <- parameter$Z
  Q <- parameter$Q
  R <- parameter$R
  augQ <- parameter$augQ
  # augR <- parameter$augR
  x00 <- parameter$x00
  V00 <- parameter$V00
  ####
  
  #xtt_1 is the predicted value
  #x1_0,x2_1,...,xTT_TT-1
  #xtt is the filtered value
  #x00,x11,...,xTT_TT
  #xtT is the smoothed value
  #xTT_TT,xTT-1_TT,...,x0_TT
  #Note length of xtt and xtT are TT+1 to include x_0
  xtt_1 <- matrix(0,nrow = dim.x,ncol = TT)
  xtt <- xtT <- matrix(0,nrow = dim.x,ncol = TT+1)
  # Vtt is analogous to S&S Ptt, var[xt,xt|y(1:t)];
  # Vtt_1 is analogous to S&S Ptt1, cov[xt,xt-1|y(1:t)]
  # VtT and Vtt1T are for the backwards smoother;
  # VtT is var[xt,xt|y(1:T)]  and Vtt_1T is cov[xt,xt-1|y(1:T)]
  # J :JTT-1,...,J0
  #Note the length of Vtt and VtT are TT+1 to include V00,V0T

  Vtt_1 <- Vtt_1T <- J <- array(0,dim=c(dim.x,dim.x,TT))
  Vtt <- VtT <- array(0,dim=c(dim.x,dim.x,TT+1))

  #vt is the observation innovation
  #Ft is its covariance
  #Kt is the Kalman gain
  vt <- matrix(0,dim.y,TT)
  Ft <- array(0,dim=c(dim.y,dim.y,TT))
  Kt <- array(0, dim=c(dim.x,dim.y,TT))

  #FORWARD PASS (Kalman Filter)
  xtt[,1] <- x00
  Vtt[,,1] <- V00
  for(t in 1:TT)
  {
    xtt_1[,t] <- B %*% xtt[,t,drop = FALSE]
    Vtt_1[,,t] <- B %*% Vtt[,,t] %*% trans(B) + augQ
    Vtt_1[,,t] = symm(Vtt_1[,,t])

    vt[,t] <- y[,t,drop=FALSE] - Z%*%xtt_1[,t,drop=FALSE]
    Ft[,,t] <- Z%*%Vtt_1[,,t]%*%trans(Z) + R
    siginv <- pcholinv(Ft[,,t])
    Kt[,,t] <- Vtt_1[,,t]%*%trans(Z)%*%siginv
    xtt[,t+1] <- xtt_1[,t,drop=FALSE] + Kt[,,t]%*%vt[,t,drop=FALSE]
    Vtt[,,t+1] <- Vtt_1[,,t]-Kt[,,t]%*%Z%*%Vtt_1[,,t]
  }

  #BACKWARD PASS (Kalman Smoother)
  xtT[,TT+1] <- xtt[,TT+1,drop=FALSE]
  VtT[,,TT+1] <- Vtt[,,TT+1]
  for(t in TT:1)
  {
    J[,,t] <- Vtt[,,t] %*% trans(B) %*% pcholinv(Vtt_1[,,t])
    xtT[,t] <- xtt[,t,drop = FALSE] + J[,,t] %*% (xtT[,t+1] - xtt_1[,t])
    VtT[,,t] <- Vtt[,,t] + J[,,t] %*%
      (VtT[,,t+1] - Vtt_1[,,t]) %*% trans(J[,,t])
  }

  #Lag-One Covariance Smoother
  Vtt_1T[,,TT] <- B %*% Vtt[,,TT] - Kt[,,TT] %*% Z %*% B %*% Vtt[,,TT]
  for (t in (TT-1):1)
  {
    Vtt_1T[,,t] <- Vtt[,,t+1] %*% trans(J[,,t]) +
      J[,,t+1]%*%(Vtt_1T[,,t+1]-B%*%Vtt[,,t+1])%*%trans(J[,,t])
  }
  return(list(xtt_1 = xtt_1, xtt = xtt, xtT = xtT, 
              Vtt_1 = Vtt_1, Vtt_1T = Vtt_1T, J = J,
              Vtt = Vtt, VtT = VtT, vt = vt,
              Ft= Ft, Kt = Kt, parameter = parameter))
}
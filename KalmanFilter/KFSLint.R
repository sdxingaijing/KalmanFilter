KFSLint <- function(y,dim.s,dim.r,dim.x,dim.y,TT,p,lambda)
{
  #Use VAR and Ridge Regression to initialize parameters
  #######
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
  #s stands for audio signal
  #r stands for neuron data
  #Note y is in wide form
  #need to transform
  s <- trans(y[1:dim.s,])
  r <- trans(y[(dim.s+1):dim.y,])
  s <- as.matrix(s)
  r <- as.matrix(r)
  G <- diag(nrow = (dim.s)*p,ncol = dim.s)
  #H <- rbind(matrix(0,nrow=dim.s,ncol=dim.r),diag(nrow=dim.r))
  x_0 <- matrix(0,dim.s*p,1)#since we already centered the data
  #use MTS to set the initial value of B and Q
  library(MTS)
  temp <- VAR(s,p=p,include.mean = FALSE,output = FALSE)
  B <- rbind(temp$Phi,diag(nrow=(dim.s)*(p-1),ncol=(dim.s)*p))
  Q <- temp$Sigma
  #Q <- G %*% Q %*% trans(G)
  rm(temp)
  #use Ridge to set the initial values of Z and R
  temps <- matrix(0,nrow=TT,ncol=p*dim.s)
  temps[,1:dim.s] <- s
  if(p != 1)
  {
    for(i in 2:p)
    {
      temps[i:TT,((i-1)*(dim.s)+1):(i*dim.s)] <- s[1:(TT-i+1),]
    }
  }
  V_0 <- crossprod(temps)
  lam <- lambda
  tempbeta <- solve(t(temps)%*%temps/(dim(r)[1]) + lam*diag(1,nrow=(dim(temps)[2])))%*%t(temps)%*%r/(dim(r)[1])
  temppre <- temps %*% tempbeta
  tempres <- r - temppre
  # temp <- lm(r~temps-1)
  # R <- cov(temp$residuals)
  R <- cov(tempres)
  #R <- H %*% R %*% trans(H)
  #Z <- rbind(diag(nrow=dim.s,ncol=p*dim.s),t(tempbeta))
  Z <- t(tempbeta)
  rm(temps)
  rm(temppre)
  rm(tempbeta)
  rm(tempres)
  #End of initialization
  # augR = H %*% R %*% trans(H)
  # H=H,
  return(list(B=B,Z=Z,Q=Q,R=R,x00=x_0,V00 = V_0,
              G=G,TT=TT,augQ = G %*% Q %*% trans(G)))
}
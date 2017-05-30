Eytxt <- function(y,TT,dim.x,dim.y,FSL)
{
  #Return expectation of ytytT
  # list(xtt_1 = xtt_1, xtt = xtt, xtT = xtT, 
  #      Vtt_1 = Vtt_1, Vtt_1T = Vtt_1T, J = J,
  #      Vtt = Vtt, VtT = VtT, vt = vt,
  #      Ft= Ft, Kt = Kt, parameter = parameter)
  xtT <- FSL$xtT
  Vtt_1 <- FSL$Vtt_1
  Vtt_1T <- FSL$Vtt_1T
  VtT <- FSL$VtT
  #
  Eytxt <- array(0,dim = c(dim.y,dim.x,TT))
  Eytyt <- array(0,dim = c(dim.y,dim.y,TT))
  Extxt_1 <- array(0,dim = c(dim.x,dim.x,TT))
  Extxt <- array(0,dim = c(dim.x,dim.x,TT+1))
  Extxt[,,1] <- tcrossprod(xtT[,1]) + VtT[,,1]
  #
  sumytyt <- matrix(0,nrow = dim.y,ncol = dim.y)
  sumytxt <- matrix(0, nrow = dim.y, ncol = dim.x)
  sumxt_1xt_1 <- sumxtxt <- matrix(0, nrow = dim.x, ncol = dim.x)
  sumxtxt_1 <- matrix(0, nrow = dim.x, ncol = dim.x)
  for(t in 1:TT)
  {
    Eytyt[,,t] <- tcrossprod(y[,t])
    Eytxt[,,t] <- tcrossprod(y[,t],xtT[,t+1])
    Extxt[,,t+1] <- tcrossprod(xtT[,t+1]) + VtT[,,t+1]
    Extxt_1[,,t] <- tcrossprod(xtT[,t+1],xtT[,t]) + Vtt_1T[,,t]
    sumytyt <- sumytyt + Eytyt[,,t]
    sumytxt <- sumytxt + Eytxt[,,t]
    sumxtxt <- sumxtxt + Extxt[,,t+1]
    sumxt_1xt_1 <- sumxt_1xt_1 + Extxt[,,t]
    sumxtxt_1 <- sumxtxt_1 + Extxt_1[,,t]
  }
  return(list(Eytyt = Eytyt,Eytxt = Eytxt,Extxt = Extxt, Extxt_1 = Extxt_1,
              sumytyt = sumytyt, sumytxt = sumytxt, sumxtxt = sumxtxt,
              sumxt_1xt_1 = sumxt_1xt_1, sumxtxt_1 = sumxtxt_1))
}
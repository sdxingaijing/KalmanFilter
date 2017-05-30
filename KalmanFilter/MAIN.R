rm(list=ls())
setwd('D:\\RICE UNIVERSITY\\640\\COMPETITION\\DATA')
source("USEFUN.R")
source("KFSL.R")
source("KFSLint.R")
source("KMoEM.R")
source("Eytxt.R")
#load("KF2.RData")
#Read and split the data set
nsent <- 100
Ecog <- datasplit(sent = nsent)
#Variable selection
priordrop <- 140
bi <- FALSE
lag <- 2
lambda <- seq(from= 0.5, to = 1, by= 0.1)
score <- selevar(Ecog$Y,Ecog$X,Ecog$training,
                 priordrop = priordrop,lag = lag,
                 scalex = TRUE,nfolds =5,lambda =lambda, bi = bi)
plot(score,xlab = names(score))
# score1 <- selevar(Ecog$Y,Ecog$X,Ecog$training,
#                  priordrop = 140,lag = 2,scalex = TRUE,lambda =1, bi = TRUE)

#
Y_ecog <- Ecog$Y
X_ecog <- Ecog$X
X_ecog <- X_ecog[,-c(1:priordrop)]
if(bi == TRUE)
{
  #80% of the time
  X_ecog <- X_ecog[,which(score>=nsent * 0.8)]
}else{
  #one third of the time: on average get selected by at least one lag
  X_ecog <- X_ecog[,which(score>= nsent*ncol(Y_ecog))]
}
PC <- TRUE
if(PC)
{
  #Since we find that X_ecog data is still highly correlated 
  #Determinent of its covariance matrix is 0
  #We use PCA to de-correlate the data
  npc <- PCAsele(X_ecog, Ecog, nsent,scale = FALSE, percent = 0.9)$npc
}else{
  npc <- ncol(X_ecog)
}

#
p <- 1
#Set the form of parameters
#parlist <- paramset(p = p,dim.s = ncol(Y_ecog),dim.r = ncol(X_ecog))
parlist <- paramset(p = p,dim.s = ncol(Y_ecog),dim.r = npc)
parlist$Z <- parlist$Z[(ncol(Y_ecog)+1):(ncol(Y_ecog)+npc),]
#Get their fixed parts and free parts
#They should be time-invariant so we can set them beforehand
parfD <- list()
for(i in 1:length(parlist))
{
  parfD[[i]] <- convert.model.mat(parlist[[i]]) 
  print(i)
}
names(parfD) <- c("B","Z","Q","R","augQ","augR")
#
#Saved the data

#Note Y_ecog X_ecog are not centered or scaled

#Begin the loop
#Read the data sentence by sentence
eta <- 0.3 
alpha <- 0.05 
maxit <- 1000
tol <- 1e-5
#Note Center and scale the data first
for(sent in 1:nsent)
{
  s <- as.matrix(Y_ecog[(Ecog$training[1,sent]):(Ecog$training[2,sent]),])
  r <- as.matrix(X_ecog[(Ecog$training[1,sent]):(Ecog$training[2,sent]),])
  s <- scale(s,center = TRUE,scale = FALSE)
  r <- scale(r,center = TRUE,scale = FALSE)
  if(PC)
  {
    pcatemp <- svd(r)$v
    r <- r %*% pcatemp[,1:npc]
  }
  p <- p
  dim.s <- ncol(s)#num of audio frequencies
  dim.r <- ncol(r)#num of neuron variables
  
  TT <- nrow(s)
  
  #Stack two data sets to get observation process (wide form): columns are obs
  dim.x <- p*dim.s
  #Kalman Filter Smoother and Lag-one Covariance Matrix
  if(sent == 1)#Should initialize parameters
  {
    y <- rbind(t(s),t(r))
    dim.y <- dim.s + dim.r
    parameter <- KFSLint(y = y,dim.s = dim.s,dim.r = dim.r,dim.x = dim.x,
                         dim.y = dim.y,TT,p, lambda =1)
    y <- t(r)
    dim.y <- dim.r
    FSL <- KFSL(y = y,dim.s = dim.s,dim.r = dim.r,dim.x = dim.x,
                dim.y = dim.y,TT =TT,p=p, parameter=parameter)
  }else{
    y <- t(r)
    dim.y <- dim.r
    FSL <- KFSL(y = y,dim.s = dim.s,dim.r = dim.r,dim.x = dim.x,
                dim.y = dim.y,TT =TT,p=p, parameter= EM)
  }
  Expect <- Eytxt(y = y,TT = TT,dim.x = dim.x, dim.y = dim.y ,FSL = FSL)
  #EM
  if(sent == 1)
  {
    EM <- KMoEM(dim.s = dim.s,dim.r = dim.r,dim.x = dim.x ,dim.y = dim.y,TT = TT,p = p,
                initpar = FSL$parameter,FSL = FSL,Expect = Expect, 
                parfD = parfD, eta = 0 ,alpha = alpha ,maxit = maxit, tol = tol)
    
  }else{
    EM <- KMoEM(dim.s = dim.s,dim.r = dim.r,dim.x = dim.x ,dim.y = dim.y,TT = TT,p = p,
                initpar = FSL$parameter,FSL = FSL,Expect = Expect, 
                parfD = parfD, eta = eta ,alpha = alpha ,maxit = maxit, tol = tol)
  }
  cat("Sentence ",sent, "is done!\n")
}
parameter <- EM
#End of Kalman Filter

#Get the smoothed valued of training set
for(sent in 1:nsent)
{
  if(sent == 1)
  {
    trueY <- s <- as.matrix(Y_ecog[(Ecog$training[1,sent]):(Ecog$training[2,sent]),])
  }else{
    s <- as.matrix(Y_ecog[(Ecog$training[1,sent]):(Ecog$training[2,sent]),])
    trueY <- rbind(trueY,s)
  }
  r <- as.matrix(X_ecog[(Ecog$training[1,sent]):(Ecog$training[2,sent]),])
  s <- scale(s,center = TRUE,scale = FALSE)
  r <- scale(r,center = TRUE,scale = FALSE)
  if(PC)
  {
    pcatemp <- svd(r)$v
    r <- r %*% pcatemp[,1:npc]
  }
  p <- p
  dim.s <- ncol(s)#num of audio frequencies
  dim.r <- ncol(r)#num of neuron variables
  
  TT <- nrow(s)
  
  #Stack two data sets to get observation process (wide form): columns are obs
  dim.x <- p*dim.s
  
  y <- t(r)
  dim.y <- dim.r
  FSL <- KFSL(y = y,dim.s = dim.s,dim.r = dim.r,dim.x = dim.x,
              dim.y = dim.y,TT =TT,p=p, parameter= parameter)
  if(sent == 1)
  {
    estY <- trans(FSL$xtT[,-1])
  }else{
    estY <- rbind(estY, trans(FSL$xtT[,-1]))
  }
  print(sent)
}

#Get the smoothed valued of validation set
for(sent in 1:(length(Ecog$breakpoint) - nsent))
{
  if(sent == 1)
  {
    trueYv <- s <- as.matrix(Y_ecog[(Ecog$validation[1,sent]):(Ecog$validation[2,sent]),])
  }else{
    s <- as.matrix(Y_ecog[(Ecog$validation[1,sent]):(Ecog$validation[2,sent]),])
    trueYv <- rbind(trueYv,s)
  }
  r <- as.matrix(X_ecog[(Ecog$validation[1,sent]):(Ecog$validation[2,sent]),])
  s <- scale(s,center = TRUE,scale = FALSE)
  r <- scale(r,center = TRUE,scale = FALSE)
  if(PC)
  {
    pcatemp <- svd(r)$v
    r <- r %*% pcatemp[,1:npc]
  }
  p <- p
  dim.s <- ncol(s)#num of audio frequencies
  dim.r <- ncol(r)#num of neuron variables
  TT <- nrow(s)
  p <- p
  #Stack two data sets to get observation process (wide form): columns are obs
  dim.x <- p*dim.s
  y <- t(r)
  dim.y <- dim.r
  FSL <- KFSL(y = y,dim.s = dim.s,dim.r = dim.r,dim.x = dim.x,
              dim.y = dim.y,TT =TT,p=p, parameter= parameter)
  if(sent == 1)
  {
    estYv <- trans(FSL$xtT[,-1])
  }else{
    estYv <- rbind(estYv, trans(FSL$xtT[,-1]))
  }
  print(sent)
}

###














# library(glmnet)
# trial <- glmnet(estY , trueY, family = "mgaussian",alpha= 0, lambda =1,intercept = FALSE)
# BETA1 <- trial$beta
# temp <- as.matrix(cbind(BETA1$V1,BETA1$V2,BETA1$V3,BETA1$V4,BETA1$V5,BETA1$V6,BETA1$V7,
#                         BETA1$V8,BETA1$V9,BETA1$V10,BETA1$V11,BETA1$V12,BETA1$V13,BETA1$V14,
#                         BETA1$V15,BETA1$V16,BETA1$V17,BETA1$V18,BETA1$V19,BETA1$V20,
#                         BETA1$V21,BETA1$V22,BETA1$V23,BETA1$V24,BETA1$V25,BETA1$V26,BETA1$V27,
#                         BETA1$V28,BETA1$V29,BETA1$V30,BETA1$V31,BETA1$V32))
# 
# pre <- estYv %*% temp
# 
# print(sqrt(sum((pre-scale(trueYv,center=TRUE,scale = FALSE))^2)/(prod(dim(trueYv)))))


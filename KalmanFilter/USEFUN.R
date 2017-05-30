#Some Useful functions
#Data Reading and Data Splitting
datasplit <- function(sent = 100)
{
  setwd('D:\\RICE UNIVERSITY\\640\\COMPETITION\\DATA')
  Y_ecog <- as.matrix(read.csv('train_Y_ecog.csv',header=FALSE))
  X_ecog <- as.matrix(read.csv('train_X_ecog.csv',header=FALSE))
  #Split the data set into training and validation
  Breakend <- as.numeric(as.matrix(read.table('train_breakpoint.txt',header=FALSE)))
  Breakstr <- c(0,Breakend[-length(Breakend)])+1
  Inde <- sample(1:length(Breakend),size=length(Breakend))
  #Select sent sentences
  TimeEndtr <- Breakend[Inde[1:sent]]
  TimeStrtr <- Breakstr[Inde[1:sent]]
  TimeEndv <- Breakend[Inde[(sent+1):length(Inde)]]
  TimeStrv <- Breakstr[Inde[(sent+1):length(Inde)]]
  return(list(Y = Y_ecog, X = X_ecog,breakpoint = Breakend ,training = rbind(TimeStrtr,TimeEndtr), validation = rbind(TimeStrv,TimeEndv)))
}

#Select Variables
selevar <- function(Y, X, breakpoint, priordrop = 140,
                    lag = 2,nfolds =5 ,lambda = seq(from= 0.5, to = 1, by= 0.1),
                    scalex = TRUE,bi = FALSE)
{
  #if bi =TRUE only care about whether a var appeared or not
  Y <- as.matrix(Y)
  X <- as.matrix(X[,-c(1:priordrop)])
  breakpoint <- as.matrix(breakpoint)
  #Perform variable selection using ``lagged" lasso
  library(glmnet)
  score <- numeric(ncol(X))
  for(sent in 1:(dim(breakpoint)[2]))
  {
    tempX <- X[(breakpoint[1,sent]):(breakpoint[2,sent]),]
    tempY <- Y[(breakpoint[1,sent]):(breakpoint[2,sent]),]
    tempX <- scale(tempX, center = TRUE, scale = scalex)
    tempY <- scale(tempY, center = TRUE, scale= FALSE)
    if(lag!=0)
    {
      tempX1 <- tempX
      nc <- ncol(tempX); nr <- nrow(tempX)
      for(i in 1:lag)
      {
        temp <- rbind(tempX[-c(1:i),],matrix(0, nrow = i, ncol = nc))    
        tempX1 <- cbind(tempX1, temp)
      }
      tempX <- tempX1[-c(nr:(nr-lag+1)),]
      tempY <- tempY[-c(nr:(nr-lag+1)),]
      rm(tempX1); rm(temp)
    }
    res <- cv.glmnet(tempX,tempY,family="mgaussian",
                     alpha=1,nfolds = nfolds,
                     lambda = lambda,,standardize = FALSE,intercept = FALSE)
    BETA1 <- coef(res,s="lambda.1se")
    temp <- as.matrix(cbind(BETA1$V1,BETA1$V2,BETA1$V3,BETA1$V4,BETA1$V5,BETA1$V6,BETA1$V7,
                            BETA1$V8,BETA1$V9,BETA1$V10,BETA1$V11,BETA1$V12,BETA1$V13,BETA1$V14,
                            BETA1$V15,BETA1$V16,BETA1$V17,BETA1$V18,BETA1$V19,BETA1$V20,
                            BETA1$V21,BETA1$V22,BETA1$V23,BETA1$V24,BETA1$V25,BETA1$V26,BETA1$V27,
                            BETA1$V28,BETA1$V29,BETA1$V30,BETA1$V31,BETA1$V32))
    temp <- temp[-1,]
    varnames <- rownames(temp)[1:ncol(X)]
    temp <- ifelse(temp!=0,1,0)
    temp <- apply(temp,1,sum)
    temp <- matrix(temp,ncol = lag +1)
    if(bi)
    {
      temp <- apply(temp,1,sum)
      temp <- ifelse(temp!=0,1,0)
      score <- score + temp
    }else{
      score <- score + apply(temp,1,sum)
    }
    print(sent)
  }
  names(score) <- varnames
  return(score)
}

#Faster Transpose
trans <- function(x)return(matrix(x, dim(x)[2], dim(x)[1], byrow = TRUE))
#
#Make symmetric matrix
symm <- function (x) 
{
  t.x = trans(x)
  x = (x + t.x)/2
  x
}
#
#Faster Inverse
pcholinv <- function (x) 
{
  dim.x = dim(x)[1]
  diag.x = x[1 + 0:(dim.x - 1) * (dim.x + 1)]
  if (any(diag.x == 0)) {
    if (any(diag.x != 0)) {
      b = chol2inv(chol(x[diag.x != 0, diag.x != 0]))
      OMG.x = diag(1, dim.x)[diag.x != 0, , drop = FALSE]
      inv.x = t(OMG.x) %*% b %*% OMG.x
    }
    else {
      inv.x = matrix(0, dim.x, dim.x)
    }
  }
  else {
    inv.x = chol2inv(chol(x))
  }
  return(inv.x)
}
#
#Vec Operator
vec <- function (x) 
{
  if (!is.array(x)) 
    stop("vec:arg must be a 2D or 3D matrix")
  len.dim.x = length(dim(x))
  if (!(len.dim.x == 2 | len.dim.x == 3)) 
    stop("vec: arg must be a 2D or 3D matrix")
  if (len.dim.x == 2) {
    attr(x, "dim") = c(length(x), 1)
    return(x)
  }
  return(array(x, dim = c(length(x[, , 1]), 1, dim(x)[3])))
}
#Find the fixed part and free part of a matrix
convert.model.mat <- function (param.matrix) 
{
  if (!is.array(param.matrix)) 
    stop("convert.model.mat: function requires a 2D or 3D matrix")
  if (!(length(dim(param.matrix)) %in% c(2, 3))) 
    stop("convert.model.mat: arg must be a 2D or 3D matrix")
  Tmax = 1
  if (length(dim(param.matrix)) == 3) 
    Tmax = dim(param.matrix)[3]
  dim.f1 = dim(param.matrix)[1] * dim(param.matrix)[2]
  fixed = array(0, dim = c(dim.f1, 1, Tmax))
  c = param.matrix
  varnames = c()
  d = array(sapply(c, is.character), dim = dim(c))
  f = array(list(0), dim = dim(c))
  f[!d] = c[!d]
  f = vec(f)
  f = array(unlist(f), dim = c(dim.f1, 1, Tmax))
  is.char = c()
  if (any(d)) {
    is.char = which(d)
    if (any(grepl("[*]", c) | grepl("[+]", c))) {
      for (i in is.char) {
        e = mystrsplit(c[[i]])
        firstel = suppressWarnings(as.numeric(e[1]))
        if (length(e) == 1) {
          e = c("0", "+", "1", "*", e)
        }
        else {
          if (is.na(firstel) || !is.na(firstel) & e[2] == 
              "*") 
            e = c("0", "+", e)
        }
        pluses = which(e == "+")
        afterplus = suppressWarnings(as.numeric(e[pluses + 
                                                    1]))
        if (any(is.na(afterplus))) {
          k = 1
          for (j in pluses[is.na(afterplus)]) {
            e = append(e, c("1", "*"), after = j + (k - 
                                                      1) * 2)
            k = k + 1
          }
        }
        stars = which(e == "*")
        pluses = which(e == "+")
        if (length(stars) != length(pluses)) 
          stop("convert.model.mat: must use eqn form a+b1*p1+b2*p2...; extra p's can be left off")
        c[[i]] = paste(e, collapse = "")
        f[[i]] = as.numeric(e[1])
        varnames = c(varnames, e[stars + 1])
      }
      varnames = unique(varnames)
    }
    else {
      varnames = unique(unlist(c[is.char]))
    }
  }
  free = array(0, dim = c(dim.f1, length(varnames), Tmax))
  if (any(d)) {
    if (any(grepl("[*]", c) | grepl("[+]", c))) {
      for (i in is.char) {
        e = mystrsplit(c[[i]])
        stars = which(e == "*")
        e.vars = e[stars + 1]
        for (p in varnames) {
          drow = i%%dim.f1
          if (drow == 0) 
            drow = dim.f1
          if (p %in% e.vars) 
            free[drow, which(p == varnames), ceiling(i/dim.f1)] = sum(as.numeric(e[(stars - 
                                                                                      1)[e[stars + 1] == p]]))
        }
      }
    }
    else {
      for (p in varnames) {
        i = which(sapply(c, function(x) {
          identical(x, p)
        }))
        drow = i%%dim.f1
        drow[drow == 0] = dim.f1
        free[drow, which(p == varnames), ceiling(i/dim.f1)] = 1
      }
    }
  }
  fixed = f
  colnames(free) = varnames
  return(list(fixed = fixed, free = free))
}
#Set parameters
paramset <- function(p = 2,dim.s = 32,dim.r = 420)
{
  B <- list()
  count <- 1
  for(i in 1:p)
  {
    for(j in 1:dim.s)
    {
      for(k in 1:dim.s)
      {
        B[[count]] <- paste("Phi",i,"_",k,"_",j,sep='')
        count <- count + 1
      }
    }
  }
  B <- matrix(B,nrow=dim.s,ncol=p*dim.s)
  b <- diag(nrow=(p-1)*dim.s,ncol=p*dim.s)
  B <- rbind(B,b)
  rm(b)
  #
  Z <- list()
  count <- 1
  for(i in 1:p)
  {
    for(j in 1:dim.s)
    {
      for(k in 1:dim.r)
      {
        Z[[count]] <- paste("Psi",i,"_",k,"_",j,sep='')
        count <- count + 1
      }
    }
  }
  Z <- matrix(Z,nrow=dim.r,ncol=p*dim.s)
  z <- diag(nrow = dim.s,ncol = p*dim.s)
  Z <- rbind(z,Z)
  rm(z)
  #
  Q <- matrix(list(),nrow = dim.s, ncol = dim.s)
  for(i in 1:dim.s)
  {
    for(j in i:dim.s)
    {
      Q[i,j] <- paste("Q",i,"_",j,sep='')
      if(i != j)
      {
        Q[j,i] <- paste("Q",i,"_",j,sep='')
      }
    }
  }
  q1 <- matrix(0,nrow = dim.s,ncol = (p-1)*dim.s)
  q2 <- matrix(0, nrow = (p-1)*dim.s, ncol = p*dim.s )
  augQ <- cbind(Q,q1)
  augQ <- rbind(augQ,q2)
  rm(q1)
  rm(q2)
  #
  R <- matrix(list(),nrow = dim.r, ncol = dim.r)
  for(i in 1:dim.r)
  {
    for(j in i:dim.r)
    {
      R[i,j] <- paste("R",i,"_",j,sep='')
      if(i != j)
      {
        R[j,i] <- paste("R",i,"_",j,sep='')
      }
    }
  }
  r1 <- matrix(0,nrow = dim.r,ncol = dim.s)
  r2 <- matrix(0, nrow = dim.s, ncol = dim.r + dim.s)
  augR <- cbind(r1,R)
  augR <- rbind(r2,augR)
  rm(r1)
  rm(r2)
  #
  
  
  return(list(B=B,Z=Z,Q=Q,R=R,augQ=augQ,augR=augR))
}
#PCA select
PCAsele <- function(X_ecog,Ecog,nsent, scale = FALSE, percent = 0.8)
{
  n <- numeric(nsent)
  for(sent in 1:nsent)
  {
    r <- as.matrix(X_ecog[(Ecog$training[1,sent]):(Ecog$training[2,sent]),])
    r <- scale(r,center = TRUE,scale = FALSE)
    a <- (svd(r)$d)^2
    n[sent] <- which((cumsum(a)/sum(a))>0.8)[1]
  }
  m <- which.max(table(n))
  return(list(npc = m, table = n))
}




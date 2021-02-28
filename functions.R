genrDatwb_2g <- function(N,lambda1,lambda2,k1,k2,expLambda,true_beta,dx){
  x1_5 <- sapply(rep(0,5),rnorm,n=N,sd=1)
  xc <- sapply(rep(0,5),rnorm,n=N,sd=1)
  x <- cbind(rep(1,N),x1_5,xc)
  
  linpred <- x%*%true_beta
  prob <- exp(linpred) / (1+exp(linpred))
  u <- runif(N)
  B <- ifelse(u<prob,1,0)+1
  # generate time to failure for each subgroup
  n1 <- sum(B==1)
  n2 <- sum(B==2)
  t <- rep(0,N)
  u <- runif(n1)
  t[B==1] <- lambda1*(-log(u))^(1/k1)
  u <- runif(n2)
  t[B==2] <- lambda2*(-log(u))^(1/k2)
  # generate censoring 
  u <- runif(N)
  c <- rexp(n=N,rate=expLambda)
  y <- pmin(t,c)
  delta <- (t<=c)*1
  dat <- data.frame(x=x[,-1],t=t,c=c,y=y,delta=delta,B=B)
  names(dat) <- c(paste("x",1:10,sep = '_'),'t','c','y','delta','B')
  return(dat)
}

#### calculate split criterion
## @param: var_name: name of variable that needs to find potential split
## dat: dataset at current node
## q.curr: a matrix of n*K, each element q_ik is the initial weight/posterior prob for subject i in group k
## beta.curr: a 2*1 matrix of current beta\
splitCal_ini <- function(var_name, dat, beta.curr,k1.curr,k2.curr,lambda1.curr,lambda2.curr){
  x <- dat[,var_name]
  delta <- dat$delta
  y <- dat$y
  n <- nrow(dat)
  x_sp <- quantile(x,seq(0.2,0.8,0.1))
  splits <- sort(unique(x_sp))
  nsplits <- length(splits)
  K <- length(beta.curr)/2 + 1 
  
  S1 <- exp(-(y/lambda1.curr)^k1.curr)
  S2 <- exp(-(y/lambda2.curr)^k2.curr)
  f1 <- k1.curr*lambda1.curr^(-k1.curr)*y^(k1.curr-1)*exp(-(y/lambda1.curr)^k1.curr)
  f2 <- k2.curr*lambda2.curr^(-k2.curr)*y^(k2.curr-1)*exp(-(y/lambda2.curr)^k2.curr)
  
  objFun <- matrix(0,nrow=nsplits,ncol=2)
  colnames(objFun) <- c("split","objFun")
  beta.new <- matrix(0,nrow=2,ncol=nsplits)
  Q1.new <- matrix(0,nrow=n,ncol=nsplits)
  Q2.new <- matrix(0,nrow=n,ncol=nsplits)
  for(i in 1:nsplits){
    sp <- splits[i]
    Z <- cbind(rep(1,n),ifelse(x<sp,1,0)) 
    
    E1 <- matrix(1,nrow=n,ncol=1)
    E2 <- exp(Z%*%beta.curr)
    Denom <- E1 + E2
    
    Q1 <- (delta*f1+(1-delta)*S1)*E1
    Q2 <- (delta*f2+(1-delta)*S2)*E2
    SQ <- Q1+Q2
    Q1 <- Q1/SQ
    Q2 <- Q2/SQ
    
    ### use weighted logistic regression
    b <- c(rep(0,n),rep(1,n))
    q <- c(Q1,Q2)
    dat.log <- data.frame(z1=Z[,1],z2=Z[,2],b=b,q=q)
    fit.log <- glm(b~z2,weight=q,data=dat.log,family=quasibinomial)
    beta.new[,i] <- as.numeric(coef(fit.log))
    
    E1 <- matrix(1,nrow=n,ncol=1)
    E2 <- exp(Z%*%beta.new[,i])
    Denom <- E1 + E2
    Q1.new[,i] <- Q1
    Q2.new[,i] <- Q2
    
    obj <- sum(Q1.new[,i] * (E1/Denom) + Q2.new[,i] * (E2/Denom))
    objFun[i,] <- c(sp,obj)
  }
  objFun <- na.omit(objFun)
  ind.max <- which(objFun[,"objFun"]==max(objFun[,"objFun"]))
  split_at <- objFun[ind.max,"split"]
  max_objFun <- objFun[ind.max,"objFun"]
  beta.opt <- beta.new[,ind.max]
  Q1.opt <- Q1.new[,ind.max]
  Q2.opt <- Q2.new[,ind.max]
  return(list(split_at = split_at, max_objFun = max_objFun, beta.new = beta.opt, 
              split_var = var_name, objFun = objFun, Q1.new = Q1.opt, Q2.new = Q2.opt))
}

###################################
### generate subset of data based on the feature and split that gives the largest value of objective function
#param: dat: data set from last split (on current node)
#       info: a data frame, columns are var_name.v: a vector of variable names
#                                       split_at.v: a vector of optimal split for each variable
#                                       max_objFun.v: a vector of maximum of objective function corresponding to the variable and split
find_split <- function(info){
  ind.max <- which(info$max_objFun.v==max(info$max_objFun.v))
  split_at <- info$split_at.v[ind.max]
  max_objFun <- info$max_objFun.v[ind.max]
  var_name <- info$var_name.v[ind.max]
  return(list(var_name=var_name,split_at=split_at,max_objFun=max_objFun))
}

split_dat <- function(dat,find_split,Q.curr,idx){
  Q.curr <- Q.curr[, c(1,(ncol(Q.curr)/2)+1)]
  colnames(Q.curr) <- c("Q1","Q2")
  dat$Q1 <- Q.curr$Q1
  dat$Q2 <- Q.curr$Q2
  splitVar <- as.character(find_split$var_name[idx])
  splitValue <- find_split$split_at[idx]
  dat_aftSp <- list()
  dat_aftSp[["le_sp"]] <- dat[(dat[,splitVar]<splitValue),]
  dat_aftSp[["ge_sp"]] <- dat[(dat[,splitVar]>=splitValue),]
  return(dat_aftSp)
}

#### combine data from each nodes for survival update
## treeInfo contains all current in treeInfo[[S.s]][[dat]]
## returns the whole dataset with current weights
comb_dat <- function(treeInfo){
  dat_all <- c()
  s_curr <- length(treeInfo)
  nNodes <- length(treeInfo[[paste("S",s_curr,sep=".")]][["dat"]])
  for(node in 1:nNodes){
    dat_all <- rbind(dat_all,treeInfo[[paste("S",s_curr,sep=".")]][["dat"]][[node]])
  }
  s_curr <- s_curr - 1
  while(s_curr >= 1){
    nNodes <- length(treeInfo[[paste("S",s_curr,sep=".")]][["dat"]])
    dat_alltemp <- c()
    for(node in 1:nNodes){
      dat_alltemp <- rbind(dat_alltemp,treeInfo[[paste("S",s_curr,sep=".")]][["dat"]][[node]])
    }
    row_num_diff <- setdiff(row.names(dat_alltemp),row.names(dat_all))
    dat_all <- rbind(dat_all,dat_alltemp[row.names(dat_alltemp)%in%row_num_diff,])
    s_curr <- s_curr - 1
  }
  return(dat_all)
}

#### estimate the survival distribution
EstSurv <- function(dat_all,k1.curr,lambda1.curr,k2.curr,lambda2.curr,tau=0.0001){
  n <- nrow(dat_all)
  Y <- matrix(dat_all$y,nrow=n,ncol=1)
  delta <- matrix(dat_all$delta,nrow=n,ncol=1)
  Q1 <- dat_all$Q1
  Q2 <- dat_all$Q2
  
  diffv <- c()
  diff <- 1
  iter <- 0
  theta.curr <- matrix(c(k1.curr,lambda1.curr,k2.curr,lambda2.curr),4,1)
  while(diff>tau & diff < 10){
    S1 <- exp(-(Y/lambda1.curr)^k1.curr)
    S2 <- exp(-(Y/lambda2.curr)^k2.curr)
    f1 <- k1.curr*lambda1.curr^(-k1.curr)*Y^(k1.curr-1)*exp(-(Y/lambda1.curr)^k1.curr)
    f2 <- k2.curr*lambda2.curr^(-k2.curr)*Y^(k2.curr-1)*exp(-(Y/lambda2.curr)^k2.curr)
    
    parDerv1_k1 <- sum(Q1*(delta*(1/k1.curr + log(Y/lambda1.curr)) - (Y/lambda1.curr)^k1.curr * log(Y/lambda1.curr)))
    parDerv1_lambda1 <- sum(Q1*(delta*(-k1.curr/lambda1.curr) + k1.curr*Y^k1.curr*lambda1.curr^(-k1.curr-1)))
    parDerv1_k2 <- sum(Q2*(delta*(1/k2.curr + log(Y/lambda2.curr)) - (Y/lambda2.curr)^k2.curr * log(Y/lambda2.curr)))
    parDerv1_lambda2 <- sum(Q2*(delta*(-k2.curr/lambda2.curr) + k2.curr*Y^k2.curr*lambda2.curr^(-k2.curr-1)))
    parDerv1_surv <- matrix(c(parDerv1_k1,parDerv1_lambda1,parDerv1_k2,parDerv1_lambda2),4,1)
    
    parDerv2_k1k1 <- sum(Q1*(delta*(-1/k1.curr^2) - (Y/lambda1.curr)^k1.curr*(log(Y/lambda1.curr))^2))
    parDerv2_k1lambda1 <- sum(Q1*(delta*(-1/lambda1.curr) + Y^k1.curr*lambda1.curr^(-k1.curr-1)*(k1.curr*log(Y/lambda1.curr)+1)))
    parDerv2_lambda1lambda1 <- sum(Q1*(delta*(k1.curr/lambda1.curr^2) - k1.curr*(k1.curr+1)*Y^k1.curr*lambda1.curr^(-k1.curr-2)))
    parDerv2_k2k2 <- sum(Q2*(delta*(-1/k2.curr^2) - (Y/lambda2.curr)^k2.curr*(log(Y/lambda2.curr))^2))
    parDerv2_k2lambda2 <- sum(Q2*(delta*(-1/lambda2.curr)+Y^k2.curr*lambda2.curr^(-k2.curr-1)*(k2.curr*log(Y/lambda2.curr)+1)))
    parDerv2_lambda2lambda2 <- sum(Q2*(delta*(k2.curr/lambda2.curr^2) - k2.curr*(k2.curr+1)*Y^k2.curr*lambda2.curr^(-k2.curr-2)))
    parDerv2_surv <- matrix(c(parDerv2_k1k1,parDerv2_k1lambda1,0,0,
                              parDerv2_k1lambda1,parDerv2_lambda1lambda1,0,0,
                              0,0,parDerv2_k2k2,parDerv2_k2lambda2,
                              0,0,parDerv2_k2lambda2,parDerv2_lambda2lambda2),4,4)
    
    size <- 0
    theta.new <- theta.curr - (2^size*solve(parDerv2_surv))%*%parDerv1_surv
    while((sum(theta.new<0)>0)|(sum(abs(theta.new)>6.5)>0)){
      size <- size-1
      theta.new <- theta.curr - (2^size*solve(parDerv2_surv))%*%parDerv1_surv
    }
    diff <- sqrt(t(theta.new-theta.curr)%*%(theta.new-theta.curr))
    diffv <- c(diffv,diff)
    theta.curr <- theta.new
    k1.curr <- theta.new[1,];k2.curr <- theta.new[3,]
    lambda1.curr <- theta.new[2,];lambda2.curr <- theta.new[4,]
    iter <- iter+1
  }
  if (diff>10) {converge=0} 
  if (diff<tau) {converge=1}
  return(list(k1=k1.curr,k2=k2.curr,lambda1=lambda1.curr,
              lambda2=lambda2.curr, iter=iter,diffv=diffv, converge=converge))
}

#### function to update weight q_ik
weightUpdate <- function(tree_Info,dat,survEst.curr){
  k1.curr <- survEst.curr$k1;k2.curr <- survEst.curr$k2
  lambda1.curr <- survEst.curr$lambda1;lambda2.curr <- survEst.curr$lambda2
  n <- nrow(dat)
  y <- dat$y
  delta <- dat$delta
  splitVar <- as.character(tree_Info$split_var)
  splitValue <- tree_Info$split_at
  beta.curr <- matrix(c(tree_Info$beta0,tree_Info$beta1),nrow=2)
  
  S1 <- exp(-(y/lambda1.curr)^k1.curr)
  S2 <- exp(-(y/lambda2.curr)^k2.curr)
  f1 <- k1.curr*lambda1.curr^(-k1.curr)*y^(k1.curr-1)*exp(-(y/lambda1.curr)^k1.curr)
  f2 <- k2.curr*lambda2.curr^(-k2.curr)*y^(k2.curr-1)*exp(-(y/lambda2.curr)^k2.curr)
  
  Z <- cbind(rep(1,n),ifelse(dat[,splitVar]<splitValue,1,0))
  
  E1 <- matrix(1,nrow=n,ncol=1)
  E2 <- exp(Z%*%beta.curr)
  Denom <- E1 + E2
  
  Q1 <- (delta*f1+(1-delta)*S1)*E1
  Q2 <- (delta*f2+(1-delta)*S2)*E2
  SQ <- Q1+Q2
  Q1 <- Q1/SQ
  Q2 <- Q2/SQ
  
  return(list(Q1=Q1, Q2=Q2))
}

splitCal <- function(var_name, dat, beta.curr,k1.curr,k2.curr,lambda1.curr,lambda2.curr,Q.curr){
  x <- dat[,var_name]
  delta <- dat$delta
  y <- dat$y
  n <- nrow(dat)
  x_sp <- quantile(x,seq(0.2,0.8,0.1))
  splits <- sort(unique(x_sp))
  nsplits <- length(splits)
  K <- length(beta.curr)/2 + 1
  
  S1 <- exp(-(y/lambda1.curr)^k1.curr)
  S2 <- exp(-(y/lambda2.curr)^k2.curr)
  f1 <- k1.curr*lambda1.curr^(-k1.curr)*y^(k1.curr-1)*exp(-(y/lambda1.curr)^k1.curr)
  f2 <- k2.curr*lambda2.curr^(-k2.curr)*y^(k2.curr-1)*exp(-(y/lambda2.curr)^k2.curr)
  
  objFun <- matrix(0,nrow=nsplits,ncol=2)
  colnames(objFun) <- c("split","objFun")
  beta.new <- matrix(0,nrow=2,ncol=nsplits)
  Q1.new <- matrix(0,nrow=n,ncol=nsplits)
  Q2.new <- matrix(0,nrow=n,ncol=nsplits)
  for(i in 1:nsplits){
    sp <- splits[i]
    Z <- cbind(rep(1,n),ifelse(x<sp,1,0))
    
    E1 <- matrix(1,nrow=n,ncol=1)
    E2 <- exp(Z%*%beta.curr)
    Denom <- E1 + E2
    Q1 <- dat$Q1
    Q2 <- dat$Q2
    
    b <- c(rep(0,n),rep(1,n))
    q <- c(Q1,Q2)
    dat.log <- data.frame(z1=rep(Z[,1],K),z2=rep(Z[,2],K),b=b,q=q)
    fit.log <- glm(b~z2,weights=q,data=dat.log,family=quasibinomial)
    beta.new[,i] <- as.numeric(coef(fit.log))
    
    E1 <- matrix(1,nrow=n,ncol=1)
    E2 <- exp(Z%*%beta.new[,i])
    Denom <- E1 + E2
    Q1.new[,i] <- Q1
    Q2.new[,i] <- Q2
    
    obj <- round(sum(Q1.new[,i] * (E1/Denom) + Q2.new[,i] * (E2/Denom)),32)
    objFun[i,] <- c(sp,obj)
  }
  objFun <- na.omit(objFun)
  ind.max <- which(objFun[,"objFun"]==max(objFun[,"objFun"]))
  split_at <- objFun[ind.max,"split"]
  max_objFun <- objFun[ind.max,"objFun"]
  beta.opt <- beta.new[,ind.max]
  Q1.opt <- Q1.new[,ind.max]
  Q2.opt <- Q2.new[,ind.max]
  return(list(split_at = split_at, max_objFun = max_objFun, beta.new = beta.opt, 
              split_var = var_name, objFun = objFun, 
              Q1.new = Q1.opt, Q2.new = Q2.opt))
}

### Function to grow a decision tree
## @param: X: n*p data matrix of covariates
## y: observed survival time
## delta: censoring status, 1= event
## beta.ini: starting value of beta coefficients
## k1, k2, lambda1, lambda2: known parameters of survival function for two subgroups
growTree <- function(X, y, delta, beta.ini,ini_k1,ini_k2,ini_lambda1,ini_lambda2){
  treeInfo <- list()
  splitInfo <- list()
  varNames <- colnames(X)
  nVars <- length(varNames)
  n <- nrow(X)
  nstop <- ceiling(sqrt(n))
  dat <- as.data.frame(cbind(X,y,delta))
  
  s=1
  node=1
  # initial spliting
  split_at.v <- c() 
  max_objFun.v <- c()
  splitInfo[[paste("S",s,sep=".")]][[node]] <- list()
  for(nvar in 1:nVars){
    var_name <- varNames[nvar]
    splitInfo[[paste("S",s,sep=".")]][[node]][[var_name]] <- splitCal_ini(var_name=var_name, dat=dat, beta.curr=beta.ini,
                                                                          k1.curr=ini_k1,k2.curr=ini_k2,lambda1.curr=ini_lambda1,
                                                                          lambda2.curr=ini_lambda2)
    if(length(splitInfo[[paste("S",s,sep=".")]][[node]][[var_name]]$split_at)>1){
      ind <- sample(1:length(splitInfo[[paste("S",s,sep=".")]][[node]][[var_name]]$split_at),1)
    }else{ind <- 1}
    split_at.v <- c(split_at.v,splitInfo[[paste("S",s,sep=".")]][[node]][[var_name]]$split_at[ind])
    max_objFun.v <- c(max_objFun.v,splitInfo[[paste("S",s,sep=".")]][[node]][[var_name]]$max_objFun[ind])
  }
  info <- data.frame(var_name.v = varNames, split_at.v = split_at.v, max_objFun.v = max_objFun.v)
  S_res <- find_split(info=info)
  if(length(S_res$var_name)>1){
    n_optVar <- length(S_res$var_name)
    idx <- sample(1:n_optVar,1)
  } else { idx <- 1 }
  Q.curr <- data.frame(Q1=splitInfo[[paste("S",s,sep=".")]][[node]][[S_res$var_name[idx]]]$Q1.new,Q2=splitInfo[[paste("S",s,sep=".")]][[node]][[S_res$var_name[idx]]]$Q2.new)# q_ik not updated yet
  next_dat <- split_dat(dat=dat,find_split=S_res,Q.curr=Q.curr,idx=idx)
  
  treeInfo[[paste("S",s,sep=".")]][["Info"]] <- data.frame(split_var=rep(S_res$var_name[idx],2), split_at=rep(S_res$split_at[idx],2),
                                                           max_objFun.v=rep(S_res$max_objFun[idx],2),beta0=rep(splitInfo[[paste("S",s,sep=".")]][[node]][[S_res$var_name[idx]]]$beta.new[1],2), 
                                                           beta1=rep(splitInfo[[paste("S",s,sep=".")]][[node]][[S_res$var_name[idx]]]$beta.new[2],2),
                                                           sign=c("le","ge"),n=c(nrow(next_dat$le_sp),nrow(next_dat$ge_sp)),
                                                           stop=c((nrow(next_dat$le_sp)<=nstop),(nrow(next_dat$ge_sp)<=nstop)))
  treeInfo[[paste("S",s,sep=".")]][["dat"]][[1]] <-  next_dat$le_sp
  treeInfo[[paste("S",s,sep=".")]][["dat"]][[2]] <- next_dat$ge_sp
  
  dat_all <- comb_dat(treeInfo=treeInfo)
  survEst.curr <- EstSurv(dat_all=dat_all,k1.curr=ini_k1,lambda1.curr=ini_lambda1,k2.curr=ini_k2,lambda2.curr=ini_lambda2,tau=0.0001)
  treeInfo[[paste("S",s,sep=".")]][["survEst"]] <- survEst.curr
  
  while(sum(treeInfo[[paste("S",s,sep=".")]][["Info"]]$stop) < 2^s){
    s=s+1
    nNodes <- nrow(treeInfo[[paste("S",s-1,sep=".")]][["Info"]])
    treeInfo[[paste("S",s,sep=".")]][["Info"]] <- c()
    splitInfo[[paste("S",s,sep=".")]] <- list()
    for(node in 1:nNodes){
      lg <- treeInfo[[paste("S",s-1,sep=".")]][["Info"]][node,"sign"]
      if(!treeInfo[[paste("S",s-1,sep=".")]][["Info"]][node,"stop"]){
        split_at.v <- c()
        max_objFun.v <- c()
        nextDat <- treeInfo[[paste("S",s-1,sep=".")]][["dat"]][[node]]
        newWeight <- weightUpdate(tree_Info=treeInfo[[paste("S",s-1,sep=".")]][["Info"]][node,],dat=nextDat,survEst.curr=survEst.curr)
        nextDat$Q1 <- newWeight$Q1
        nextDat$Q2 <- newWeight$Q2
        splitInfo[[paste("S",s,sep=".")]][[node]] <- list()
        for(nvar in 1:nVars){
          var_name <- varNames[nvar]
          splitInfo[[paste("S",s,sep=".")]][[node]][[var_name]] <- splitCal(var_name=var_name, dat=nextDat, beta.curr=matrix(c(treeInfo[[paste("S",s-1,sep=".")]][["Info"]]$beta0[node],0),2,1),
                                                                            k1.curr=survEst.curr$k1,k2.curr=survEst.curr$k2,lambda1.curr=survEst.curr$lambda1,lambda2.curr=survEst.curr$lambda2)
          if(length(splitInfo[[paste("S",s,sep=".")]][[node]][[var_name]]$split_at)>1){
            ind <- sample(1:length(splitInfo[[paste("S",s,sep=".")]][[node]][[var_name]]$split_at),1)
          }else{ind <- 1}
          split_at.v <- c(split_at.v,splitInfo[[paste("S",s,sep=".")]][[node]][[var_name]]$split_at[ind])
          max_objFun.v <- c(max_objFun.v,splitInfo[[paste("S",s,sep=".")]][[node]][[var_name]]$max_objFun[ind])
        }
        info <- data.frame(var_name.v = varNames, split_at.v = split_at.v, max_objFun.v = max_objFun.v)
        S_res <- find_split(info=info)
        if(length(S_res$var_name)>1){
          n_optVar <- length(S_res$var_name)
          idx <- sample(1:n_optVar,1)
        } else { idx <- 1 }
        
        Q.curr <- data.frame(Q1=splitInfo[[paste("S",s,sep=".")]][[node]][[S_res$var_name[idx]]]$Q1.new,Q2=splitInfo[[paste("S",s,sep=".")]][[node]][[S_res$var_name[idx]]]$Q2.new)
        next_dat.new <- split_dat(dat=nextDat,find_split=S_res,Q.curr=Q.curr,idx=idx)
        
        tree_info <- data.frame(split_var=rep(S_res$var_name[idx],2), split_at=rep(S_res$split_at[idx],2),
                                max_objFun.v=rep(S_res$max_objFun[idx],2),beta0=rep(splitInfo[[paste("S",s,sep=".")]][[node]][[S_res$var_name[idx]]]$beta.new[1],2), 
                                beta1=rep(splitInfo[[paste("S",s,sep=".")]][[node]][[S_res$var_name[idx]]]$beta.new[2],2),
                                sign=c("le","ge"),n=c(nrow(next_dat.new$le_sp),nrow(next_dat.new$ge_sp)),
                                stop=c((nrow(next_dat.new$le_sp)<nstop),(nrow(next_dat.new$ge_sp)<nstop)),node_prev=rep(node,2))
        treeInfo[[paste("S",s,sep=".")]][["Info"]] <- rbind(treeInfo[[paste("S",s,sep=".")]][["Info"]],tree_info)
        treeInfo[[paste("S",s,sep=".")]][["dat"]][[2*node-1]] <-  next_dat.new$le_sp
        treeInfo[[paste("S",s,sep=".")]][["dat"]][[2*node]] <- next_dat.new$ge_sp
        
      } else { 
        tree_info <- data.frame(split_var=rep(NA,2), split_at=rep(NA,2),
                                max_objFun.v=rep(NA,2),beta0=rep(NA,2), 
                                beta1=rep(NA,2),
                                sign=c("le","ge"),n=c(NA,NA),
                                stop=c(TRUE,TRUE),node_prev=rep(node,2))
        treeInfo[[paste("S",s,sep=".")]][["Info"]] <- rbind(treeInfo[[paste("S",s,sep=".")]][["Info"]],tree_info)
        treeInfo[[paste("S",s,sep=".")]][["dat"]][[2*node-1]] <-  NULL
        treeInfo[[paste("S",s,sep=".")]][["dat"]][[2*node]] <- NULL
      }
    }
    dat_all <- comb_dat(treeInfo=treeInfo)
    survEst.curr <- EstSurv(dat_all=dat_all,k1.curr=survEst.curr$k1,lambda1.curr=survEst.curr$lambda1,k2.curr=survEst.curr$k2,lambda2.curr=survEst.curr$lambda2,tau=0.0001)
    treeInfo[[paste("S",s,sep=".")]][["survEst"]] <- survEst.curr  
  }
  
  return(list("treeInfo"=treeInfo,"splitInfo"=splitInfo))
}

BICTree <- function(treeInfo, lambda1, lambda2, k1, k2, nbranch,N){
  obs.loglik <- c()
  n.obs <- c()
  
  for(branch in nbranch:1){
    if(branch==nbranch){
      ind <- which(!is.na(treeInfo[[branch]][["Info"]]$n))
    }
    if(branch<nbranch){
      ind <- which((!is.na(treeInfo[[branch]][["Info"]]$n) & (treeInfo[[branch]][["Info"]]$stop==TRUE)))
    }
    node.prev <- treeInfo[[branch]][["Info"]]$node_prev[ind]
    n.subdat <- length(ind)
    for(subdat.ind in ind){
      subdat <- treeInfo[[branch]][["dat"]][[subdat.ind]]
      split.var <- treeInfo[[branch]][["Info"]]$split_var[subdat.ind]
      split.at <- treeInfo[[branch]][["Info"]]$split_at[subdat.ind]
      beta0 <- treeInfo[[branch]][["Info"]]$beta0[subdat.ind]
      beta1 <- treeInfo[[branch]][["Info"]]$beta1[subdat.ind]
      if(treeInfo[[branch]][["Info"]]$sign[subdat.ind]=="le"){
        pi1 <- 1 / (1+exp(beta0+beta1))
        pi2 <- exp(beta0+beta1) / (1+exp(beta0+beta1))
      }
      if(treeInfo[[branch]][["Info"]]$sign[subdat.ind]=="ge"){
        pi1 <- 1 / (1+exp(beta0))
        pi2 <- exp(beta0) / (1+exp(beta0))
      }
      y <- subdat$y; delta <- subdat$delta
      S1 <- exp(-(y/lambda1)^k1)
      S2 <- exp(-(y/lambda2)^k2)
      f1 <- k1*lambda1^(-k1)*y^(k1-1)*exp(-(y/lambda1)^k1)
      f2 <- k2*lambda2^(-k2)*y^(k2-1)*exp(-(y/lambda2)^k2)
      obs.loglik.curr <- sum(delta * log(f1*pi1 + f2*pi2) + (1-delta) * log(S1*pi1 + S2*pi2))
      obs.loglik <- c(obs.loglik,obs.loglik.curr)
      n.obs <- c(n.obs,nrow(subdat))
    }
  }
  if(sum(n.obs)!=N){print("Sample Size is incorrect!")}
  n.leaves <- length(obs.loglik)
  BIC <- -2*sum(obs.loglik) + log(N)*n.leaves
  
  return(list(BIC=BIC,n.leaves=n.leaves))
}

PruneTree <- function(treeInfo,N){
  nbranch <- length(treeInfo)
  BIC <- BICTree(treeInfo=treeInfo,lambda1=treeInfo[[nbranch]][["survEst"]]$lambda1,
                 lambda2=treeInfo[[nbranch]][["survEst"]]$lambda2,
                 k1=treeInfo[[nbranch]][["survEst"]]$k1,
                 k2=treeInfo[[nbranch]][["survEst"]]$k2,nbranch=nbranch,N=N)$BIC
  BIC.new <- BICTree(treeInfo=treeInfo,lambda1=treeInfo[[(nbranch-1)]][["survEst"]]$lambda1,
                     lambda2=treeInfo[[(nbranch-1)]][["survEst"]]$lambda2,
                     k1=treeInfo[[(nbranch-1)]][["survEst"]]$k1,
                     k2=treeInfo[[(nbranch-1)]][["survEst"]]$k2,nbranch=nbranch-1,N=N)$BIC
  while(BIC>BIC.new){
    BIC <- BIC.new
    nbranch <- nbranch-1
    if(nbranch==1){
      BIC.new <- BICTree(treeInfo=treeInfo,lambda1=treeInfo[[nbranch]][["survEst"]]$lambda1,
                         lambda2=treeInfo[[nbranch]][["survEst"]]$lambda2,
                         k1=treeInfo[[nbranch]][["survEst"]]$k1,
                         k2=treeInfo[[nbranch]][["survEst"]]$k2,nbranch=nbranch,N=N)$BIC
      BIC <- BIC.new
    }
    if(nbranch>1){
      BIC.new <- BICTree(treeInfo=treeInfo,lambda1=treeInfo[[(nbranch-1)]][["survEst"]]$lambda1,
                         lambda2=treeInfo[[(nbranch-1)]][["survEst"]]$lambda2,
                         k1=treeInfo[[(nbranch-1)]][["survEst"]]$k1,
                         k2=treeInfo[[(nbranch-1)]][["survEst"]]$k2,nbranch=nbranch-1,N=N)$BIC
    }
  }
  return(list(branch.opt=nbranch,BIC=BIC))
}

##### function to test tree, calculate predicted class
TestTree <- function(treeInfo,N,testDat,prune){
  if(prune){
    branch.opt <- PruneTree(treeInfo=treeInfo,N=N)$branch.opt
    n.leaves <- BICTree(treeInfo=treeInfo,lambda1=treeInfo[[branch.opt]][["survEst"]]$lambda1,
                        lambda2=treeInfo[[branch.opt]][["survEst"]]$lambda2,
                        k1=treeInfo[[branch.opt]][["survEst"]]$k1,
                        k2=treeInfo[[branch.opt]][["survEst"]]$k2,nbranch=branch.opt,N=N)$n.leaves
  }
  if(!prune){
    branch.opt <- length(treeInfo)
    n.leaves <- BICTree(treeInfo=treeInfo,lambda1=treeInfo[[branch.opt]][["survEst"]]$lambda1,
                        lambda2=treeInfo[[branch.opt]][["survEst"]]$lambda2,
                        k1=treeInfo[[branch.opt]][["survEst"]]$k1,
                        k2=treeInfo[[branch.opt]][["survEst"]]$k2,nbranch=branch.opt,N=N)$n.leaves
  }
  
  # calculate predicted class
  for(branch in branch.opt:1){
    if(branch==branch.opt){
      ind <- which(!is.na(treeInfo[[branch]][["Info"]]$n))
    }
    if(branch<branch.opt){
      ind <- which((!is.na(treeInfo[[branch]][["Info"]]$n)) & (treeInfo[[branch]][["Info"]]$stop==TRUE))
    }
    treeInfo[[branch]][["Info"]]$B_pred <- NA
    if(length(ind)>0){
      for(subdat.ind in ind){
        subdat <- treeInfo[[branch]][["dat"]][[subdat.ind]]
        beta0 <- treeInfo[[branch]][["Info"]]$beta0[subdat.ind]
        beta1 <- treeInfo[[branch]][["Info"]]$beta1[subdat.ind]
        if(treeInfo[[branch]][["Info"]]$sign[subdat.ind]=="le"){
          subdat$pi1 <- 1 / (1+exp(beta0+ifelse(is.na(beta1),0,beta1)))
          subdat$pi2 <- exp(beta0+ifelse(is.na(beta1),0,beta1)) / (1+exp(beta0+ifelse(is.na(beta1),0,beta1)))
        }
        if(treeInfo[[branch]][["Info"]]$sign[subdat.ind]=="ge"){
          subdat$pi1 <- 1 / (1+exp(beta0))
          subdat$pi2 <- exp(beta0) / (1+exp(beta0))
        }
        if(sum(subdat$pi1==subdat$pi2)>0){
          subdat$B.pred[which(subdat$pi1==subdat$pi2)] <- sample(1:2,1)
        } else {
          subdat$B.pred <- ifelse(subdat$pi1>subdat$pi2, 1,2)
        }
        
        if(length(unique(subdat$B.pred))==1){
          Bhat <- unique(subdat$B.pred)
        }else{
          Bhat <- which(table(subdat$B.pred)==max(table(subdat$B.pred)))
        }
        if(length(Bhat)>1){
          Bhat <- sample(1:2,1)
        }
        treeInfo[[branch]][["Info"]]$B_pred[subdat.ind] <- Bhat
      }
    }
  }
  
  # predict on testDat
  tDat.tree <- list()
  tDat <- NULL
  for(branch.t in 1:branch.opt){
    n.node <- nrow(treeInfo[[branch.t]][["Info"]])
    ind <- which(!is.na(treeInfo[[branch.t]][["Info"]]$B_pred))
    tDat.tree[[branch.t]] <- list()
    for(node in 1:n.node){
      split.var <- as.character(treeInfo[[branch.t]][["Info"]]$split_var[node])
      split.at <- treeInfo[[branch.t]][["Info"]]$split_at[node]
      if(is.na(split.at)){
        tDat.tree[[branch.t]][[node]] <- NULL
      }else{
        if(treeInfo[[branch.t]][["Info"]]$sign[node]=="le"){
          if(branch.t==1){tDat.tree[[branch.t]][[node]] <- testDat[testDat[,split.var]<split.at,]}
          if(branch.t>1){tDat.tree[[branch.t]][[node]] <- tDat.tree[[branch.t-1]][[ceiling(node/2)]][tDat.tree[[branch.t-1]][[ceiling(node/2)]][,split.var]<split.at,]}
        }
        if(treeInfo[[branch.t]][["Info"]]$sign[node]=="ge"){
          if(branch.t==1){tDat.tree[[branch.t]][[node]] <- testDat[testDat[,split.var]>=split.at,]}
          if(branch.t>1){tDat.tree[[branch.t]][[node]] <- tDat.tree[[branch.t-1]][[ceiling(node/2)]][tDat.tree[[branch.t-1]][[ceiling(node/2)]][,split.var]>=split.at,]}
        }
      }
    }
    
    if(length(ind)>0){
      for(id in ind){
        Bpred <- treeInfo[[branch.t]][["Info"]]$B_pred[id]
        if(nrow(tDat.tree[[branch.t]][[id]])>0){
          tDat <- rbind(tDat,cbind(tDat.tree[[branch.t]][[id]],Bpred))
        }
      }
    }
  }
  accuracy <- sum(tDat$B==tDat$Bpred)/nrow(tDat)
  return(list(accuracy=accuracy,branch.opt=branch.opt, n.leaves=n.leaves,Bpred=tDat$Bpred,Btrue=tDat$B,tDat=tDat))
}

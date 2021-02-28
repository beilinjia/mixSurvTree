########################################################
##### This file is an example for simulation study #####
##### including data generation                    #####
########################################################
##### Author: Beilin Jia
source("./functions.R")

N=1000
lambda1=exp(0);lambda2=exp(1.5)
k1=1;k2=3
ini_lambda1=0.5;ini_lambda2=1.5;ini_k1=0.5;ini_k2=1.5
expLambda=0.15
true_beta=c(0.4,0.2,-0.6,-0.3,rep(0,7))
initial_beta=matrix(rep(0,11),11,1)
n.test=10000;dx=10

#### generate training data 
set.seed(22)
train <- genrDatwb_2g(N=N,dx=dx,lambda1=lambda1,lambda2=lambda2,k1=k1,k2=k2,
                      expLambda=expLambda,true_beta=true_beta)

#### grow a tree
tree <- growTree(X=train[,1:10],y=train$y,delta=train$delta,beta.ini=matrix(0,2,1),
                      ini_lambda1=ini_lambda1,ini_lambda2=ini_lambda2,ini_k1=ini_k1,ini_k2=ini_k2)

#### testing
test <- genrDatwb_2g(N=10000,dx=dx,lambda1=lambda1,lambda2=lambda2,k1=k1,k2=k2,
                     expLambda=expLambda,true_beta=true_beta)

testtree <- TestTree(treeInfo=tree$treeInfo,N=N,testDat=test,prune=TRUE)
testtree$accuracy
table(testtree$tDat$B,testtree$tDat$Bpred)

testtree_notprune <- TestTree(treeInfo=tree$treeInfo,N=N,testDat=test,prune=FALSE)
testtree_notprune$accuracy
table(testtree_notprune$tDat$B,testtree_notprune$tDat$Bpred)




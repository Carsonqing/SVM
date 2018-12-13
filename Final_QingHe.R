data=read.table("E:/Academic/6106 Stat Comp/hw/Project3/pb2.txt")
colnames(data)=c("y","x1","x2","x3","x4")
Data = as.matrix(data, ncol=5)
head(Data)
x.train =Data[Data[,"y"]==1,2:5];head(x.train)
Y = Data[Data[,"y"]==1,1];head(Y)
N<-nrow(x.train);N

rbf_kernel <- function(x1,x2,gamma){
  K<-exp(-(1/gamma^2)*t(x1-x2)%*%(x1-x2))
  return(K)
}
gamma=3
C=3
############define the function to solve alpha
svddtrain <- function(X,C=Inf, gamma=1.5,esp=1e-10){
  X<-as.matrix(X)
  N=nrow(X)
  Dm<-matrix(0,N,N)
  
  require('quadprog')
  for(i in 1:N){
    for(j in 1:N){
      Dm[i,j]<- rbf_kernel(X[i,],X[j,],gamma)
    }
  }
  Dm<-2*Dm+diag(N)*1e-12 # adding a very small number to the diag, some trick
  
  dv<-c(rep(1,N))###dv is a constant vector with value 1
  Am<-t(rbind(rep(1,N),diag(N),-1*diag(N)))
  bv<-c(1,rep(0,N)) # the 1 is for the sum(alpha)==0, others for each alpha_i >= 0
  if(C!=Inf){
    bv<-c(bv,rep(-C,N))
  }
  
  alpha_org<-solve.QP(Dm,dv,Am,meq=1,bvec=bv)$solution
  alphaindx<-which(alpha_org>esp,arr.ind=TRUE)
  alpha<-alpha_org[alphaindx]
  nSV<-length(alphaindx)
  
  #below we compute the radius
  Xv<-X[alphaindx,]###subset the support vector
  K<-matrix(0,nSV,nSV)###compute the K matrix again because this time we only use the sv
  for(i in 1:nSV){
    for(j in 1:nSV){
      K[i,j]<-rbf_kernel(Xv[i,],Xv[j,],gamma)
    }
  }
  
  temp=t(alpha)%*%K%*%alpha
  d=rep(1+temp,nSV)
  for (i in 1:nSV) {
    for (j in 1:nSV){
      d[i] = d[i] -2*alpha[j]*rbf_kernel(Xv[i,],Xv[j,],gamma)
    }
  }
  R.2=mean(d)
  list(alpha=alpha, R.2=R.2, nSV=nSV, Xv=Xv, gamma=gamma)
  }
############end of definition

svddtrain(x.train,C=3,gamma=3)



#######################prediction#############

distance = function(x,Xv,nSV,alpha,gamma){
  K<-matrix(0,nSV,nSV)
  for(i in 1:nSV){
    for(j in 1:nSV){
      K[i,j]<-rbf_kernel(Xv[i,],Xv[j,],gamma)
    }
  }
  
  temp=t(alpha)%*%K%*%alpha;temp
  d=1+temp
  nSV=nrow(Xv)
  for (j in 1:nSV){
    d = d -2*alpha[j]*rbf_kernel(x,Xv[j,],gamma)
  }
  return(d)
}

pred = function(z,model){
  R2=model$R.2###calculate the radius
  nSV=model$nSV
  Xv=model$Xv
  alpha=model$alpha
  gamma=model$gamma
  
  dz=distance(z,Xv,nSV,alpha,gamma)
  if (dz<R2) {y.hat=1
  }else{
    y.hat=-1
  }
  list(distance=dz,class=y.hat)
}


gamma=3;C=3
model=svddtrain(x.train,C,gamma)


###training error
pred(x.train[1,],model)
apply(x.train,1,function(x) pred(x,model)$class)
sum(apply(x.train,1,function(x) pred(x,model)$class)==-1)/nrow(x.train)

###testing error
x.test =Data[Data[,"y"]==2,2:5];head(x.test)
#apply(x.test,1,function(x) distance(x,Xv,alpha,gamma))
apply(x.test,1,function(x) pred(x,model)$class)
sum(apply(x.test,1,function(x) pred(x,model)$class)==1)/nrow(x.test)

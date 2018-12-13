data=read.table("E:/Academic/6106 Stat Comp/hw/Project3/pb2.txt")
colnames(data)=c("y","x1","x2","x3","x4")
Data = as.matrix(data, ncol=5)
x.train = scale(Data[Data[,"y"]==1,2:5])
y.train = as.matrix(Data[Data[,"y"]==1,1],ncol=1)
var(x.train)
n = nrow(x.train);n

rbf_kernel <- function(x1,x2,gamma){
  K<-exp(-(1/gamma^2)*t(x1-x2)%*%(x1-x2))
  return(K)
}

###############function to solve a and b#################
solve.ab= function(X,gamma,c){
  n=dim(X)[1]
  ### compute K matrix ###
  K<-matrix(0,n,n)
  
  for(i in 1:n){
    for(j in 1:n){
      K[i,j]<-rbf_kernel(X[i,],X[j,],gamma)
    }
  }
  temp1=cbind(K+diag(n)/c,rep(1,n))
  temp2=cbind(c(rep(1,n),0))
  A=rbind(temp1,t(temp2))
  temp3=t(cbind(t(rep(0,n)),1))
  value=solve(A,temp3)
  a=value[1:n]
  b=value[n+1]
  return(list(A,a,b))
}


#############prediction#############
classify=function(xt,x,a,b,gamma,c){
  yt=b
  n=nrow(x)
  for (j in 1:n){
    yt = yt+a[j]*rbf_kernel(xt,x[j,],gamma)
  }
  return(yt)
}

##############take gamma=c=2 for example############
gamma=2
c=2
results=solve.ab(x.train,gamma,c)
A=results[[1]]
a=results[[2]];a
b=results[[3]];b
A.inv=solve(A)

yt=rep(0,n)
for (i in 1:n) {
  yt[i]=classify(x.train[i,],x.train,a,b,gamma,c)
}

yt
yt.new=1-2*abs(yt)/(c*(sign(yt)*(max(yt)+min(yt))+(max(yt)-min(yt))));yt.new
sum(yt.new<0)

###write the function to find the best c############
gamma=2
c=100
step=2/3
error=0
while(error==0){
  results=solve.ab(x.train,gamma,c)
  A=results[[1]]
  a=results[[2]]
  b=results[[3]]
  A.inv=solve(A)

  y.loo=rep(0,31)
  y.loo.new=rep(0,31)
  for (i in 1:31) {
    y.loo[i]=-a[i]/A.inv[i,i]
    
    results.r=solve.ab(x.train[-i,],gamma,c)
    A.r=results.r[[1]]
    a.r=results.r[[2]]
    b.r=results.r[[3]]
    A.r.inv=solve(A)
    
    yt=apply(x.train[-i,], 1,function(x) classify(x,x.train[-i,],a.r,b.r,gamma,c))
    ymax=max(yt)
    ymin=min(yt)
    y.loo.new[i]=1-2*abs(y.loo[i])/(c*(sign(y.loo[i])*(ymax+ymin)+(ymax-ymin)))
  }
 
  error=sum(y.loo.new<0)/nrow(x.train)
  
  print(c(c,error))
  if (error==0) {
    c=step*c  
  }
}






#############################problem 2###############################
###yt is a p-dimension observation#######3
install.packages("Matrix")
library(Matrix)


y = Data[Data[,"y"]==1,2:5];head(y)
s=var(y);s
s.inv=solve(s)
(x=chol(solve(s)))

###################write a for loop to calculate lambda#############
n=nrow(y);n
lambda=rep(0,n)
for (i in 1:n) {
  yt=y[i,]
  if (length(unique(yt))==1) {
    yt[1]=yt[1]+0.001
    cat("warning: identital values at",i)
  }
  yt
  zt=x%*%yt
  dt=data.frame(cbind(yt,x))
  null=lm(yt~1,data=dt)
  full=lm(yt~.,data=dt)
  lm.fw=step(null, scope=list(lower=null, upper=full), direction="forward",steps = 1,trace = FALSE)
  var=names(lm.fw$coefficients[2])###the non-zero variable
  coef=lm.fw$coefficients[2]      ###the non-zero coefficient
  if (is.na(var))  {
    ut=c(0,0,0,0)
  }else if (var=="x1"){
    ut=c(coef,0,0,0)
    }else if (var=="x2") {
    ut=c(0,coef,0,0)
    }else if (var=="x3") {
    ut=c(0,0,coef,0)
  }else if (var=="x4") {
  ut=c(0,0,0,coef)
  }
  lambda[i]=2*t(yt)%*%s.inv%*%ut-t(ut)%*%s.inv%*%ut
}

lambda

###################end of for loop#####################


#####################problem 3###############


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
####useless here#######
for (i in 1:n){
  if (Y[i] > 1){
    Y[i]<--1
  }
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



###################problem 3.b################
data=as.matrix(data)
train=sample(1:62,43,replace = F)###0.7
data.train=data[train,]
x.train=data[train,2:5]
y.train=data[train,1]
for (i in 1:length(y.train)){
  if (y.train[i] > 1){
    y.train[i]<--1
  }
}
x1.train=data.train[data.train[,"y"]==1,2:5]
x2.train=data.train[data.train[,"y"]==2,2:5]
#x.test=data[-train,2:5];x.test
#y.test=data[-train,1]
#for (i in 1:length(y.test)){
#  if (y.test[i] > 1){
#    y.test[i]<--1
#  }
#}
###the 2 classes svdd model###

gamma=C=3
model.1=svddtrain(x1.train,C=C,gamma=gamma)
model.2=svddtrain(x2.train,C=C,gamma=gamma)

#center

a.1=as.vector(t(model.1$alpha)%*%model.1$Xv);a.1
a.2=as.vector(t(model.2$alpha)%*%model.2$Xv);a.2

R2.1=model.1$R.2
R2.2=model.2$R.2

nSV.1=model.1$nSV
nSV.2=model.2$nSV

Xv.1=model.1$Xv
Xv.2=model.2$Xv

alpha.1=model.1$alpha
alpha.2=model.2$alpha

#confidence measure
#x=x.train[2,];x

#confidence.1=exp(-t(x-a.1)%*%(x-a.1)/R2.1);confidence.1
#confidence.2=exp(-t(x-a.2)%*%(x-a.2)/R2.2);confidence.2
#theta=exp(-100)

y.hat=numeric(43)
####1st layer####
###write a function to loop all the training data points

x2=matrix(NA, nrow=43,ncol=2)
for (i in 1:43) {
  y.1=c(pred(x.train[i,],model.1)$class)
  y.2=c(pred(x.train[i,],model.2)$class)
  
  u.1=c(R2.1-pred(x.train[i,],model.1)$distance)
  u.2=c(R2.2-pred(x.train[i,],model.2)$distance)
  
  x2[i,]=c(u.1,u.2)
  ##make the decision
  if (u.1*u.2 <0 ) {###not suspicious
    y.hat[i]=y.1
  }
  
  
  }
y.hat
sum(y.train==y.hat)/17
##end of the 1st layer
##now pass the new data set to the 2nd layer
x2=x2[which(y.hat==0),]
y2=y.train[which(y.hat==0)]

x2.1=x2[which(y2==1),]
x2.2=x2[which(y2==2),]

###################2nd layer##############
gamma=2;C=4
model2.1=svddtrain(x2.1,C=C,gamma=gamma)
model2.2=svddtrain(x2.2,C=C,gamma=gamma)

#center

a2.1=as.vector(t(model2.1$alpha)%*%model2.1$Xv);a2.1
a2.2=as.vector(t(model2.2$alpha)%*%model2.2$Xv);a2.2

R2.1=model2.1$R.2
R2.2=model2.2$R.2

nSV2.1=model2.1$nSV
nSV2.2=model2.2$nSV

Xv2.1=model2.1$Xv
Xv2.2=model2.2$Xv

alpha2.1=model2.1$alpha
alpha2.2=model2.2$alpha

#confidence measure
#confidence.1=exp(-t(x-a.1)%*%(x-a.1)/R2.1);confidence.1
#confidence.2=exp(-t(x-a.2)%*%(x-a.2)/R2.2);confidence.2
#theta=exp(-100)

n=nrow(x2)
y2.hat=numeric(n)
####1st layer####
###write a function to loop all the training data points

x3=matrix(NA, nrow=n,ncol=2)
for (i in 1:n) {
  y.1=c(pred(x2[i,],model2.1)$class)
  y.2=c(pred(x2[i,],model2.2)$class)
  
  u.1=c(R2.1-pred(x2[i,],model2.1)$distance)
  u.2=c(R2.2-pred(x2[i,],model2.2)$distance)
  
  x3[i,]=c(u.1,u.2)
  ##make the decision
  if (u.1*u.2 <0 ) {###not suspicious
    y2.hat[i]=y.1
  }
  
  
}
y2.hat
y2.train=y.train[which(y.hat==0)]
sum(y2.train==y2.hat)/22
##end of the 2nd layer

##now pass the new data set to the 3rd layer
x3=x3[which(y2.hat==0),]
y3=y2[which(y2.hat==0)]

x3.1=x3[which(y3==1),];x3.1
x3.2=x3[which(y3==2),];x3.2

###end of 2nd layer


###################3rd layer##############
gamma=2;C=5
model3.1=svddtrain(x3.1,C=C,gamma=gamma)
model3.2=svddtrain(x3.2,C=C,gamma=gamma)

#center

a3.1=as.vector(t(model3.1$alpha)%*%model3.1$Xv);a3.1
a3.2=as.vector(t(model3.2$alpha)%*%model3.2$Xv);a3.2

R3.1=model3.1$R.2
R3.2=model3.2$R.2

nSV3.1=model3.1$nSV
nSV3.2=model3.2$nSV

Xv3.1=model3.1$Xv
Xv3.2=model3.2$Xv

alpha3.1=model3.1$alpha
alpha3.2=model3.2$alpha

#confidence measure
#confidence.1=exp(-t(x-a.1)%*%(x-a.1)/R2.1);confidence.1
#confidence.2=exp(-t(x-a.2)%*%(x-a.2)/R2.2);confidence.2
#theta=exp(-100)

n=nrow(x3)
y3.hat=numeric(n)
####1st layer####
###write a function to loop all the training data points

x4=matrix(NA, nrow=n,ncol=2)
for (i in 1:n) {
  y.1=c(pred(x3[i,],model3.1)$class)
  y.2=c(pred(x3[i,],model3.2)$class)
  
  u.1=c(R3.1-pred(x3[i,],model3.1)$distance)
  u.2=c(R3.2-pred(x3[i,],model3.2)$distance)
  
  x4[i,]=c(u.1,u.2)
  ##make the decision
  if (u.1*u.2 <0 ) {###not suspicious
    y3.hat[i]=y.1
  }
  
  
}
y3.hat
##end of the 3rd layer


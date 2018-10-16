import numpy as np
import csv
import random
from numpy.linalg import inv
import matplotlib.pyplot as plt

#define the kernel function
def f(x1,x2,sigma):
    u=x1-x2        
    return(np.exp(-np.dot(u,u)/(sigma**2)))

#define the inverse function and return b and alpha
def inverse(x,y,c,sigma):         #d1 and d2 are x and y 
    size=x.shape[0]        #size is the sample size
    omega=np.empty((size,size))     #define the omega matrix
    for i in range(len(x)):
        for j in range(len(x)):
            omega[i,j]=y[i]*y[j]*f(x[i,],x[j,],sigma)
    
    i=np.identity(size)     #identity matrix
    omega_i=omega+i/c
    
    m=np.append([0],y).reshape((size+1,1)) #m is the left  matrix
    yt=y.reshape((1,size))
    n=np.concatenate((yt,omega_i),axis=0)   #n is the right  matrix
    final=np.column_stack((m,n))            #stack m and n by column then we get the final matrix
    final_inv=inv(final)
    v=np.append([0],[1]*size).reshape((size+1),1)   #v is the constant matrix at the right of the equation
    answer=np.matmul(final_inv,v)
    #get b and alpha
    b=answer[0]
    alpha=answer[1:]
    return(b,alpha)
    
#predict the test set
def predict(x,x0,y0,b,a,sigma):
    """x is the data point you want to predict
    x0 is the traning data
    b and a is the parameters in formula 5
    sigma is the tuning variance
    """
    #n is the test sample size
    #n0 is the training sample size
    n = len(x) 
    n0 = len(x0)
    y_h = np.zeros(n)
    for i in range(n):
        y_h[i]=b
        for j in range(n0):
            y_h[i] += a[j] * y0[j]*f(x[i],x0[j],sigma)  
    return(y_h)    
    
#data preparation
with open("E:/Academic/6106 Stat Comp/hw/charlie.csv")as file:
    reader = csv.reader(file)
    data_raw = [row for row in reader]
 
data = np.array(data_raw)[1:,]

#create response vector
y=data[:,0]
def category(y):
    if y == "Original":
        return 1
    if y == "New":
        return -1 
y = np.array([category(i) for i in y])

#create predictors vector
x=data[:,2:6]
x=x.astype(float)

#random sampling 
random.seed(2)
idx=np.random.choice(30, 23, replace=False)
y_tr=y[idx]
x_tr=x[idx]

idx2=[i for i in range(30) if i not in idx]
y_te=y[idx2]
x_te=x[idx2]

#use the training data to solve b and alpha
b_tr,alpha_tr=inverse(x_tr,y_tr,3,2)#b=0.15
print("b_tr:",b_tr,"\nalpha_tr:",alpha_tr)

#test error rate
y_h=predict(x_te,x_tr,y_tr,b_tr,alpha_tr,2)
y_h[y_h>0]=1
y_h[y_h<0]=-1
misscla=np.sum(y_h!=y_te)/7

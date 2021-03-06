import numpy as np
import csv
import random
from numpy.linalg import inv
import matplotlib.pyplot as plt
def H(x, y, sigma,c):
    length = len(x)
    K = np.empty((length, length))
    I = np.identity(length)
    for i in range(length):
        for j in range(length):
            K[i, j] = f(x[i], x[j], sigma)
    H = K+I/c
    return H
    
def solve(x,c,sigma):
    length=len(x)
    e=np.ones(length).reshape((length,1))
    et=np.ones(length)
    H_matrix=H(x,y,sigma,c)
    H_inv=inv(H_matrix)
    temp=np.matmul(et,H_inv)
    temp2=np.matmul(temp,e)
    rho=-1/temp2#negative sign
    alpha=np.matmul(H_inv,e)/temp2
    return(rho,alpha)

def predict2(x,x0,b,a,sigma):
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
            y_h[i] += a[j] *f(x[i],x0[j],sigma)
    
    return(y_h)
    
#prepare the training and testing data    
z=data[data[:,0]=="Original",-2:].astype(np.float)
test=data[data[:,0]=="New",-2:].astype(np.float)

#evaluate the classification function performance
p1,q1=solve(z,3.5,3)#sigma=2
y_h=predict2(test,z,p1,q1,3)#input sigma

#get the testing outlier rate
outlier=np.sum(y_h<0)/len(test)
outlier     
print("p1:",p1,"\nq1",q1)
y_h0=predict2(z,z,p1,q1,3)#input sigma

#get the training outlier rate
outlier0 = np.sum(y_h<0)/len(z)
outlier0   

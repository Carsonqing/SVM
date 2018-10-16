#read the data into python
with open("E:/Academic/6106 Stat Comp/charlie.csv")as file:
    reader = csv.reader(file)
    data_raw = [row for row in reader]
 
data = np.array(data_raw)[1:,]
x=data[data[:,0]=="Original",-2:].astype(np.float)

#define the gaussian kernel
def k(x1,x2,sigma):
    u=x1-x2        
    return(np.exp(-np.dot(u,u)/(sigma**2)))
def solve(x,sigma,c):
    n=x.shape[0]        
    #calculate vector k
    k_v=np.zeros((n,1))
    for i in range(len(x)):
        k_v[i]=k(x[i,],x[i,],sigma)
    
    #calculate matrix K
    K=np.empty((n,n))    #matrix K   
    for i in range(len(x)):
        for j in range(len(x)):
            K[i,j]=k(x[i,],x[j,],sigma)
    
    #calculate matrix H
    H=K+np.identity(n)/(2*c)
    H_inv=inv(H)
    e=np.ones((n,1))
    et=np.ones(n)
    temp1=np.matmul(np.matmul(et,H_inv),k_v)    #e'*H_inv*k
    temp2=np.matmul(np.matmul(et,H_inv),e)      #e'*H_inv*e
    alpha=1/2*np.matmul(H_inv,k_v+(2-temp1)/temp2*e)
    return(alpha,K)

(alpha,K)=solve(x,2,3)
print("alpha:\n",alpha)

#use matrix form to calculate sum(a_j*a_l*k(x_j,x_l))
alpha_t=np.reshape(alpha,(1,len(alpha)))
part3=np.matmul(np.matmul(alpha_t,K),alpha)
print("part3:",part3)


#define the distance function
def distance(z,x,alpha,sigma):
    """
    z is testing data
    x is training data
    """
    m=z.shape[0]#test sample size
    n=x.shape[0]#training sample size
    d=np.zeros(m)
    for i in range(m):
        d[i]=k(z[i,],z[i,],sigma)
        for j in range(n):
            d[i] += -2*alpha[j]*k(z[i,],x[j,],sigma)+part3
    return d
#calculate the R^2
R2=np.mean(distance(x,x,alpha,2))
print("R2:",R2)

test=data[data[:,0]=="New",-2:].astype(np.float)

#define the classification function
def classify(x):
    if x  <= R2:
        return 1
    if x > R2:
        return 0

#calculate the distance of test data
d_test=distance(test,x,alpha,2)
print("test data distance:",d_test)
#calculate the distance of training data
d_x=distance(x,x,alpha,2)
print("training data distance:",d_x)
#test error
y_test = np.array([classify(i) for i in d_test])
print("test error:",np.sum(y_test!=0)/len(test))
#training error
y_x = np.array([classify(i) for i in d_x])
print("training error:",np.sum(y_x==0)/len(x))

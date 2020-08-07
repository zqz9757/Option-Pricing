######COS期权定价
#将分布函数f(x)转化为cos expansion（以N(0,1)为例）
#并比较正态pdf和cos expansion


def cos_ak(x):
    def normal_cf(v):
        return np.exp(-0.5*(v**2))
    def getAk_normal(k,a,b):
        ak=(normal_cf(k*np.pi/(b-a))*np.exp(-1j*k*np.pi*a/(b-a))*2)/(b-a)
        return ak.real
    a=-10
    b=10
    N=64
    temp=0
    for k in range(N):
        if k==0:
            temp+=getAk_normal(k,a,b)*np.cos(k*(x-a)*np.pi/(b-a))/2
        else:
            temp+=getAk_normal(k,a,b)*np.cos(k*(x-a)*np.pi/(b-a))
    return temp
for i in range(-5,5,1):
    print('x=',i,'\n')
    print('get by pdf:',norm.pdf(i))
    print('get by COS:',cos_ak(i))
#COS数值解
S0=100.0
K=80.0
T=0.1
r=0.1
sigma=0.25
def cf_log_x_div_k(v,S0=S0,T=T,r=r,sigma=sigma):            #log(ST/K)的特征函数
    cf_value = np.exp(((np.log(S0/K))+ (r - 0.5 * sigma ** 2)*T )* 1j * v - 0.5 * sigma ** 2 * v ** 2*T)
    return cf_value
def getAk(k,a,b,cf=cf_log_x_div_k):
    ak=2/(b-a)*cf(k*np.pi/(b-a))*np.exp(-1j*k*np.pi*a/(b-a))
    return ak.real
def chi(a,b,k,c,d):
    A=np.cos(k*np.pi*(d-a)/(b-a))*np.exp(d)
    B=-np.cos(k*np.pi*(c-a)/(b-a))*np.exp(c)
    C=k*np.pi/(b-a)*np.sin(k*np.pi*(d-a)/(b-a))*np.exp(d)
    D=-k*np.pi/(b-a)*np.sin(k*np.pi*(c-a)/(b-a))*np.exp(c)
    return (A+B+C+D)/(1+(k*np.pi/(b-a))**2)
def psi(a,b,k,c,d):
    if k==0:
        return d-c
    else:
        A=np.sin(k*np.pi*(d-a)/(b-a))
        B=np.sin(k*np.pi*(c-a)/(b-a))
        return (A-B)*(b-a)/(k*np.pi)
def vk_call(a,b,K,k):
    res=2/(b-a)*K*(chi(a,b,k,0,b)-psi(a,b,k,0,b))
    return res
def Cos_call(S0,K,T,r,sigma):
    N=256
    a=r*T-10*np.sqrt((sigma**2)*T)
    b=r*T+10*np.sqrt((sigma**2)*T)
    sum=0
    x=np.full(N,0.0)
    for k in range(N):
        if k==0:
            x[k]=getAk(k,a,b,cf=cf_log_x_div_k)*vk_call(a,b,K,k)/2
        else:
            x[k]=getAk(k,a,b,cf=cf_log_x_div_k)*vk_call(a,b,K,k)
    for k in range(N):
        sum+=x[k]
    Ak=np.asarray([getAk(k,a,b) for k in range(N)])
    Vk=np.asarray([vk_call(a,b,K,k) for k in range(N)])
    sum_mult=np.full(N,1.0)
    sum_mult[0]=0.5
    #return sum,np.sum(Ak*Vk*sum_mult)
    return np.exp(-r*T)*sum*(b-a)/2,np.sum(np.exp(-r*T)*Ak*Vk*sum_mult)*(b-a)/2
    #return np.exp(-r*T)*sum*(b-a)/2,

print(Cos_call(S0,K,T,r,sigma))

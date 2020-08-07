## European call FFT


S0=100.0
K=90.0
T=1.0
r=0.05
sigma=0.3
#---------------------------------------------------

def BSM_call_value_FFT(S0,K,T,r,sigma):
    k=np.log(K)
    x0=np.log(S0)
    N=2**10
    alpha=1.5
    eta=0.15
    lambda_=2*np.pi/(N*eta)
    beta=x0-0.5*lambda_*N
    k_m=np.asarray([beta+i*lambda_ for i in range(N)])
    W=SimpsonW(N,eta)
    v=np.asarray([i*eta for i in range(N)])
    Psi=np.asarray([BSM_call_characteristic_function(vj,alpha,x0,T,r,sigma) for vj in v])
    FFTFunc= Psi*np.exp(-1j*beta*v)*W
    y=fft(FFTFunc).real
    cT=np.exp(-alpha*k_m)*y/np.pi
    return np.exp(k_m),cT        ##返回不同k_m下不同cT
k,c = BSM_call_value_FFT(S0,K,T,r,sigma)
print(np.interp(K,k,c))  #插补法

#画图
def BSM_call_value_FFT(S0, K, T, r, sigma):
    k = np.log(K)
    x0 = np.log(S0)
    N = 2 ** 10
    alpha = 1.5
    eta = 0.15
    lambda_ = 2 * np.pi / (N * eta)
    beta = x0 - 0.5 * lambda_ * N
    k_m = np.asarray([beta + i * lambda_ for i in range(N)])
    W = SimpsonW(N, eta)
    v = np.asarray([i * eta for i in range(N)])
    Psi = np.asarray([BSM_call_characteristic_function(vj, alpha, x0, T, r, sigma) for vj in v])
    FFTFunc = Psi * np.exp(-1j * beta * v) * W
    y = fft(FFTFunc).real
    cT = np.exp(-alpha * k_m) * y / np.pi
    return np.exp(k_m), cT
k,c = BSM_call_value_FFT(S0,K,T,r,sigma)
print(np.interp(K,k,c))
plt.figure(figsize=(20,15))
x=np.array(range(80,150,5))
y=np.interp(x,k,c)
plt.scatter(x,y)
plt.show()

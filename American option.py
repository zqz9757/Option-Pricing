### 美式期权定价
#Gauss-Seidel
import sys
class FDCnAm(FDCnEu):
    def __init__(self,S0, K, r, T, sigma, Smax, M, N, omega, tol, is_call=True): #Omega是超松弛参数，omega越大收敛的越快
        super(FDCnAm,self).__init__(S0,K,r,T ,sigma, Smax, M, N, is_call)
        self.omega = omega
        self.tol=tol
        self.i_values=np.arange(self.M+1)
        self.j_values=np.arange(self.N+1)
    def _setup_boundary_conditions_(self):
        if self.is_call:
            self.payoffs=np.maximum(self.boundary_conds[1:self.M]-self.K,0)
        else:
            self.payoffs=np.maximum(self.K-self.boundary_conds[1:self.M],0)
        self.past_values=self.payoffs
        self.boundary_values=self.K*np.exp(-self.r*self.dt*(self.N-self.j_values)) #K折现
    def _traverse_grid_(self):
        aux=np.zeros(self.M-1)
        new_values=np.zeros(self.M-1)

        for j in reversed(range(self.N)):
            aux[0]=self.alpha[1]*\
                   (self.boundary_conds[j]+self.boundary_conds[j+1])
            rhs=np.dot(self.M2,self.past_values)+aux
            old_values=np.copy(self.past_values)
            error=sys.float_info.max

            while self.tol<error:
                new_values[0]=max(self.payoffs[0],old_values[0]+\
                                  self.omega/(1-self.beta[1])*\
                                  (rhs[0]-(1-self.beta[1])*old_values[0]+(self.gamma[1]*old_values[1])))
                for k in range(self.M-2)[1:]:
                    new_values[k]=max(self.payoffs[k],old_values[k]+\
                                      self.omega/(1-self.beta[k+1])*\
                                      (rhs[k]+self.alpha[k + 1] * new_values[k-1]-\
                                       (1-self.beta[k+1])*old_values[k]+self.gamma[k+1]*old_values[k+1]))
                new_values[-1]=max(self.payoffs[-1],old_values[-1]+self.omega/(1-self.beta[-2])* \
                                   (rhs[-1] + self.alpha[-2] * new_values[-2] - \
                                    (1 - self.beta[-2]) * old_values[-1]))
                error=np.linalg.norm(new_values-old_values)
                old_values=np.copy(new_values)
                self.past_values=np.copy(new_values)
                self.values=np.concatenate(([self.boundary_conds[0]],new_values,[0]))
     def _interpolate_(self):
     return np.interp(self.S0,self.boundary_conds,self.values)
option_6=FDCnAm(50,50,0.1,5./12, 0.4, 100, 100, 42, 1.2, 0.001)
print(option_6.price())




# 不定义类
#显式差分格式

import time
def AmFDCEm(S,E,D0,r,sigma,T,M,N,is_call):
    M=int(M)
    N=int(N)
    grid = np.zeros(shape=(M + 1, N + 1))  # 定义网格
    boundary_conds = np.linspace(0, 1, M + 1)  # 做插值
    i_values = np.arange(M)  # arange函数返回一个数据
    j_values = np.arange(N)

    dksi = 1 / float(M)
    dt = T / float(N)

    ksi = S/(S+E)
    start=time.time()
    if is_call==True:
        grid[:, 0] = np.maximum(2* boundary_conds-1, 0)
    else:
        grid[:, 0] = np.maximum(1 - 2 * boundary_conds, 0)

    grid[0, :-1] = grid[0, 0] * np.exp(-r * j_values * dt)
    grid[0,-1]=np.exp(-r*T)
    grid[-1, :-1] = grid[-1, 0] * np.exp(-D0 * j_values * dt)
    grid[-1,-1]=np.exp(-D0*T)

    alpha = dt / (dksi ** 2)
    ksi_values = i_values * dksi
    a = 0.5 * (sigma ** 2 * ksi_values ** 2 * (1 - ksi_values) ** 2 -
               (r - D0) * ksi_values * (1 - ksi_values) * dksi) * alpha
    b = 1 - sigma ** 2 * ksi_values ** 2 * (1 - ksi_values) ** 2 * alpha -\
        (r * (1 - ksi_values) + D0 * ksi_values) * dt
    c = 0.5 * (sigma ** 2 * ksi_values ** 2 * (1 - ksi_values) ** 2 +\
        (r - D0) * ksi_values * (1 - ksi_values) * dksi) * alpha
    for j in j_values:  # j=0,1,...,M-1
        for i in range(M)[1:]:  # i=1,2,...,M-1
            grid[i, j + 1] = a[i] * grid[i - 1, j] + \
                             b[i] * grid[i, j] + \
                             c[i] * grid[i + 1, j]
            if is_call==True:
                grid[i, j+1]=np.maximum(grid[i, j+1], np.maximum(2*ksi_values[i]-1 , 0))
            else:
                grid[i, j + 1] = np.maximum(grid[i, j + 1], np.maximum(1 - 2 * ksi_values[i], 0))
    V_ = np.interp(ksi, boundary_conds,grid[:, -1])
    V = (S + E) * V_
    end=time.time()
    cpu=end-start
    C=9.94092345
    P=5.92827717
    Errors=C-V
    Errors2=P-V
    return V ,cpu,Errors,Errors2

#隐式格式,直接法
def AmAmFDCIm_1(S,E,D0,r,sigma,T,M,N,is_call):
    M = int(M)
    N = int(N)
    grid = np.zeros(shape=(M + 1, N + 1))  # 定义网格
    boundary_conds = np.linspace(0, 1, M + 1)  # 做插值
    i_values = np.arange(M+1)  # arange函数返回一个数据
    j_values = np.arange(N+1)

    dksi = 1 / float(M)
    dt = T / float(N)

    ksi = S / (S + E)
    start=time.time()
    if is_call == True:
        grid[:, 0] = np.maximum(2 * boundary_conds - 1, 0)
    else:
        grid[:, 0] = np.maximum(1 - 2 * boundary_conds, 0)


    ksi_values = i_values * dksi

    a=dt/4*((r-D0)*i_values*(1-ksi_values)-sigma**2 * i_values **2 * (1-ksi_values)**2)
    b=1+dt/2*(sigma**2 * i_values **2 * (1-ksi_values)**2 + r * (1-ksi_values) + D0 * ksi_values)
    c=dt/4*(-(r-D0) * i_values * (1-ksi_values)- sigma**2 * i_values **2 * (1-ksi_values)**2 )

    coeffs = np.diag(a[1:], -1) + \
             np.diag(b) + \
             np.diag(c[:M], 1)

    aux=np.diag(-a[1:], -1) + np.diag(2-b) + np.diag(-c[:M], 1)

    import scipy
    for j in j_values[:N]:
        qj=np.dot(aux,grid[:,j])

        P, L, U = scipy.linalg.lu(coeffs)
        x1=scipy.linalg.solve(L,qj)
        grid[:,j+1]=scipy.linalg.solve(U,x1)
        for i in range(M+1):
            if is_call == True:
                grid[i,j+1] = np.maximum(grid[i, j+1], np.maximum(2*ksi_values[i]-1 , 0))
            else:
                grid[i,j+1] = np.maximum(grid[i, j + 1], np.maximum(1-2 * ksi_values[i], 0))
    V_ = np.interp(ksi, boundary_conds, grid[:, -1])
    V = (S + E) * V_
    end=time.time()
    C=9.94092345
    P=5.92827717
    Error=C-V
    Error2=P-V
    cpu=end-start
    return V, Error, Error2, cpu

#隐式格式，迭代法
def AmAmFDCIm_2(S,E,D0,r,sigma,T,M,N,omega,tol,is_call=True):
    M = int(M)
    N = int(N)
    grid = np.zeros(shape=(M + 1, N + 1))  # 定义网格
    boundary_conds = np.linspace(0, 1, M + 1)  # 做插值
    i_values = np.arange(M+1)  # arange函数返回一个数据
    j_values = np.arange(N+1)

    dksi = 1 / float(M)
    dt = T / float(N)

    ksi = S / (S + E)

    start=time.time()
    if is_call == True:
        grid[:, 0] = np.maximum(2 * boundary_conds - 1, 0)
    else:
        grid[:, 0] = np.maximum(1 - 2 * boundary_conds, 0)

    grid[0, :] = grid[0, 0] * np.exp(-r * j_values * dt)

    grid[-1, :] = grid[-1, 0] * np.exp(-D0 * j_values * dt)

    ksi_values = i_values * dksi

    a=dt/4*((r-D0)*i_values*(1-ksi_values)-sigma**2 * i_values **2 * (1-ksi_values)**2)
    b=1+dt/2*(sigma**2 * i_values **2 * (1-ksi_values)**2 + r * (1-ksi_values) + D0 * ksi_values)
    c=dt/4*(-(r-D0) * i_values * (1-ksi_values)- sigma**2 * i_values **2 * (1-ksi_values)**2 )

    # coeffs = np.diag(a[1:], -1) + \
             #np.diag(b) + \
             #np.diag(c[:M], 1)

    aux=np.diag(-a[1:], -1) + np.diag(2-b) + np.diag(-c[:M], 1)

    for j in j_values[:N]:

        qj=np.dot(aux,grid[:,j])

        old_values=grid[:,j]
        new_values=np.copy(old_values)
        error=sys.float_info.max

        while error>tol:
            if is_call==True:
                new_values[1]=np.maximum((1-omega)*old_values[0]+ omega / b[0] * (qj[0]-a[0]*new_values[0]-c[0]*old_values[1]),
                                          np.maximum(2 * ksi_values[0] - 1, 0))
                for i in range(1,M):
                new_values[i]=np.maximum((1-omega)*old_values[i]+ omega / b[i] * (qj[i]-a[i]*new_values[i-1]-c[i]*old_values[i+1]),
                                          np.maximum(2 * ksi_values[i] - 1, 0))
                new_values[M]= grid[-1,j+1]
                error=np.linalg.norm(new_values-old_values)
                old_values=np.copy(new_values)
            else:
                new_values[1] = np.maximum(
                    (1 - omega) * old_values[0] + omega / b[0] * (qj[0] - a[0] * new_values[0] - c[0] * old_values[1]),
                    np.maximum(1-2 *ksi_values[0] , 0))
                for i in range(1, M):
                    new_values[i] = np.maximum((1 - omega) * old_values[i] + omega / b[i] * (
                                qj[i] - a[i] * new_values[i - 1] - c[i] * old_values[i + 1]),
                                               np.maximum(1-2 * ksi_values[i] , 0))
                new_values[M] = grid[-1, j + 1]
                error = np.linalg.norm(new_values - old_values)
                old_values = np.copy(new_values)
        grid[:,j+1]=new_values

    V_ = np.interp(ksi, boundary_conds, grid[:, -1])
    V = (S + E) * V_

    end = time.time()
    C = 9.94092345
    P = 5.92827717
    Error = C - V
    Error2 = P - V
    cpu = end - start

    return V, Error,Error2,cpu



#蒙特卡洛
import scipy.stats as ss
import numpy as np

def LSM(S0, K, T, r, sig, D0, payoff, N, paths, order):

    start=time.time()
    dt = T/(N-1)          # time interval
    df = np.exp(-r * dt)  # discount factor per time time interval

    X0 = np.zeros((paths,1))
    increments = ss.norm.rvs(loc=(r- D0-0.5*sig**2)*dt, scale=np.sqrt(dt)*sig, size=(paths,N-1))
    X = np.concatenate((X0,increments), axis=1).cumsum(1)
    S = S0 * np.exp(X)
    if payoff == "put":
        H = np.maximum(K - S, 0)   # intrinsic values for put option
    if payoff == "call":
        H = np.maximum(S - K, 0)   # intrinsic values for call option
    V = np.zeros_like(H)            # value matrix
    V[:,-1] = H[:,-1]

    # Valuation by LS Method
    for t in range(N-2, 0, -1):
        good_paths = H[:,t] > 0  #
        # polynomial regression：将EV>0的部分挑出来回归
        rg = np.polyfit( S[good_paths, t], V[good_paths, t+1] * df, order)
        # 估计E(HV)
        C = np.polyval( rg, S[good_paths,t] )
        # 如果E(HV)<EV，那么行权
        exercise = np.zeros( len(good_paths), dtype=bool)
        exercise[good_paths] = H[good_paths,t] > C
        V[exercise,t] = H[exercise,t]
        V[exercise,t+1:] = 0
        # 剩下的是持续持有价格，折现即可
        discount_path = (V[:,t] == 0)
        V[discount_path,t] = V[discount_path,t+1] * df
    V0 = np.mean(V[:,1]) * df  #
    end = time.time()
    C = 9.94092345
    P = 5.92827717
    Error = C - V0
    Error2 = P - V0
    cpu = end - start
    return V0, Error,Error2,cpu

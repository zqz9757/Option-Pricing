def dN(x):
    return math.exp(-0.5 * x **2)/math.sqrt(2 * math.pi)
def N(d):
    return quad(lambda x: dN(x),-20,d,limit=50)[0]
def d1f(St,K,t,T,r,sigma):
    d1 = (math.log(St/K) + (r + 0.5*sigma**2)*(T - t)) /(sigma*math.sqrt(T - t))
    return d1
def BSM_call_value(St,K,t,T,r,sigma):
    d1 = d1f(St,K,t,T,r,sigma)
    d2 = d1 - sigma*math.sqrt(T - t)
    call_value = St*N(d1)-math.exp(-r*(T - t) )* K * N(d2)
    return  call_value
St=100.0
K=90.0
T=1.0
r=0.05
sigma=0.3
t=0
print(BSM_call_value(St,K,t,T,r,sigma))

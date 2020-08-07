### 障碍期权定价
##C-N格式 下降敲出期权
##只需修改边界条件即可
class FDCnDOC(FDCnEu):

    def __init__(self, S0, K, r, T, sigma, Sbarrier, Smax, M, N,is_call=True):
        super(FDCnEu,self).__init__(
            S0, K, r, T,sigma, Smax,M,N,is_call)
        self.dS=(Smax-Sbarrier)/float(self.M)
        self.boundary_conds=np.linspace(Sbarrier,Smax,self.M+1)
        self.i_values=self.boundary_conds/self.dS
option_4=FDCnDOC(50, 50, 0.1, 5./12, 0.4, 40, 100, 120, 500)
print(option_4.price())

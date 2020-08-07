##欧式看涨（跌）期权显式有限差分
class FiniteDifferences(object):          #定义一个有限差分类
    def __init__(self, S0, K, r, T, sigma, Smax, M, N,is_call=True):
        self.S0 = S0
        self.K = K
        self.r = r
        self.T = T
        self.sigma = sigma
        self.Smax = Smax
        self.M, self.N = int(M), int(N)  # Ensure M&N are integers
        self.is_call = is_call

        self.dS = Smax / float(self.M)
        self.dt = T / float(self.N)
        self.i_values = np.arange(self.M)    #arange函数返回一个数据
        self.j_values = np.arange(self.N)
        self.grid = np.zeros(shape=(self.M+1, self.N+1))
        self.boundary_conds = np.linspace(0, Smax, self.M+1)  #生成0到Smax之间的等间隔序列
    def _setup_boundary_conditions_(self):
        pass
    def _setup_coefficients_(self):
        pass
    def _traverse_grid_(self):
        """  Iterate the grid backwards in time """
        pass
    def _interpolate_(self):   #线性插补，找到离S0最近的初始列的格点
        """
        Use piecewise linear interpolation on the initial
        grid column to get the closest price at S0.
        """
        return np.interp(self.S0,
                         self.boundary_conds,
                         self.grid[:, 0])
    def price(self):
        self._setup_boundary_conditions_()
        self._setup_coefficients_()
        self._traverse_grid_()
        return self._interpolate_()

class FDExplicitEu(FiniteDifferences):
    def _setup_boundary_conditions_(self):  #定义边界条件
        if self.is_call:  #call
            self.grid[:, -1] = np.maximum(self.boundary_conds - self.K, 0)  #期末价值就是max(S-K,0)
            self.grid[-1, :-1] = (self.Smax - self.K) *np.exp(-self.r *self.dt *(self.N-self.j_values))
        else:   #put
            self.grid[:, -1] = np.maximum(self.K-self.boundary_conds, 0)
            self.grid[0, :-1] = (self.K - self.Smax) * np.exp(-self.r * self.dt *(self.N-self.j_values))
    def _setup_coefficients_(self):
        self.a = 0.5*self.dt*((self.sigma**2) * (self.i_values**2) - self.r*self.i_values)
        self.b = 1 - self.dt*((self.sigma**2) * (self.i_values**2) + self.r)
        self.c = 0.5*self.dt*((self.sigma**2) * (self.i_values**2) + self.r*self.i_values)
    def _traverse_grid_(self):
        for j in reversed(self.j_values):
            for i in range(self.M)[2:]:
                self.grid[i,j] = self.a[i]*self.grid[i-1,j+1] +\
                                 self.b[i]*self.grid[i,j+1] + \
                                 self.c[i]*self.grid[i+1,j+1]
option_1= FDExplicitEu(100, 90, 0.05, 1., 0.3, 250, 100, 1000, True)
print("European_put: ", option_1.price())
option_1 = FDExplicitEu(100, 90, 0.05, 1., 0.3, 250, 100, 100, True)
print("European_put: ", option_1.price())

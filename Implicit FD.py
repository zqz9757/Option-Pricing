#欧式看涨（跌）期权隐式有限差分

import scipy.linalg as linalg
class FDImplicitEu(FDExplicitEu):
    def _setup_coefficients_(self):
        self.a = 0.5*(self.r*self.dt*self.i_values -(self.sigma**2)*self.dt*(self.i_values**2))
        self.b = 1 + (self.sigma**2)*self.dt*(self.i_values**2) + self.r*self.dt
        self.c = -0.5*(self.r * self.dt*self.i_values +(self.sigma**2)*self.dt*(self.i_values**2))
        self.coeffs = np.diag(self.a[2:self.M], -1) + \
                      np.diag(self.b[1:self.M]) + \
                      np.diag(self.c[1:self.M-1], 1)           #系数矩阵A
    def _traverse_grid_(self):
        """ Solve using linear systems of equations """
        P, L, U = linalg.lu(self.coeffs)
        aux = np.zeros(self.M-1)
        for j in reversed(range(self.N)):
            aux[0] = np.dot(-self.a[1], self.grid[0, j])
            x1 = linalg.solve(L, self.grid[1:self.M, j+1]+aux)
            x2 = linalg.solve(U, x1)
            self.grid[1:self.M, j] = x2
option_2 = FDImplicitEu(100, 90, 0.05, 1., 0.3, 250, 100, 1000, True)
print(option_2.price())
option_2 = FDImplicitEu(100, 90, 0.05, 1., 0.3, 250, 100, 100, True)
print(option_2.price())

#Crank-Nicolson解法

class FDCnEu(FDExplicitEu):

    def _setup_coefficients_(self):
        self.alpha = 0.25*self.dt*((self.sigma**2)*(self.i_values**2) - self.r*self.i_values)
        self.beta = -self.dt*0.5*((self.sigma**2)*(self.i_values**2) + self.r)
        self.gamma = 0.25*self.dt*((self.sigma**2)*(self.i_values**2) + self.r*self.i_values)
        self.M1 = -np.diag(self.alpha[2:self.M], -1) + \
                  np.diag(1-self.beta[1:self.M]) - \
                  np.diag(self.gamma[1:self.M-1], 1)
        self.M2 = np.diag(self.alpha[2:self.M], -1) + \
                  np.diag(1+self.beta[1:self.M]) + \
                  np.diag(self.gamma[1:self.M-1], 1)

    def _traverse_grid_(self):
        """ Solve using linear systems of equations """
        P, L, U = linalg.lu(self.M1)

        for j in reversed(range(self.N)):
            x1 = linalg.solve(L,np.dot(self.M2,self.grid[1:self.M, j+1]))
            x2 = linalg.solve(U, x1)
            self.grid[1:self.M, j] = x2
option_3 = FDCnEu(100, 90, 0.05, 1., 0.3, 250, 100, 1000, True)
print(option_3.price())
option_3 = FDCnEu(100, 90, 0.05, 1., 0.3, 250, 100, 100, True)
print(option_3.price())

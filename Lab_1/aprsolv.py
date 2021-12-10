import numpy as np
import scipy.integrate
import matplotlib.pyplot as plt

class ApproxSolver:
    def __init__(self):
        self.a = 1
        self.b = 5
        self.b1 = 2
        self.b2 = 1
        self.b3 = 1
        self.k1 = 2
        self.k2 = 1
        self.c1 = 1
        self.c2 = 0
        self.c3 = 1
        self.p1 = 3
        self.p2 = 0
        self.d1 = 2
        self.d2 = 1
        self.d3 = 3
        self.q1 = 2
        self.q2 = 1
        self.a1 = 1
        self.a2 = 2
        self.a3 = 3
        self.a4 = 4
        self.n1 = 7
        self.n2 = 4
        self.n3 = 2
        self.alpha = self.a1*self.a**self.n1 + self.a2*self.a**self.n2 + self.a3*self.a**self.n3 + self.a4
        self.beta =  self.a1*self.n1*self.a**(self.n1 - 1) + self.a2*self.n2*self.a**(self.n2 - 1) + self.a3*self.n3*self.a**(self.n3 - 1)
        self.gamma = self.a1*self.b**self.n1 + self.a2*self.b**self.n2 + self.a3*self.b**self.n3 + self.a4
        self.sigma = -(self.a1*self.n1*self.b**(self.n1 - 1) + self.a2*self.n2*self.b**(self.n2 - 1) + self.a3*self.n3*self.b**(self.n3 - 1))
        self.X = np.poly1d([1, 0])
        self.A = self.b + self.gamma*(self.b-self.a)/(2*self.gamma+self.sigma*(self.b-self.a))
        self.B = self.a + self.alpha*(self.b-self.a)/(2*self.alpha-self.beta*(self.b-self.a))

    def k(self):
        return self.b1*self.X**self.k1 + self.b2*self.X**self.k2 + self.b3

    def p(self):
        return self.c1*self.X**self.p1 + self.c2*self.X**self.p2 + self.c3

    def q(self):
        return self.d1*self.X**self.q1 + self.d2*self.X**self.q2 + self.d3

    def power(self, poly, p):
        if p >= 0:
            return poly**p
        else:
            return np.poly1d([1])/(poly**abs(p))

    def f(self):
        return -ApproxSolver.k(self)*(self.a1*self.n1*(self.n1-1)*ApproxSolver.power(self, self.X, self.n1-2)+self.a2*self.n2*(self.n2-1)*ApproxSolver.power(self, self.X, self.n2-2)+
                                        +self.a3*self.n3*(self.n3-1)*ApproxSolver.power(self, self.X, self.n3-2))-(self.b1*self.k1*ApproxSolver.power(self, self.X, self.k1-1)+
                                        + self.b2*self.k2*ApproxSolver.power(self, self.X, self.k2-1))*(self.a1*self.n1*ApproxSolver.power(self, self.X, self.n1-1)+self.a2*self.n2*ApproxSolver.power(self, self.X, self.n2-1)+
                                        +self.a3*self.n3*ApproxSolver.power(self, self.X, self.n3-1))+ApproxSolver.p(self)*(self.a1*self.n1*ApproxSolver.power(self, self.X, self.n1-1)+self.a2*self.n2*ApproxSolver.power(self, self.X, self.n2-1)+
                                        +self.a3*self.n3*ApproxSolver.power(self, self.X, self.n3-1))+ApproxSolver.q(self)*(self.a1*ApproxSolver.power(self, self.X, self.n1)+self.a2*ApproxSolver.power(self, self.X, self.n2)+self.a3*ApproxSolver.power(self, self.X, self.n3) + self.a4)

    def Phi(self, n):
        if n < 2: return
        phi = []
        phi.append(np.poly1d([1,-self.a])**2*np.poly1d([1,-self.A]))
        phi.append(np.poly1d([1,-self.B])*np.poly1d([-1, self.b])**2)
        for i in range(2, n):
            phi.append(np.poly1d([1, -self.a])**(i-1)*np.poly1d([-1, self.b])**2)
        return phi

    def A(self, func):
        return -np.polyder(ApproxSolver.k(self)*np.polyder(func)) + ApproxSolver.p(self)*np.polyder(func) + ApproxSolver.q(self)*func

    def real(self):
        return self.a1*ApproxSolver.power(self, self.X, self.n1)+\
               +self.a2*ApproxSolver.power(self, self.X, self.n2)+\
               + self.a3*ApproxSolver.power(self, self.X, self.n3) + self.a4

    def BubnGal(self, n):
        phi_i = ApproxSolver.Phi(self, n)
        right = []
        for j in range(0, n):
            right.append(scipy.integrate.quad(ApproxSolver.f(self)*phi_i[j], self.a, self.b)[0])

        left = []
        for j in range(0, n):
            temp = []
            for i in range(0,n):
                temp.append(scipy.integrate.quad(ApproxSolver.A(self, phi_i[i])*phi_i[j], self.a, self.b)[0])
            left.append(temp)
        C = np.linalg.solve(left, right)

        res = np.poly1d([1])
        for i in range(0, n):
            res += C[i]*phi_i[i]
        delta = abs(np.sqrt(scipy.integrate.quad((1/(self.b-self.a))*(res - ApproxSolver.real(self))**2, self.a, self.b)[0]))
        x = np.arange(self.a, self.b, 0.001)
        y = np.polyval(res, x)
        y_res = np.polyval(ApproxSolver.real(self), x)
        plt.plot(x, y, color='red', label='Accurate u')
        plt.plot(x, y_res, color='green', label='u')
        plt.legend()
        plt.title(f'Galerkin method, n = {n}, $\\vartriangle$={delta}')
        plt.show()


    def Ritz(self, n):
        self.c1 = 0
        self.c2 = 0
        self.c3 = 0
        phi_i = ApproxSolver.Phi(self, n)
        right = []
        for j in range(0, n):
            right.append(scipy.integrate.quad(ApproxSolver.f(self) * phi_i[j], self.a, self.b)[0])

        left = []
        for j in range(0, n):
            temp = []
            for i in range(0, n):
                temp.append(scipy.integrate.quad(ApproxSolver.A(self, phi_i[i]) * phi_i[j], self.a, self.b)[0])
            left.append(temp)
        C = np.linalg.solve(left, right)

        res = np.poly1d([1])
        for i in range(0, n):
            res += C[i] * phi_i[i]
        delta = abs(np.sqrt(scipy.integrate.quad((1/(self.b-self.a))*(res - ApproxSolver.real(self)) ** 2, self.a, self.b)[0]))
        x = np.arange(self.a, self.b, 0.001)
        y = np.polyval(res, x)
        y_res = np.polyval(ApproxSolver.real(self), x)
        plt.plot(x, y, color='red')
        plt.plot(x, y_res, color='blue')
        # plt.legend()
        plt.title(f'Метод Рітца, n = {n}, deviation = {delta}')
        plt.show()












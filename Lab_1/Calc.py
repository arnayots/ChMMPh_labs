import math
import numpy as np
import matplotlib.pyplot as plt

class Calc:

    def __init__(self):
        print('hi')

        self.k = 3

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

        self.alpha = self.a1 * self.a ** self.n1 \
                     + self.a2 * self.a ** self.n2 \
                     + self.a3 * self.a ** self.n3 \
                     + self.a4

        self.beta = self.a1 * self.n1 * self.a ** (self.n1 - 1) \
                    + self.a2 * self.n2 * self.a ** (self.n2 - 1) \
                    + self.a3 * self.n3 * self.a ** (self.n3 - 1)

        self.gamma = self.a1 * self.b ** self.n1 \
                     + self.a2 * self.b ** self.n2 \
                     + self.a3 * self.b ** self.n3 \
                     + self.a4

        self.delta = -(self.a1 * self.n1 * self.b ** (self.n1 - 1)
                     + self.a2 * self.n2 * self.b ** (self.n2 - 1)
                     + self.a3 * self.n3 * self.b ** (self.n3 - 1))

        self.k = self.get_k_str()
        self.p = self.get_p_str()
        self.q = self.get_q_str()

        self.phi_type = 1
        self.n_phi = 3

        self.colloc_pnts = []
        # print(self.k)

    def get_k_str(self):
        tmp = f'{self.b1} * x ^ {self.k1} + {self.b2} * x ^ {self.k2} + {self.b3}'
        return tmp

    def get_p_str(self):
        tmp = f'{self.c1} * x ^ {self.p1} + {self.c2} * x ^ {self.p2} + {self.c3}'
        return tmp

    def get_q_str(self):
        tmp = f'{self.d1} * x ^ {self.q1} + {self.d2} * x ^ {self.q2} + {self.d3}'
        return tmp

    def print_k(self):
        print(self.k)

    def k_x(self, x):
        return self.b1 * x ** self.k1 + self.b2 * x ** self.k2 + self.b3

    def dk_x(self, x):
        return self.b1 * self.k1 * x ** (self.k1 - 1) + self.b2 * self.k2 * x ** (self.k2 - 1)

    def p_x(self, x):
        return self.c1 * x ** self.p1 + self.c2 * x ** self.p2 + self.c3

    def q_x(self, x):
        return self.d1 * x ** self.q1 + self.d2 * x ** self.q2 + self.d3



    def phi_i_x(self, i, x):
        if self.phi_type == 1:
            return ((x - self.a) ** i) * (self.b - x)
        elif self.phi_type == 2:
            return (x - self.a) * ((self.b - x) ** i)
        elif self.phi_type == 3:
            return math.sin((x - self.a) / (self.b - self.a) * i * math.pi)
        else:
            return 0

    def d_phi_i_x(self, i, x):
        if self.phi_type == 1:
            return (x - self.a) ** (i - 1) * (self.a + i * (self.b - x) - x)
        elif self.phi_type == 2:
            return (self.b - x) ** (i - 1) * (self.a * i + self.b - (i + 1) * x)
        elif self.phi_type == 3:
            return math.pi * i * math.cos(math.pi * i * (self.a - x) / (self.a - self.b)) / (self.b - self.a)
        else:
            return 0

    def d2_phi_i_x(self, i, x):
        if self.phi_type == 1:
            return -i * (x - self.a) ** (i - 2) * (-2 * self.a - self.b * i + self.b + i * x + x)
        elif self.phi_type == 2:
            return i * (self.b - x) ** (i - 2) * (self.a * (i - 1) + 2 * self.b - (i + 1) * x)
        elif self.phi_type == 3:
            return i ** 2 * math.pi ** 2 * math.sin(math.pi * i * (self.a - x) - (self.a - self.b)) / ((self.b - self.a) ** 2)
        else:
            return 0

    def get_const_i_multiplier_x(self, i, x):
        return -self.k_x(x) * self.d2_phi_i_x(i, x) + (self.p_x(x) - self.dk_x(x)) * self.d_phi_i_x(i, x) + self.q_x(x) * self.phi_i_x(i, x)

    def get_f_x(self, x):
        tmp = -self.k_x(x) * (self.a1 * self.n1 * (self.n1 - 1) * x ** (self.n1 - 2)
                              + self.a2 * self.n2 * (self.n2 - 1) * x ** (self.n2 - 2)
                              + self.a3 * self.n3 * (self.n3 - 1) * x ** (self.n3 - 2))
        tmp += - (self.b1 * self.k1 * x ** (self.k1 - 1) + self.b2 * self.k2 * x ** (self.k2 - 1)) \
               * (self.a1 * self.n1 * x ** (self.n1 - 1) + self.a2 * self.n2 * x ** (self.n2 - 1) + self.a3 * self.n3 * x ** (self.n3 - 1))
        tmp += self.p_x(x) * (self.a1 * self.n1 * x ** (self.n1 - 1) + self.a2 * self.n2 * x ** (self.n2 - 1) + self.a3 * self.n3 * x ** (self.n3 - 1))
        tmp += self.q_x(x) * (self.a1 * x ** self.n1 + self.a2 * x ** self.n2 + self.a3 * x ** self.n3 + self.a4)
        return tmp

    def solve(self):
        N = len(self.colloc_pnts)
        if N == 0:
            return False

        A = np.zeros((N, N))
        F = np.zeros((N, 1))

        for i in range(1, N + 1):
            cur_pnt = self.colloc_pnts[i - 1]
            F[i - 1][0] = self.get_f_x(cur_pnt)
            for j in range(1, self.n_phi + 1):
                A[i - 1][j - 1] = self.get_const_i_multiplier_x(j, cur_pnt)

        # print(A)
        # print(F)

        res = np.linalg.solve(A, F)
        self.const_res_vector = res
        # print(res)
        return res

    def solution_real_x(self, x):
        return self.a1 * x ** self.n1 + self.a2 * x ** self.n2 + self.a3 * x ** self.n3 + self.a4

    def solution_modeled_x(self, x):
        res = 0
        for i in range(0, len(self.const_res_vector)):
            res += self.const_res_vector[i] * self.phi_i_x(i + 1, x)
        return res


    def print_comparation_plot(self):
        x = np.linspace(self.a, self.b, 100)
        y_1 = self.solution_modeled_x(x)
        y_2 = self.solution_real_x(x)

        plt.plot(x, y_1, 'red',  x, y_2, 'blue')
        plt.show()



    def set_phi_type(self, phi):
        self.phi_type = phi

    def get_phi_type(self):
        return(self.phi_type)

    def push_colloc_pnt(self, pnt):
        self.colloc_pnts.append(pnt)



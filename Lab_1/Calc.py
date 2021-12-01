import math
import numpy as np
import matplotlib.pyplot as plt

def near(x, y):
    return abs(x - y) < 1E-4

def calc_integral(f, a, b, n_pts = 500):
    delta = (b - a) / n_pts
    res = 0
    for i in range(1, n_pts):
        res += f(a + delta * i) * delta
    res += 0.5 * f(a) * delta
    res += 0.5 * f(b) * delta
    return res

class Calc:
    TASK_TYPE_COLLOC = 0
    TASK_TYPE_RITZ = 1

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
        self.a2 = -3
        self.a3 = 3
        self.a4 = 4
        self.n1 = 1
        self.n2 = 2
        self.n3 = 3

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

        self.mu_1 = 0
        self.mu_2 = 0
        self.A_param = self.b + (self.gamma * (self.b - self.a)) / (2 * self.gamma + self.delta * (self.b - self.a))
        self.B_param = self.a + (self.a * (self.b - self.a)) / (2 * self.alpha - self.beta * (self.b - self.a))

        self.k = self.get_k_str()
        self.p = self.get_p_str()
        self.q = self.get_q_str()

        self.phi_type = 1
        self.n_phi = 0
        self.task_type = 0

        self.colloc_pnts = []
        # print(self.k)

    def set_task_type(self, type):
        self.task_type = type

    def get_task_type(self):
        return self.task_type

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
        return self.b1 * (x ** self.k1) + self.b2 * (x ** self.k2) + self.b3

    def dk_x(self, x):
        return self.b1 * self.k1 * (x ** (self.k1 - 1)) + self.b2 * self.k2 * (x ** (self.k2 - 1))

    def p_x(self, x):
        if near(x, self.a):
            return - self.beta
        if near(x, self.b):
            return self.delta
        return self.c1 * (x ** self.p1) + self.c2 * (x ** self.p2) + self.c3

    def q_x(self, x):
        return self.d1 * (x ** self.q1) + self.d2 * (x ** self.q2) + self.d3

    def phi_i_x(self, i, x):
        if i == 1:
            return (x - self.a) ** 2 * (x - self.A_param)
        if i == 2:
            return (x - self.B_param) * (self.b - x) ** 2

        if self.phi_type == 1:
            return (x - self.a) ** 2 * (self.b - x) ** (i - 1)
            # return ((x - self.a) ** i) * (self.b - x)
        elif self.phi_type == 2:
            return ((x - self.a) ** (i - 1)) * ((self.b - x) ** 2)
            # return (x - self.a) * ((self.b - x) ** i)
        elif self.phi_type == 3:
            return False
            # return math.sin((x - self.a) / (self.b - self.a) * i * math.pi)
        else:
            return 0

    def d_phi_i_x(self, i, x):
        if i == 1:
            return (x - self.a) * (3 * x - self.a - 2 * self.A_param)
        if i == 2:
            return (x + 2 * self.B_param - 3 * x) * (self.b - x)

        if self.phi_type == 1:
            return 2 * (x - self.a) * (self.b - x) ** (i - 1)
            # return (x - self.a) ** (i - 1) * (self.a + i * (self.b - x) - x)
        elif self.phi_type == 2:
            return (self.b - x) * (x - self.a) ** (i - 2) * (self.a * 2 + self.b * (i - 1) - (i + 1) * x)
            # return (i - 1) * (self.b - x) ** 2 * (x - self.a) ** (i - 2) - 2 * (self.b - x) - (x - self.a) ** (i - 1)
            # return (self.b - x) ** (i - 1) * (self.a * i + self.b - (i + 1) * x)
        elif self.phi_type == 3:
            return False
            # return math.pi * i * math.cos(math.pi * i * (self.a - x) / (self.a - self.b)) / (self.b - self.a)
        else:
            return 0

    def d2_phi_i_x(self, i, x):
        if i == 1:
            return -2 * (2 * self.a + self.A_param - 3 * x)
        if i == 2:
            return -2 * (2 * self.b + self.B_param - 3 * x)

        if self.phi_type == 1:
            return (i - 2) * (i - 1) * (x - self.a) ** 2 * (self.b - x) ** (i - 3) - 4 * (i - 1) * (x - self.a) * (self.b - x) ** (i - 2) + 2 * (self.b - x) ** (i - 1)
            # return -i * (x - self.a) ** (i - 2) * (-2 * self.a - self.b * i + self.b + i * x + x)
        elif self.phi_type == 2:
            return (i - 2) * (i - 1) * (self.b - x) ** 2 * (x - self.a) ** (i - 3) - 4 * (i - 1) * (self.b - x) * (x - self.a) ** (i - 2) + 2 * (x - self.a) * (i - 1)
        #     return i * (self.b - x) ** (i - 2) * (self.a * (i - 1) + 2 * self.b - (i + 1) * x)
        elif self.phi_type == 3:
            return False
            # return i ** 2 * math.pi ** 2 * math.sin(math.pi * i * (self.a - x) - (self.a - self.b)) / ((self.b - self.a) ** 2)
        else:
            return 0

    def u_x(self, x):
        sum = 0
        for i in range(0, len(self.colloc_pnts)):
            sum += self.const_res_colloc_vector[i] * self.phi_i_x(i + 1, x)
        return sum
    def du_x(self, x):
        sum = 0
        for i in range(0, len(self.colloc_pnts)):
            sum += self.const_res_colloc_vector[i] * self.d_phi_i_x(i + 1, x)
        return sum
    def d2_u_x(self, x):
        sum = 0
        for i in range(0, len(self.colloc_pnts)):
            sum += self.const_res_colloc_vector[i] * self.d2_phi_i_x(i + 1, x)
        return sum

    def get_const_i_multiplier_x(self, i, x):
        return -self.k_x(x) * self.d2_phi_i_x(i, x) + (self.p_x(x)
                                                       - self.dk_x(x)
                                                       ) * self.d_phi_i_x(i, x) + self.q_x(x) * self.phi_i_x(i, x)

    def get_f_x(self, x):
        if near(x, self.a):
            return self.mu_1
        if near(x, self.b):
            return self.mu_2
        tmp = -self.k_x(x) * (self.a1 * self.n1 * (self.n1 - 1) * (x ** (self.n1 - 2))
                              + self.a2 * self.n2 * (self.n2 - 1) * (x ** (self.n2 - 2))
                              + self.a3 * self.n3 * (self.n3 - 1) * (x ** (self.n3 - 2)))
        tmp += - (self.b1 * self.k1 * (x ** (self.k1 - 1)) + self.b2 * self.k2 * (x ** (self.k2 - 1))) \
               * (self.a1 * self.n1 * (x ** (self.n1 - 1)) + self.a2 * self.n2 * x ** (self.n2 - 1) + self.a3 * self.n3 * (x ** (self.n3 - 1)))
        tmp += self.p_x(x) * (self.a1 * self.n1 * (x ** (self.n1 - 1)) + self.a2 * self.n2 * (x ** (self.n2 - 1)) + self.a3 * self.n3 * (x ** (self.n3 - 1)))
        tmp += self.q_x(x) * (self.a1 * (x ** self.n1) + self.a2 * (x ** self.n2) + self.a3 * (x ** self.n3) + self.a4)
        return tmp

    def calc_int_a_phi_phi(self, i, j):
        tmp_func = lambda x : (-self.k_x(x) * self.d2_phi_i_x(i, x) - self.dk_x(x) * self.d_phi_i_x(i, x) + self.q_x(x) * self.phi_i_x(i, x)) * self.phi_i_x(j, x)
        return calc_integral(tmp_func, self.a, self.b)

    def calc_int_f_phi(self, j):
        tmp_func = lambda x: self.get_f_x(x) * self.phi_i_x(j, x)
        return calc_integral(tmp_func, self.a, self.b)

    def solve_colloc(self):
    # solves a task with method of collocations
    # result vector of constants saves to self.const_res_colloc_vector
    # if result is correct, function also returns it. In other way, returns False
        self.task_type = 0
        N = len(self.colloc_pnts)
        if N == 0:
            return False

        A = np.zeros((N, N))
        F = np.zeros((N, 1))

        for i in range(1, N + 1):
            cur_pnt = self.colloc_pnts[i - 1]
            F[i - 1][0] = self.get_f_x(cur_pnt)
            for j in range(1, len(self.colloc_pnts) + 1):
                A[i - 1][j - 1] = self.get_const_i_multiplier_x(j, cur_pnt)
        # print(A)
        # print(F)
        res = np.linalg.solve(A, F)
        self.const_res_colloc_vector = res
        # print(A[0][0] * res[0] + A[0][1] * res[1] +A[0][2] * res[2])
        # x_c = self.colloc_pnts[0]
        # print(self.get_f_x(x_c))
        # print('test res')
        # print(self.get_const_i_multiplier_x(1, x_c) * res[0] + self.get_const_i_multiplier_x(2, x_c) * res[1] + self.get_const_i_multiplier_x(3, x_c) * res[2])
        #
        # print(self.phi_i_x(1, x_c) * res[0] + self.phi_i_x(2, x_c) * res[1] + self.phi_i_x(3, x_c) * res[2])
        # print(self.solution_real_x(x_c))
        # print(res)
        return res

    def solve_ritz(self, N):
        self.task_type = 1

        A = np.zeros((N, N))
        F = np.zeros((N, 1))

        for j in range(1, N + 1):
            F[j - 1][0] = self.calc_int_f_phi(j)
            for i in range(1, N + 1):
                A[j - 1][i - 1] = self.calc_int_a_phi_phi(i, j)
        res = np.linalg.solve(A, F)
        self.const_res_ritz_vector = res
        return res


    def solution_real_x(self, x):
        return self.a1 * (x ** self.n1) + self.a2 * (x ** self.n2) + self.a3 * (x ** self.n3) + self.a4

    def solution_modeled_x(self, x):
        res = 0
        if self.task_type == 0:
            for i in range(0, len(self.const_res_colloc_vector)):
                res += self.const_res_colloc_vector[i] * self.phi_i_x(i + 1, x)
        elif self.task_type == 1:
            for i in range(0, len(self.const_res_ritz_vector)):
                res += self.const_res_ritz_vector[i] * self.phi_i_x(i + 1, x)
        return res

    def print_comparation_plot(self):
        x = np.linspace(self.a, self.b, 100)
        y_1 = self.solution_modeled_x(x)
        y_2 = self.solution_real_x(x)

        plt.plot(x, y_1, 'red',  x, y_2, 'blue')
        if self.task_type == 0:
            plt.title(f'Метод колокацій, n = {len(self.colloc_pnts)}, deviation = {self.calc_deviation()}')
        elif self.task_type == 1:
            plt.title(f'Метод Рітца, n = {len(self.const_res_ritz_vector)}, deviation = {self.calc_deviation()}')
        plt.show()

    def calc_discrepancy_x(self, x):
        return -self.k_x(x) * self.d2_u_x(x) + (self.p_x(x) - self.dk_x(x)) * self.du_x(x) + self.q_x(x) * self.u_x(x) - self.get_f_x(x)

    def print_discrepancy_table(self):
        if self.task_type == 0:
            print('--- Нев\'язка ---')
            for i in range(0, len(self.colloc_pnts)):
                x = self.colloc_pnts[i]
                print('| ' + "%4f" % x + ' | ' + "%12f" % self.calc_discrepancy_x(x) + ' | ')

    def calc_deviation(self):
        tmp_func = lambda x: (self.solution_real_x(x) - self.solution_modeled_x(x)) ** 2
        return math.sqrt(calc_integral(tmp_func, self.a, self.b))


    def set_phi_type(self, phi):
        self.phi_type = phi

    def get_phi_type(self):
        return(self.phi_type)

    def push_colloc_pnt(self, pnt):
        self.colloc_pnts.append(pnt)
        self.n_phi = len(self.colloc_pnts)


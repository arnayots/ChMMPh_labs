import numpy as np
import matplotlib.pyplot as plt

debug = True
debug = False

def near(x, y):
    return abs(x - y) < 1E-4

class calc:

    def __init__(self):
        print('init')

        # self.k = 3

        self.a = 1
        self.b = 3
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


        self.n1 = 1
        self.n2 = 2
        self.n3 = 3

        self.a1 = 1
        self.a2 = -3
        self.a3 = 2
        self.a4 = 5

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

    def k_x(self, x):
        return self.b1 * (x ** self.k1) + self.b2 * (x ** self.k2) + self.b3

    def k_i(self, i):
        return self.k_x(self.y_i[i])

    def p_x(self, x):
        return 0
        if near(x, self.a):
            return - self.beta
        if near(x, self.b):
            return self.delta
        return self.c1 * (x ** self.p1) + self.c2 * (x ** self.p2) + self.c3

    def p_i(self, i):
        return self.p_x(self.y_i[i])

    def q_x(self, x):
        return self.d1 * (x ** self.q1) + self.d2 * (x ** self.q2) + self.d3

    def q_i(self, i):
        return self.q_x(self.y_i[i])

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

    def __apply_new_grid__(self):
        self.h = (self.b - self.a) / (self.N + 1)
        self.y_i = []
        for i in range(0, self.N + 1):
            self.y_i.append(self.a + i * self.h)
        self.y_i.append(self.b)

    def __apply_new_a_vec__(self):
        self.vec_a = []
        for i in range(0, self.N + 2):
            self.vec_a.append(self.__get_a_i__(i))

    def __get_A_i__(self, i):
        if self.diff_type == 0:
            return (2 * self.vec_a[i] + self.vec_a[i + 1]) / (self.h ** 2)
        elif self.diff_type == 1:
            return self.vec_a[i] / (self.h ** 2) + 1 / self.h
        elif self.diff_type == 2:
            return self.vec_a[i] / (self.h ** 2)
        elif self.diff_type == 3:
            return self.vec_a[i] / (self.h ** 2)
        return 0

    def __get_B_i__(self, i):
        if self.diff_type == 0:
            return self.vec_a[i] / (self.h ** 2) - 1 / self.h
        elif self.diff_type == 1:
            return self.vec_a[i + 1] / (self.h ** 2)
        elif self.diff_type == 2:
            return self.vec_a[i + 1] / (self.h ** 2)
        elif self.diff_type == 3:
            return self.vec_a[i + 1] / (self.h ** 2) - self.p_i(i) / self.h
        return 0

    def __get_C_i__(self, i):
        if self.diff_type == 0:
            return - ((self.vec_a[i + 1] - 3 * self.vec_a[i]) / (self.h ** 2) + 1 / self.h)
        elif self.diff_type == 1:
            return self.vec_a[i + 1] / (self.h ** 2) + self.vec_a[i] / (self.h ** 2) + 1 / self.h
        elif self.diff_type == 2:
            return self.vec_a[i + 1] / (self.h ** 2) + self.vec_a[i] / (self.h ** 2) + self.q_i(i)
        elif self.diff_type == 3:
            return self.vec_a[i + 1] / (self.h ** 2) + self.vec_a[i] / (self.h ** 2) + self.q_i(i) - self.p_i(i) / self.h
        return 0

    def __get_F_i__(self, i):
        return self.get_f_x(self.y_i[i])

    def __get_x1__(self):
        if self.diff_type == 0:
            return self.vec_a[1] * self.y_i[1] / (self.vec_a[1] + self.h * self.sigma1_new)
        elif self.diff_type == 1 or self.diff_type == 2 or self.diff_type == 3:
            return self.vec_a[1] / (self.h * self.sigma1_new + self.vec_a[1])
        return 0

    def __get_x2__(self):
        if self.diff_type == 0:
            return self.vec_a[self.N] * self.y_i[self.N - 1] / (self.vec_a[self.N] + self.h * self.sigma2_new)
        elif self.diff_type == 1 or self.diff_type == 2 or self.diff_type == 3:
            return - self.vec_a[self.N] / (self.sigma2_new * self.h - self.vec_a[self.N])
        return 0

    def __get_mu1_new__(self):
        if self.diff_type == 0:
            return (self.mu_1 + 0.5 * self.h * self.__get_F_i__(0)) / (self.vec_a[0] / self.h + self.sigma1_new)
        elif self.diff_type == 1 or self.diff_type == 2 or self.diff_type == 3:
            return self.mu_1_middle / (self.sigma1_new + self.vec_a[1] / self.h)
        return 0

    def __get_mu2_new__(self):
        if self.diff_type == 0:
            return (self.mu_2 + 0.5 * self.h * self.__get_F_i__(self.N)) / (self.vec_a[self.N] / self.h + self.sigma2_new)
        elif self.diff_type == 1 or self.diff_type == 2 or self.diff_type == 3:
            return self.mu_2_middle / (self.sigma2_new - self.vec_a[self.N] / self.h)
        return 0

    def __get_sigma1_new__(self):
        if self.diff_type == 3:
            sigma1_old = self.beta * self.k_i(0) / self.alpha
            return self.h * 0.5 * self.q_i(0) + sigma1_old * (1 + self.h * 0.5 * self.p_i(0) / self.k_i(0))
            return sigma1_old + self.h * 0.5 * self.q_i(0) + sigma1_old * self.h * 0.5 * self.p_i(0) / self.k_i(0)
        return self.beta * self.k_x(self.y_i[0]) / self.alpha + 0.5 * self.h * self.q_x(self.y_i[0])

    def __get_sigma2_new__(self):
        if self.diff_type == 3:
            sigma2_old = self.delta * self.k_i(self.N) / self.gamma
            return self.h * 0.5 * self.q_i(self.N) + sigma2_old * (1 + self.h * 0.5 * self.p_i(self.N) / self.k_i(self.N))
            return sigma2_old + self.h * 0.5 * self.q_i(self.N) + sigma2_old * self.h * 0.5 * self.p_i(self.N) / self.k_i(self.N)
        return self.delta * self.k_x(self.y_i[self.N]) / self.gamma + 0.5 * self.h * self.q_x(self.y_i[self.N])

    def __get_mu_1_middle__(self):
        if self.diff_type == 3:
            return self.mu_1 * (1 + self.h * 0.5 * self.p_i(0) / self.k_i(0)) + self.h * 0.5 * self.__get_F_i__(0)
            return self.mu_1 + self.h * self.__get_F_i__(0) * 0.5 + self.mu_1 * self.h * 0.5 * self.p_i(0) / self.k_i(0)
        return - self.mu_1 * self.k_x(self.y_i[0]) + 0.5 * self.h * self.get_f_x(self.y_i[0])

    def __get_mu_2_middle__(self):
        if self.diff_type == 3:
            return self.mu_2 * (1 + self.h * 0.5 * self.p_i(self.N) / self.k_i(self.N)) + self.h * 0.5 * self.__get_F_i__(self.N)
            return self.mu_2 + self.h * self.__get_F_i__(self.N) * 0.5 + self.mu_2 * self.h * 0.5 * self.p_i(self.N) / self.k_i(self.N)
        return self.mu_2 * self.k_x(self.y_i[self.N]) + 0.5 * self.h * self.get_f_x(self.y_i[self.N])

    def __get_a_i__(self, i):
        return self.k_x(self.y_i[i] - 0.5 * self.h)
        # return self.k_x( - 0.5 * self.h)

    def __fill_matr_left__(self):
        self.left = np.zeros((self. N + 1, self. N + 1))
        self.left[0][0] = 1
        self.left[0][1] = - self.x_1
        self.left[self.N][self.N] = 1
        self.left[self.N][self.N - 1] = - self.x_2
        for i in range(1, self.N):
            self.left[i][i - 1] = - self.__get_A_i__(i)
            self.left[i][i] = - self.__get_C_i__(i)
            self.left[i][i + 1] = self.__get_B_i__(i)

    def __fill_vec_right__(self):
        self.right = np.zeros((self. N + 1, 1))
        self.right[0][0] = self.mu_1_new
        self.right[self.N][0] = self.mu_2_new
        for i in range(1, self.N):
            self.right[i][0] = - self.__get_F_i__(i)

    def solution_real_x(self, x):
        return self.a1 * (x ** self.n1) + self.a2 * (x ** self.n2) + self.a3 * (x ** self.n3) + self.a4

    def solve_MAIT(self, n, diff_type=0):
        self.N = n
        self.diff_type = diff_type
        self.__apply_new_grid__()
        self.__apply_new_a_vec__()
        self.sigma1_new = self.__get_sigma1_new__()
        self.sigma2_new = self.__get_sigma2_new__()
        self.mu_1_middle = self.__get_mu_1_middle__()
        self.mu_2_middle = self.__get_mu_2_middle__()
        self.x_1 = self.__get_x1__()
        self.x_2 = self.__get_x2__()
        self.mu_1_new = self.__get_mu1_new__()
        self.mu_2_new = self.__get_mu2_new__()
        # self.__fill_matr_left__()
        # self.__fill_vec_right__()

        self.left = np.zeros((self.N + 1, self.N + 1))
        # self.left[0][0] = - self.vec_a[1] / self.h - (self.beta * self.k_i(0) / self.alpha + 0.5 * self.h * self.q_i(0))
        # self.left[0][1] = self.vec_a[1] / self.h
        # self.left[self.N][self.N] = self.vec_a[self.N] / self.h - (self.delta * self.k_i(self.N) / self.gamma + 0.5 * self.h * self.q_i(self.N))
        # self.left[self.N][self.N - 1] = - self.vec_a[self.N] / self.h
        mu_1_new = self.h * self.__get_F_i__(0) * 0.5
        mu_2_new = self.h * self.__get_F_i__(self.N) * 0.5
        sigma_1_old = self.beta * self.k_i(0) / self.alpha
        sigma_2_old = self.delta * self.k_i(self.N) / self.gamma
        sigma_1_new = sigma_1_old * (1 + self.h * 0.5 * self.p_i(0) / self.k_i(0)) + self.h * 0.5 * self.q_i(0)
        sigma_2_new = sigma_2_old * (1 + self.h * 0.5 * self.p_i(self.N) / self.k_i(self.N)) + self.h * 0.5 * self.q_i(self.N)


        self.left[0][0] = sigma_1_new + self.k_x(self.y_i[0] + self.h * 0.5) / self.h
        self.left[0][1] = - self.k_x(self.y_i[0] + self.h * 0.5) / self.h
        self.left[self.N][self.N] = sigma_2_new + self.k_x(self.y_i[self.N] - self.h * 0.5) / self.h
        self.left[self.N][self.N - 1] = - self.k_x(self.y_i[self.N] - self.h * 0.5) / self.h
        for i in range(1, self.N):
            self.left[i][i - 1] = - self.vec_a[i] / (self.h ** 2) - self.p_x(self.y_i[i] - self.h * 0.5) * 0.5 / self.h
            self.left[i][i] = self.q_i(i) + (self.vec_a[i + 1] + self.vec_a[i]) / (self.h ** 2) + (self.p_x(self.y_i[i] - self.h * 0.5) - self.p_x(self.y_i[i] + self.h * 0.5)) * 0.5 / self.h
            self.left[i][i + 1] = - self.vec_a[i + 1] / (self.h ** 2) + self.p_x(self.y_i[i] + self.h * 0.5) * 0.5 / (self.h ** 2)
            # self.left[i][i - 1] = self.k_i(i) / (self.h ** 2)
            # self.left[i][i] = - self.k_i(i + 1) / (self.h ** 2) - self.k_i(i) / (self.h ** 2) + (self.p_i(i) + self.q_i(i)) / self.h
            # self.left[i][i + 1] = self.k_i(i + 1) / (self.h ** 2) - (self.p_i(i) + self.q_i(i)) / self.h

        self.right = np.zeros((self.N + 1, 1))
        self.right[0][0] = mu_1_new
        self.right[self.N][0] = mu_2_new
        for i in range(1, self. N):
            self.right[i][0] = self.__get_F_i__(i)

        # print(self.left)
        # print(self.left.shape)
        # print(self.right)
        # print(self.right.shape)
        res = np.linalg.solve(self.left, self.right)
        if debug:
            print(res)
            print(res.shape)

        x_arr = []
        u_arr = []
        y_arr = []

        for i in range(0, self.N + 1):
            if debug:
                print(f'i: {i}, x_i: {self.y_i[i]}, u_i: {self.solution_real_x(self.y_i[i])}, y_i: {res[i]}')
            x_arr.append(self.y_i[i])
            u_arr.append(self.solution_real_x(self.y_i[i]))
            y_arr.append(res[i][0])

        # fig = plt.figure()
        plt.plot(x_arr, u_arr, x_arr, y_arr)
        plt.plot(x_arr, u_arr, x_arr, y_arr)
        plt.title(f"Number of points = {self.N + 1}")
        plt.show()

        if debug:
            print('division:')
            tmp_left = self.vec_a[1] * (res[1] - res[0]) / self.h
            tmp_right = self.sigma1_new * res[0] - self.mu_1_middle
            print(f'left  : {tmp_left}   {tmp_right}   {tmp_right - tmp_left}')
            tmp_left = self.vec_a[self.N] * (res[self.N] - res[self.N - 1]) / self.h
            tmp_right = self.sigma2_new * res[self.N] - self.mu_2_middle
            print(f'right : {tmp_left}   {tmp_right}   {tmp_right - tmp_left}')


            if self.diff_type == 0:
                for i in range(1, self.N):
                    tmp_left = (self.vec_a[i+1] - self.vec_a[i]) * (self.y_i[i] - self.y_i[i-1]) / (self.h ** 2) + \
                               self.vec_a[i] * (self.y_i[i-1] - 2 * self.y_i[i] + self.y_i[i+1]) / (self.h ** 2) - \
                               self.y_i[i+1] / self.h + self.y_i[i] / self.h
                    tmp_right = self.get_f_x(self.y_i[i])
                    print(f' {i} -- {tmp_left}   --   {tmp_right}   ---   {tmp_left - tmp_right}')
            elif self.diff_type == 1:
                for i in range(1, self.N):
                    # tmp_left = (self.vec_a[i+1] - self.vec_a[i]) * (self.y_i[i+1] - self.y_i[i]) / (self.h ** 2) - (self.y_i[i+1] - self.y_i[i]) / self.h
                    # tmp_right = self.get_f_x(self.y_i[i])
                    tmp_left = res[i - 1] * self.__get_A_i__(i) - res[i] * self.__get_C_i__(i) + res[i + 1] * self.__get_B_i__(i)
                    tmp_right = - self.__get_F_i__(i)
                    print(f' {i}   {tmp_left}      {tmp_right}      {tmp_left - tmp_right}')
                print('=============================================================================================')
                for i in range(1, self.N):
                    tmp_left = self.vec_a[i] / (self.h ** 2) + 1 / self.h
                    tmp_right = self.__get_A_i__(i)
                    print(f' {i}   A   {tmp_left}      {tmp_right}      {tmp_left - tmp_right}')
                    tmp_left = self.vec_a[i + 1] / (self.h ** 2) + self.vec_a[i] / (self.h ** 2) + 1 / self.h
                    tmp_right = self.__get_C_i__(i)
                    print(f' {i}   C   {tmp_left}      {tmp_right}      {tmp_left - tmp_right}')
                    tmp_left = self.vec_a[i + 1] / (self.h ** 2)
                    tmp_right = self.__get_B_i__(i)
                    print(f' {i}   B   {tmp_left}      {tmp_right}      {tmp_left - tmp_right}')


                print('=============================================================================================')
                for i in range(1, self.N):
                    tmp_left = self.vec_a[i + 1] * (res[i + 1] - res[i]) / (self.h ** 2) - self.vec_a[i] * (
                                res[i] - res[i - 1]) / (self.h ** 2) - (res[i] - res[i - 1]) / self.h
                    tmp_right = - self.__get_F_i__(i)
                    print(f' {i}   {tmp_left}      {tmp_right}      {tmp_left - tmp_right}')

                print('=============================================================================================')
                for i in range(1, self.N):
                    tmp_left = (self.vec_a[i + 1] * (res[i + 1] - res[i]) / self.h - self.vec_a[i] * (res[i] - res[i - 1]) / self.h) / self.h - (res[i] - res[i - 1]) / self.h
                    tmp_right = - self.__get_F_i__(i)
                    print(f' {i}   {tmp_left}      {tmp_right}      {tmp_left - tmp_right}')

            elif self.diff_type == 3:
                for i in range(1, self.N):
                    tmp_left = res[i - 1] * self.__get_A_i__(i) - res[i] * self.__get_C_i__(i) + res[i + 1] * self.__get_B_i__(i)
                    tmp_right = - self.__get_F_i__(i)
                    print(f' {i}   {tmp_left}      {tmp_right}      {tmp_left - tmp_right}')
                print('=============================================================================================')
            for i in range(1, self.N):
                tmp_left = self.vec_a[i] / (self.h ** 2)
                tmp_right = self.__get_A_i__(i)
                # print(f' {i}   A   {tmp_left}      {tmp_right}      {tmp_left - tmp_right}')
                tmp_left = self.vec_a[i + 1] / (self.h ** 2) + self.vec_a[i] / (self.h ** 2) - self.p_i(i) / self.h + self.q_i(i)
                tmp_right = self.__get_C_i__(i)
                # print(f' {i}   C   {tmp_left}      {tmp_right}      {tmp_left - tmp_right}')
                tmp_left = self.vec_a[i + 1] / (self.h ** 2) - self.p_i(i) / self.h
                tmp_right = self.__get_B_i__(i)
                # print(f' {i}   B   {tmp_left}      {tmp_right}      {tmp_left - tmp_right}')
            print('=============================================================================================')
            for i in range(1, self.N):
                tmp_left = res[i - 1] * self.vec_a[i] / (self.h ** 2) + \
                           res[i] * (-self.vec_a[i + 1] / (self.h ** 2) - self.vec_a[i] / (self.h ** 2) - self.q_i(i)) + \
                           res[i + 1] * (self.vec_a[i + 1] / (self.h ** 2)) - \
                           self.p_i(i) * res[i + 1] / self.h + \
                           self.p_i(i) * res[i] / self.h
                tmp_right = - self.__get_F_i__(i)
                print(f' {i}   {tmp_left}      {tmp_right}      {tmp_left - tmp_right}')


        print(f"N = {self.N + 1}")
        print("i  |  x_i  |  u_i  |  y_i  |  y_i - u_i")
        for i in range(0, self.N + 1):
            u_i = self.solution_real_x(self.y_i[i])
            print(f'{i}  |  {self.y_i[i]:.8f}   |  {u_i:.8f}  |  {res[i][0]:.8f}  |  {res[i][0] - u_i:.8f}')























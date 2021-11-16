from Calc import Calc


# class Calc:
#     def d12(self):
#         print(1)

# a = 1
#     b = 5
#     b1 = 2
#     b2 = 1
#     b3 = 1
#     k1 = 2
#     k2 = 1
#     c1 = 1
#     c2 = 0
#     c3 = 1
#     p1 = 3
#     p2 = 0
#     d1 = 2
#     d2 = 1
#     d3 = 3
#     q1 = 2
#     q2 = 1
#     a1 = 1
#     a2 = 2
#     a3 = 3
#     a4 = 4
#     n1 = 7
#     n2 = 4
#     n3 = 2
#
#     def __init__(self):
#         self.alpha = self.a1 * self.a ** self.n1 \
#                      + self.a2 * self.a ** self.n2 \
#                      + self.a3 * self.a ** self.n3 \
#                      + self.a4
#
#         self.beta = self.a1 * self.n1 * self.a ** (self.n1 - 1) \
#                     + self.a2 * self.n2 * self.a ** (self.n2 - 1) \
#                     + self.a3 * self.n3 * self.a ** (self.n3 - 1)
#
#         self.gamma = self.a1 * self.b ** self.n1 \
#                      + self.a2 * self.b ** self.n2 \
#                      + self.a3 * self.b ** self.n3 \
#                      + self.a4
#
#         self.delta = -(self.a1 * self.n1 * self.b ** (self.n1 - 1)
#                      + self.a2 * self.n2 * self.b ** (self.n2 - 1)
#                      + self.a3 * self.n3 * self.b ** (self.n3 - 1))
#
#         self.k = self.get_k_str()
#         self.p = self.get_p_str()
#         self.q = self.get_q_str()
#
#         # return True
#
#         print('hi')
#
#     def get_k_str(self):
#         tmp = f'{self.b1} * x ^ {self.k1} + {self.b2} * x ^ {self.k2} + {self.b3}'
#         return tmp
#
#     def get_p_str(self):
#         tmp = f'{self.c1} * x ^ {self.p1} + {self.c2} * x ^ {self.p2} + {self.c3}'
#         return tmp
#
#     def get_q_str(self):
#         tmp = f'{self.d1} * x ^ {self.q1} + {self.d2} * x ^ {self.q2} + {self.d3}'
#         return tmp




if __name__ == '__main__':
    print('start')
    tmp = Calc()
    tmp.set_phi_type(1)

    # print(tmp.phi_i_x(1, 3))
    # tmp.push_colloc_pnt(1.1)
    # tmp.push_colloc_pnt(2)
    # tmp.push_colloc_pnt(3)
    # tmp.push_colloc_pnt(4)
    # tmp.push_colloc_pnt(4.2)
    # tmp.push_colloc_pnt(4.5)
    # tmp.push_colloc_pnt(4.7)
    # tmp.push_colloc_pnt(4.9)

    for i in range (11, 49):
        tmp.push_colloc_pnt(float(i) / 10)

    tmp.solve_colloc()
    tmp.print_discrepancy_table()
    tmp.print_comparation_plot()

    print('end')


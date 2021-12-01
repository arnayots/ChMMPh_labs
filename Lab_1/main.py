from Calc import Calc
# import aprsolv

if __name__ == '__main__':
    print('start')
    tmp = Calc()
    tmp.set_phi_type(2)

    # print(tmp.phi_i_x(1, 3))
    tmp.push_colloc_pnt(1.1)
    # tmp.push_colloc_pnt(1.3)
    tmp.push_colloc_pnt(1.6)
    # tmp.push_colloc_pnt(1.9)
    tmp.push_colloc_pnt(2.2)
    # tmp.push_colloc_pnt(2.5)
    tmp.push_colloc_pnt(2.8)
    # tmp.push_colloc_pnt(3.1)
    tmp.push_colloc_pnt(3.4)
    # tmp.push_colloc_pnt(3.7)
    tmp.push_colloc_pnt(4.0)
    # tmp.push_colloc_pnt(4.3)
    tmp.push_colloc_pnt(4.6)
    # tmp.push_colloc_pnt(4.7)
    # tmp.push_colloc_pnt(4.9)

    print(tmp.colloc_pnts)

    # for i in range (11, 49, 1):
    #     tmp.push_colloc_pnt(float(i) / 10)

    # tmp.solve_colloc()
    # tmp.solve_ritz(10)

    # tmp2 = aprsolv.ApproxSolver()
    # tmp2.Ritz(4)

    print('end')

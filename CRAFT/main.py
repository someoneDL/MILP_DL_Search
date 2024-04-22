
from copy import deepcopy
from asyncore import read
import re
from gurobipy import *
from gurobipy import GRB
import math

def var_x(x, r, z):
    return [[x + "_" + str(r) + "_" + str(i) + "_" + str(j) for j in range(2)] for i in range(z)]

def var_y(y, r, z):
    return [[y + "_" + str(r) + "_" + str(i) + "_" + str(j) for j in range(4)] for i in range(z)]

def var_p(r):
    return [["p_" + str(r) + "_" + str(i) + "_" + str(j) for j in range(2)] for i in range(16)]

def flat(state, num):
    a = []
    for i in range(num):
        a += state[i]
    return a

def Model_I(s_in, f,z):
    MLI_coe = []
    file_name = 'min_model_in_factor.txt'
    fl = open(file_name, 'r')
    k = fl.readlines()
    for a in k:
        a = a.strip()
        a = a.strip(',')
        a = a.strip('[')
        a = a.strip(']')
        a = a.split(',')
        MLI_coe.append(a)
    fl.close()
    for j in range(len(MLI_coe)):
        MLI_coe[j] = list(map(lambda x: int(x), MLI_coe[j]))
    for i in range(z):
        for j in range(len(MLI_coe)):
            for k in range(0, 2):
                f.write(str(" %+d" % int(str(MLI_coe[j][k]))) + " " + s_in['x'][i][k] + " ")
            for k in range(2, 6):
                f.write(str(" %+d" % int(str(MLI_coe[j][k]))) + " " + s_in['y'][i][k - 2] + " ")
            f.write(">= " + str(-MLI_coe[j][-1]) + "\n")
    return 0

def Model_s(s_in, s_out, p, f):
    MLS_s_coe = []
    MLS_s_coe2 = []
    file_name = 'min_model_s_factor1.txt'
    fl = open(file_name, 'r')
    k = fl.readlines()
    for a in k:
        a = a.strip()
        a = a.strip(',')
        a = a.strip('[')
        a = a.strip(']')
        a = a.split(',')
        MLS_s_coe.append(a)
    fl.close()
    for j in range(len(MLS_s_coe)):
        MLS_s_coe[j] = list(map(lambda x: int(x), MLS_s_coe[j]))

    file_name = 'min_model_s_factor2.txt'
    fl = open(file_name, 'r')
    k = fl.readlines()
    for a in k:
        a = a.strip()
        a = a.strip(',')
        a = a.strip('[')
        a = a.strip(']')
        a = a.split(',')
        MLS_s_coe2.append(a)
    fl.close()
    for j in range(len(MLS_s_coe2)):
        MLS_s_coe2[j] = list(map(lambda x: int(x), MLS_s_coe2[j]))

    for i in range(16):
        for j in range(len(MLS_s_coe)):
            for k in range(0, 2):
                f.write(str(" %+d" % int(str(MLS_s_coe[j][k]))) + " " + s_in['x'][i][k] + " ")
            for k in range(2, 4):
                f.write(str(" %+d" % int(str(MLS_s_coe[j][k]))) + " " + s_out['x'][i][k - 2] + " ")
            f.write(">= " + str(- MLS_s_coe[j][-1]) + "\n")
    for i in range(16):
        for j in range(len(MLS_s_coe2)):
            for k in range(0, 4):
                f.write(str(" %+d" % int(str(MLS_s_coe2[j][k]))) + " " + s_in['y'][i][k] + " ")
            for k in range(4, 8):
                f.write(str(" %+d" % int(str(MLS_s_coe2[j][k]))) + " " + s_out['y'][i][k - 4] + " ")
            for k in range(8,10):
                f.write(str(" %+d" % int(str(MLS_s_coe2[j][k]))) + " " + p[i][k - 8] + " ")
            f.write(">= " + str(- MLS_s_coe2[j][-1]) + "\n")
    return 0

def Model_XOR(s_in1, s_in2, s_out, XOR_coe, f):
    for i in range(4):
        for j in range(len(XOR_coe)):
            f.write(str(" %+d" % int(str(XOR_coe[j][0]))) + " " + s_in1['x'][0] + " ")
            f.write(str(" %+d" % int(str(XOR_coe[j][1]))) + " " + s_in2['x'][0] + " ")
            f.write(str(" %+d" % int(str(XOR_coe[j][2]))) + " " + s_out['x'][0] + " ")
            f.write(str(" %+d" % int(str(XOR_coe[j][3]))) + " " + s_in1['y'][i] + " ")
            f.write(str(" %+d" % int(str(XOR_coe[j][4]))) + " " + s_in2['y'][i] + " ")
            f.write(str(" %+d" % int(str(XOR_coe[j][5]))) + " " + s_out['y'][i] + " ")
            f.write(">= " + str(-XOR_coe[j][-1]) + "\n")
    # f.write("\n")
    return 0

def Model_MDS(s_in, s_mi, s_out, f):
    xor_coe = []
    file_name = 'min_model_xor_factor.txt'
    fl = open(file_name, 'r')
    k = fl.readlines()
    for a in k:
        a = a.strip()
        a = a.strip(',')
        a = a.strip('[')
        a = a.strip(']')
        a = a.split(',')
        xor_coe.append(a)
    fl.close()
    for j in range(len(xor_coe)):
        xor_coe[j] = list(map(lambda x: int(x), xor_coe[j]))
    for i in range(4):
        # (1,0,1,1)
        c_in0 = {'x': s_in['x'][i], 'y': s_in['y'][i]}
        c_in1 = {'x': s_in['x'][8 + i], 'y': s_in['y'][8 + i]}
        c_in2 = {'x': s_in['x'][12 + i], 'y': s_in['y'][12 + i]}
        c_m = {'x': s_mi['x'][i], 'y': s_mi['y'][i]}
        c_out = {'x': s_out['x'][i], 'y': s_out['y'][i]}
        Model_XOR(c_in0, c_in1, c_m, xor_coe, f)
        Model_XOR(c_m, c_in2, c_out, xor_coe, f)
        #(0,1,0,1)
        c_in0 = {'x': s_in['x'][4 + i], 'y': s_in['y'][4 + i]}
        c_in1 = {'x': s_in['x'][12 + i], 'y': s_in['y'][12 + i]}
        c_out = {'x': s_out['x'][4 + i], 'y': s_out['y'][4 + i]}
        Model_XOR(c_in0, c_in1, c_out, xor_coe, f)
        # (0,0,1,0)
        c_in = {'x': s_in['x'][8+i], 'y': s_in['y'][8+i]}
        c_out = {'x': s_out['x'][8 + i], 'y': s_out['y'][8 + i]}
        for k in range(2):
            f.write(c_in['x'][k] + " - " + c_out['x'][k] + " = 0\n")
        for k in range(2, 6):
            f.write(c_in['y'][k - 2] + " - " + c_out['y'][k - 2] + " = 0\n")
        #(0,0,0,1)
        c_in = {'x': s_in['x'][12 + i], 'y': s_in['y'][12+i]}
        c_out = {'x': s_out['x'][12 + i], 'y': s_out['y'][12 + i]}
        for k in range(2):
            f.write(c_in['x'][k] + " - " + c_out['x'][k] + " = 0\n")
        for k in range(2, 6):
            f.write(c_in['y'][k - 2] + " - " + c_out['y'][k - 2] + " = 0\n")

    return 0

def Model_P(s_out, s_n, f):
    Ph = [15, 12, 13, 14, 10, 9, 8, 11, 6, 5, 4, 7, 1, 2, 3, 0]
    for i in range(16):
        for k in range(2):
            f.write(s_out['x'][i][k] + " - " + s_n['x'][Ph[i]][k] + " = 0\n")
        for k in range(0,4):
            f.write(s_out['y'][i][k] + " - " + s_n['y'][Ph[i]][k] + " = 0\n")
    return 0


def product_lp(var_set, name, R, z_in=[], z_out=16, d_in=[], c_out=[], flag=0,AP_a=[],AP_b=[]):
    f = open(name, "w+")
    P = []
    for r in range(R):
        P.append(var_p(r))
        var_set.extend(flat(P[r], 16))
    obj = " + "
    objc = ""
    for ri in range(R):
        for i in range(16):
            objc += P[ri][i][0] + obj + "2 " + P[ri][i][1] + obj
    objc = objc[:-2]
    f.write("Minimize\n" + objc + "\n")
    f.write("Subject To\n")
    if flag != 0:
        f.write(objc + " = " + str(flag) + "\n")

    s_in = {'x': var_x('x', 0, 16), 'y': var_y('y', 0, 16)}
    Model_I(s_in, f, 16)
    var_set.extend(flat(s_in['x'], 16))
    var_set.extend(flat(s_in['y'], 16))

    obj = ''
    for i in range(16):
        obj += str(s_in['x'][i][0]) + " + "
    obj = obj[:-3]
    f.write(obj + " = 0\n")
    obj = ''
    for i in range(16):
        obj += str(s_in['x'][i][1]) + " + "
    obj = obj[:-3]
    f.write(obj + " >= 1\n")
    if z_in != []:
        for z in z_in:
            f.write(s_in['x'][z][1] + " = 1\n")
    if d_in != []:
        for z in z_in:
            for i in range(4):
                f.write(s_in['y'][z][i] + " = " + str(d_in[i]) + "\n")

    for r in range(R):
        s_m = {'x': var_x('xm', r, 16), 'y': var_y('ym', r, 16)}
        Model_I(s_m, f, 16)
        var_set.extend(flat(s_m['x'], 16))
        var_set.extend(flat(s_m['y'], 16))

        s_n = {'x': var_x('xn', r, 16), 'y': var_y('yn', r, 16)}
        Model_I(s_n, f, 16)
        var_set.extend(flat(s_n['x'], 16))
        var_set.extend(flat(s_n['y'], 16))

        s_mi = {'x': var_x('xmi', r, 4), 'y': var_y('ymi', r, 4)}
        Model_I(s_mi, f, 4)
        var_set.extend(flat(s_mi['x'], 4))
        var_set.extend(flat(s_mi['y'], 4))

        s_out = {'x': var_x('x', r + 1, 16), 'y': var_y('y', r + 1, 16)}
        Model_I(s_out, f, 16)
        var_set.extend(flat(s_out['x'], 16))
        var_set.extend(flat(s_out['y'], 16))
        f.write("\n")
        Model_s(s_in, s_m, P[r], f)
        Model_MDS(s_m, s_mi, s_n, f)
        Model_P(s_out,s_n,f)
        # f.write("\n")
        if r + 1 != R:
            s_in = s_out
    obj = ''
    for i in range(16):
        obj += str(s_out['x'][i][1]) + " + "
    obj = obj[:-3]
    f.write(obj + " >= 1\n")
    if z_out != 16:
        f.write(s_out['x'][z_out][0] + " = 0\n")
        f.write(s_out['x'][z_out][1] + " = 1\n")
    if c_out != []:
        for i in range(4):
            f.write(s_out['y'][z_out][i] + " = " + str(c_out[i]) + "\n")
    if AP_a != []:
        liss = ''
        for a in range(len(AP_a)):
            liss += AP_a[a] + " + "
        liss = liss[:-3]
        f.write(liss + " = " + str(len(AP_a)) + '\n')
        obj = ''
        for r in range(R):
            for i in range(16):
                obj += 'x_' + str(r) + '_' + str(i) + '_1' + " + "
        obj = obj[:-3]
        f.write(obj + " = " + str(len(AP_a)) + "\n")
    if AP_b != []:
        liss = ''
        for b in range(len(AP_b)):
            liss += AP_b[b] + " + "
        liss = liss[:-3]
        f.write(liss + " = " + str(len(AP_b)) + '\n')
        obj = ''
        for r in range(R):
            for i in range(16):
                obj += 'xm_' + str(r) + '_' + str(i) + '_1' + " + "
        obj = obj[:-3]
        f.write(obj + " = " + str(len(AP_b)) + "\n")

    f.write("Binary\n")
    for v in var_set:
        f.write(v + "\n")
    f.write('End' + '\n')
    f.close()
    return 0

def Research_first(name, R):
    '''Identifying the optimal differential pattern trail'''
    z_in = []
    z_out = []
    var_set = []
    product_lp(var_set, name, R)
    m = read(name)
    m.Params.Method = 3
    m.Params.PoolSearchMode = 2
    m.optimize()
    if m.status == 2 or m.status == 11:
        m.write("result_look_first.sol")
    AP_a = []
    AP_b = []
    for r in range(R):
        for i in range(16):
            nx = 'x_' + str(r) + '_' + str(i) + '_1'
            nxm = 'xm_' + str(r) + '_' + str(i) + '_1'
            if abs(m.getVarByName(name=nx).X - 1) < 0.1:
                AP_a.append(nx)
            if abs(m.getVarByName(name=nxm).X - 1) < 0.1:
                AP_b.append(nxm)
    print('-------------print trails-----\n', AP_a, AP_b)

    aim = m.ObjVal
    for i in range(16):
        nx = 'x_0_' + str(i) + '_1'
        if abs(m.getVarByName(name=nx).X - 1) < 0.1:
            z_in.append(i)
    for i in range(16):
        nx = 'x_' + str(R) + '_' + str(i) + '_1'
        if abs(m.getVarByName(name=nx).X - 1) < 0.1:
            W=i
            break
    return aim, z_in, W, AP_a, AP_b

def Research_second(name, R):
    '''Selecting the concrete input difference'''
    aim, z_in, W, AP_a,AP_b=Research_first(name, R)
    print(aim,z_in,W, AP_a,AP_b)
    print('--------', aim, z_in, W, AP_a,AP_b, '--------')
    F_c=[0 for i in range(16)]
    for c in range(16):
        cl = [((c >> i) & 1) for i in range(4)][::-1]
        var_set = []
        product_lp(var_set, name, R,z_in = z_in,z_out = W, d_in= cl ,flag=aim,AP_a=AP_a,AP_b=AP_b)
        m = read(name)
        m.Params.Method = 3
        m.Params.PoolSearchMode = 2
        m.Params.PoolSolutions = 2**20
        m.optimize()
        if m.status == 2:
            nSolutions = m.SolCount
            F_c[c]=nSolutions
    print(F_c)
    max_c=F_c.index(max(F_c))
    d_in=[((max_c >> i) & 1) for i in range(4)][::-1]
    print('----Completed the comparison calculation of the input difference and found the best input difference value----',d_in)
    return aim, z_in, W, d_in, AP_a, AP_b

def operation(name, R, z_in, W, d_in, c_out, aim, AP_a, AP_b):
    '''Enhancing truncated differentials'''
    pr = 0
    var_set = []

    for i in range(aim, aim + 20):
        print(i, aim)
        product_lp(var_set, name, R, z_in, W, d_in, c_out, flag=i, AP_a=AP_a, AP_b=AP_b)
        m = read(name)
        m.Params.Method = 3
        m.Params.PoolSearchMode = 2
        # m.Params.IntegralityFocus = 0
        m.Params.PoolSolutions = 781250
        m.Params.OutputFlag = 0
        m.optimize()
        print(m.status)
        print('m.status',m.status)
        if m.status == 2:
            nSolutions = m.SolCount
            print(nSolutions)
            pr = pr + nSolutions * (2 ** (-i))

    return pr


def M_Truncated_differntial(name, R):
    '''Generating the DL distinguihser'''
    aim, z_in, W, d_in, AP_a,AP_b = Research_second(name, R)
    print(aim, z_in, W, d_in)
    print("pattern trail")
    print("AP_a:  ", AP_a)
    print("AP_b:  ", AP_b)
    aim = math.ceil(aim)
    Bia = 0  # Bia =  pr(c->c0)*dlct0 + pr(c->c1)*dlct1 +  pr(c->c2)*dlct2 + ...
    MASK = {'1': [8, 4, 0, 0, 0, 0, -4, -4, 0, 0, -4, -4, 4, 0, 0, 0],
            '2': [8, 0, 8, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -8, 0, -8],
            '3': [8, -4, 0, 0, 0, 0, 4, -4, 0, 0, -4, 4, -4, 0, 0, 0],
            '4': [8, 0, 0, -4, 4, 0, 0, -4, 0, 4, -4, 0, 0, 0, -4, 0],
            '5': [8, 0, -8, 0, 0, 0, 0, 0, -8, 0, 8, 0, 0, 0, 0, 0],
            '6': [8, 0, 0, 4, -4, 0, 0, -4, 0, -4, -4, 0, 0, 0, 4, 0],
            '7': [8, 0, -8, 0, 0, -8, 0, 8, 0, 0, 0, 0, 0, 0, 0, 0],
            '8': [8, 4, 0, 0, 4, 0, 0, 0, 0, 0, -4, -4, 0, 0, -4, -4],
            '9': [8, 0, 0, 4, 0, -4, -4, 0, 0, 0, 0, 4, 0, -4, -4, 0],
            '10': [8, 0, 0, -4, 0, 0, -4, 0, 0, -4, 4, 0, -4, 0, 0, 4],
            '11': [8, -4, 0, 0, -4, 4, 0, 0, 0, -4, 0, 0, -4, 4, 0, 0],
            '12': [8, 0, 0, -4, 0, -4, 4, 0, 0, 0, 0, -4, 0, -4, 4, 0],
            '13': [8, 0, 0, -4, 0, 0, -4, 0, 0, 4, -4, 0, 4, 0, 0, -4],
            '14': [8, -4, 0, 0, -4, 4, 0, 0, 0, -4, 0, 0, -4, 4, 0, 0],
            '15':[8, -4, 0, 0, -4, 0, 0, 0, 0, 0, 4, -4, 0, 0, -4, 4]}
    PR = []
    Prm = [[] for i in range(16)]
    BIA = []
    Mask = 0
    k = 0
    F = 0
    print('-------------Start calculating truncated differentials------------')
    for c in range(16):
        d_out = [((c >> i) & 1) for i in range(4)][::-1]
        pr = operation(name, R, z_in, W, d_in, d_out, aim, AP_a, AP_b)
        lis = [hex(c), pr]
        PR.append(lis)
    print(PR)
    print('---------Start selecting output masks---------')

    for key in MASK.keys():
        bia = 0
        for c in range(16):
            if MASK[key][c] != 0:
                bia += PR[c][1] * (MASK[key][c] / 16)
                Prm[k].append(PR[c])
        BIA.append([key, bia])
        if abs(bia) > abs(Bia):
            Mask = key
            Bia = bia
            F = k
        k += 1

    print("objective function value: ", aim)
    print("the position of input active words: ", z_in)
    print("the position of output active words:", W)
    print("input difference: ", d_in)
    print("output mask: ", Mask)
    print("The available truncation differential corresponding to the current mask, classified by the output difference:",Prm[F])
    print("Bia: ", Bia, math.log2(abs(Bia)))
    print("The bias corresponding to different masks:", BIA)
    return 0


if __name__ == "__main__":
    R = 7
    name = "CRAFT64_td_" + str(R) + ".lp"
    M_Truncated_differntial(name,R)
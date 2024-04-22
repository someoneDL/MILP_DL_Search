def W(a):
    #Hamming Weight
    w=0
    for i in range(5):
        w = w + (a>>i) & 1
    return w
def branch(S):
    flag=10
    n=2**5
    for i in range(n):
        c=W(i)+W(S[i])
        if c<flag:
            flag=c

def Inner_product(a,b):
    c=0
    for i in range(5):
        c ^= ((a>>i)&1)&((b>>i)&1)
    return c

def Pr_Lat(S):
    n=2**4
    LAT=[[-2**3 for i in range(n)]for j in range(n)]
    for i in range(n):
        for j in range(n):
            for k in range(n):
                if Inner_product(i,k)==Inner_product(j,S[k]):
                    LAT[i][j] += 1
    return LAT
def Pr_DDT(S):
    n=2**4
    DDT=[[0 for i in range(n)]for j in range(n)]
    for i in range(n):
        for j in range(n):
            DDT[i][S[j]^S[j^i]]+=1
    return DDT
def Pr_DLCT(S):
    n=2**4
    DLCT=[[-2**3 for i in range(n)]for j in range(n)]
    for i in range(n):
        for j in range(n):
            for k in range(n):
                if Inner_product((S[k]^S[k^i]),j)==0:
                    DLCT[i][j] += 1
    return DLCT
def D_characteristic(DDT):
    n=2**4
    D_C=[]
    for i in range(n):
        for j in range(n):
            # list=[i >>d & 1 for d in range(4)][::-1]
            # list1=[j >>d & 1 for d in range(4)][::-1]
            # list.extend(list1)
            list=[i,j]
            if DDT[i][j]==0:
                continue
            elif DDT[i][j]==2:
                # list.extend([1,1])
                list.append(3)
            elif DDT[i][j]==4:
                #概率 2^-2
                # list.extend([1,0])
                list.append(2)
            else:
                # list.extend([0,0])
                list.append(0)
            D_C.append(list)
    return D_C
def L_characteristic(LAT):
    n=2**4
    L_C=[]
    for i in range(n):
        for j in range(n):
            # list=[i,j]
            list=[i >>d & 1 for d in range(4)][::-1]
            list1=[j >>d & 1 for d in range(4)][::-1]
            list.extend(list1)
            if LAT[i][j]==0:
                continue  
            elif LAT[i][j]==4 or LAT[i][j]==-4:
                #2p0+2p1
                list.extend([1,1])
                # list.append(3)
            elif LAT[i][j]==8 or LAT[i][j]==-8:
                list.extend([1,0])
                # list.append(2)
            else:#16
                list.extend([0,0])
                # list.append(0)
            L_C.append(list)
    return L_C
def DL_characteristic(DLCT):
    n=2**4
    DL_C=[]
    for i in range(1,n):
        for j in range(1,n):
            list=[]
            list=[i >>d & 1 for d in range(4)][::-1]
            list1=[j >>d & 1 for d in range(4)][::-1]
            list.extend(list1)
            if DLCT[i][j]==0:
                continue
            elif DLCT[i][j]==16 or DLCT[i][j]==-16:
                list.extend([0])
            DL_C.append(list)
    return DL_C


if __name__=="__main__":
    # S=[4,11,31,20,26,21,9,2,27,5,8,18,29,3,6,28,30,19,7,14,0,13,17,24,16,12,1,25,22,10,15,23] #ascon
    # S=[0xc,0xa,0xd,0x3,0xe,0xb,0xf,0x7,0x8,0x9,0x1,0x5,0x0,0x2,0x4,0x6] #Midori64
    # S=[0xc,0x6,0x9,0x0,0x1,0xa,0x2,0xb,0x3,0x8,0x5,0xd,0x4,0xe,0x7,0xf] #skinny64-64
    S=[0xc,0xa,0xd,0x3,0xe,0xb,0xf,0x7,0x8,0x9,0x1,0x5,0x0,0x2,0x4,0x6] #CRAFT64-128
    DDT=Pr_DDT(S)
    LAT=Pr_Lat(S)
    DLCT=Pr_DLCT(S)
    print("DDT")
    for i in range(16):
        print(DDT[i])
    # print("LAT")
    # for i in range(16):
    #     print(LAT[i])
    print("DLCT")
    for i in range(16):
        print(DLCT[i])
    # D_c=D_characteristic(DDT)
    # for i in range(len(D_c)):
    #     print(D_c[i])
    # L_c=L_characteristic(LAT)
    # for i in range(len(L_c)):
    #     print(L_c[i])
    # DL_c=DL_characteristic(DLCT)
    # for i in range(len(DL_c)):
    #     print(DL_c[i])
    f=open('tmp.txt','w')
    print("DLCT")
    for i in range(16):
        for j in range(16):
            f.write(str(DLCT[j][i])+ ", ")
        f.write("\n")
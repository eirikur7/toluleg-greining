import numpy as np
from numpy import linalg as LA
import matplotlib.pyplot as plt
from mpl_toolkits import mplot3d
import time

def createAMatrix(n, m, LSpan, Lx, Ly, H, K, delta):
    A = np.zeros((n*m, m*n))
    hx = Lx/(m - 1)
    hy = Ly/(n -1)
    testDic = {}
    for j in range(n):
        for i in range(m):
            eqNr = i + j*m
            if((j==(n-1)) and (i != 0) and (i != (m-1))):     #Top
                A[eqNr, eqNr]       = ((2*hy*H)/K) - 3
                A[eqNr, eqNr - m]   = 4
                A[eqNr, eqNr - 2*m] = -1
            elif(( j==0 ) and (i != 0) and (i != (m-1))):  #Bottom
                A[eqNr, eqNr]       = ((2*hy*H)/K) - 3
                A[eqNr, eqNr + m]   = 4
                A[eqNr, eqNr + 2*m] = -1
            elif(( i==(m-1) )):                                 #Right
                A[eqNr, eqNr] = ((-2*hx*H)/K) + 3
                A[eqNr, eqNr - 1] = -4
                A[eqNr, eqNr - 2] = 1
            elif((i==0) and ((j*hy < LSpan[0]) or (j*hy > LSpan[1]))):             #Left # (Lx-L - (hx*i))>0 
                A[eqNr, eqNr] = ((2*hx*H)/K) - 3
                A[eqNr, eqNr + 1] = 4
                A[eqNr, eqNr + 2] = -1
            elif(i==0):                                         #Power
                A[eqNr, eqNr] = 3
                A[eqNr, eqNr + 1] = -4
                A[eqNr, eqNr + 2] = 1
            else:                                               #Base Plate
                A[eqNr, eqNr] = -2*((H/(K*delta)) + (1/(hx**2)) + (1/(hy**2)))#-(2*(hy**2)) - (2*(hx**2)) - ( ((hx**2)*(hy**2)*2*H)/(K*delta) )
                A[eqNr, eqNr-1] = 1/(hx**2)
                A[eqNr, eqNr+1] = 1/(hx**2)
                A[eqNr, eqNr+m] = 1/(hy**2)
                A[eqNr, eqNr-m] = 1/(hy**2)
    return A

def createBMatrix(n, m, L,Lx, P, delta, K): #Needs to be fixed for the future
    B = np.zeros((m*n, 1))
    # for i in range(m*n):
    #     B[i, 0] = 20
    hx = Lx/(m-1)
    for j in range(n):
        B[j*m] = (2*hx*P)/(L*delta*K)

    return B

def prob3():
    Lx = 2
    Ly = 2
    delta = 0.1
    P = 5
    L = 2
    LSpan = [0, 2]
    K = 1.68
    H = 0.005
    m = 10
    n = 10
    A = createAMatrix(n, m, LSpan, Lx, Ly, H, K, delta)
    B = createBMatrix(n, m, L, Lx, P, delta, K)
    V = LA.solve(A, B) + 20
    print(V[0])

    fig = plt.figure()
    ax = plt.axes(projection='3d')
    x = np.linspace(0, Lx, m)
    y = np.linspace(0, Ly, n)
    X1, Y1 = np.meshgrid(x, y)
    ax.plot_surface(X1, Y1, V.reshape((n,m)), cmap='coolwarm', edgecolor='none')
    ax.set_xlabel(xlabel='x[cm]')
    ax.set_ylabel(ylabel='y[cm]')
    ax.set_zlabel(zlabel='temperature[°C]')
    plt.show()

    # fig = plt.figure()
    # ax = plt.axes(projection='3d')
    # x = np.linspace(0, Lx, m)
    # xAll = np.zeros((n*m, 1))
    # for i in range(n):
    #     xAll[i*m:(m*i+m), 0] = x

    # print(xAll)
    # y = np.linspace(0, Ly, n)
    # X1, Y1 = np.meshgrid(x, y)
    # hy = 1/(n-1)
    # yAll = np.zeros((n*m, 1))
    # for i in range(m):
    #     yAll[i*n:(n*i+n), 0] = [hy*i]*n
    # print(xAll.shape)
    # print(yAll.shape)
    # print(V.shape)
    # ax.scatter3D(xAll, yAll, V)
    # ax.plot_surface(xAll, yAll, V)
    # plt.show()

    # for j in range(n):
    #         for i in range(m):
    #             eq = i + j*m
    #             eq2 = eq + m*(m-j-1)
    #             for ii in range(m*n):
    #                 if A[eq2,ii] >= 0:
    #                     print(" ", end="")
    #                 print(f"{A[eq2,ii]:.0f}", end=" ")
    #             print()
    # for key in tDic:
    #     x,y = tDic[key]
    #     eqNr = x + m*y
    #     print("-------------------------{}-------------------".format(key))
    #     print("x={},y={}".format(x,y))
    #     for j in range(n):
    #         for i in range(m):
    #             print(A[eqNr, i + j*m],end="| ")
    #         print()
    #     print('----------------------------------------------')


    # print(A.shape)
    # print(B)
    # for j in range(m*n):
    #     tempStr = ""
    #     for i in range(m*n):
    #         tempStr += str(A[j, i]) + " "
    #     print("{}. {}".format(j, tempStr))
    #     print()


def prob4():
    Lx = 2
    Ly = 2
    delta = 0.1
    P = 5
    L = 2
    LSpan = [0, 2]
    K = 1.68
    H = 0.005
    compareIndex = 0
    totalValues = 9
    nmMatrix = np.zeros((totalValues, totalValues))
    V_ref = LA.solve(createAMatrix(100, 100, LSpan, Lx, Ly, H, K, delta), createBMatrix(100, 100, L, Lx, P, delta, K))
    V_ref_compare = V_ref[0]
    # nmMatrix[10,10] = V100[compareIndex]

    for i in range(totalValues):
        for ii in range(totalValues):
            time_to_compare = time.time()
            n = 10*(i + 1)
            m = 10*(ii + 1)
            V = LA.solve(createAMatrix(n, m, LSpan, Lx, Ly, H, K, delta), createBMatrix(n, m, L, Lx, P, delta, K))
            nmMatrix[i, ii] = V[compareIndex] - V_ref_compare
            if nmMatrix[i, ii] < 0.01 and time.time()-time_to_compare <= 0.5:
                print("m={}, n={}: V[0]={:.2f}, diff={:.2f},time={:.2f}".format(m, n, V[compareIndex, 0], nmMatrix[i, ii], time.time()-time_to_compare))
    
    # error_array = np.zeros((totalValues**2,1))
    # for i in range(totalValues):
    #     for ii in range(totalValues):
    #         eq = ii + (i*totalValues)
    #         error_array[eq,0] = nmMatrix[i, ii] - V_ref_compare


    fig = plt.figure()
    ax = plt.axes(projection='3d')
    allM = np.linspace(0, totalValues*10, totalValues)
    allN = np.linspace(0, totalValues*10, totalValues)
    M1, N1 = np.meshgrid(allM, allN)
    print(M1.shape)
    print(N1.shape)
    ax.plot_surface(M1, N1, nmMatrix, cmap='coolwarm', edgecolor='none')
    ax.set_xlabel(xlabel='M')
    ax.set_ylabel(ylabel='N')
    ax.set_zlabel(zlabel='Deviation in (0,0)[°C]')
    plt.show()


def prob5():
    Lx = 4
    Ly = 4
    L = 2
    LSpan = [0, 2]
    delta = 0.1
    P = 5
    K = 1.68
    H = 0.005
    n = 10
    m = 10
    A = createAMatrix(n, m, LSpan, Lx, Ly, H, K, delta)
    B = createBMatrix(n, m, L, Lx, P, delta, K)
    V = LA.solve(A, B) + 20
    print(V[0])

    fig = plt.figure()
    ax = plt.axes(projection='3d')
    x = np.linspace(0, Lx, m)
    y = np.linspace(0, Ly, n)
    X1, Y1 = np.meshgrid(x, y)
    ax.plot_surface(X1, Y1, V.reshape((n,m)), cmap='coolwarm', edgecolor='none')
    ax.set_xlabel(xlabel='x[cm]')
    ax.set_ylabel(ylabel='y[cm]')
    ax.set_zlabel(zlabel='temperature[°C]')
    plt.show()

     
    
def prob6():
    Lx = 4
    Ly = 4
    L = 2
    delta = 0.1
    P = 5
    K = 1.68
    H = 0.005
    n = 10
    m = 10
    LSpansTries = 20
    LSpansEnd = np.linspace(L, Ly, LSpansTries)
    for i in range(LSpansTries):
        LSpan = [LSpansEnd[i]-L, LSpansEnd[i]]
        print(LSpan)
        A = createAMatrix(n, m, LSpan, Lx, Ly, H, K, delta)
        B = createBMatrix(n, m, L, Lx, P, delta, K)
        V = LA.solve(A, B) + 20
        fig = plt.figure()
        ax = plt.axes(projection='3d')
        x = np.linspace(0, Lx, m)
        y = np.linspace(0, Ly, n)
        X1, Y1 = np.meshgrid(x, y)
        ax.plot_surface(X1, Y1, V.reshape((n,m)), cmap='coolwarm', edgecolor='none')
        ax.set_xlabel(xlabel='x[cm]')
        ax.set_ylabel(ylabel='y[cm]')
        ax.set_zlabel(zlabel='temperature[°C]')
        ax.set_title("Span:[{:.3f}, {:.3f}]".format(LSpan[0], LSpan[1]))
        fig.savefig("hrafnljoti/proj3Fig/prob6_[{:.3f}, {:.3f}].png".format(i+1, LSpan[0], LSpan[1]))
        # plt.show()




if __name__ == "__main__":
    # prob3()
    # prob4()
    # prob5()
    prob6()




# def createAMatrix2(n, m, L, Lx, Ly, H, K, delta):
#     A = np.zeros((n*m+1, m*m+1))
#     hx = Lx/(m- 1)
#     hy = Ly/(n -1 )
#     for j in range(1, n+1):
#         for i in range(1, m+1):
#             if((j==1) and (i != 1) and (i != m)):     #Top
#                 A[i+(j-1)*m, i]       = ((2*hy*H)/K) + 3
#                 A[i+(j-1)*m, i + m]   = -4
#                 A[i+(j-1)*m, i + 2*m] = 1
#             elif(( j==n ) and (i != 1) and (i != m)):  #Bottom
#                 A[i+(j-1)*m, i]       = ((2*hy*H)/K) - 3
#                 A[i+(j-1)*m, i + m]   = 4
#                 A[i+(j-1)*m, i + 2*m] = -1
#             elif(( i==m )):                            #Right
#                 A[i+(j-1)*m, m + (j-1)*m] = ((2*hy*H)/K) + 3
#                 if((m + 0 + (j-1)*m) < (m*n -1)):
#                     A[i+(j-1)*m, m + 1 + (j-1)*m] = -4
#                 if((m + 1 + (j-1)*m) < (m*n -1)):
#                     A[i+(j-1)*m, m + 2 + (j-1)*m] = 1
#             elif((i==1) and (hx*i < L)):               #Left
#                 A[i+(j-1)*m, 1 + (j-1)*m] = ((2*hy*H)/K) - 3
#                 A[i+(j-1)*m, 2 + (j-1)*m] = 4
#                 A[i+(j-1)*m, 3 + (j-1)*m] = -1
#             elif(i==1):                                #Power
#                 A[i+(j-1)*m, 1 + (j-1)*m] = 3
#                 A[i+(j-1)*m, 2 + (j-1)*m] = -4
#                 A[i+(j-1)*m, 3 + (j-1)*m] = 1
#             else:                                      #Base Plate
#                 A[i+(j-1)*m, i+(j-1)*m] = -(2*(hy**2)) - (2*(hx**2)) - (2*(hx**2)*(hy**2)*((2*H)/(K*delta)))
#                 A[i+(j-1)*m, i-1+(j-1)*m] = hy**2
#                 A[i+(j-1)*m, i+1+(j-1)*m] = hy**2
#                 A[i+(j-1)*m, i+(j-2)*m]   = hx**2
#                 A[i+(j-1)*m, i+j*m]       = hx**2
#     # print(A)
#     # print(A.shape)

#     tempStr = ""
#     j = 2 + (1-1)*m
#     for ii in range(1, m+1):
#         for iii in range(1, n+1):
#             i = ii+(iii-1)*m
#             tempStr += str(A[j, i]) + " "
#         tempStr += 'XXX\n'
#     print("{}".format(tempStr))
#     print()

#     A = A[1:n*m+1, 1:n*m+1]
#     # print(A)
#     # print(A.shape)

#     return A
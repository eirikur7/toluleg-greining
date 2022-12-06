import numpy as np
from numpy import linalg as LA
import matplotlib.pyplot as plt

def createAMatrix2(n, m, L, Lx, Ly, H, K, delta):
    A = np.zeros((n*m, m*m))
    hx = Lx/(m - 1)
    hy = Ly/(n -1)
    testDic = {}
    for j in range(n):
        for i in range(m):
            eqNr = i + j*m
            if((j==0) and (i != 0) and (i != (m-1))):     #Top
                A[eqNr, eqNr]       = ((2*hy*H)/K) + 3
                A[eqNr, eqNr + m]   = -4
                A[eqNr, eqNr + 2*m] = 1
                testDic["Top"] = (i, j)
            elif(( j==(n-1) ) and (i != 0) and (i != (m-1))):  #Bottom
                A[eqNr, eqNr]       = ((-2*hy*H)/K) - 3
                # if((eqNr + m ) <= (m*n-1)):
                A[eqNr, eqNr - m]   = 4
                # if((eqNr + 2*m ) <= (m*n-1)):
                A[eqNr, eqNr - 2*m] = -1
                testDic["Bottom"] = (i, j)
            elif(( i==(m-1) )):                            #Right
                A[eqNr, eqNr] = ((-2*hy*H)/K) + 3
                # if((eqNr + 1) < (m*n -1)):
                A[eqNr, eqNr - 1] = -4
                # if((eqNr + 2) < (m*n - 1)):
                A[eqNr, eqNr - 2] = 1
                testDic["Right"] = (i, j)
            elif((i==0) and ( (Lx-L - (hx*i))>0 )):               #Left
                A[eqNr, eqNr] = ((2*hy*H)/K) - 3
                A[eqNr, eqNr + 1] = 4
                A[eqNr, eqNr + 2] = -1
            elif(i==0):                                #Power
                A[eqNr, eqNr] = 3
                A[eqNr, eqNr + 1] = -4
                A[eqNr, eqNr + 2] = 1
                testDic["Power"] = (i, j)
            else:                                      #Base Plate
                A[eqNr, eqNr] = -(2*(hy**2)) - (2*(hx**2)) - ( ((hx**2)*(hy**2)*2*H)/(K*delta) )
                A[eqNr, eqNr-1] = hy**2
                A[eqNr, eqNr+1] = hy**2
                A[eqNr, eqNr+m] = hx**2
                A[eqNr, eqNr-m] = hx**2
                testDic["Plate"] = (i, j)
    # print(A)
    return A, testDic



def createAMatrix(n, m, L, Lx, Ly, H, K, delta):
    A = np.zeros((n*m+1, m*m+1))
    hx = Lx/(m- 1)
    hy = Ly/(n -1 )
    for j in range(1, n+1):
        for i in range(1, m+1):
            if((j==1) and (i != 1) and (i != m)):     #Top
                A[i+(j-1)*m, i]       = ((2*hy*H)/K) + 3
                A[i+(j-1)*m, i + m]   = -4
                A[i+(j-1)*m, i + 2*m] = 1
            elif(( j==n ) and (i != 1) and (i != m)):  #Bottom
                A[i+(j-1)*m, i]       = ((2*hy*H)/K) - 3
                A[i+(j-1)*m, i + m]   = 4
                A[i+(j-1)*m, i + 2*m] = -1
            elif(( i==m )):                            #Right
                A[i+(j-1)*m, m + (j-1)*m] = ((2*hy*H)/K) + 3
                if((m + 0 + (j-1)*m) < (m*n -1)):
                    A[i+(j-1)*m, m + 1 + (j-1)*m] = -4
                if((m + 1 + (j-1)*m) < (m*n -1)):
                    A[i+(j-1)*m, m + 2 + (j-1)*m] = 1
            elif((i==1) and (hx*i < L)):               #Left
                A[i+(j-1)*m, 1 + (j-1)*m] = ((2*hy*H)/K) - 3
                A[i+(j-1)*m, 2 + (j-1)*m] = 4
                A[i+(j-1)*m, 3 + (j-1)*m] = -1
            elif(i==1):                                #Power
                A[i+(j-1)*m, 1 + (j-1)*m] = 3
                A[i+(j-1)*m, 2 + (j-1)*m] = -4
                A[i+(j-1)*m, 3 + (j-1)*m] = 1
            else:                                      #Base Plate
                A[i+(j-1)*m, i+(j-1)*m] = -(2*(hy**2)) - (2*(hx**2)) - (2*(hx**2)*(hy**2)*((2*H)/(K*delta)))
                A[i+(j-1)*m, i-1+(j-1)*m] = hy**2
                A[i+(j-1)*m, i+1+(j-1)*m] = hy**2
                A[i+(j-1)*m, i+(j-2)*m]   = hx**2
                A[i+(j-1)*m, i+j*m]       = hx**2
    # print(A)
    # print(A.shape)

    tempStr = ""
    j = 2 + (1-1)*m
    for ii in range(1, m+1):
        for iii in range(1, n+1):
            i = ii+(iii-1)*m
            tempStr += str(A[j, i]) + " "
        tempStr += 'XXX\n'
    print("{}".format(tempStr))
    print()

    A = A[1:n*m+1, 1:n*m+1]
    # print(A)
    # print(A.shape)

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
    K = 1.68
    H = 0.005
    m = 10
    n = 10
    A,tDic = createAMatrix2(n, m, L, Lx, Ly, H, K, delta)
    B = createBMatrix(n, m, L, Lx, P, delta, K)

    for key in tDic:
        x,y = tDic[key]
        eqNr = x + m*y
        print("-------------------------{}-------------------".format(key))
        print("x={},y={}".format(x,y))
        for j in range(n):
            for i in range(m):
                print(A[eqNr, i + j*m],end="| ")
            print()
        print('----------------------------------------------')


    # print(A.shape)
    # print(B)
    # for j in range(m*n):
    #     tempStr = ""
    #     for i in range(m*n):
    #         tempStr += str(A[j, i]) + " "
    #     print("{}. {}".format(j, tempStr))
    #     print()

    V = LA.solve(A, B)
    # print(V)
    print(V[0])





if __name__ == "__main__":
    prob3()
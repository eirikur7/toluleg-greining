import numpy as np
from numpy import linalg as LA
import matplotlib.pyplot as plt
import matplotlib
import matplotlib.cm as cm
import time
from mpl_toolkits.mplot3d import Axes3D
from mpl_toolkits import mplot3d
import matplotlib.pyplot as plt
from matplotlib import cm
from matplotlib.ticker import LinearLocator
import numpy as np

ANIMATION_TEXT = "Generating animation..."
PLOT_TEXT = "Generating plot..."
FIG_PATH = "src/figures_project3/"
ANIMATION_PATH = FIG_PATH + "{}_animation.gif"
PLOT_PATH = FIG_PATH + "{}_plot.png"

# ---------------------------- Helper Functions ---------------------------- #
def createAMatrix2(n, m, L, Lx, Ly, H, K, delta):
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
            elif((i==0) and ( (Lx-L - (hx*i))>0 )):             #Left
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

def solveSys(Lx, Ly, delta, P, Lspan, K, H, m, n):
    # Lx = 2
    # Ly = 2
    # delta = 0.1
    # P = 5
    # L = 2
    # K = 1.68
    # H = 0.005
    A = createAMatrix(n, m, Lspan, Lx, Ly, H, K, delta)
    B = createBMatrix(n, m, Lspan, Lx, P, delta, K)
    V = LA.solve(A, B) + 20
    return A, B, V

#-----------------------Problems-----------------------#
def prob1():
    print('--- Problem 1 ---')
    print('Calculations can be found in the report')

def prob2():
    print('--- Problem 2 ---')
    print("calculations can be found in the report")


def prob3():
    Lx, Ly, m, n = 2, 2, 10, 10
    A, B, V = solveSys(Lx=2, Ly=2, delta=0.1, P=5, Lspan=[0,2], K=1.68, H=0.005, m=10, n=10)
    print(A.shape)
    print(B.shape)
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
    print('--- Problem 4 ---')

    # reference solution
    A_ref, B_ref, V_ref = solveSys(Lx=2, Ly=2, delta=0.1, P=5, L=2, K=1.68, H=0.005, m=100, n=100)

    # dictionary to store deviations and execution times
    deviations = {(10,10): {'dev':float('inf'), 'time':float('inf')}}
    nArr = range(10,91,10)  # 10, 20, 30, ..., 90
    mArr = range(10,91,10)  # 10, 20, 30, ..., 90

    # Loop over all combinations of m and n
    fig = plt.figure()
    ax = plt.axes()
    for n in nArr:
        for m in mArr:
            start_time = time.time()    # start timer
            A, B, V = solveSys(Lx=2, Ly=2, delta=0.1, P=5, L=2, K=1.68, H=0.005, m=m, n=n)    # solve system of equations for m and n

            # Compute the deviation from the reference solution
            dev = abs(V_ref[0] - V[0])
            exec_time = time.time() - start_time # stop timer

            deviations[(m,n)] = {'dev':dev, 'time':exec_time}   # save deviation and execution time

            # only plot if deviation is less than 0.01 and time is less than 0.5
            if dev < 0.01 and exec_time < 0.5:
                # add dot to the plot
                sc = ax.scatter(m, n, c=dev, vmin=0, vmax=0.01, cmap='coolwarm')

    plt.colorbar(sc, label='deviation[°C]')
    ax.set_xlabel(xlabel='M')
    ax.set_ylabel(ylabel='N')
    # plt.show()

    # Further analysis, find the best combination of m and n
    fig = plt.figure()
    ax = plt.axes(projection='3d')
    for key in deviations:
        m,n = key
        dev = deviations[key]['dev']
        exec_time = deviations[key]['time']
        if exec_time < 0.3 and m > 20:
            sc = ax.scatter3D(m, n, dev, c=exec_time, vmin=0, vmax=0.5, cmap='coolwarm')
            # if m == 80 and n == 40:
            #     ax.scatter3D(m, n, dev, c='green', marker='o', s=100)
            if m == 90 and n == 30:
                ax.scatter3D(m, n, dev, c='green', marker='o', s=100)
            # elif m == 80 and n == 40:
            #     ax.scatter3D(m, n, dev, c='red', marker='x', s=100)
    plt.colorbar(sc, label='time[s]')
    ax.set_xlabel(xlabel='M')
    ax.set_ylabel(ylabel='N')
    ax.set_zlabel(zlabel='deviation[°C]')

    # calculate the slope of the deviation from m=40 to m=90
    dev40 = deviations[(40,40)]['dev']
    dev90 = deviations[(90,40)]['dev']
    slope = (dev90 - dev40) / (90 - 40)
    print(f'slope m=[40,90]: {slope}')

    # calculate the slope of the deviation from n=20 to n=90
    dev20 = deviations[(40,20)]['dev']
    dev90 = deviations[(40,90)]['dev']
    slope = (dev90 - dev20) / (90 - 20)
    print(f'slope n=[20,90]: {slope}')


    # print the valeus of the best combination of m and n
    print('Best combination of m and n:')
    print(f'm=80, n=40 : deviation {deviations[(80,40)]} ')
    print(f'm=90, n=30 : deviation {deviations[(90,30)]} ')


    # which matters more m or n?
    fig = plt.figure()
    ax = plt.axes(projection='3d')
    M1, N1 = np.meshgrid(mArr, nArr)
    dev = np.matrix([deviations[key]['dev'] for key in deviations]).reshape((len(nArr),len(mArr)))
    ax.plot_surface(M1, N1, dev, cmap='coolwarm', edgecolor='none')
    ax.set_xlabel(xlabel='M')
    ax.set_ylabel(ylabel='N')
    ax.set_zlabel(zlabel='Deviation in (0,0)[°C]')
    plt.show()


def prob5():
    print('--- Problem 5 ---')
    print('Not implemented')

def prob6():
    print('--- Problem 6 ---')
    print('Not implemented')

def prob7():
    print('--- Problem 7 ---')
    print('Not implemented')

def prob8():
    print('--- Problem 8 ---')
    print('Not implemented')

def prob9():
    print('--- Problem 9 ---')
    print('Not implemented')



if __name__ == "__main__":
    # prob3()
    prob4()

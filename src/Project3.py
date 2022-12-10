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
PLOT_SPAN = FIG_PATH + "problem{}_({},{}).png"

bestN = 30
bestM = 90

# ---------------------------- Helper Functions ---------------------------- #
def createAMatrix(n, m, LSpan, Lx, Ly, H, K, delta):
    A = np.zeros((n*m, m*n))
    hx = Lx/(m - 1)
    hy = Ly/(n -1)
    for j in range(n):
        for i in range(m):
            eqNr = i + j*m
            if((j==(n-1)) and (i != 0) and (i != (m-1))):               #Top
                A[eqNr, eqNr]       = ((2*hy*H)/K) - 3
                A[eqNr, eqNr - m]   = 4
                A[eqNr, eqNr - 2*m] = -1
            elif(( j==0 ) and (i != 0) and (i != (m-1))):               #Bottom
                A[eqNr, eqNr]       = ((2*hy*H)/K) - 3
                A[eqNr, eqNr + m]   = 4
                A[eqNr, eqNr + 2*m] = -1
            elif(( i==(m-1) )):                                         #Right
                A[eqNr, eqNr]       = ((2*hx*H)/K) - 3
                A[eqNr, eqNr - 1]   = 4
                A[eqNr, eqNr - 2]   = -1
            elif((i==0) and ((j*hy < LSpan[0]) or (j*hy > LSpan[1]))):  #Left
                A[eqNr, eqNr]       = ((2*hx*H)/K) - 3
                A[eqNr, eqNr + 1]   = 4
                A[eqNr, eqNr + 2]   = -1
            elif(i==0):                                                 #Power
                A[eqNr, eqNr]       = 3/(2*hx)
                A[eqNr, eqNr + 1]   = -4/(2*hx)
                A[eqNr, eqNr + 2]   = 1/(2*hx)
            else:                                                       #Base Plate
                A[eqNr, eqNr]       = -2*( (H/(K*delta)) + (1/(hx**2)) + (1/(hy**2)) )
                A[eqNr, eqNr-1]     = 1/(hx**2)
                A[eqNr, eqNr+1]     = 1/(hx**2)
                A[eqNr, eqNr+m]     = 1/(hy**2)
                A[eqNr, eqNr-m]     = 1/(hy**2)
    return A

def createAMatrixCut(cutN, cutM, n, m, LSpan, Lx, Ly, H, K, delta):
    A = np.zeros((n*m, m*n))
    hx = Lx/(m - 1)
    hy = Ly/(n -1)
    cutN += 1
    cutM += 1
    for j in range(n):
        for i in range(m):
            eqNr = i + j*m
            setSpecialCase = False
            if ((j >= (n-cutN)) or (i >= (m-cutM))): #Not regular
                setSpecialCase = True
                if((j > (n-cutN)) and (i > (m-cutM))): #air
                    A[eqNr, eqNr]       = 1
                elif( (j == (n-cutN)) and (i == (m-cutM))): #Base Plate
                    A[eqNr, eqNr]       = -2*( (H/(K*delta)) + (1/(hx**2)) + (1/(hy**2)) )
                    A[eqNr, eqNr-1]     = 1/(hx**2)
                    A[eqNr, eqNr+1]     = 1/(hx**2)
                    A[eqNr, eqNr+m]     = 1/(hy**2)
                    A[eqNr, eqNr-m]     = 1/(hy**2)
                elif( (j >= (n-cutN)) and (i == (m-cutM))): #Right
                    A[eqNr, eqNr]       = ((2*hx*H)/K) - 3
                    A[eqNr, eqNr - 1]   = 4
                    A[eqNr, eqNr - 2]   = -1
                elif( (j == (n-cutN)) and (i == (m-1))):    #Right (/top)
                    A[eqNr, eqNr]       = ((2*hx*H)/K) - 3
                    A[eqNr, eqNr - 1]   = 4
                    A[eqNr, eqNr - 2]   = -1
                elif( (j == (n-cutN)) and (i >= (m-cutM))): #Top
                    A[eqNr, eqNr]       = ((2*hy*H)/K) - 3
                    A[eqNr, eqNr - m]   = 4
                    A[eqNr, eqNr - 2*m] = -1
                else:
                    setSpecialCase = False
            if(not setSpecialCase):
                if((j==(n-1)) and (i != 0) and (i != (m-1))):               #Top
                    A[eqNr, eqNr]       = ((2*hy*H)/K) - 3
                    A[eqNr, eqNr - m]   = 4
                    A[eqNr, eqNr - 2*m] = -1
                elif(( j==0 ) and (i != 0) and (i != (m-1))):               #Bottom
                    A[eqNr, eqNr]       = ((2*hy*H)/K) - 3
                    A[eqNr, eqNr + m]   = 4
                    A[eqNr, eqNr + 2*m] = -1
                elif(( i==(m-1) )):                                         #Right
                    A[eqNr, eqNr]       = ((2*hx*H)/K) - 3
                    A[eqNr, eqNr - 1]   = 4
                    A[eqNr, eqNr - 2]   = -1
                elif((i==0) and ((j*hy < LSpan[0]) or (j*hy > LSpan[1]))):  #Left
                    A[eqNr, eqNr]       = ((2*hx*H)/K) - 3
                    A[eqNr, eqNr + 1]   = 4
                    A[eqNr, eqNr + 2]   = -1
                elif(i==0):                                                 #Power
                    A[eqNr, eqNr]       = 3/(2*hx)
                    A[eqNr, eqNr + 1]   = -4/(2*hx)
                    A[eqNr, eqNr + 2]   = 1/(2*hx)
                else:                                                       #Base Plate
                    A[eqNr, eqNr]       = -2*( (H/(K*delta)) + (1/(hx**2)) + (1/(hy**2)) )
                    A[eqNr, eqNr-1]     = 1/(hx**2)
                    A[eqNr, eqNr+1]     = 1/(hx**2)
                    A[eqNr, eqNr+m]     = 1/(hy**2)
                    A[eqNr, eqNr-m]     = 1/(hy**2)
    return A

def createBMatrix(n, m, LSpan, Ly, P, delta, K):
    B = np.zeros((m*n, 1))
    hy = Ly/(n-1)
    L = LSpan[1]-LSpan[0]
    for j in range(n):
        if (((j*hy >= LSpan[0]) and (j*hy <= LSpan[1]))):
            B[j*m] = (P)/(L*delta*K)

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

def findHighestTemp(P, LSpan, K):
    Lx, Ly, delta, H = 4, 4, 0.1, 0.005
    A, B, V = solveSys(Lx, Ly, delta, P, LSpan, K, H, bestM, bestN)
    return np.max(V) - 100


def bisection(f,a,b,tol,LSpan, K):
    if f(a, LSpan, K)*f(b, LSpan, K) >= 0:
        print("Bisection method fails.")
        return None
    else:
        fa=f(a,LSpan, K)
        while( ((b-a)/2>tol) or (f((a+b)/2, LSpan, K) > 0) ):
            c=(a+b)/2
            fc=f(c,LSpan, K)
            if fc==0:break
            if fc*fa<0:
                b=c
            else:
                a=c
                fa=fc
    return((a+b)/2)

def findHighestPowerBeforeLimit(tempLimit, tol, LSpan, K):
    Lx, Ly, delta, H = 4, 4, 0.1, 0.005
    #Preform exponential increase of P until max temperature becomes more than limit. To get upper limit(b)
    P_b = 2
    maxTemp = 0
    while maxTemp <= tempLimit:
        P_b = P_b*2
        V = LA.solve(createAMatrix(bestN, bestM, LSpan, Lx, Ly, H, K, delta), createBMatrix(bestN, bestM, LSpan, Ly, P_b, delta, K))
        maxTemp = np.max(V)
    #Preform exponential decrease of P until max temperature becomes less than limit. To get lower limit(a)
    P_a = 4
    maxTemp = tempLimit
    while maxTemp >= tempLimit:
        P_a = P_a/2
        V = LA.solve(createAMatrix(bestN, bestM, LSpan, Lx, Ly, H, K, delta), createBMatrix(bestN, bestM, LSpan, Ly, P_a, delta, K))
        maxTemp = np.max(V)
        # print("P={}, T={}".format(P_a, maxTemp))
    
    maxP = bisection(findHighestTemp, P_a, P_b, tol, LSpan, K)
    return maxP

def getOnlyPlate(V, n, m, cutN, cutM, Lx, Ly):
    hx = Lx/(m-1)
    hy = Ly/(n-1)
    V_plate = np.zeros((m*n-(cutM*cutN), 1))
    x_plate = np.zeros((m*n-(cutM*cutN), 1))
    y_plate = np.zeros((m*n-(cutM*cutN), 1))

    V_notch = np.zeros(((cutM*cutN), 1))
    x_notch = np.zeros(((cutM*cutN), 1))
    y_notch = np.zeros(((cutM*cutN), 1))
    y_cur = 0
    cntNot = 0
    for y in range(n):
        x_cur = 0
        for x in range(m):
            eq = x + y*m
            if (((y<(n-cutN)) or (x <(m-cutM)))):
                V_plate[eq-cntNot,0] = V[eq, 0]
                x_plate[eq-cntNot,0] = x_cur
                y_plate[eq-cntNot,0] = y_cur
            else:
                V_notch[cntNot,0] = V[eq, 0]
                x_notch[cntNot,0] = x_cur
                y_notch[cntNot,0] = y_cur
                cntNot +=1
            x_cur += hx
        y_cur += hy
    return [x_plate, y_plate, V_plate], [x_notch, y_notch, V_notch]

def findHighestTempFor9(P, LSpan, cutN, cutM):
    Lx, Ly, delta, K, H = 4, 4, 0.1, 1.68, 0.005
    L = 2
    A = createAMatrixCut(cutN, cutM, bestN, bestM, LSpan, Lx, Ly, H, K, delta)
    B = createBMatrix(bestN, bestM, LSpan, Ly, P, delta, K)
    V = LA.solve(A, B) + 20
    plate, _ = getOnlyPlate(V, bestN, bestM, cutN, cutM, Lx, Ly)
    return np.max(plate[2]) - 100


def bisectionFor9(f, a, b, tol, LSpan, cutN, cutM):
    if f(a, LSpan, cutN, cutM)*f(b, LSpan, cutN, cutM) >= 0:
        print("Bisection method fails.")
        return None
    else:
        fa=f(a,LSpan, cutN, cutM)
        while (b-a)/2>tol or (f((a+b)/2, LSpan, cutN, cutM) > 0):
            c=(a+b)/2
            fc=f(c,LSpan, cutN, cutM)
            if fc==0:break
            if fc*fa<0:
                b=c
            else:
                a=c
                fa=fc
            print(f((a+b)/2, LSpan, cutN, cutM))
    
    return((a+b)/2)


def findHighestPowerBeforeLimitFor9(tempLimit, tol, LSpan, cutN, cutM):
    Lx, Ly, delta, K, H = 4, 4, 0.1, 1.68, 0.005 
    #Preform exponential increase of P until max temperature becomes more than limit. To get upper limit(b)
    P_b = 2
    maxTemp = 0
    while maxTemp <= tempLimit:
        P_b = P_b*2
        V = LA.solve(createAMatrixCut(cutN, cutM, bestN, bestM, LSpan, Lx, Ly, H, K, delta), createBMatrix(bestN, bestM, LSpan, Ly, P_b, delta, K))
        plate, _ = getOnlyPlate(V, bestN, bestM, cutN, cutM, Lx, Ly)
        maxTemp = np.max(plate[2])
        # if(maxTemp == tempLimit):
        #     return P_b
        # print("P_b={}, temp={}".format(P_b, maxTemp))
    
    #Preform exponential decrease of P until max temperature becomes less than limit. To get lower limit(a)
    P_a = 2
    maxTemp = tempLimit
    while maxTemp >= tempLimit:
        P_a = P_a/2
        V = LA.solve(createAMatrixCut(cutN, cutM, bestN, bestM, LSpan, Lx, Ly, H, K, delta), createBMatrix(bestN, bestM, LSpan, Ly, P_a, delta, K))
        plate, _ = getOnlyPlate(V, bestN, bestM, cutN, cutM, Lx, Ly)
        maxTemp = np.max(plate[2])
        # if(maxTemp == tempLimit):
        #     return P_a
        # print("P_a={}, temp={}".format(P_a, maxTemp))
    
    maxP = bisectionFor9(findHighestTempFor9, P_a, P_b, tol, LSpan, cutN, cutM)
    return maxP


def prob6For9(cutN, cutM, plot=True):
    L = 2
    Lx, Ly, delta, P, K, H= 4, 4, 0.1, 5, 1.68, 0.005
    
    LSpansTries = int(bestN/2)
    LSpansEnd = np.linspace(L, Ly, LSpansTries)
    maxTemp = np.zeros((LSpansTries, 1))
    for i in range(LSpansTries):
        LSpan = [LSpansEnd[i]-L, LSpansEnd[i]]
        A = createAMatrixCut(cutN, cutM, bestN, bestM, LSpan, Lx, Ly, H, K, delta)
        B = createBMatrix(bestN, bestM, LSpan, Ly, P, delta, K)
        V = LA.solve(A, B) + 20
        plate, _ = getOnlyPlate(V, bestN, bestM, cutN, cutM, Lx, Ly)
        if plot:
            fig = plt.figure()
            ax = plt.axes(projection='3d')
            for ii in range(plate[2].size):
                ax.scatter(plate[0][ii], plate[1][ii], plate[2][ii], vmin=np.min(plate[2]), vmax=np.max(plate[2]), cmap='coolwarm', c=plate[2][ii])
            ax.set_xlabel(xlabel='x[cm]')
            ax.set_ylabel(ylabel='y[cm]')
            ax.set_zlabel(zlabel='temperature[°C]')
            ax.set_title("Span:[{:.3f}, {:.3f}]".format(LSpan[0], LSpan[1]))
            # fig.savefig("hrafnljoti/proj3Fig/prob9.6_[{:.3f}, {:.3f}].png".format(LSpan[0], LSpan[1]))
            fig.savefig(PLOT_SPAN.format(9, LSpan[0], LSpan[1]))
            plt.close(fig)

        print("Span={}, temp={}".format(LSpan, np.max(plate[2])))
        maxTemp[i] = np.max(plate[2])
        # plt.show()

    bestIndex = np.argmin(maxTemp)
    LSpan = [LSpansEnd[bestIndex]-L, LSpansEnd[bestIndex]]
    A = createAMatrixCut(cutN, cutM, bestN, bestM, LSpan, Lx, Ly, H, K, delta)
    B = createBMatrix(bestN, bestM, LSpan, Ly, P, delta, K)
    V = LA.solve(A, B) + 20
    plate, _ = getOnlyPlate(V, bestN, bestM, cutN, cutM, Lx, Ly)

    fig2, ax2 = plt.subplots(1,1)
    ax2.plot(LSpansEnd, maxTemp, 'b-')
    ax2.plot(LSpansEnd, maxTemp, 'r*')
    ax2.set_xlabel("End Span[cm]")
    ax2.set_ylabel("Max Temperature[◦C]")
    plt.show()

    # fig = plt.figure()
    # ax = plt.axes(projection='3d')
    # for ii in range(plate[2].size):
    #     ax.scatter(plate[0][ii], plate[1][ii], plate[2][ii], vmin=np.min(plate[2]), vmax=np.max(plate[2]), cmap='coolwarm', c=plate[2][ii])
    # ax.set_xlabel(xlabel='x[cm]')
    # ax.set_ylabel(ylabel='y[cm]')
    # ax.set_zlabel(zlabel='temperature[°C]')
    # ax.set_title("Span:[{:.3f}, {:.3f}]".format(LSpan[0], LSpan[1]))
    # fig.savefig(PLOT_SPAN.format("9_best", LSpan[0], LSpan[1]))
    # print("Best span: [{}, {}]: maxTemp = {}".format(LSpansEnd[bestIndex]-L, LSpansEnd[bestIndex], maxTemp[bestIndex]))
    # plt.show()
    
    return LSpansEnd[bestIndex]

def prob7For9(bestEnd6 = 3, cutN = 0, cutM = 0):
    L = 2
    Lx, Ly, LSpan, delta, K, H = 4, 4, [bestEnd6 - L, bestEnd6], 0.1, 1.68, 0.005

    maxP = findHighestPowerBeforeLimitFor9(100, 10**(-2), LSpan, cutN, cutM)
    V = LA.solve(createAMatrixCut(cutN, cutM, bestN, bestM, LSpan, Lx, Ly, H, K, delta), createBMatrix(bestN, bestM, LSpan, Ly, maxP, delta, K)) + 20
    plate, _ = getOnlyPlate(V, bestN, bestM, cutN, cutM, Lx, Ly)
    maxTemp = np.max(plate[2])
    print("Max Power={}, maxTemp={}".format(maxP, maxTemp))
    fig = plt.figure()
    ax = plt.axes(projection='3d')
    for i in range(plate[2].size):
        ax.scatter(plate[0][i], plate[1][i], plate[2][i], vmin=np.min(plate[2]), vmax=np.max(plate[2]), cmap='coolwarm', c=plate[2][i])
    ax.set_xlabel(xlabel='x[cm]')
    ax.set_ylabel(ylabel='y[cm]')
    ax.set_zlabel(zlabel='temperature[°C]')
    ax.set_title("Span:[{:.3f}, {:.3f}]".format(LSpan[0], LSpan[1]))
    plt.show()
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
    print("Temp_(0, 0) = {}, Temp_(0, Ly)".format(V[0], V[n*m]))
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
    A_ref, B_ref, V_ref = solveSys(Lx=2, Ly=2, delta=0.1, P=5, Lspan=[0, 2], K=1.68, H=0.005, m=100, n=100)

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
            A, B, V = solveSys(Lx=2, Ly=2, delta=0.1, P=5, Lspan=[0, 2], K=1.68, H=0.005, m=m, n=n)    # solve system of equations for m and n

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

    # calculate the slope of the deviation from m=10 to m=30
    dev10 = deviations[(10,10)]['dev']
    dev40 = deviations[(40,10)]['dev']
    slopeM_10_40 = (dev40 - dev10) / (40 - 10)
    print(f'slope m=[10,40]: {slopeM_10_40}')

    # calculate the slope of the deviation from m=10 to n=30
    dev10 = deviations[(10,10)]['dev']
    dev20 = deviations[(10,20)]['dev']
    slopeN_10_20 = (dev20 - dev10) / (20 - 10)
    print(f'slope n=[10,20]: {slopeN_10_20}')

    # calculate the slope of the deviation from m=30 to m=90
    dev40 = deviations[(40,40)]['dev']
    dev90 = deviations[(90,40)]['dev']
    slopeM_40_90 = (dev90 - dev40) / (90 - 40)
    print(f'slope m=[40,90]: {slopeM_40_90}')

    # calculate the slope of the deviation from n=30 to n=90
    dev20 = deviations[(40,20)]['dev']
    dev90 = deviations[(40,90)]['dev']
    slopeN_20_90 = (dev90 - dev20) / (90 - 20)
    print(f'slope n=[20,90]: {slopeN_20_90}')

    # print the ratio of the slopes from m=10 to m=30 and n=10 to n=30
    print(f'ratio m=[10,40]/n=[10,20]: {slopeM_10_40/slopeN_10_20}')

    # print the ratio of the slopes from m=30 to m=90 and n=30 to n=90
    print(f'ratio m=[40,90]/n=[20,90]: {slopeM_40_90/slopeN_20_90}')

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
    Lx, Ly, LSpan, delta, P, K, H = 4, 4, [0, 2], 0.1, 5, 1.68, 0.005
    L = 2
    A, B, V = solveSys(Lx, Ly, delta, P, LSpan, K, H, bestM, bestN)

    fig = plt.figure()
    ax = plt.axes(projection='3d')
    x = np.linspace(0, Lx, bestM)
    y = np.linspace(0, Ly, bestN)
    X1, Y1 = np.meshgrid(x, y)
    ax.plot_surface(X1, Y1, V.reshape((bestN,bestM)), cmap='coolwarm', edgecolor='none')
    ax.set_xlabel(xlabel='x[cm]')
    ax.set_ylabel(ylabel='y[cm]')
    ax.set_zlabel(zlabel='temperature[°C]')
    plt.show()

def prob6(plot=True):
    print('--- Problem 6 ---')
    Lx, Ly, L, delta, P, K, H = 4, 4, 2, 0.1, 5, 1.68, 0.005
    LSpansTries = int(bestN/2)
    LSpansEnd = np.linspace(L, Ly, LSpansTries)
    maxTemp = np.zeros((LSpansTries, 1))
    for i in range(LSpansTries):
        LSpan = [LSpansEnd[i]-L, LSpansEnd[i]]
        A, B, V = solveSys(Lx, Ly, delta, P, LSpan, K, H, bestM, bestN)
        if (plot):
            fig = plt.figure()
            ax = plt.axes(projection='3d')
            x = np.linspace(0, Lx, bestM)
            y = np.linspace(0, Ly, bestN)
            X1, Y1 = np.meshgrid(x, y)
            ax.plot_surface(X1, Y1, V.reshape((bestN,bestM)), cmap='coolwarm', edgecolor='none')
            ax.set_xlabel(xlabel='x[cm]')
            ax.set_ylabel(ylabel='y[cm]')
            ax.set_zlabel(zlabel='temperature[°C]')
            ax.set_title("Span:[{:.3f}, {:.3f}]".format(LSpan[0], LSpan[1]))
            # fig.savefig("hrafnljoti/proj3Fig/prob6_[{:.3f}, {:.3f}].png".format(LSpan[0], LSpan[1]))
            fig.savefig(PLOT_SPAN.format(6, LSpan[0], LSpan[1]))
            plt.close(fig)

        maxTemp[i] = np.max(V)
        print("Span={}, temp={}".format(LSpan, np.max(V)))
        # plt.show()

    bestIndex = np.argmin(maxTemp)
    LSpan = [LSpansEnd[bestIndex]-L, LSpansEnd[bestIndex]]
    A, B, V = solveSys(Lx, Ly, delta, P, LSpan, K, H, bestM, bestN)
    fig = plt.figure()
    ax = plt.axes(projection='3d')
    x = np.linspace(0, Lx, bestM)
    y = np.linspace(0, Ly, bestN)
    X1, Y1 = np.meshgrid(x, y)
    ax.plot_surface(X1, Y1, V.reshape((bestN,bestM)), cmap='coolwarm', edgecolor='none')
    ax.set_xlabel(xlabel='x[cm]')
    ax.set_ylabel(ylabel='y[cm]')
    ax.set_zlabel(zlabel='temperature[°C]')
    ax.set_title("Span:[{:.3f}, {:.3f}]".format(LSpan[0], LSpan[1]))
    # fig.savefig("hrafnljoti/proj3Fig/prob6_[{:.3f}, {:.3f}].png".format(LSpan[0], LSpan[1]))
    fig.savefig(PLOT_SPAN.format("6_best", LSpan[0], LSpan[1]))
    print("Best span: [{}, {}]: maxTemp = {}".format(LSpansEnd[bestIndex]-L, LSpansEnd[bestIndex], maxTemp[bestIndex]))
    
    fig2, ax2 = plt.subplots(1,1)
    print()
    ax2.plot(LSpansEnd, maxTemp, 'b-')
    ax2.plot(LSpansEnd, maxTemp, 'r*')
    ax2.set_xlabel("End Span[cm]")
    ax2.set_ylabel("Max Temperature[◦C]")
    plt.show()
    return LSpansEnd[bestIndex] #bestIndex, LSpansEnd

def prob7(bestEnd6 = 3):
    print('--- Problem 7 ---')
    L = 2
    Lx, Ly, LSpan, delta, P, K, H = 4, 4, [bestEnd6 - L, bestEnd6], 0.1, 5, 1.68, 0.005

    maxP = findHighestPowerBeforeLimit(100, 10**(-2), LSpan, K)
    print(maxP)
    V = LA.solve(createAMatrix(bestN, bestM, LSpan, Lx, Ly, H, K, delta), createBMatrix(bestN, bestM, LSpan, Ly, maxP, delta, K)) + 20
    maxTemp = np.max(V)
    print("Max Power={}, maxTemp={}".format(maxP, maxTemp))
    fig = plt.figure()
    ax = plt.axes(projection='3d')
    x = np.linspace(0, Lx, bestM)
    y = np.linspace(0, Ly, bestN)
    X1, Y1 = np.meshgrid(x, y)
    ax.plot_surface(X1, Y1, V.reshape((bestN,bestM)), cmap='coolwarm', edgecolor='none')
    ax.set_xlabel(xlabel='x[cm]')
    ax.set_ylabel(ylabel='y[cm]')
    ax.set_zlabel(zlabel='temperature[°C]')
    ax.set_title("Span:[{:.3f}, {:.3f}]".format(LSpan[0], LSpan[1]))
    plt.show()

def prob8(bestEnd6 = 3):
    print('--- Problem 8 ---')
    L = 2
    LSpan = [bestEnd6 - L, bestEnd6]
    kRuns = 20
    K = np.linspace(1, 5, kRuns)
    highestPower = np.zeros((kRuns, 1))
    for i in range(kRuns):
        highestPower[i,0] = findHighestPowerBeforeLimit(100, 10**(-2), LSpan, K=K[i])
        print("K={}, Max Power={:.2f}".format(K[i], highestPower[i,0]))

    fig, ax = plt.subplots(1,1)
    ax.plot(K, highestPower, 'b-')
    ax.plot(K, highestPower, 'r*')
    ax.set_xlabel("K[W/cm◦C]")
    ax.set_ylabel("Max Power Before exceeding 100◦C[W]")
    plt.show()

def prob9():
    print('--- Problem 9 ---')
    cutM = 18
    cutN = 6

    bestEnd6 = prob6For9(cutN, cutM, plot=False)
    prob7For9(bestEnd6=bestEnd6, cutN=cutN, cutM=cutN)



if __name__ == "__main__":
    # prob3()
    prob4()
    # bestEnd6 = prob6(plot=False)
    # prob7(bestEnd6=bestEnd6)
    # prob8()
    # prob9()
import numpy as np
from numpy import linalg as LA
import matplotlib.pyplot as plt
from mpl_toolkits import mplot3d
import time

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
    hy = Ly/(n)
    L = LSpan[1]-LSpan[0]
    for j in range(n):
        if (((j*hy >= LSpan[0]) and (j*hy <= LSpan[1]))):
            B[j*m] = (P)/(L*delta*K)

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
    B = createBMatrix(n, m, LSpan, Ly, P, delta, K)
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
    V_ref = LA.solve(createAMatrix(100, 100, LSpan, Lx, Ly, H, K, delta), createBMatrix(100, 100, LSpan, Ly, P, delta, K))
    V_ref_compare = V_ref[0]

    for i in range(totalValues):
        for ii in range(totalValues):
            time_to_compare = time.time()
            n = 10*(i + 1)
            m = 10*(ii + 1)
            V = LA.solve(createAMatrix(n, m, LSpan, Lx, Ly, H, K, delta), createBMatrix(n, m, LSpan, Ly, P, delta, K))
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
    B = createBMatrix(n, m, LSpan, Ly, P, delta, K)
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
    maxTemp = np.zeros((LSpansTries, 1))
    for i in range(LSpansTries):
        LSpan = [LSpansEnd[i]-L, LSpansEnd[i]]
        A = createAMatrix(n, m, LSpan, Lx, Ly, H, K, delta)
        B = createBMatrix(n, m, LSpan, Ly, P, delta, K)
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
        fig.savefig("hrafnljoti/proj3Fig/prob6_[{:.3f}, {:.3f}].png".format(LSpan[0], LSpan[1]))
        maxTemp[i] = np.max(V)
        # plt.show()
        plt.close(fig)

    bestIndex = np.argmin(maxTemp)
    print("Best span: [{}, {}]: maxTemp = {}".format(LSpansEnd[bestIndex]-L, LSpansEnd[bestIndex], maxTemp[bestIndex]))
    return bestIndex, LSpansEnd




def findHighestTemp(P, LSpan, K):
    Lx = 4
    Ly = 4
    L = 2
    delta = 0.1
    H = 0.005
    n = 10
    m = 10
    A = createAMatrix(n, m, LSpan, Lx, Ly, H, K, delta)
    B = createBMatrix(n, m, LSpan, Ly, P, delta, K)
    V = LA.solve(A, B) + 20
    return np.max(V) - 100


def bisection(f,a,b,tol,LSpan, K):
    if f(a, LSpan, K)*f(b, LSpan, K) >= 0:
        print("Bisection method fails.")
        return None
    else:
        fa=f(a,LSpan, K)
        while (b-a)/2>tol:
            c=(a+b)/2
            fc=f(c,LSpan, K)
            if fc==0:break
            if fc*fa<0:
                b=c
            else:
                a=c
                fa=fc
    return((a+b)/2)

# def bisection2(P, a, b, tol):



def findHighestPowerBeforeLimit(tempLimit, tol, LSpan, K):
    Lx = 4
    Ly = 4
    delta = 0.1
    H = 0.005
    n = 10
    m = 10
    #Preform exponential increase of P until max temperature becomes more than limit. To get upper limit(b)
    P_b = 2
    maxTemp = 0
    while maxTemp <= tempLimit:
        P_b = P_b*2
        V = LA.solve(createAMatrix(n, m, LSpan, Lx, Ly, H, K, delta), createBMatrix(n, m, LSpan, Ly, P_b, delta, K))
        maxTemp = np.max(V)

    #Preform exponential decrease of P until max temperature becomes less than limit. To get lower limit(a)
    P_a = 2
    maxTemp = tempLimit
    while maxTemp > tempLimit:
        P_a = P_a/2
        V = LA.solve(createAMatrix(n, m, LSpan, Lx, Ly, H, K, delta), createBMatrix(n, m, LSpan, Ly, P_b, delta, K))
        maxTemp = np.max(V)

    maxP = bisection(findHighestTemp, P_a, P_b, tol, LSpan, K)
    return maxP


def prob7(bestIndex6, LSpansEnd6):
    Lx = 4
    Ly = 4
    L = 2
    LSpan = [LSpansEnd6[bestIndex6] - L, LSpansEnd6[bestIndex6]]
    delta = 0.1
    K = 1.68
    H = 0.005
    n = 10
    m = 10

    maxP = findHighestPowerBeforeLimit(100, 10**(-2), LSpan, K)
    V = LA.solve(createAMatrix(n, m, LSpan, Lx, Ly, H, K, delta), createBMatrix(n, m, LSpan, Ly, maxP, delta, K)) + 20
    maxTemp = np.max(V)
    print("Max Power={}, maxTemp={}".format(maxP, maxTemp))
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
    plt.show()

def prob8(bestIndex6, LSpansEnd6):
    L = 2
    LSpan = [LSpansEnd6[bestIndex6] - L, LSpansEnd6[bestIndex6]]
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

def prob6For9(cutN, cutM):
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
    maxTemp = np.zeros((LSpansTries, 1))
    for i in range(LSpansTries):
        LSpan = [LSpansEnd[i]-L, LSpansEnd[i]]
        A = createAMatrixCut(cutN, cutM, n, m, LSpan, Lx, Ly, H, K, delta)
        B = createBMatrix(n, m, LSpan, Ly, P, delta, K)
        V = LA.solve(A, B) + 20
        plate, _ = getOnlyPlate(V, n, m, cutN, cutM, Lx, Ly)
        fig = plt.figure()
        ax = plt.axes(projection='3d')
        for ii in range(plate[2].size):
            ax.scatter(plate[0][ii], plate[1][ii], plate[2][ii], vmin=np.min(plate[2]), vmax=np.max(plate[2]), cmap='coolwarm', c=plate[2][ii])
        ax.set_xlabel(xlabel='x[cm]')
        ax.set_ylabel(ylabel='y[cm]')
        ax.set_zlabel(zlabel='temperature[°C]')
        ax.set_title("Span:[{:.3f}, {:.3f}]".format(LSpan[0], LSpan[1]))
        fig.savefig("hrafnljoti/proj3Fig/prob9.6_[{:.3f}, {:.3f}].png".format(LSpan[0], LSpan[1]))
        maxTemp[i] = np.max(plate[2])
        # plt.show()
        plt.close(fig)

    bestIndex = np.argmin(maxTemp)
    print("Best span: [{}, {}]: maxTemp = {}".format(LSpansEnd[bestIndex]-L, LSpansEnd[bestIndex], maxTemp[bestIndex]))
    return bestIndex, LSpansEnd

def findHighestTempFor9(P, LSpan, cutN, cutM):
    Lx = 4
    Ly = 4
    L = 2
    delta = 0.1
    K = 1.68
    H = 0.005
    n = 10
    m = 10
    A = createAMatrixCut(cutN, cutM,n, m, LSpan, Lx, Ly, H, K, delta)
    B = createBMatrix(n, m, LSpan, Ly, P, delta, K)
    V = LA.solve(A, B) + 20
    plate, _ = getOnlyPlate(V, n, m, cutN, cutM, Lx, Ly)
    return np.max(plate[2]) - 100


def bisectionFor9(f,a,b,tol,LSpan, cutN, cutM):
    if f(a, LSpan, cutN, cutM)*f(b, LSpan, cutN, cutM) >= 0:
        print("Bisection method fails.")
        return None
    else:
        fa=f(a,LSpan, cutN, cutM)
        while (b-a)/2>tol:
            c=(a+b)/2
            fc=f(c,LSpan, cutN, cutM)
            if fc==0:break
            if fc*fa<0:
                b=c
            else:
                a=c
                fa=fc
    return((a+b)/2)


def findHighestPowerBeforeLimitFor9(tempLimit, tol, LSpan, cutN, cutM):
    Lx = 4
    Ly = 4
    delta = 0.1
    K = 1.68
    H = 0.005
    n = 10
    m = 10
    #Preform exponential increase of P until max temperature becomes more than limit. To get upper limit(b)
    P_b = 2
    maxTemp = 0
    while maxTemp <= tempLimit:
        P_b = P_b*2
        V = LA.solve(createAMatrixCut(cutN, cutM, n, m, LSpan, Lx, Ly, H, K, delta), createBMatrix(n, m, LSpan, Ly, P_b, delta, K))
        plate, _ = getOnlyPlate(V, n, m, cutN, cutM, Lx, Ly)
        maxTemp = np.max(plate[2])

    #Preform exponential decrease of P until max temperature becomes less than limit. To get lower limit(a)
    P_a = 2
    maxTemp = tempLimit
    while maxTemp > tempLimit:
        P_a = P_a/2
        V = LA.solve(createAMatrixCut(cutN, cutM, n, m, LSpan, Lx, Ly, H, K, delta), createBMatrix(n, m, LSpan, Ly, P_b, delta, K))
        plate, _ = getOnlyPlate(V, n, m, cutN, cutM, Lx, Ly)
        maxTemp = np.max(plate[2])

    maxP = bisectionFor9(findHighestTempFor9, P_a, P_b, tol, LSpan, cutN, cutM)
    return maxP




def prob7For9(bestIndex, LSpansEnd, cutN, cutM):
    Lx = 4
    Ly = 4
    L = 2
    LSpan = [LSpansEnd[bestIndex] - L, LSpansEnd[bestIndex]]
    delta = 0.1
    K = 1.68
    H = 0.005
    n = 10
    m = 10

    maxP = findHighestPowerBeforeLimit(100, 10**(-2), LSpan, K)
    V = LA.solve(createAMatrix(n, m, LSpan, Lx, Ly, H, K, delta), createBMatrix(n, m, LSpan, Ly, maxP, delta, K)) + 20
    maxTemp = np.max(V)
    print("Max Power={}, maxTemp={}".format(maxP, maxTemp))
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
    plt.show()


def prob9():
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
    cutM = 3
    cutN = 2
    # A = createAMatrixCut(cutN, cutM, n, m, LSpan, Lx, Ly, H, K, delta)
    # B = createBMatrix(n, m, LSpan, Ly, P, delta, K)
    # V = LA.solve(A, B) + 20
    # print(V.shape)
    # eq = 5 + 7*m
    # for i in range(n):
    #     for ii in range(m):
    #         eq2 = ii + i*m
    #         print(A[eq, eq2], end=' ')
    #     print()


    bestIndex, LSpansEnd = prob6For9(cutN, cutM)
    prob7For9(bestIndex, LSpansEnd, cutN, cutM)



    # plate, notch = getOnlyPlate(V, n, m, cutN, cutM, Lx, Ly)
    # #both
    # fig0 = plt.figure()
    # x = np.linspace(0, Lx, m)
    # y = np.linspace(0, Ly, n)
    # X1, Y1 = np.meshgrid(x, y)
    # ax0 = plt.axes(projection='3d')
    # ax0.plot_surface(X1, Y1, V.reshape((n,m)), cmap='coolwarm', edgecolor='none')
    # ax0.set_xlabel(xlabel='x[cm]')
    # ax0.set_ylabel(ylabel='y[cm]')
    # ax0.set_zlabel(zlabel='temperature[°C]')
    # #plate
    # fig = plt.figure()
    # ax = plt.axes(projection='3d')
    # for i in range(plate[2].size):
    #     ax.scatter(plate[0][i], plate[1][i], plate[2][i], vmin=np.min(plate[2]), vmax=np.max(plate[2]), cmap='coolwarm', c=plate[2][i])
    # ax.set_xlabel(xlabel='x[cm]')
    # ax.set_ylabel(ylabel='y[cm]')
    # ax.set_zlabel(zlabel='temperature[°C]')
    # #Notch
    # fig2  = plt.figure()
    # ax2 = plt.axes(projection='3d')
    # ax2.scatter(notch[0], notch[1], notch[2])
    # # ax.plot_surface(X1, Y1, V.reshape((n,m)), cmap='coolwarm', edgecolor='none')
    # ax2.set_xlabel(xlabel='x[cm]')
    # ax2.set_ylabel(ylabel='y[cm]')
    # ax2.set_zlabel(zlabel='temperature[°C]')
    # plt.show()



def test():
    bestIndex6, LSpansEnd6 = prob6()
    Lx = 4
    Ly = 4
    L = 2
    LSpan = [LSpansEnd6[bestIndex6] - L, LSpansEnd6[bestIndex6]]
    delta = 0.1
    H = 0.005
    n = 10
    m = 10
    P = 10
    P_runs = 200
    P_range = np.linspace(-1000, 1000, P_runs)
    K_range = np.linspace(1, 5, 10)

    fig, ax = plt.subplots(1,1)
    for j in range(10):
        K = K_range[j]
        tempMatrix = np.zeros((P_runs, 1))
        for i in range(P_range.size):
            A = createAMatrix(n, m, LSpan, Lx, Ly, H, K, delta)
            B = createBMatrix(n, m, LSpan, Ly, P_range[i], delta, K)
            V = LA.solve(A, B) + 20
            tempMatrix[i] = np.max(V)

        # A = createAMatrix(n, m, LSpan, Lx, Ly, H, K, delta)
        # B = createBMatrix(n, m, L, Lx, P, delta, K)
        # V = LA.solve(A, B) + 20
        # print(np.max(V))

        ax.plot(P_range, tempMatrix, label="K={}".format(K))
        ax.set_xlabel("P")
        ax.set_ylabel("Temp")




    ax.axhline(y=0, color='k', linestyle="--")
    ax.legend()
    plt.show()

    # fig = plt.figure()
    # ax = plt.axes(projection='3d')
    # x = np.linspace(0, Lx, m)
    # y = np.linspace(0, Ly, n)
    # X1, Y1 = np.meshgrid(x, y)
    # ax.plot_surface(X1, Y1, V.reshape((n,m)), cmap='coolwarm', edgecolor='none')
    # ax.set_xlabel(xlabel='x[cm]')
    # ax.set_ylabel(ylabel='y[cm]')
    # ax.set_zlabel(zlabel='temperature[°C]')
    # plt.show()



if __name__ == "__main__":
    # prob3()
    prob4()
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
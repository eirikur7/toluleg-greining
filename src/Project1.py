import numpy as np
from numpy import linalg as LA
import matplotlib.pyplot as plt
import math

initial = np.array([0,0,6370,0])
correct_position = [0, 0, 6370]

A = np.array([15600, 18760, 17610, 19170])
B = np.array([7540, 2750, 14630, 610])
C = np.array([20140, 18610, 13480, 18390])
t = np.array([0.07074,0.07220,0.07690,0.07242])

c = 299792.458
c2 = c*c
RHO = 26570

ABCT = np.matrix([A, B, C, t])

def F(inArr, ABCT_array):
    x = inArr[0]
    y = inArr[1]
    z = inArr[2]
    d = inArr[3]
    row = ABCT_array.shape[1]
    ret_matrix = np.zeros((row, 1))
    for i in range(row):
        ret_matrix[i] = (x-ABCT_array[0,i])**2 + (y-ABCT_array[1,i])**2 + (z-ABCT_array[2,i])**2 - c2*((ABCT_array[3,i]-d)**2)
    return ret_matrix

def DF(inArr, ABCT_array):
    x= inArr[0]
    y= inArr[1]
    z= inArr[2]
    d= inArr[3]
    row = ABCT_array.shape[1]
    col = ABCT_array.shape[0]
    ret_matrix = np.zeros((row, col))
    for i in range(row):
        ret_matrix[i,0] = 2*(x-ABCT_array[0,i])
        ret_matrix[i,1] = 2*(y-ABCT_array[1,i])
        ret_matrix[i,2] = 2*(z-ABCT_array[2,i])
        ret_matrix[i,3] = 2*c2*(ABCT_array[3,i]-d)
    return ret_matrix

def calc_pos_error(wrongPos, dimentional=2, correctPos=correct_position):
    if dimentional == 2:
        return np.sqrt((correctPos[0] - wrongPos[0,0])**2 + (correctPos[1] - wrongPos[1,0])**2 + (correctPos[2] - wrongPos[2,0])**2)
    elif dimentional == 1:
        return np.sqrt((correctPos[0] - wrongPos[0])**2 + (correctPos[1] - wrongPos[1])**2 + (correctPos[2] - wrongPos[2])**2)
#----------------------ANALYSIS METHODS----------------------#
def newtonMethod(x0, ABCT_arr, tol, cap=None):
    x0 = np.reshape(x0, (x0.size, 1))
    x=x0
    oldx =x + 2*tol
    iterations = 0
    while LA.norm(x-oldx, np.inf) > tol:
        if((cap != None) and (cap <= iterations)):
            return np.empty(shape=(0))
        oldx=x
        s=LA.solve(DF(x, ABCT_arr),F(x, ABCT_arr))
        x=x-s
        iterations += 1
    return(x)

def bisection(f,a,b,tol):
    '''gert ráð fyrir að búið se að skilgreina f(x) fyrir utan t.d.
    def f(x):
        return(x**2-2)
    '''
    if f(a)*f(b) >= 0:
        print("Bisection method fails.")
        return None
    else:
        fa=f(a)
        while (b-a)/2>tol:
            c=(a+b)/2
            fc=f(c)
            if fc==0:break
            if fc*fa<0:
                b=c
            else:
                a=c
                fa=fc
    return((a+b)/2)

def gaussNewton(x0, ABCT_arr, tol, cap=None):
    x = np.matrix.reshape(x0, (4,1))
    oldx = x + 2*tol
    iterations = 0
    while LA.norm(x-oldx, np.inf) > tol:
        if((cap != None) and (cap <= iterations)):
            return np.empty(shape=(0))
        oldx = x
        jacobi = DF(x, ABCT_arr)
        jacobiTrans = np.matrix.transpose(jacobi)
        f = F(x, ABCT_arr)
        # fTrans = np.matrix.transpose(f)
        DfTransF = np.matmul(jacobiTrans, f)
        DfTransDf = np.matmul(jacobiTrans, jacobi)
        s = LA.solve(DfTransDf, DfTransF)
        x = x-s
        iterations += 1
    return x


#----------------------PROBLEMS----------------------#

def prob1(x0, ABCT_arr, tol):
    x = newtonMethod(x0, ABCT_arr, tol)
    print("Problem 1")
    print("x = {:.6f}, y = {:.6f}, z = {:.6f}, d = {:.6e}".format(x[0, 0], x[1, 0], x[2, 0], x[3, 0]))
    print("-"*55)
    return x


def find_abc(phi, theta):
    '''This really is the requirements for problem 2'''
    new_a, new_b, new_c, new_t = np.zeros(phi.size), np.zeros(phi.size), np.zeros(phi.size), np.zeros(phi.size)
    for i in range(phi.size):
        new_a[i] = RHO*np.sin(phi[i])*np.cos(theta[i])
        new_b[i] = RHO*np.sin(phi[i])*np.sin(theta[i])
        new_c[i] = RHO*np.cos(phi[i])
        new_t[i] = (np.sqrt((new_a[i]-initial[0])**2 + (new_b[i]-initial[1])**2+(new_c[i]-initial[2])**2))/c
    return new_a, new_b, new_c, new_t

def prob3():
    phi_orig = np.array([(np.pi)/8,(np.pi)/6,(3*(np.pi))/8,(np.pi)/4])
    theta_orig = np.array([-(np.pi)/4,(np.pi)/2,(2*(np.pi))/3,((np.pi))/6])
    phi_w_errors = np.array([((np.pi)/8)+(10**(-8)),((np.pi)/6)+(10**(-8)),((3*(np.pi))/8)-(10**(-8)),((np.pi)/4)-(10**(-8))])
    corr_a, corr_b, corr_c, corr_t = find_abc(phi_orig, theta_orig)
    incorr_a, incorr_b, incorr_c, incorr_t = find_abc(phi_w_errors, theta_orig)
    incorr_x = newtonMethod(initial, np.array([incorr_a, incorr_b, incorr_c, corr_t]), 10**(-8))
    distance_error = calc_pos_error(wrongPos=incorr_x)
    print("Problem 3")
    print("WrongPos = ({},{},{}) , error = {:.4e} km".format(incorr_x[1,0],incorr_x[1,0],incorr_x[2,0], distance_error))
    print("-"*55)


def prob4():
    phi_orig = np.array([(np.pi)/8,(np.pi)/6,(3*(np.pi))/8,(np.pi)/4])
    theta_orig = np.array([-(np.pi)/4,(np.pi)/2,(2*(np.pi))/3,((np.pi))/6])
    new_a, new_b, new_c, correct_t = find_abc(phi_orig, theta_orig)
    ABCT_arr = np.matrix([new_a, new_b, new_c, correct_t])
    # correct_values = newtonMethod(initial, ABCT_arr, 10**(-6))
    incorrect_values = np.zeros((16, 4))
    check = 15
    for i in range(16):
        phi_temp = phi_orig.copy()
        for ii in range(phi_orig.size):
            if((check>>ii)&1) == 1:
                phi_temp[ii] += 10**(-8)
            else:
                phi_temp[ii] += -10**(-8)

        check -= 1
        new_a, new_b, new_c, new_t = find_abc(phi_temp, theta_orig)
        ABCT_arr = np.matrix([new_a, new_b, new_c, correct_t])
        incorrect_values[i] = newtonMethod(initial, ABCT_arr, 10**(-6))[:, 0]

    # correct_distance = np.sqrt((correct_values[0]**2) + (correct_values[1]**2) + (correct_values[2]**2))
    error_cur = calc_pos_error(wrongPos=incorrect_values[0], dimentional=1)
    error_cur_index = 0
    # np.sqrt( (incorrect_values[0,0]**2) + (incorrect_values[0,1]**2) + (incorrect_values[0,2]**2) )
    for i in range(1, 16):
        error = calc_pos_error(wrongPos=incorrect_values[i], dimentional=1)
        if(error > error_cur ):
            error_cur = error
            error_cur_index = i

    print("Problem 4:")
    print("WrongPos=({:.2f},{:.2f},{:.2f}), Error={:.4e} km".format(
        incorrect_values[error_cur_index][0], incorrect_values[error_cur_index][1], incorrect_values[error_cur_index][2], error_cur))
    print("-"*55)

def prob5():

    total_diff_array = np.zeros(9)
    error_array = np.zeros(9)
    for i in range(9):
        diffence = 10**(-i)
        diffence2 = 2*10**(-8)
        theta_close = np.array([-((np.pi)/4)-diffence,-((np.pi)/4)+diffence,-((np.pi)/4)-diffence,-((np.pi)/4)+diffence])
        phi_close = np.array([((np.pi)/8)-diffence-diffence2,((np.pi)/8)-diffence+diffence2,((np.pi)/8)+diffence-diffence2,((np.pi)/8)+diffence+diffence2])
        close_a,close_b,close_c,close_t = find_abc(phi_close, theta_close)
        close_x = newtonMethod(initial, np.array([close_a, close_b, close_c, close_t]), 10**(-8))

        total_diff_array[i] = diffence
        error_array[i] = calc_pos_error(close_x)


    # total_diff_array = np.log10(total_diff_array)
    plt.plot(np.flip(np.log10(total_diff_array)*-1), np.flip(error_array), 'ro')
    plt.ylabel('error[km]')
    plt.xlabel('angle difference[10^(-x)]')
    plt.title('Problem 5')
    
    plt.show()


def randomAnglesError(sets, error, nrSatilites):
    i = 0
    measuring_array = np.zeros((sets))
    angles_arry = np.zeros((sets, nrSatilites*2))
    while i < sets:
        phi_arr = np.random.uniform(low=0, high=np.pi/2, size=(nrSatilites))
        theta_arr = np.random.uniform(low=0, high=2*np.pi, size=(nrSatilites))
        new_a, new_b, new_c, original_t = find_abc(phi_arr, theta_arr)
        ABCT_arr = np.matrix([new_a, new_b, new_c, original_t])

        new_a, new_b, new_c, new_t = find_abc(phi_arr + error, theta_arr + error)
        ABCT_arr = np.matrix([new_a, new_b, new_c, original_t])
        incorrect_xyzd = gaussNewton(initial, ABCT_arr, 10e-8, 50)
        if(incorrect_xyzd.size > 0):
            measuring_array[i] = calc_pos_error(incorrect_xyzd)
            angles_arry[i][0:nrSatilites] = phi_arr
            angles_arry[i][nrSatilites:] = theta_arr
            i += 1

    return measuring_array, angles_arry

def prob6():
    measuring_array, angles_arry = randomAnglesError(10000, 10e-8, 4)
    # if running_once:
    print("PROBLEM 6:")
    print("Error: max={:.4e}, min={:.4e}, mean={:.4e}, median={:.4e}, std={:.4e}".format(np.max(measuring_array), np.min(measuring_array), np.average(measuring_array), np.median(measuring_array), np.std(measuring_array)))
    max_index = np.argmax(measuring_array)
    # print(angles_arry[max_index][0:4] - (np.pi/2))
    # print(angles_arry[max_index][4:]  - (2*np.pi))
    plt.hist(measuring_array, 500, range=(0, 0.01))
    plt.show()
    print("-"*55)
    # return np.max(measuring_array)

def prob7():
    print(bisection(f, 1e-14, 1e-8, 10**(-6)))
    print("Problem 7")
    print("-"*55)

def f(y):
    measuring_error, angles = randomAnglesError(sets=100, error=y, nrSatilites=4)
    return np.max(measuring_error) - 0.0001

def printLocation(ABCT_arr):
    str = "x: {:.2f}km\ny: {:.2f}km\nz: {:.2f}km\nd: {:.2e}km"
    try:
        print(str.format(ABCT_arr[0], ABCT_arr[1], ABCT_arr[2], ABCT_arr[3]))
    except TypeError:
        print(str.format(ABCT_arr[0,0], ABCT_arr[1,0], ABCT_arr[2,0], ABCT_arr[3,0]))

def prob8():
    measuring_array, angles_arry = randomAnglesError(10000, 10e-8, 5)
    print("PROBLEM 8:")
    print("Error: max={:.4e}, min={:.4e}, mean={:.4e}, median={:.4e}, std={:.4e}".format(np.max(measuring_array), np.min(measuring_array), np.average(measuring_array), np.median(measuring_array), np.std(measuring_array)))
    max_index = np.argmax(measuring_array)
    # print(angles_arry[max_index][0:4] - (np.pi/2))
    # print(angles_arry[max_index][4:]  - (2*np.pi))
    plt.hist(measuring_array, 500, range=(0, 0.01))
    plt.show()
    print("-"*55)






# def satelliteSim(noOfIterations):
#     noOfSatellites = 4
#     theta = pi
#     satPhi = [0 for i in range(noOfSatellites)]

#     for i in range(noOfSatellites, noOfSatellites + noOfIterations*2, 2):


    
    
def prob9():
    nrSatilites = [6, 7, 8, 9]
    colors = ["r", "b", "g", "k"]
    sets = 1000
    bins = 100
    allMeasure = np.zeros((len(nrSatilites), sets))
    for i in range(len(nrSatilites)):
        measuring_array, angles_arry = randomAnglesError(sets, 10e-8, nrSatilites[i])
        allMeasure[i] = measuring_array

    hist_range = [0, 0.005]
    fig_rows = 2
    figure, ax = plt.subplots(2, math.ceil(len(nrSatilites)/fig_rows))
    allPlotsFig, allAx = plt.subplots()
    for i in range(len(nrSatilites)):
        ax[int(math.floor(i/fig_rows)), i%fig_rows].hist(allMeasure[i], bins, range=(hist_range[0], hist_range[1]), color=colors[i])
        ax[int(math.floor(i/fig_rows)), i%fig_rows].title.set_text("K="+str(nrSatilites[i]))
        allAx.hist(allMeasure[i], bins, range=(hist_range[0], hist_range[1]), label="k="+str(nrSatilites[i]), color=colors[i])

    print("PROBLEM 9:")
    for i in range(len(nrSatilites)):
        print("Error for k={}: max={:.4e}, min={:.4e}, mean={:.4e}, median={:.4e}, std={:.4e}".format(nrSatilites[i], 
            np.max(allMeasure[i]), np.min(allMeasure[i]), np.average(allMeasure[i]), np.median(allMeasure[i]), np.std(allMeasure[i])))
    print("-"*55)
    allPlotsFig.legend(loc="upper right")
    allPlotsFig.show()
    figure.show()
    plt.show()

def calculateError(pos1, pos2):
    return np.sqrt((pos1[0] - pos2[0])**2 + (pos1[1] - pos2[1])**2 + (pos1[2] - pos2[2])**2)
 
def test_sporbaug():
    theta_arr = np.array([0, 0, math.pi])
    
    phi_arr1 = np.array([math.pi/2, 0, math.pi/2]) 
    new_a, new_b, new_c, original_t = find_abc(phi_arr1, theta_arr)
    ABCT_arr1 = np.matrix([new_a, new_b, new_c, original_t])

    phi_arr2 = np.array([math.pi/4, 0, math.pi/4])
    new_a, new_b, new_c, original_t = find_abc(phi_arr2, theta_arr)
    ABCT_arr2 = np.matrix([new_a, new_b, new_c, original_t])
   
    pos1 = gaussNewton(initial, ABCT_arr1, 10e-8, 50)
    pos2 = gaussNewton(initial, ABCT_arr2, 10e-8, 50)

    print(pos1)
    print()
    print(pos2)
    print()
    print(calculateError(initial, pos1))
    print()
    print(calculateError(initial, pos2))

if __name__ == "__main__":


    # print(F(initial, ABCT))
    # print(DF(initial, ABCT))

    # prob8()
    test_sporbaug()

    # theta_1 = np.array([(np.pi)/8,(np.pi)/6,(3*(np.pi))/8,(np.pi)/4])
    # phi_1 = np.array([-(np.pi)/4,(np.pi)/2,(2*(np.pi))/3,((np.pi))/6])
    # prob1(initial, ABCT, 10e-8)
    # print("Problem 2, finished")
    # print("-"*55)
    # # prob3()
    # # prob4()
    # # prob5()
    # # prob6()
    # # prob7()
    # # prob8()
    # # prob9()
    # theta_2 = np.array([((np.pi)/8)+(10**(-8)),((np.pi)/6)+(10**(-8)),((3*(np.pi))/8)-(10**(-8)),((np.pi)/4)-(10**(-8))])
    # new_a, new_b, new_c, new_t = find_abc(theta_1, phi_1)

    # new_a2,new_b2,new_c2,new_t2 = find_abc(theta_2, phi_1)
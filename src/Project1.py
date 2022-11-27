import numpy as np
from numpy import linalg as LA
import matplotlib.pyplot as plt
import math

initial = np.array([0,0,6370,0]) #x,y,z,d -> km, km, km, km
correct_position = [0, 0, 6370] #x, y, z -> km, km, km

A = np.array([15600, 18760, 17610, 19170])
B = np.array([7540, 2750, 14630, 610])
C = np.array([20140, 18610, 13480, 18390])
t = np.array([0.07074,0.07220,0.07690,0.07242])

c = 299792.458 #km/s
c2 = c*c # c squared
RHO = 26570 #km

ABCT = np.matrix([A, B, C, t])

def F(inArr, ABCT_array): 
    '''F(x,y,z,d) = [f1,f2,f3,f4]'''
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
    ''' DF is the Jacobian matrix of F(x,y,z,d)'''
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
    ''' This function takes in two arguments, the correct position and the wrong position. 
    It then calculates the distance between the two positions using the Pythagorean theorem.
    It then returns the distance between the two positions.'''

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
    _, _, _, corr_t = find_abc(phi_orig, theta_orig)
    incorr_a, incorr_b, incorr_c, _ = find_abc(phi_w_errors, theta_orig)
    incorr_x = newtonMethod(initial, np.array([incorr_a, incorr_b, incorr_c, corr_t]), 10**(-8))
    distance_error = calc_pos_error(wrongPos=incorr_x)
    print("Problem 3")
    print("WrongPos = ({},{},{}) , error = {:.4e} km".format(incorr_x[1,0],incorr_x[1,0],incorr_x[2,0], distance_error))
    print("-"*55)


def prob4():
    phi_orig = np.array([(np.pi)/8,(np.pi)/6,(3*(np.pi))/8,(np.pi)/4])
    theta_orig = np.array([-(np.pi)/4,(np.pi)/2,(2*(np.pi))/3,((np.pi))/6])
    _, _, _, correct_t = find_abc(phi_orig, theta_orig)
    incorrect_values = np.zeros((16, 4))
    for i in range(15, -1, -1):
        phi_temp = phi_orig.copy()
        for ii in range(phi_orig.size):
            if((i>>ii)&1) == 1:
                phi_temp[ii] = phi_temp[ii] + ((10**(-8)))
            else:
                phi_temp[ii] = phi_temp[ii] - ((10**(-8)))
        new_a, new_b, new_c, _ = find_abc(phi_temp, theta_orig)
        ABCT_arr = np.matrix([new_a, new_b, new_c, correct_t])
        incorrect_values[i] = newtonMethod(initial, ABCT_arr, 10**(-6))[:, 0]

    error_cur = calc_pos_error(wrongPos=incorrect_values[0], dimentional=1)
    error_cur_index = 0
    for i in range(1, 16):
        error = calc_pos_error(wrongPos=incorrect_values[i], dimentional=1)
        if(error > error_cur):
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
    plt.plot(np.flip(np.log10(total_diff_array)*-1), np.flip(error_array), 'ro')
    plt.ylabel('error[km]')
    plt.xlabel('angle difference[10^(-x)]')
    plt.title('Problem 5')
    
    plt.show()


def randomAnglesError(sets, error, nrSatilites, UseNewtomMethod = False):
    i = 0
    measuring_array = np.zeros((sets))
    angles_arry = np.zeros((sets, nrSatilites*2))
    while i < sets:
        phi_arr = np.random.uniform(low=0, high=np.pi/2, size=(nrSatilites))
        theta_arr = np.random.uniform(low=0, high=2*np.pi, size=(nrSatilites))
        _, _, _, original_t = find_abc(phi_arr, theta_arr)
        new_a, new_b, new_c, _ = find_abc(phi_arr+error, theta_arr)
        ABCT_arr = np.matrix([new_a, new_b, new_c, original_t])
        if UseNewtomMethod:
            incorrect_xyzd = newtonMethod(initial, ABCT_arr, 10**(-8), 50)
        else:
            incorrect_xyzd = gaussNewton(initial, ABCT_arr, 10**(-8), 50)
        if(incorrect_xyzd.size > 0):
            measuring_array[i] = calc_pos_error(incorrect_xyzd)
            angles_arry[i][0:nrSatilites] = phi_arr
            angles_arry[i][nrSatilites:] = theta_arr
            i += 1

    return measuring_array, angles_arry

def prob6():
    measuring_array, angles_arry = randomAnglesError(10000, (10**(-8)), 4)
    print("PROBLEM 6:")
    mean = np.mean(measuring_array)
    std = np.std(measuring_array)
    meadian = np.median(measuring_array)
    print("Error: max={:.4e}, min={:.4e}, mean={:.4e}, median={:.4e}, std={:.4e}".format(np.max(measuring_array), np.min(measuring_array), mean, meadian, std))
    fig, ax = plt.subplots(2, 1)
    ax[0].hist(measuring_array, 500,range=(0, 0.001))
    ax[0].set_xlabel("Positional Error[km]")
    ax[0].set_ylabel("Counts")
    ax[0].set_title("99% of Data")
    ax[1].hist(measuring_array, 500, range=(0, (10*meadian)))
    ax[1].set_xlabel("Positional Error[km]")
    ax[1].set_ylabel("Counts")
    ax[1].set_title("Range 0 to 10*median")
    plt.show()

def prob7():
    print(bisection(f, (10**(-14)), (10**(-8)), 10**(-6)))
    print("Problem 7")
    print("-"*55)

def f(y):
    measuring_error, _ = randomAnglesError(sets=100, error=y, nrSatilites=4)
    return np.max(measuring_error) - 0.0001

def printLocation(ABCT_arr):
    str = "x: {:.2f}km\ny: {:.2f}km\nz: {:.2f}km\nd: {:.2e}km"
    try:
        print(str.format(ABCT_arr[0], ABCT_arr[1], ABCT_arr[2], ABCT_arr[3]))
    except TypeError:
        print(str.format(ABCT_arr[0,0], ABCT_arr[1,0], ABCT_arr[2,0], ABCT_arr[3,0]))

def prob8():
    measuring_array, angles_arry = randomAnglesError(10000, (10**(-8)), 5)
    print("PROBLEM 8:")
    print("Error: max={:.4e}, min={:.4e}, mean={:.4e}, median={:.4e}, std={:.4e}".format(np.max(measuring_array), np.min(measuring_array), np.average(measuring_array), np.median(measuring_array), np.std(measuring_array)))
    max_index = np.argmax(measuring_array)
    print(angles_arry[max_index][0:4] - (np.pi/2))
    print(angles_arry[max_index][4:]  - (2*np.pi))
    plt.hist(measuring_array, 500, range=(0, 0.005))
    plt.show()
    print("-"*55)


    
def prob9(sets = 10000):
    nrSatilites = [6, 7, 8, 9]
    colors = ["r", "b", "g", "k"]
    bins = 500
    allMeasure = np.zeros((len(nrSatilites), sets))
    for i in range(len(nrSatilites)):
        measuring_array, _ = randomAnglesError(sets, (10**(-8)), nrSatilites[i])
        allMeasure[i] = measuring_array

    hist_range = [0, 0.0005]
    fig_rows = 2
    figure, ax = plt.subplots(fig_rows, math.ceil(len(nrSatilites)/fig_rows))
    allPlotsFig, allAx = plt.subplots()
    for i in range(len(nrSatilites)):
        ax[int(math.floor(i/fig_rows)), i%fig_rows].hist(allMeasure[i], bins, range=(hist_range[0], hist_range[1]), color=colors[i])
        ax[int(math.floor(i/fig_rows)), i%fig_rows].title.set_text("K="+str(nrSatilites[i]))
        ax[int(math.floor(i/fig_rows)), i%fig_rows].set_ylabel("Counts")
        allAx.hist(allMeasure[i], bins, range=(hist_range[0], hist_range[1]), label="k="+str(nrSatilites[i]), color=colors[i])
    ax[1, 0].set_xlabel("Positional Error[km]")
    ax[1, 1].set_xlabel("Positional Error[km]")
    print("PROBLEM 9:")
    for i in range(len(nrSatilites)):
        print("Error for k={}: max={:.4e}, min={:.4e}, mean={:.4e}, median={:.4e}, std={:.4e}".format(nrSatilites[i],
            np.max(allMeasure[i]), np.min(allMeasure[i]), np.average(allMeasure[i]), np.median(allMeasure[i]), np.std(allMeasure[i])))
    print("-"*55)
    allAx.set_ylabel("Counts")
    allAx.set_xlabel("Positional Error[km]")
    allPlotsFig.legend(loc='center right', bbox_to_anchor=(0.8, 0.5))
    allPlotsFig.show()
    figure.show()
    plt.show()

def calculateError(pos1, pos2):
    return np.sqrt((pos1[0] - pos2[0])**2 + (pos1[1] - pos2[1])**2 + (pos1[2] - pos2[2])**2)
 
def approximateNoOfSatellites():
    print("Simulating 5 to 40 satellites, in each iteration we take X amount of measurements. For the report a 1000 were done for each iteration, however currently it's set to a 100 to speed up the sim. Change the variable 'sets' in the function 'randomAnglesError' to 1000 to get similar results as in the report.")
    print("printing out the current number of satellites to track progress.")
    iterations = list(range(5,40,1))
    sets = 100
    satellite_error = [10**(-8), 10**(-9), 10**(-10), 10**(-11)]
    measurements = np.zeros((len(iterations), 2*len(satellite_error)))
    counter = 0
    for i in iterations:
        temp_arr = []
        for error in satellite_error:
            measure_arr, _ = randomAnglesError(sets, error, i)
            temp_arr += [np.average(measure_arr), np.std(measure_arr)]
        measurements[counter] = temp_arr
        counter += 1
        print(i, end=" ")
    print(i)
    
    plt.plot(iterations, measurements[:,0]*10**5, '.', label=satellite_error[0])
    plt.plot(iterations, measurements[:,2]*10**5, '.', label=satellite_error[1])
    plt.plot(iterations, measurements[:,4]*10**5, '.', label=satellite_error[2])
    plt.plot(iterations, measurements[:,6]*10**5, '.', label=satellite_error[3])
    plt.legend(loc="upper right")
    plt.xlabel("Number of satellites")
    plt.ylabel("Error [cm]")
    plt.show()

    plt.plot(iterations, measurements[:,1]*10**5, '.', label=satellite_error[0])
    plt.plot(iterations, measurements[:,3]*10**5, '.', label=satellite_error[1])
    plt.plot(iterations, measurements[:,5]*10**5, '.', label=satellite_error[2])
    plt.plot(iterations, measurements[:,7]*10**5, '.', label=satellite_error[3])
    plt.legend(loc="upper right")
    plt.xlabel("Number of satellites")
    plt.ylabel("Error [cm]")
    plt.show()


def calcError(original_t,phi,theta,nr_sattelites=4):
    check = 2**nr_sattelites - 1

    incorrect_values_phi = np.zeros((check+1, 4))

    incorrect_values_theta = np.zeros((check+1, 4))

    incorrect_values_both = np.zeros((check+1, 4))
    
    for i in range(check+1):
        phi_temp = phi.copy()
        theta_temp = theta.copy()
        for ii in range(phi.size):
            if((check>>ii)&1) == 1:
                phi_temp[ii] += 10**(-8)
                theta_temp[ii] += 10**(-8)
            else:
                phi_temp[ii] += -10**(-8)
                theta_temp[ii] += -10**(-8)

        check -= 1
        #Phi
        new_a, new_b, new_c, _ = find_abc(phi_temp, theta)
        ABCT_arr = np.matrix([new_a, new_b, new_c, original_t])
        incorrect_values_phi[i] = newtonMethod(initial, ABCT_arr, 10**(-6))[:, 0]

        #Theta
        new_a, new_b, new_c, _ = find_abc(phi, theta_temp)
        ABCT_arr = np.matrix([new_a, new_b, new_c, original_t])
        incorrect_values_theta[i] = newtonMethod(initial, ABCT_arr, 10**(-6))[:, 0]

        #both
        new_a, new_b, new_c, _ = find_abc(phi_temp, theta_temp)
        ABCT_arr = np.matrix([new_a, new_b, new_c, original_t])
        incorrect_values_both[i] = newtonMethod(initial, ABCT_arr, 10**(-6))[:, 0]


    #phi
    error_cur_phi = calc_pos_error(wrongPos=incorrect_values_phi[0], dimentional=1)
    error_cur_index_phi = 0
    #Theta
    error_cur_theta = calc_pos_error(wrongPos=incorrect_values_theta[0], dimentional=1)
    error_cur_index_theta = 0

    #both
    error_cur_both = calc_pos_error(wrongPos=incorrect_values_both[0], dimentional=1)
    error_cur_index_both = 0

    for i in range(1, incorrect_values_phi.shape[0]):

        error_phi = calc_pos_error(wrongPos=incorrect_values_phi[i], dimentional=1)
        error_theta = calc_pos_error(wrongPos=incorrect_values_theta[i], dimentional=1)
        error_both = calc_pos_error(wrongPos=incorrect_values_both[i], dimentional=1)

        if(error_phi > error_cur_phi ):
            error_cur_phi = error_phi
            error_cur_index_phi = i

        if(error_theta > error_cur_theta ):
            error_cur_theta = error_theta
            error_cur_index_theta = i

        if(error_both > error_cur_both ):
            error_cur_both = error_both
            error_cur_index_both = i

    return incorrect_values_phi[error_cur_index_phi], incorrect_values_theta[error_cur_index_theta], incorrect_values_both[error_cur_index_both]

def prob10_2(sets=1000,nrSatilites=4):
    i = 0
    error_arr_phi = np.zeros((sets, 1))
    error_arr_theta = np.zeros((sets, 1))
    error_arr_both = np.zeros((sets, 1))

    while i < sets:
        phi_arr = np.random.uniform(low=0, high=np.pi/2, size=(nrSatilites))
        theta_arr = np.random.uniform(low=0, high=2*np.pi, size=(nrSatilites))
        
        _, _, _, original_t = find_abc(phi_arr, theta_arr)
        
        errored_phi,errored_theta,errored_both = calcError(original_t,phi_arr,theta_arr,nr_sattelites=nrSatilites)

        error_phi = calc_pos_error(wrongPos=errored_phi, dimentional=1)
        error_theta = calc_pos_error(wrongPos=errored_theta, dimentional=1)
        error_both = calc_pos_error(wrongPos=errored_both, dimentional=1)

        error_arr_phi[i] = error_phi
        error_arr_theta[i] = error_theta
        error_arr_both[i] = error_both

        i += 1
        #subplot
     
    print("Mean error for phi: ",error_arr_phi.mean())
    print("Mean error for theta:", error_arr_theta.mean())
    print("Mean error for both:", error_arr_both.mean())
    plt.subplot(1, 3, 1)
    plt.title("Phi")
    plt.hist(error_arr_phi, bins=100, range=(0, 0.005))
    plt.subplot(1, 3, 2)
    plt.title("Theta")
    plt.hist(error_arr_theta, bins=100,range=(0, 0.0000000005))
    plt.subplot(1, 3, 3)
    plt.title("Both")
    plt.hist(error_arr_both, bins=100,range=(0, 0.005))
    plt.show()
    


if __name__ == "__main__":
    user_input = input("What problem do you want to run? (1-10) ")
    allowed = ["1","2","3","4","5","6","7","8","9","10"]
    while user_input in allowed:
        if user_input == "1":
            prob1(initial, ABCT, 10**(-6))
        elif user_input == "2":
            print("Problem 2, finished (see problem 3)")
        elif user_input == "3":
            prob3()
        elif user_input == "4":
            prob4()
        elif user_input == "5":
            prob5()
        elif user_input == "6":
            prob6()
        elif user_input == "7":
            prob7()
        elif user_input == "8":
            prob8()
        elif user_input == "9":
            prob9(sets=1000)
        elif user_input == "10":
            user_input2 = input("Do you want to run problem 10.1 or 10.2? (1-2) ")
            if user_input2 == "1":
                approximateNoOfSatellites()
            elif user_input2 == "2":
                prob10_2()
        user_input = input("What problem do you want to run? (1-10) ")
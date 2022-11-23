import numpy as np
from numpy import linalg as LA
import matplotlib.pyplot as plt

initial = np.array([0,0,6370,0])

A = np.array([15600, 18760, 17610, 19170])
B = np.array([7540, 2750, 14630, 610])
C = np.array([20140, 18610, 13480, 18390])
t = np.array([0.07074,0.07220,0.07690,0.07242])
c = 299792.458
c2 = c*c
ABCT = np.matrix([A, B, C, t])

def F(inArr, ABCT_array):
    x = inArr[0]
    y = inArr[1]
    z = inArr[2]
    d = inArr[3]
    f1 = (x-ABCT_array[0,0])**2 + (y-ABCT_array[1,0])**2 + (z-ABCT_array[2,0])**2 - c2*((ABCT_array[3,0]-d)**2)
    f2 = (x-ABCT_array[0,1])**2 + (y-ABCT_array[1,1])**2 + (z-ABCT_array[2,1])**2 - c2*((ABCT_array[3,1]-d)**2)
    f3 = (x-ABCT_array[0,2])**2 + (y-ABCT_array[1,2])**2 + (z-ABCT_array[2,2])**2 - c2*((ABCT_array[3,2]-d)**2)
    f4 = (x-ABCT_array[0,3])**2 + (y-ABCT_array[1,3])**2 + (z-ABCT_array[2,3])**2 - c2*((ABCT_array[3,3]-d)**2)
    return np.array([f1,f2,f3,f4])

def DF(inArr, ABCT_array):
    x= inArr[0]
    y= inArr[1]
    z= inArr[2]
    d= inArr[3]
    return np.matrix([[2*(x-ABCT_array[0,0]), 2*(y-ABCT_array[1,0]), 2*(z-ABCT_array[2,0]), 2*c2*(ABCT_array[3,0]-d)],
                      [2*(x-ABCT_array[0,1]), 2*(y-ABCT_array[1,1]), 2*(z-ABCT_array[2,1]), 2*c2*(ABCT_array[3,1]-d)],
                      [2*(x-ABCT_array[0,2]), 2*(y-ABCT_array[1,2]), 2*(z-ABCT_array[2,2]), 2*c2*(ABCT_array[3,2]-d)], 
                      [2*(x-ABCT_array[0,3]), 2*(y-ABCT_array[1,3]), 2*(z-ABCT_array[2,3]), 2*c2*(ABCT_array[3,3]-d)]])

def newtonMethod(x0, ABCT_arr, tol, cap=None):
    x=x0
    oldx =x + 2*tol
    cnt = 0
    while LA.norm(x-oldx, np.inf) > tol:
        if((cap != None) and (cap <= cnt)):
            return np.empty(shape=(0))
        oldx=x
        s=LA.solve(DF(x, ABCT_arr),F(x, ABCT_arr))
        x=x-s
        cnt += 1
    return(x)

def prob1(x0, ABCT_arr, tol):
    x = newtonMethod(x0, ABCT_arr, tol)
    print("x={}, y={}, z={}, d={}".format(x[0], x[1], x[2], x[3]))
    return x
    

def find_abc(phi, theta): 
    rho = 26570
    new_a, new_b, new_c, new_t = np.zeros(phi.size), np.zeros(phi.size), np.zeros(phi.size), np.zeros(phi.size)
    for i in range(phi.size):
        new_a[i] = rho*np.sin(phi[i])*np.cos(theta[i])
        new_b[i] = rho*np.sin(phi[i])*np.sin(theta[i])
        new_c[i] = rho*np.cos(phi[i])
        new_t[i] = (np.sqrt((new_a[i]-initial[0])**2 + (new_b[i]-initial[1])**2+(new_c[i]-initial[2])**2))/c
    return new_a, new_b, new_c, new_t



def prob4():
    phi_orig = np.array([(np.pi)/8,(np.pi)/6,(3*(np.pi))/8,(np.pi)/4])
    theta_orig = np.array([-(np.pi)/4,(np.pi)/2,(2*(np.pi))/3,((np.pi))/6])
    new_a, new_b, new_c, correct_t = find_abc(phi_orig, theta_orig)
    ABCT_arr = np.matrix([new_a, new_b, new_c, correct_t])
    correct_values = newtonMethod(initial, ABCT_arr, 10**(-8))
    # print(correct_values)
    incorrect_values = np.zeros((16, 4))
    for i in range(16): 
        phi_temp = np.copy(phi_orig)
        for ii in range(phi_orig.size):
            print("i={}, ii={}, &={}".format(i, ii, i&(1<<ii)))
            if((i&(1<<ii)) > 0):
                phi_temp[ii] += 10**(-8)
            else:
                phi_temp[ii] -= 10**(-8)
        
        print(phi_temp)
        new_a, new_b, new_c, new_t = find_abc(phi_temp, theta_orig)
        ABCT_arr = np.matrix([new_a, new_b, new_c, correct_t])
        incorrect_values[i] = newtonMethod(initial, ABCT_arr, 10**(-8))
        #incorrect_values[i] = new_t
    print("correct")
    print(correct_values)
    print("incorect")
    print(incorrect_values)
    correct_distance = np.sqrt((correct_values[0]**2) + (correct_values[1]**2) + (correct_values[2]**2))
    max_row = 0
    incorect_distance_cur = np.sqrt( (incorrect_values[0,0]**2) + (incorrect_values[0,1]**2) + (incorrect_values[0,2]**2) )
    # max_column = 0
    for i in range(1, 16):
        # for ii in range(correct_values.size):
            # if( (abs(correct_values[ii] - incorrect_values[i,ii])) > (abs(correct_values[ii] - incorrect_values[max_row,max_column])) ):
            #    max_row = i
            #    max_column = ii 
        incorect_distance = np.sqrt( (incorrect_values[i,0]**2) + (incorrect_values[i,1]**2) + (incorrect_values[i,2]**2) )
        if( abs(incorect_distance - correct_distance ) > abs(incorect_distance_cur - correct_distance) ):
            max_row = i
            incorect_distance_cur = incorect_distance

    print("PROBLEM4:")
    # print("-------------------------------------------------------")
    print("-"*55)
    prob4SolError = abs(correct_distance - incorect_distance_cur)
    print("row={}, Real={}, Wrong={} Error={}, Percentage={}%".format(max_row, correct_distance, incorect_distance_cur, prob4SolError, 100*(prob4SolError/correct_distance)))
    # print("----------------------------------------------------")
    print("-"*55)
    
    return incorrect_values, correct_values


def prob6(sets, bins, error):
    i = 0
    measuring_array = np.zeros((sets))
    angles_arry = np.zeros((sets, 8))
    while i < sets:
        phi_arr = np.random.uniform(low=0, high=np.pi/2, size=(4))
        theta_arr = np.random.uniform(low=0, high=2*np.pi, size=(4))
        new_a, new_b, new_c, original_t = find_abc(phi_arr, theta_arr)
        ABCT_arr = np.matrix([new_a, new_b, new_c, original_t])
        correct_xyzd = newtonMethod(initial, ABCT_arr, 10e-8, 50)
        
        if(correct_xyzd.size > 0):
            new_a, new_b, new_c, new_t = find_abc(phi_arr + error, theta_arr + error)
            ABCT_arr = np.matrix([new_a, new_b, new_c, original_t])
            incorrect_xyzd = newtonMethod(initial, ABCT_arr, 10e-8, 50)
            if(incorrect_xyzd.size > 0):
                correct_distance   = np.sqrt((correct_xyzd[0]**2) + (correct_xyzd[1]**2) + (correct_xyzd[2]**2))
                incorrect_distance = np.sqrt((incorrect_xyzd[0]**2) + (incorrect_xyzd[1]**2) + (incorrect_xyzd[2]**2))
                measuring_array[i] = abs(correct_distance - incorrect_distance)
                angles_arry[i][0:4] = phi_arr
                angles_arry[i][4:] = theta_arr
                i += 1
    
    
    print("PROBLEM 6:")
    print("-"*55)
    # print(measuring_array)
    # print(angles_arry)
    
    print("Error: max={}, min={}, average={}, std={}".format(np.max(measuring_array), np.min(measuring_array), np.average(measuring_array), np.std(measuring_array)))
    max_index = np.argmax(measuring_array)
    print(angles_arry[max_index][0:4] - (np.pi/2))
    print(angles_arry[max_index][4:]  - (2*np.pi))
    plt.hist(measuring_array, bins, range=(0, 0.01))
    plt.show()
    print("-"*55)


def prob6_2(sets, bins, phi_original, theta_original):
    new_a, new_b, new_c, original_t = find_abc(phi_original, theta_original)
    ABCT_arr = np.matrix([new_a, new_b, new_c, original_t])
    correct_xyzd = newtonMethod(initial, ABCT_arr, 10e-8)
    correct_distance = np.sqrt((correct_xyzd[0]**2) + (correct_xyzd[1]**2) + (correct_xyzd[2]**2))
    measuring_array = np.zeros((sets))
    # theta_arr_all = np.random.uniform(low=0, high=2*np.pi, size=(4*sets))
    # phi_arr_all = np.random.uniform(low=0, high=np.pi/2, size=(4*sets))
    # print("Phi={}, theta={}".format(theta_arr_all, phi_arr_all))
    i = 0
    print(correct_distance)
    while i < sets:
        theta_arr = np.random.uniform(low=0, high=2*np.pi, size=(4))#np.copy(theta_arr_all[i:i+4])
        # print(theta_arr_all[i:i+4])
        phi_arr = np.random.uniform(low=0, high=2*np.pi, size=(4)) #np.array(phi_arr_all[i:i+4])
        # print("Phi={}, theta={}".format(phi_arr, theta_arr))
        new_a, new_b, new_c, new_t = find_abc(phi_arr, theta_arr) #randomize just phi
        ABCT_arr = np.matrix([new_a, new_b, new_c, original_t])
        incorect_xyzd = newtonMethod(initial, ABCT_arr, 10e-8, 100)
        if(incorect_xyzd.size > 0):
            print(incorect_xyzd)
            incorect_distance = np.sqrt((incorect_xyzd[0]**2) + (incorect_xyzd[1]**2) + (incorect_xyzd[2]**2))
            measuring_array[i] = abs(incorect_distance - correct_distance)
            i+=1

    # counts, binsFound = np.histogram(measuring_array)
    # figure = plt.figure()
    # ax = figure.add_subplot(111)
    # ax.hist(binsFound[:-1], binsFound, weights=counts)
    # plt.show()
    # counts, binsFound = np.histogram(measuring_array)
    # plt.hist(binsFound[:-1], binsFound, weights=counts)
    # print(measuring_array)
    print("PROBLEM 6:")
    print("-"*55)
    print(measuring_array)
    print("Error: max={}, min={}, average={}, std={}".format(np.max(measuring_array), np.min(measuring_array), np.average(measuring_array), np.std(measuring_array)))
    print("-"*55)



    plt.hist(measuring_array, bins, range=(np.min(measuring_array), np.max(measuring_array)))
    plt.show()
    
def prob6_1(sets, phi_original, theta_original, minError, maxError):
    #Calculates with Original Phi and Theta and find the correct distance
    new_a, new_b, new_c, original_t = find_abc(phi_original, theta_original)
    ABCT_arr = np.matrix([new_a, new_b, new_c, original_t])
    correct_xyzd = newtonMethod(initial, ABCT_arr, 10e-8)
    correct_distance = np.sqrt((correct_xyzd[0]**2) + (correct_xyzd[1]**2) + (correct_xyzd[2]**2))

    measuring_array = np.zeros((2, sets))
    step = (maxError-minError)/sets
    for i in range(sets):
        phi_arr = np.copy(phi_original) + (step*i + minError)#np.random.uniform(low=0, high=np.pi/2, size=(4))
        # theta_arr = np.random.uniform(low=0, high=2*np.pi, size=(4))
        new_a, new_b, new_c, new_t = find_abc(phi_arr, theta_original)
        ABCT_arr = np.matrix([new_a, new_b, new_c, original_t])
        incorect_xyzd = newtonMethod(initial, ABCT_arr, 10e-8)
        incorect_distance = np.sqrt((incorect_xyzd[0]**2) + (incorect_xyzd[1]**2) + (incorect_xyzd[2]**2))
        measuring_array[0][i] = (step*i + minError)
        measuring_array[1][i] = incorect_distance
    
    # print(measuring_array)
    # x = measuring_array[0]
    # y = measuring_array[1]
    # print(y)
    # print(np.min(y))
    fig = plt.figure()
    ax = fig.add_subplot(111)
    ax.plot(measuring_array[0], abs(measuring_array[1]-correct_distance))
    ax.set_xlim(minError, maxError)
    # ax.set_ylim(np.min(measuring_array[1])*0.95, np.max(measuring_array[1])*1.05)
    plt.show()
if __name__ == "__main__":
    # print(F(initial, ABCT))
    # print(DF(initial, ABCT))
    # prob1(initial, ABCT, 10e-8)
    theta = np.array([(np.pi)/8,(np.pi)/6,(3*(np.pi))/8,(np.pi)/4])
    phi = np.array([-(np.pi)/4,(np.pi)/2,(2*(np.pi))/3,((np.pi))/6])
    prob6(10000, 500, 10e-8)
    # prob6_2(100, 10, phi, theta)
    #prob6_1(1000, phi, theta, -10e-3, 10e-3)
    
    # # print("PROBLEM4:")
    # # print("----------------------------------------------------")
    # # prob4SolArr, prob4SolIndex, prob4CorrectT = prob4()
    # prob4()
    # # print(prob4SolArr[prob4SolIndex])
    # # print(prob4SolIndex)
    # # prob4SolError = abs(prob4CorrectT - prob4SolArr[prob4SolIndex][3])
    # # print("Real={}, Wrong={} Error={}, Percentage={}%".format(prob4CorrectT, prob4SolArr[prob4SolIndex][3], prob4SolError, prob4SolError/prob4CorrectT))
    # # print("----------------------------------------------------")
    
    # theta_2 = np.array([((np.pi)/8)+(10**(-8)),((np.pi)/6)+(10**(-8)),((3*(np.pi))/8)-(10**(-8)),((np.pi)/4)-(10**(-8))])
    # new_a, new_b, new_c, new_t = find_abc(theta_1, phi_1)

    # new_a2,new_b2,new_c2,new_t2 = find_abc(theta_2, phi_1)

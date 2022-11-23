import numpy as np
from numpy import linalg as LA

initial = np.array([0,0,6370,0])

A = np.array([15600, 18760, 17610, 19170])
B = np.array([7540, 2750, 14630, 610])
C = np.array([20140, 18610, 13480, 18390])
t = np.array([0.07074,0.07220,0.07690,0.07242])
c = 299792.458
c2 = c*c
ABCT = np.matrix([A, B, C, t])
RHO = 26570

def F(inArr, ABCT_Matrix):
    x = inArr[0]
    y = inArr[1]
    z = inArr[2]
    d = inArr[3]
    f1 = (x-ABCT_Matrix[0,0])**2 + (y-ABCT_Matrix[1,0])**2 + (z-ABCT_Matrix[2,0])**2 - c2*((ABCT_Matrix[3,0]-d)**2)
    f2 = (x-ABCT_Matrix[0,1])**2 + (y-ABCT_Matrix[1,1])**2 + (z-ABCT_Matrix[2,1])**2 - c2*((ABCT_Matrix[3,1]-d)**2)
    f3 = (x-ABCT_Matrix[0,2])**2 + (y-ABCT_Matrix[1,2])**2 + (z-ABCT_Matrix[2,2])**2 - c2*((ABCT_Matrix[3,2]-d)**2)
    f4 = (x-ABCT_Matrix[0,3])**2 + (y-ABCT_Matrix[1,3])**2 + (z-ABCT_Matrix[2,3])**2 - c2*((ABCT_Matrix[3,3]-d)**2)
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

def newtonMethod(x0, ABCT_arr, tol):
    x=x0
    oldx =x + 2*tol
    while LA.norm(x-oldx, np.inf) > tol:
        oldx=x
        s=LA.solve(DF(x, ABCT_arr),F(x, ABCT_arr))
        x=x-s
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

def prob1(x0, ABCT_arr, tol):
    x = newtonMethod(x0, ABCT_arr, tol)
    print("Problem 1")
    print("x = {:.2f}, y = {:.2f}, z = {:.2f}, d = {:.2e}".format(x[0], x[1], x[2], x[3]))
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
    corr_x = newtonMethod(initial,  np.array([corr_a, corr_b, corr_c, corr_t]),10**(-8))
    incorr_x = newtonMethod(initial, np.array([incorr_a, incorr_b, incorr_c, corr_t]), 10**(-8))
    distance_incorr_x = np.sqrt(incorr_x[0]**2 + incorr_x[1]**2 + incorr_x[2]**2)
    distance_corr_x = np.sqrt(corr_x[0]**2 + corr_x[1]**2 + corr_x[2]**2)
    print("Problem 3")
    print("Real = {:f} km, Wrong = {:f} , error = {:f} km".format(distance_corr_x,distance_incorr_x,abs(distance_incorr_x-distance_corr_x)))
    print("-"*55)


def prob4():
    phi_orig = np.array([(np.pi)/8,(np.pi)/6,(3*(np.pi))/8,(np.pi)/4])
    theta_orig = np.array([-(np.pi)/4,(np.pi)/2,(2*(np.pi))/3,((np.pi))/6])
    new_a, new_b, new_c, correct_t = find_abc(phi_orig, theta_orig)
    ABCT_arr = np.matrix([new_a, new_b, new_c, correct_t])
    correct_values = newtonMethod(initial, ABCT_arr, 10**(-6))
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
        incorrect_values[i] = newtonMethod(initial, ABCT_arr, 10**(-6))

    correct_distance = np.sqrt((correct_values[0]**2) + (correct_values[1]**2) + (correct_values[2]**2))
    incorrect_distance_cur = np.sqrt( (incorrect_values[0,0]**2) + (incorrect_values[0,1]**2) + (incorrect_values[0,2]**2) )
    incorr_row = 0
    for i in range(1, 16):
        incorrect_distance = np.sqrt( (incorrect_values[i,0]**2) + (incorrect_values[i,1]**2) + (((incorrect_values[i,2]))**2) )
        if( abs(incorrect_distance - correct_distance ) > abs(incorrect_distance_cur - correct_distance) ):
            incorrect_distance_cur = incorrect_distance
            incorr_row = i

    print("Problem 4:")
    prob4SolError = abs(correct_distance - incorrect_distance_cur)
    x  = incorrect_values[incorr_row]
    print("correct: x = {:f}, y = {:f}, z = {:f}, d = {:f}".format(correct_values[0], correct_values[1], correct_values[2], correct_values[3]))
    print("incorrect: x = {:f}, y = {:f}, z = {:f}, d = {:f}".format(x[0], x[1], x[2], x[3]))
    print("Real={:f} km, Wrong={:f} km, Error={:f} km".format(correct_distance, incorrect_distance_cur, prob4SolError))
    print("-"*55)
    return incorrect_values, correct_values

def prob5():
    diffence = 10**(-8)
    diffence2 = 2*10**(-8)
    phi_orig = np.array([(np.pi)/8,(np.pi)/6,(3*(np.pi))/8,(np.pi)/4])
    theta_orig = np.array([-(np.pi)/4,(np.pi)/2,(2*(np.pi))/3,((np.pi))/6])
    theta = np.array([-((np.pi)/4)-diffence,-((np.pi)/4)+diffence,-((np.pi)/4)-diffence,-((np.pi)/4)+diffence])
    phi = np.array([((np.pi)/8)-diffence-diffence2,((np.pi)/8)-diffence+diffence2,((np.pi)/8)+diffence-diffence2,((np.pi)/8)+diffence+diffence2])

    close_a,close_b,close_c,close_t = find_abc(phi, theta)
    a,b,c,t = find_abc(phi_orig, theta_orig)

    close_x = newtonMethod(initial, np.array([close_a, close_b, close_c, close_t]), 10**(-8))
    x = newtonMethod(initial, np.array([a, b, c, t]), 10**(-8))
    distance_x = np.sqrt(x[0]**2 + x[1]**2 + x[2]**2)
    distance_corr_x = np.sqrt(close_x[0]**2 + close_x[1]**2 + close_x[2]**2)
    print("Problem 5")
    print("error = {:f} km".format(abs(distance_x-distance_corr_x)))
    print("-"*55)

def prob7():
    pass

def prob8():
    pass

def prob9():
    pass

def prob10():
    pass

if __name__ == "__main__":
    
 
    # print(F(initial, ABCT))
    # print(DF(initial, ABCT))
    prob1(initial, ABCT, 10e-8)


    theta_1 = np.array([(np.pi)/8,(np.pi)/6,(3*(np.pi))/8,(np.pi)/4])
    phi_1 = np.array([-(np.pi)/4,(np.pi)/2,(2*(np.pi))/3,((np.pi))/6])

    prob3()
    prob4()
    prob5()
    #prob6()
    theta_2 = np.array([((np.pi)/8)+(10**(-8)),((np.pi)/6)+(10**(-8)),((3*(np.pi))/8)-(10**(-8)),((np.pi)/4)-(10**(-8))])
    new_a, new_b, new_c, new_t = find_abc(theta_1, phi_1)

    new_a2,new_b2,new_c2,new_t2 = find_abc(theta_2, phi_1)

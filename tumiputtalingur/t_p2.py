import numpy as np
import matplotlib.pyplot as plt
from matplotlib.animation import PillowWriter
from matplotlib import animation
g = 9.81
L1 = L2 = 2 #m
m1 = m2 = 1 #kg
ANIMATION_TEXT = "Generating animation..."
PLOT_TEXT = "Generating plot..."
FIG_PATH = "src/figures/"
ANIMATION_PATH = FIG_PATH + "{}_animation.gif"
PLOT_PATH = FIG_PATH + "{}_plot.png"


#/ ------------------------Not related to problems----------------------------------------/
def getInput(print_str, available_inputs):
    """Gets input from user and checks if it is valid"""
    while True:
        user_input = input(print_str)
        if user_input.lower() not in available_inputs:
            print("Invalid input, try again")
        else:
            return user_input

def printFigText(name):
    print(f"I just saved a figure as {name}")

#/-------------------------Math Stuff----------------------------------------/
def F1(theta):
    """Math function for a single pendulum"""
    L = 2
    z1 = theta[1,0]
    z2 = (-g/L)*np.sin(theta[0,0])
    return np.matrix([[z1],[z2]])

def F2(theta):
    """theta[0]=theta_1,theta[1]=(theta_1)',theta[2]=theta_2,theta[3]=(theta_2)' """
    m1,m2 = 1, 1
    L1,L2 = 2, 2
    ang1 = theta[0,0]
    w1 = theta[1,0]
    ang2 = theta[2,0]
    w2 = theta[3,0]
    delta = ang2 - ang1
    z1 = w1
    z2_numinator = ( (m2*L1*(w1**2) * np.sin(delta) * np.cos(delta)) + (m2*g*np.sin(ang2)*np.cos(delta)) + (m2*L2*(w2**2)*np.sin(delta)) - ((m1+m2)*g*np.sin(ang1)) )
    z2_denominator = ( ((m1+m2)*L1) - (m2*L1*(np.cos(delta)**2)) )
    z2 = z2_numinator / z2_denominator
    z3 = w2
    z4_numinator = ( (-m2*L2*(w2**2)*np.sin(delta)*np.cos(delta)) +  ((m1+m2)*((g*np.sin(ang1)*np.cos(delta)) - (L1*(w1**2)*np.sin(delta)) - (g*np.sin(ang2)))) )
    z4_denominator = ( ((m1+m2)*L1) - (m2*L2*(np.cos(delta)**2)))
    z4 = z4_numinator / z4_denominator
    return np.matrix([[z1],[z2],[z3],[z4]])

def F3(theta,m1,m2,L1,L2,G):
    """theta[0]=theta_1,theta[1]=(theta_1)',theta[2]=theta_2,theta[3]=(theta_2)' """
    ang1 = theta[0,0]
    w1 = theta[1,0]
    ang2 = theta[2,0]
    w2 = theta[3,0]
    delta = ang2 - ang1
    z1 = w1
    z2_numinator = ( (m2*L1*(w1**2) * np.sin(delta) * np.cos(delta)) + (m2*G*np.sin(ang2)*np.cos(delta)) + (m2*L2*(w2**2)*np.sin(delta)) - ((m1+m2)*G*np.sin(ang1)) )
    z2_denominator = ( ((m1+m2)*L1) - (m2*L1*(np.cos(delta)**2)) )
    z2 = z2_numinator / z2_denominator
    z3 = w2
    z4_numinator = ( (-m2*L2*(w2**2)*np.sin(delta)*np.cos(delta)) +  ((m1+m2)*((G*np.sin(ang1)*np.cos(delta)) - (L1*(w1**2)*np.sin(delta)) - (G*np.sin(ang2)))) )
    z4_denominator = ( ((m1+m2)*L1) - (m2*L2*(np.cos(delta)**2)))
    z4 = z4_numinator / z4_denominator
    return np.matrix([[z1],[z2],[z3],[z4]])

def euler(n, T, initalValues, F):
    h = T/n
    theta = np.zeros((initalValues.size, n+1))
    theta[:,0:1] = initalValues
    for i in range(n):
        theta[:,(i+1):(i+2)] = theta[:,i:(i+1)] + F(theta[:,i:(i+1)])*h
    return theta

def RungeKutta(n, T, initalValues, F):
    h = T/n
    theta = np.zeros((initalValues.size, n+1))
    theta[:,0:1] = initalValues
    for i in range(n):
        k1 = F(theta[:,i:(i+1)])
        k2 = F(theta[:,i:(i+1)] + k1*(h/2))
        k3 = F(theta[:,i:(i+1)] + k2*(h/2))
        k4 = F(theta[:,i:(i+1)] + k3*h)
        theta[:,(i+1):(i+2)] = theta[:,i:(i+1)] + (k1 + 2*k2 + 2*k3 + k4)*(h/6)
    return theta

def RungeKuttaModified(n, T, initalValues, F,L1,L2,m1,m2,G):
    h = T/n
    theta = np.zeros((initalValues.size, n+1))
    theta[:,0:1] = initalValues
    for i in range(n):
        k1 = F(theta[:,i:(i+1)],m1,m2,L1,L2,G)
        k2 = F(theta[:,i:(i+1)] + k1*(h/2),m1,m2,L1,L2,G)
        k3 = F(theta[:,i:(i+1)] + k2*(h/2),m1,m2,L1,L2,G)
        k4 = F(theta[:,i:(i+1)] + k3*h,m1,m2,L1,L2,G)
        theta[:,(i+1):(i+2)] = theta[:,i:(i+1)] + (k1 + 2*k2 + 2*k3 + k4)*(h/6)
    return theta

def penduliFinalPosition(initial_value, n, T):
    '''Returns the position of the penduli at the end of the time interval'''
    theta = euler(n, T, initial_value, F2)
    return theta[0,-1], theta[2,-1]

#/-------------------------Animation and Plots-------------------------/
def animateOnePendulum(theta, n, T, name):
    print(ANIMATION_TEXT)
    def animate(i):
        line.set_data([0,x[i]], [0,y[i]])
        ln1.set_data([x[i]], [y[i]])
        if i != 0:
            trace.set_data(x[0:i], y[0:i])

    L = 2
    x = np.cos(theta[0]-(np.pi/2)) * L
    y = np.sin(theta[0]-(np.pi/2)) * L
    
    # print("x:{}\ny:{}".format(x, y))
    fig,ax = plt.subplots(1,1)
    ax.grid()
    ln1, = plt.plot([], [], 'ro', markersize=10)
    line, = plt.plot([], [], 'b-', linewidth=3)
    trace, = plt.plot([], [], 'r--')
    ax.set_ylim(-2.5, 2.5)
    ax.set_xlim(-2.5, 2.5)
    ani = animation.FuncAnimation(fig, animate, frames=n, interval=T)
    ani.save(ANIMATION_PATH.format(name), writer='pillow', fps=n/T)
    printFigText(name)

def animateTwoPendulums(theta, n, T, name):
    print(ANIMATION_TEXT)
    def animate(i):
        if i != 0:
            trace.set_data(x2[0:i], y2[0:i])
        line1.set_data([0,x1[i]], [0,y1[i]])
        point1.set_data([x1[i]], [y1[i]])
        line2.set_data([x1[i], x2[i]], [y1[i], y2[i]])
        point2.set_data([x2[i]], [y2[i]])

    L = 2
    x1 = np.cos(theta[0]-(np.pi/2)) * L
    y1 = np.sin(theta[0]-(np.pi/2)) * L
    x2 = (np.cos(theta[2]-(np.pi/2)) * L) + x1
    y2 = (np.sin(theta[2]-(np.pi/2)) * L) + y1
    
    fig,ax = plt.subplots(1,1)
    ax.grid()
    point1, = plt.plot([], [], 'ro', markersize=10)
    line1, = plt.plot([], [], 'b-', linewidth=3)
    point2, = plt.plot([], [], 'ko', markersize=10)
    line2, = plt.plot([], [], 'g-', linewidth=3)
    trace, = plt.plot([], [], 'k--')
    
    ax.set_ylim(-4.5, 4.5)
    ax.set_xlim(-4.5, 4.5)
    ani = animation.FuncAnimation(fig, animate, frames=n, interval=T)
    ani.save(ANIMATION_PATH.format(name), writer='pillow', fps=n/T)
    printFigText(name)

def animateAllPendulums(theta1, theta2, n, T, name):
    print(ANIMATION_TEXT)
    def animate(i):
        line1_1.set_data([0, x1_1[i]], [0, y1_1[i]])
        line1_2.set_data([0,x1_2[i]], [0,y1_2[i]])
        line2_2.set_data([x1_2[i], x2_2[i]], [y1_2[i], y2_2[i]])
        if i != 0:
            trace1.set_data(x1_1[0:i], y1_1[0:i])
            trace2.set_data(x2_2[0:i], y2_2[0:i])
    
    L = 2
    x1_1 = np.cos(theta1[0]-(np.pi/2)) * L
    y1_1 = np.sin(theta1[0]-(np.pi/2)) * L

    x1_2 = np.cos(theta2[0]-(np.pi/2)) * L
    y1_2 = np.sin(theta2[0]-(np.pi/2)) * L
    x2_2 = (np.cos(theta2[2]-(np.pi/2)) * L) + x1_2
    y2_2 = (np.sin(theta2[2]-(np.pi/2)) * L) + y1_2

    fig,ax = plt.subplots(1,1)
    ax.grid()
    line1_1, = plt.plot([], [], 'k-', linewidth=3)
    line1_2, = plt.plot([], [], 'b-', linewidth=3)
    line2_2, = plt.plot([], [], 'g-', linewidth=3)
    trace1, = plt.plot([], [], 'k--')
    trace2, = plt.plot([], [], 'g--')

    ax.set_ylim(-4.5, 4.5)
    ax.set_xlim(-4.5, 4.5)
    ani = animation.FuncAnimation(fig, animate, frames=n, interval=T)
    ani.save(ANIMATION_PATH.format(name), writer='pillow', fps=n/T)
    printFigText(name)

def animateTwoDoublePendulums(theta1, theta2, n, T, name):
    print(ANIMATION_TEXT)
    def animate(i):
        line1_1.set_data([0, x1_1[i]], [0, y1_1[i]])
        line2_1.set_data([x1_1[i], x2_1[i]], [y1_1[i], y2_1[i]])

        line1_2.set_data([0,x1_2[i]], [0,y1_2[i]])
        line2_2.set_data([x1_2[i], x2_2[i]], [y1_2[i], y2_2[i]])
        time_text.set_text(time_template % (i*T/n))
        if i != 0:
            trace1.set_data(x2_1[0:i], y2_1[0:i])
            trace2.set_data(x2_2[0:i], y2_2[0:i])
    L = 2
    x1_1 = np.cos(theta1[0]-(np.pi/2)) * L
    y1_1 = np.sin(theta1[0]-(np.pi/2)) * L
    y2_1 = (np.sin(theta1[2]-(np.pi/2)) * L) + y1_1
    x2_1 = (np.cos(theta1[2]-(np.pi/2)) * L) + x1_1

    x1_2 = np.cos(theta2[0]-(np.pi/2)) * L
    y1_2 = np.sin(theta2[0]-(np.pi/2)) * L
    x2_2 = (np.cos(theta2[2]-(np.pi/2)) * L) + x1_2
    y2_2 = (np.sin(theta2[2]-(np.pi/2)) * L) + y1_2

    fig,ax = plt.subplots(1,1)
    ax.grid()

    time_template = 'time = %.1fs'
    time_text = ax.text(0.05, 0.9, '', transform=ax.transAxes)

    line1_1, = plt.plot([], [], 'k-', linewidth=3)
    line2_1, = plt.plot([], [], 'r-', linewidth=3)
    line1_2, = plt.plot([], [], 'b-', linewidth=3)
    line2_2, = plt.plot([], [], 'g-', linewidth=3)
    trace1, = plt.plot([], [], 'r--')
    trace2, = plt.plot([], [], 'g--')
   

    ax.set_ylim(-4.5, 4.5)
    ax.set_xlim(-4.5, 4.5)
    ani = animation.FuncAnimation(fig, animate, frames=n, interval=T)
    ani.save(ANIMATION_PATH.format(name), writer='pillow', fps=n/T)
    printFigText(name)


def plotOnePendulum(theta, n, T, verbose=True, name=None, fig=None, ax=None, label=None):
    if verbose: print(PLOT_TEXT)
    x = np.linspace(0, T, n+1)
    y = theta[0]

    if fig == None and ax == None:
        fig, ax = plt.subplots(1,1)
    ax.grid(visible=True)
    ax.plot(x, y, label=label)
    ax.legend(loc="upper right")
    ax.set_xlabel("Time")
    ax.set_ylabel("Angle")

    if name != None:
        fig.savefig(PLOT_PATH.format(name))

    if verbose: printFigText(name)
    return fig, ax

def plotTwoPendulums(theta, n, T, verbose=True, name=None, fig=None, ax=None, label=None):
    if verbose: print(PLOT_TEXT)
    x = np.linspace(0, T, n+1)
    y1 = theta[0]
    y2 = theta[2]

    if fig == None and ax == None:
        fig, ax = plt.subplots(1,1)
    ax.grid(visible=True)
    if label:
        ax.plot(x, y1, label=label[0])
        ax.plot(x, y2, label=label[1])
    ax.legend(loc="upper right")
    ax.set_xlabel("Time")
    ax.set_ylabel("Angle")

    if name != None:
        fig.savefig(PLOT_PATH.format(name))

    if verbose: printFigText(name)
    return fig, ax

#/-------------------------Problems------------------------------------/
def prob1():
    print("\n----Problem 1")
    print("See report for the math.")

def prob2():
    print("\n---- Problem 2")
    prob2InitialVal = np.matrix([[1], [0]])
    theta = euler(200, 5, prob2InitialVal, F1)
    print("First 20 digits from the Euler method:")
    [print(theta[0][i].round(2), end=' ') for i in range(20)]
    print()

def prob3(plot=True, animate=True):
    print("\n---- Problem 3")
    prob3InitialVal = np.matrix([[np.pi/12], [0]])
    theta = euler(500, 20, prob3InitialVal, F1)
    animateOnePendulum(theta=theta, n=500, T=20, name="problem3")
    plotOnePendulum(theta=theta, n=500, T=20, name="problem3")

def prob4():
    print("\n---- Problem 4")
    
    # problem 3 to compare on plot
    prob3InitialVal = np.matrix([[np.pi/12], [0]])
    theta3 = euler(500, 20, prob3InitialVal, F1)

    # problem 4 different initial conditions
    prob4InitialVal = np.matrix([[np.pi/2], [0]])
    theta4 = euler(500, 20, prob4InitialVal, F1)

    animateOnePendulum(theta4, 500, 20, "problem4")
    fig, ax = plotOnePendulum(theta=theta3, n=500, T=20, verbose=False, name=None, label="Problem 3")
    plotOnePendulum(theta=theta4, n=500, T=20, verbose=True, name="problem4_combined_w_prob3", fig=fig, ax=ax, label="Problem 4")
    plotOnePendulum(theta=theta4, n=500, T=20, verbose=True, name="problem4")

def prob5():
    print("\n---- Problem 5")
    prob5InitialVal = np.matrix([[np.pi/12], [0]])
    theta5a = RungeKutta(500, 20, prob5InitialVal, F1)
    animateOnePendulum(theta5a, 500, 20, "problem5_a")
    fig, ax = plotOnePendulum(theta5a, 500, 20, name="problem5_a", label="θ(0)=π/12")

    prob5InitialVal = np.matrix([[np.pi/2], [0]])
    theta5b = RungeKutta(500, 20, prob5InitialVal, F1)
    animateOnePendulum(theta5b, 500, 20, "problem5_b")
    plotOnePendulum(theta5b, 500, 20, name="problem5_b")

    # combine plots for comparison
    plotOnePendulum(theta5b, 500, 20, verbose=False, name='problem5_combined_ab', fig=fig, ax=ax, label="θ(0)=π/2")

    # combine problems 5 with 3 and 4
    # problem 3 values
    prob3InitialVal = np.matrix([[np.pi/12], [0]])
    theta3 = euler(500, 20, prob3InitialVal, F1)

    # problem 4 values
    prob4InitialVal = np.matrix([[np.pi/2], [0]])
    theta4 = euler(500, 20, prob4InitialVal, F1)

    fig, ax = plotOnePendulum(theta=theta3, n=500, T=20, verbose=False, name=None, label="Eueler, θ(0)=π/12")
    plotOnePendulum(theta=theta5a, n=500, T=20, verbose=False, name='problem5a_combine_w_prob3', fig=fig, ax=ax, label="Runge-Kutta, θ(0)=π/12")

    fig, ax = plotOnePendulum(theta=theta4, n=500, T=20, verbose=False, name=None, label="Eueler, theta= pi/2")
    plotOnePendulum(theta=theta5b, n=500, T=20, verbose=False, name='problem5b_combine_w_prob4', fig=fig, ax=ax, label="Runge-Kutta, θ(0)=π/2")
  
def prob6():
    print("\n---- Problem 6")
    print("See report for the math.")

def prob7():
    print("\n---- Problem 7")
    prob6InitialVal = np.matrix([[np.pi/3], [0], [np.pi/6], [0]])
    theta = euler(200, 8, prob6InitialVal, F2)
    animateTwoPendulums(theta, 200, 8, "problem7")
    plotTwoPendulums(theta=theta, n=200, T=20, name="problem7", label=["θ(0)=π/3", "θ(0)=π/6"])

def prob8():
    print("\n---- Problem 8")
    n = 200
    T = 8
    prob8InitialVal1 = np.matrix([[np.pi/3], [0]])
    prob8InitialVal2 = np.matrix([[np.pi/3], [0], [np.pi/6], [0]])
    theta1 = euler(n, T, prob8InitialVal1, F1)
    theta2 = euler(n, T, prob8InitialVal2, F2)
    animateAllPendulums(theta1, theta2, n, T, "problem8")

def prob9():
    print("\n---- Problem 9")
    initial_values = np.matrix([[np.pi/3], [0], [np.pi/6], [0]])
    n_values = [100, 200, 400, 800, 1600, 3200, 6400]
    T = 20
    for initial_val in initial_values:
        for n in n_values:
            pos1, pos2 = penduliFinalPosition(initial_val, n, T)

def prob10():

    print("\n---- Problem 10")
    print("Not done")
    
    
def prob11(initial =  np.matrix([[2*np.pi/3], [0], [np.pi/6], [0]]),k = [1,2,3,4,5],n=1000,t=40,gif=True):
    if gif:
        print("\n---- Problem 11")
    prob11InitialVal1 = initial.copy()
    prob11InitialVal2 = initial.copy()
    theta1 = RungeKutta(n,t, prob11InitialVal1, F2)
    results = np.zeros((4,(n+1)*(len(k)+1)))
    #add all the values from theta1 to results
    for i in range(0,4):
        results[i,0:(n+1)] = theta1[i,:]

    for i in range(1,len(k)+1):
        prob11InitialVal2Error = prob11InitialVal2.copy() + np.matrix([[10**(-k[i-1])], [0], [10**(-k[i-1])], [0]])
        theta2 = RungeKutta(n, t, prob11InitialVal2Error, F2)
        if gif:
            animateTwoDoublePendulums(theta1,theta2, n, t, "problem11_k_{}".format(k[i-1]))
        for j in range(0,4):
            results[j,(n+1)*i:(n+1)*(i+1)] = theta2[j,:]    

    #compare the results
    for i in range(1,len(k)+1):
        print("k = {}, i.e. ε = 10^-{}".format(k[i-1],k[i-1]))
        done1 = False
        done2 = False
        for j in range(n):
            theta1margin = abs(results[0,j]*0.01)
            theta2margin = abs(results[2,j]*0.01)
            if (abs(abs(results[0,((n+1)*i)+j]) - abs(results[0,j]))) > theta1margin and not done1:
                done1 = True
                print("theta 1 is not within 1% of the control at n = {}".format(j))

            if (abs(results[2,((n+1)*i)+j]) - abs(results[2,j])) > theta2margin and not done2:
                print("theta 2 is not within 1% of the control at n = {}".format(j))
                done2 = True
            if done1 and done2:
                break
        if not done1:
            print("theta 1 is always within 1% of the control for n = {} and T = {}".format(n,t))
        if not done2:
            print("theta 2 is always within 1% of the control for n = {} and T = {}".format(n,t))
        print("\n")
        
    
def prob12():
    print("\n---- Problem 12")
    initial = np.matrix([[2*np.pi/3], [0], [np.pi/6], [0]])
    T = 40
    k = [1,2,3,4,5,6,7,8,9,10,11,12]
    n = 1000
    print("We begin by studying the error in the Runge-Kutta method for the double pendulum for different values of k. i.e. 10^(-k).")
    prob11(initial,k,n,T,False)

    print("Now we study the error in the Runge-Kutta method for the double pendulum for different values of n.")
    print("We will be using k = 12, i.e. 10^(-12).")
    new_k = [12]
    for _ in range(0,5):
        print("n = {}".format(n))
        prob11(initial,new_k,n,T,False)
        n = n*2
    print("We will now study the error in the Runge-Kutta method for the double pendulum for different values of T.")
    print("We will be using k = 12, i.e. 10^(-12) and n = 16000.")
    n = 16000
    for _ in range(0,5):
        print("T = {}".format(T))
        prob11(initial,new_k,n,T,False)
        T = T*2
    print("We will now study the error in the Runge-Kutta method for the double pendulum for different initial values.")
    print("We will be using k = 12, i.e. 10^(-12) and n = 8000 and T = 40.")
    T = 40
    n = 1000
    theta1_options = [2*np.pi/3, 2*np.pi/5, 2*np.pi/7, 2*np.pi/13]
    theta2_options = [np.pi/6, np.pi/10, np.pi/14, np.pi/26]
    for i in range(len(theta1_options)):
        for j in range(len(theta2_options)):
            print("Initial values: θ(0) = {}, θ'(0) = 0, θ(0) = {}, θ'(0) = 0".format(theta1_options[i],theta2_options[j]))
            initial = np.matrix([[theta1_options[i]], [0], [theta2_options[j]], [0]])
            prob11(initial=initial,k = new_k,n = n,t = T,gif = False)
    
    
def prob13():
    print("\n---- Problem 13")
    initial = np.matrix([[2*np.pi/3], [0], [np.pi/6], [0]])
    t = 40
    n = 10000
    L1 = 2
    L2 = 2
    m1 = 1
    m2 = 1
    g = 9.81
    options_for_L1 = [1,1.5,2,2.5,3]
    results_for_L1 = np.zeros((1,len(options_for_L1)))
    normal = RungeKutta(n,t,initial,F2)
    for i,j in enumerate(options_for_L1):
        done = False
        theta1 = RungeKuttaModified(n, t, initial, F3, L1 = j, L2 = L2, m1 = m1, m2 = m2, G = g)
        for k in range(n):
            if abs(normal[0,k] - theta1[0,k]) > 0.01:
                results_for_L1[0,i] = k
                done = True
                break
        if not done:
            results_for_L1[0,i] = n
    options_for_L2 = [1,1.5,2,2.5,3]
    results_for_L2 = np.zeros((1,len(options_for_L2)))
    for i,j in enumerate(options_for_L2):
        done = False
        theta1 = RungeKuttaModified(n, t, initial, F3, L1 = L1, L2 = j, m1 = m1, m2 = m2, G = g)
        for k in range(n):
            if abs(normal[0,k] - theta1[0,k]) > 0.01:
                results_for_L2[0,i] = k
                done = True
                break
        if not done:
            results_for_L2[0,i] = n
    options_for_m1 = [1,1.5,2,2.5,3,3.5,4,4.5,5]
    results_for_m1 = np.zeros((1,len(options_for_m1)))
    for i,j in enumerate(options_for_m1):
        done = False
        theta1 = RungeKuttaModified(n, t, initial, F3, L1 = L1, L2 = L2, m1 = j, m2 = m2, G = g)
        for k in range(n):
            if abs(normal[0,k] - theta1[0,k]) > 0.01:
                results_for_m1[0,i] = k
                done = True
                break
        if not done:
            results_for_m1[0,i] = n
    options_for_m2 = [1,1.5,2,2.5,3,3.5,4,4.5,5]
    results_for_m2 = np.zeros((1,len(options_for_m2)))
    for i,j in enumerate(options_for_m2):
        done = False
        theta1 = RungeKuttaModified(n, t, initial, F3, L1 = L1, L2 = L2, m1 = m1, m2 = j, G = g)
        for k in range(n):
            if abs(normal[0,k] - theta1[0,k]) > 0.01:
                results_for_m2[0,i] = k
                done = True
                break
        if not done:
            results_for_m2[0,i] = n
    options_for_g = [9.81,1.62,3.7,8.87,275,24.5]
    results_for_g = np.zeros((1,len(options_for_g)))
    for i,j in enumerate(options_for_g):
        done = False
        theta1 = RungeKuttaModified(n, t, initial, F3, L1 = L1, L2 = L2, m1 = m1, m2 = m2, G = j)
        for k in range(n):
            if abs(normal[0,k] - theta1[0,k]) > 0.01:
                results_for_g[0,i] = k
                done = True
                break
        if not done:
            results_for_g[0,i] = n
    
    #subplots for 4 subplots
    # fig, axs = plt.subplots(2, 2)
    # axs[0, 0].plot(options_for_L1, results_for_L1[0,:])
    # axs[0, 0].set_title('L1')
    # axs[0, 1].plot(options_for_L2, results_for_L2[0,:])
    # axs[0, 1].set_title('L2')
    # axs[1, 0].plot(options_for_m1, results_for_m1[0,:])
    # axs[1, 0].set_title('m1')
    # axs[1, 1].plot(options_for_m2, results_for_m2[0,:])
    # axs[1, 1].set_title('m2')
    # for ax in axs.flat:
    #     ax.set(xlabel='value', ylabel='time')
    # # Hide x labels and tick labels for top plots and y ticks for right plots.
    # for ax in axs.flat:
    #     ax.label_outer()
    plt.show()
    plt.plot(options_for_g, results_for_g[0,:])
    plt.title('g')
    plt.xlabel('value')
    plt.ylabel('time')
    plt.show()

    # plt.plot(options_for_L1,results_for_L1[0,:],label = "L1")
    # plt.plot(options_for_L2,results_for_L2[0,:],label = "L2")
    # plt.plot(options_for_m1,results_for_m1[0,:],label = "m1")
    # plt.plot(options_for_m2,results_for_m2[0,:],label = "m2")
    # plt.plot(options_for_g,results_for_g[0,:],label = "g")
    # plt.legend()
    # plt.xlabel("Value of parameter")
    # plt.ylabel("Time until divergence")
    # plt.title("Time until divergence for different parameters")
    # plt.show()





if __name__ == "__main__":
    # prob1()
    # prob2()
    # prob3()
    # prob4()
    # prob5()
    # prob6()
    # prob7()
    # prob8()
    # prob9()
    # prob10()
    # prob11()
    # prob12()
    prob13()
    # question_string = "Which question would you like to run (1-13, q to quit): "
    # question_available = [str(i) for i in range(1,13)] + ["q"]
    # question = getInput(question_string, question_available)
    # while True:
    #     if question == "q":
    #         break
    #     elif question == "1":
    #         prob1()
    #     elif question == "2":
    #         prob2()
    #     elif question == "3":
    #         prob3()
    #     elif question == "4":
    #         prob4()
    #     elif question == "5":
    #         prob5()
    #     elif question == "6":
    #         prob6()
    #     elif question == "7":
    #         prob7()
    #     elif question == "8":
    #         prob8()
    #     elif question == "9":
    #         prob9()
    #     elif question == "10":
    #         prob10()
    #     elif question == "11":
    #         prob11()
    #     elif question == "12":
    #         prob12()
    #     elif question == "13":
    #         prob13()
    #     else:
    #         print("Invalid input")
    #     question = getInput(question_string, question_available)

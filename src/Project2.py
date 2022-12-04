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

def euler(n, T, initalValues, F):
    h = T/n
    theta = np.zeros((initalValues.size, n+1))
    theta[:,0:1] = initalValues
    for i in range(n):
        theta[:,(i+1):(i+2)] = theta[:,i:(i+1)] + F(theta[:,i:(i+1)])*h
        print(theta)
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

def penduliFinalPosition(initial_value, n, T):
    '''Returns the position of the penduli at the end of the time interval'''
    theta = euler(n, T, initial_value, F2)
    return theta[0,-1], theta[2,-1]

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
            trace2.set_data(x1[0:i], y1[0:i])
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
    trace2, = plt.plot([], [], 'r--')
    
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
    line1_1, = plt.plot([], [], 'k-', linewidth=3)
    line2_1, = plt.plot([], [], 'r-', linewidth=3)
    line1_2, = plt.plot([], [], 'b-', linewidth=3)
    line2_2, = plt.plot([], [], 'g-', linewidth=3)
    trace1, = plt.plot([], [], 'k--')
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

def changeThetaBetween2Pi(theta, nrPends):
    '''To make the angle of theta between -pi and pi. nrPends = number of pendulums used. Returns [theta1, theta2]'''
    tempTheta = np.zeros((nrPends, theta.shape[1]))
    for i in range(theta.shape[1]):
        if((theta[0,i] > np.pi) or (theta[0,i] < -np.pi)):
            modPi = (theta[0, i])%(np.pi)
            mod2Pi = (theta[0, i])%(2*np.pi)
            tempTheta[0, i] = 2*modPi - mod2Pi
        else:
            tempTheta[0, i] = theta[0, i]

        if(nrPends == 2):
            if((theta[2,i] > np.pi) or (theta[2,i] < -np.pi)):
                modPi = (theta[2, i])%(np.pi)
                mod2Pi = (theta[2, i])%(2*np.pi)
                tempTheta[1, i] = 2*modPi - mod2Pi                
            else:
                tempTheta[1, i] = theta[2, i]
        
    return tempTheta 

def animateOnlyParametrizedCurve(theta, n, T, name):
    def animate(i):
        if i != 0:
            trace.set_data(x[0:i], y[0:i])

    tempTheta = changeThetaBetween2Pi(theta, 2)
    x = tempTheta[0]
    y = tempTheta[1]

    fig, ax = plt.subplots(1,1)
    trace, = ax.plot([], [], 'r', linewidth=0.5)
    ax.set_ylim(-(np.pi + 0.5), (np.pi + 0.5))
    ax.set_xlim(-(np.pi + 0.5), (np.pi + 0.5))
    ax.set_xlabel("Theta 1 [rad]")
    ax.set_ylabel("Theta 2 [rad]")
    ax.grid()

    ani = animation.FuncAnimation(fig, animate, frames=n, interval=T)
    ani.save(ANIMATION_PATH.format(name), writer='pillow', fps=n/T)
    printFigText(name)

def animateParametrizedCurve(theta, n, T,time, name):
    def animate(i):
        pend1.set_data([0, x1[i]], [0, y1[i]])
        pend2.set_data([x1[i], x2[i]], [y1[i], y2[i]])
        if i != 0:
            trace.set_data(x[0:i], y[0:i])
            ang1.set_data(time[0:i], x[0:i])
            ang2.set_data(time[0:i], y[0:i])

    tempTheta = changeThetaBetween2Pi(theta, 2)
    x = tempTheta[0]
    y = tempTheta[1]

    L = 2
    x1 = np.cos(theta[0]-(np.pi/2)) * L
    y1 = np.sin(theta[0]-(np.pi/2)) * L
    y2 = (np.sin(theta[2]-(np.pi/2)) * L) + y1
    x2 = (np.cos(theta[2]-(np.pi/2)) * L) + x1

    fig, ax = plt.subplots(2,2)
    ang1, = ax[0, 0].plot([], [], 'g', linewidth=1)
    ax[0, 0].set_ylim(-(np.pi + 0.5), (np.pi + 0.5))
    ax[0, 0].set_xlim(np.min(time),np.max(time))
    ax[0, 0].grid()
    ax[0, 0].set_ylabel("Theta 1 [rad]")

    ang2, = ax[1, 0].plot([], [], 'b', linewidth=1)
    ax[1, 0].set_ylim(-(np.pi + 0.5), (np.pi + 0.5))
    ax[1, 0].set_xlim(np.min(time),np.max(time))
    ax[1, 0].grid()
    ax[1, 0].set_xlabel("Time[s]")
    ax[1, 0].set_ylabel("Theta 2 [rad]")

    pend1, = ax[0,1].plot([], [], 'g', linewidth=3)
    pend2, = ax[0,1].plot([], [], 'b', linewidth=3)
    ax[0, 1].set_ylim(-(L*2 + 0.5), L*2 + 0.5)
    ax[0, 1].set_xlim(-(L*2 + 0.5), (L*2 + 0.5))
    ax[0, 1].grid()
    
    trace, = ax[1,1].plot([], [], 'r', linewidth=0.5)
    ax[1, 1].set_ylim(-(np.pi + 0.5), (np.pi + 0.5))
    ax[1, 1].set_xlim(-(np.pi + 0.5), (np.pi + 0.5))
    ax[1, 1].grid()
    ani = animation.FuncAnimation(fig, animate, frames=n, interval=T)
    ani.save(ANIMATION_PATH.format(name), writer='pillow', fps=n/T)
    printFigText(name)

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

def prob13Helper(initial, initialError, n, T, constantsArray, constantIndex, options):
    '''constantArr = [L1, L2, m1, m2, g]'''
    h = T/n
    # options_for_L1 = [1,1.5,2,2.5,3]
    results = np.zeros((1,len(options)))
    for i,j in enumerate(options):
        constantsArray[constantIndex] = j
        done = False
        normal = RungeKuttaModified(n, T, initial, F3, L1 = constantsArray[0], L2 = constantsArray[1], m1 = constantsArray[2], m2 = constantsArray[3], G = constantsArray[4])
        theta1 = RungeKuttaModified(n, T, initialError, F3, L1 = constantsArray[0], L2 = constantsArray[1], m1 = constantsArray[2], m2 = constantsArray[3], G = constantsArray[4])
        for k in range(n):
            if (abs(normal[0,k] - theta1[0,k]) > 0.01) or (abs(normal[0,k] - theta1[0,k]) > 0.01):
                results[0,i] = h*k#k
                done = True
                break
        if not done:
            results[0,i] = h*n#n

    return results

def animateBestProb13(theta1, theta2, n, T, name, L1, L2):
    print(ANIMATION_TEXT)
    def animate(i):
        time_text.set_text(time_template % (i*T/n))
        line1_1.set_data([0, x1_1[i]], [0, y1_1[i]])
        line2_1.set_data([x1_1[i], x2_1[i]], [y1_1[i], y2_1[i]])

        line1_2.set_data([0,x1_2[i]], [0,y1_2[i]])
        line2_2.set_data([x1_2[i], x2_2[i]], [y1_2[i], y2_2[i]])
        if i != 0:
            trace1.set_data(x2_1[0:i], y2_1[0:i])
            trace2.set_data(x2_2[0:i], y2_2[0:i])
    h = T/n
    x1_1 = np.cos(theta1[0]-(np.pi/2)) * L1
    y1_1 = np.sin(theta1[0]-(np.pi/2)) * L1
    y2_1 = (np.sin(theta1[2]-(np.pi/2)) * L2) + y1_1
    x2_1 = (np.cos(theta1[2]-(np.pi/2)) * L2) + x1_1

    x1_2 = np.cos(theta2[0]-(np.pi/2)) * L1
    y1_2 = np.sin(theta2[0]-(np.pi/2)) * L1
    x2_2 = (np.cos(theta2[2]-(np.pi/2)) * L2) + x1_2
    y2_2 = (np.sin(theta2[2]-(np.pi/2)) * L2) + y1_2

    fig,ax = plt.subplots(1,1)
    ax.grid()
    line1_1, = plt.plot([], [], 'k-', linewidth=3)
    line2_1, = plt.plot([], [], 'r-', linewidth=3)
    line1_2, = plt.plot([], [], 'b-', linewidth=3)
    line2_2, = plt.plot([], [], 'g-', linewidth=3)
    trace1, = plt.plot([], [], 'r--')
    trace2, = plt.plot([], [], 'g--')
    time_template = 'time = %.1fs'
    time_text = ax.text(0.05, 0.9, '', transform=ax.transAxes)

    ax.set_ylim(-(L1 + L2 + 0.5), (L1 + L2 + 0.5))
    ax.set_xlim(-(L1 + L2 + 0.5), (L1 + L2 + 0.5))
    ani = animation.FuncAnimation(fig, animate, frames=n, interval=T)
    ani.save(ANIMATION_PATH.format(name), writer='pillow', fps=n/T)
    printFigText(name)

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
    theta = RungeKutta(500, 20, prob6InitialVal, F2)
    animateTwoPendulums(theta, 500, 20, "problem7")

def prob8():
    print("\n---- Problem 8")
    n = 500
    T = 20
    prob8InitialVal2_1 = np.matrix([[np.pi/3], [0], [-np.pi/3], [0]])
    prob8InitialVal2_2 = np.matrix([[np.pi/2], [0], [np.pi/3], [0]])
    prob8InitialVal2_3 = np.matrix([[np.pi/5], [0], [np.pi/2], [0]])
    prob8InitialVal2_4 = np.matrix([[-np.pi/3], [0], [-np.pi/6], [0]])
    theta1 = RungeKutta(n, T, prob8InitialVal2_1, F2)
    theta2 = RungeKutta(n, T, prob8InitialVal2_2, F2)
    theta3 = RungeKutta(n, T, prob8InitialVal2_3, F2)
    theta4 = RungeKutta(n, T, prob8InitialVal2_4, F2)
    animateTwoPendulums(theta1, 500, 20, "problem8_1")
    animateTwoPendulums(theta2, 500, 20, "problem8_2")
    animateTwoPendulums(theta3, 500, 20, "problem8_3")
    animateTwoPendulums(theta4, 500, 20, "problem8_4")
 

def prob9():
    print("\n---- Problem 9")
    num_of_InitialVal = 20
    estimate_theta = np.zeros((num_of_InitialVal,4))
    sign = -1 # to alternate between positive and negative for theta1 and theta2
    for i in range(1, num_of_InitialVal+1):
        
        prob9InitialVal = np.matrix([[-1*sign*np.pi*i*0.02], [0], [sign*np.pi*i*0.02], [0]])
        estimate_theta[i-1] = RungeKutta(12800, 20, prob9InitialVal, F2)[:,12800]
        sign *= -1
    n = 100
    n_array = np.zeros(7)
    error_array = np.zeros((7,num_of_InitialVal))
    counter = 0
    while n <= 6400:
        temp_error_array = np.zeros(num_of_InitialVal)
        sign = -1
        for i in range(1, num_of_InitialVal+1):
            prob9InitialVal = np.matrix([[-1*sign*np.pi*i*0.02], [0], [sign*np.pi*i*0.02], [0]])
            temp_error_array[i-1] = np.linalg.norm(estimate_theta[i-1] - RungeKutta(n, 20, prob9InitialVal, F2)[:,n])
            sign *= -1
        # mean_error = np.mean(error_array)

        error_array[counter] = temp_error_array
        n_array[counter] = n
        counter += 1
        print(n)
        n = n*2

    slope = np.zeros(num_of_InitialVal)

    for i in range(0,num_of_InitialVal):
        slope[i] = np.polyfit(np.log10(n_array), np.log10(error_array[:,i]), 1)[0]

    print("mean: {} standard deviation: {}".format(np.mean(slope), np.std(slope)))
    a,b = np.polyfit(np.log10(n_array), np.log10(error_array[:,3]), 1)
    plt.plot(np.log10(n_array),a*np.log10(n_array)+b)
    plt.plot(np.log10(n_array),np.log10(error_array[:,3]), 'ro')
    plt.xlabel("log10(n)")
    plt.ylabel("log10(error)")
    plt.title("theta1 = 0.06pi, theta2 = -0.06pi")
    plt.show()

def prob10():
    n = 400
    T = 20
    time = np.zeros((n, 1))
    h = T/n
    for i in range(n):
        time[i,0] = h*i
    initialValues = np.matrix([[3*np.pi/8], [0], [np.pi/4], [0]])
    theta = RungeKutta(n, T, initialValues, F2)
    animateParametrizedCurve(theta, n, T, time, "problem10.1a")
    animateOnlyParametrizedCurve(theta, n, T, "problem10.1b")

    # initialValues = np.matrix([[2*np.pi/3], [0], [np.pi/6], [0]])
    # initialValues = np.matrix([[1.01*np.pi/4], [0], [np.pi/4], [0]])
    # theta = RungeKutta(n, T, initialValues, F2)
    # animateParametrizedCurve(theta, n, T, time, "problem10.2a")
    # animateOnlyParametrizedCurve(theta, n, T, "problem10.2b")

    # initialValues = np.matrix([[np.pi/5], [0], [np.pi/10], [0]])
    # initialValues = np.matrix([[1.01*np.pi/4], [0], [np.pi/4], [0]])
    # theta = RungeKutta(n, T, initialValues, F2)
    # animateParametrizedCurve(theta, n, T, time, "problem10.3a")
    # animateOnlyParametrizedCurve(theta, n, T, "problem10.3b")
    
    # initialValues = np.matrix([[1.01*np.pi/4], [0], [np.pi/4], [0]])
    # theta = RungeKutta(n, T, initialValues, F2)
    # animateParametrizedCurve(theta, n, T, time, "problem10.4a")
    # animateOnlyParametrizedCurve(theta, n, T, "problem10.4b")
    print("\n---- Problem 10")
    
    
def prob11(initial =  np.matrix([[2*np.pi/3], [0], [np.pi/6], [0]]),k = [1,2,3,4,5],n=1000,t=40,gif=True):
    h = t/n
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
    return_val = np.zeros((2,len(k)))
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
                return_val[0,i-1] = j*h

            if (abs(results[2,((n+1)*i)+j]) - abs(results[2,j])) > theta2margin and not done2:
                print("theta 2 is not within 1% of the control at n = {}".format(j))
                done2 = True
                return_val[1,i-1] = j*h

            if done1 and done2:
                break
        if not done1:
            print("theta 1 is always within 1% of the control for n = {} and T = {}".format(n,t))
            return_val[0,i-1] = n*h
        if not done2:
            print("theta 2 is always within 1% of the control for n = {} and T = {}".format(n,t))
            return_val[1,i-1] = n*h
        print("\n")
    return return_val
        
    
def prob12():
    print("\n---- Problem 12")
    initial = np.matrix([[2*np.pi/3], [0], [np.pi/6], [0]])
    T = 40
    k = [1,2,3,4,5,6,7,8,9,10,11,12]
    n = 10000
    print("We begin by studying the error in the Runge-Kutta method for the double pendulum for different values of k. i.e. 10^(-k).")
    plot_vals = prob11(initial,k,n,T,False)
    # plot the results
    plt.plot(k,plot_vals[0,:],label="theta 1")
    plt.plot(k,plot_vals[1,:],label="theta 2")
    plt.xlabel("k")
    plt.ylabel("time of deviation from 1%[s]")
    plt.title("time of deviation from 1% vs k for the double pendulum")
    plt.legend()
    plt.savefig("problem12Differntk.png")
    plt.show()

    print("Now we study the error in the Runge-Kutta method for the double pendulum for different values of n.")
    print("We will be using k = 12, i.e. 10^(-12).")
    new_k = [12]
    plot_vals = np.zeros((2,5))
    for i in range(0,5):
        print("n = {}".format(n))
        val = prob11(initial,new_k,n,T,False)
        n = n*2
        plot_vals[0,i] = val[0,0]
        plot_vals[1,i] = val[1,0]
    #plot the results
    plt.plot([10000,20000,40000,80000,160000],plot_vals[0,:],label="theta 1")
    plt.plot([10000,20000,40000,80000,160000],plot_vals[1,:],label="theta 2")
    plt.xlabel("n")
    plt.ylabel("time of deviation from 1%[s]")
    plt.title("time of deviation from 1% vs n for the double pendulum")
    plt.legend()
    plt.savefig("problem12Differntn.png")
    plt.show()

    print("We will now study the error in the Runge-Kutta method for the double pendulum for different values of T.")
    print("We will be using k = 12, i.e. 10^(-12) and n = 16000.")
    n = 16000
    plot_vals = np.zeros((2,5))
    for _ in range(0,5):
        print("T = {}".format(T))
        val = prob11(initial,new_k,n,T,False)
        T = T*2
        plot_vals[0,_] = val[0,0]
        plot_vals[1,_] = val[1,0]
    #plot the results
    plt.plot([40,80,160,320,640],plot_vals[0,:],label="theta 1")
    plt.plot([40,80,160,320,640],plot_vals[1,:],label="theta 2")
    plt.xlabel("T")
    plt.ylabel("time of deviation from 1%[s]")
    plt.title("time of deviation from 1% vs T for the double pendulum")
    plt.legend()
    plt.savefig("problem12DifferntT.png")
    plt.show()

    print("We will now study the error in the Runge-Kutta method for the double pendulum for different initial values.")
    print("We will be using k = 12, i.e. 10^(-12) and n = 8000 and T = 40.")
    T = 40
    n = 1000
    new_k = [12]
    plot_vals = np.zeros((2,50))
    #first we start off only changing theta 1 and keeping theta 2 constant
    for i in range(0,50):
        print("theta 1 = {}".format(i*np.pi/50))
        val = prob11(np.matrix([[i*np.pi/50],[0],[np.pi/6],[0]]),new_k,n,T,False)
        plot_vals[0,i] = val[0,0]
        plot_vals[1,i] = val[1,0]
    #plot the results
    plt.plot([i*np.pi/50 for i in range(0,50)],plot_vals[0,:],label="theta 1")
    plt.plot([i*np.pi/50 for i in range(0,50)],plot_vals[1,:],label="theta 2")
    plt.xlabel("theta 1")
    plt.ylabel("time of deviation from 1%[s]")
    plt.title("time of deviation from 1% vs theta 1 for the double pendulum")
    plt.legend()
    plt.savefig("problem12Differnttheta1.png")
    plt.show()

    #now we change theta 2 and keep theta 1 constant
    plot_vals = np.zeros((2,50))
    for i in range(0,50):
        print("theta 2 = {}".format(i*np.pi/50))
        val = prob11(np.matrix([[np.pi/6],[0],[i*np.pi/50],[0]]),new_k,n,T,False)
        plot_vals[0,i] = val[0,0]
        plot_vals[1,i] = val[1,0]
    #plot the results
    plt.plot([i*np.pi/50 for i in range(0,50)],plot_vals[0,:],label="theta 1")
    plt.plot([i*np.pi/50 for i in range(0,50)],plot_vals[1,:],label="theta 2")
    plt.xlabel("theta 2")
    plt.ylabel("time of deviation from 1%[s]")
    plt.title("time of deviation from 1% vs theta 2 for the double pendulum")
    plt.legend()
    plt.savefig("problem12Differnttheta2.png")
    plt.show()
    



def prob13():
    print("\n---- Problem 13")
    error = 10**(-3)
    initial = np.matrix([[2*np.pi/3], [0], [np.pi/6], [0]])
    initialError = np.matrix([[2*np.pi/3 + error], [0], [np.pi/6 + error], [0]])
    t = 40
    n = 10000
    L1 = 2
    L2 = 2
    m1 = 3
    m2 = 1
    g = 9.81
    constantArray = [L1, L2, m1, m2, g]
    # subplots for 4 subplots
    options_for_L1 = [1,1.5,2,2.5,3]
    results_for_L1 = prob13Helper(initial, initialError, n, t, constantArray, 0, options_for_L1)
    options_for_L2 = [1,1.5,2,2.5,3]
    results_for_L2 = prob13Helper(initial, initialError, n, t, constantArray, 1, options_for_L2)
    options_for_m1 = [1,1.5,2,2.5,3,3.5,4,4.5,5]
    results_for_m1 = prob13Helper(initial, initialError, n, t, constantArray, 2, options_for_m1)
    options_for_m2 = [1,1.5,2,2.5,3,3.5,4,4.5,5]
    results_for_m2 = prob13Helper(initial, initialError, n, t, constantArray, 3, options_for_m2)
    options_for_g = [9.81,1.62,3.7,8.87,275,24.5]
    results_for_g = prob13Helper(initial, initialError, n, t, constantArray, 4, options_for_g)
    bestL1Index = np.argmax(results_for_L1)
    bestL2Index = np.argmax(results_for_L2)
    bestM1Index = np.argmax(results_for_m1)
    bestM2Index = np.argmax(results_for_m2)
    bestGIndex = np.argmax(results_for_g)
    maxArr = np.array([np.max(results_for_L1), np.max(results_for_L2), 0, np.max(results_for_m2)])
    maxArr2 = np.array([np.max(results_for_L1), np.max(results_for_L2), np.max(results_for_m1), np.max(results_for_m2), np.max(results_for_g)])
    bestbest = np.argmax(maxArr)
    #[L1, L2, m1, m2, g]
    constantArrayBest = [options_for_L1[bestL1Index], options_for_L2[bestL2Index], options_for_m1[bestM1Index], options_for_m2[bestM2Index], options_for_g[bestGIndex]] 
    # print("Best: L1={}, L2={}, m1={}, m2={}, g={}".format(bestL1Index, bestL2Index, bestM1Index, bestM2Index, bestGIndex))
    print("Best: L1={}, L2={}, m1={}, m2={}, g={}, bestIndex={}, allBest={}".format(constantArrayBest[0], constantArrayBest[1], constantArrayBest[2], constantArrayBest[3], constantArrayBest[4],bestbest , maxArr2))
    
    t2 = 80
    n2 = 2000
    theta1_1 = RungeKuttaModified(n2, t2, initial, F3, L1 = constantArrayBest[0], L2 = constantArrayBest[1], m1 = constantArrayBest[2], m2 = constantArrayBest[3], G = g)
    theta2_1 = RungeKuttaModified(n2, t2, initialError, F3, L1 = constantArrayBest[0], L2 = constantArrayBest[1], m1 = constantArrayBest[2], m2 = constantArrayBest[3], G = g)
    animateBestProb13(theta1_1,theta2_1, n2, t2, "problem13_2a", constantArrayBest[0], constantArrayBest[1])

    constantArrayBest2 = constantArray
    constantArrayBest2[bestbest] = constantArrayBest[bestbest]
    theta1_2 = RungeKuttaModified(n2, t2, initial, F3, L1 = constantArrayBest2[0], L2 = constantArrayBest2[1], m1 = constantArrayBest2[2], m2 = constantArrayBest2[3], G = g)
    theta2_2 = RungeKuttaModified(n2, t2, initialError, F3, L1 = constantArrayBest2[0], L2 = constantArrayBest2[1], m1 = constantArrayBest2[2], m2 = constantArrayBest2[3], G = g)
    animateBestProb13(theta1_2,theta2_2, n2, t2, "problem13_2b", constantArrayBest2[0], constantArrayBest2[1])

    fig1, axs1 = plt.subplots(2, 1)
    axs1[0].plot(options_for_L1, results_for_L1[0,:])
    axs1[0].set_ylabel("Time[s]")
    axs1[0].set_xlabel("L1 Length[m]")
    axs1[0].set_title("L1")
    axs1[1].plot(options_for_L2, results_for_L2[0,:])
    axs1[1].set_title("L2")
    axs1[1].set_ylabel("Time[s]")
    axs1[1].set_xlabel("L2 Length[m]")

    fig2, axs2 = plt.subplots(3, 1)
    axs2[0].plot(options_for_m1, results_for_m1[0,:])
    axs2[0].set_ylabel("Time[s]")
    axs2[0].set_xlabel("m1[kg]")
    axs2[0].set_title("M1")
    axs2[1].plot(options_for_m2, results_for_m2[0,:])
    axs2[1].set_ylabel("Time[s]")
    axs2[1].set_xlabel("m2[kg]")
    axs2[1].set_title("M2")
    axs2[2].plot(options_for_g, results_for_g[0,:])
    axs2[2].set_ylabel("Time[s]")
    axs2[2].set_xlabel("g[m/s^2]")
    axs2[2].set_title("Gravity")
    # plt.tight_layout()
    fig1.tight_layout()
    fig2.tight_layout()
    fig1.show()
    fig2.show()
    plt.show()


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
    prob12()
    # prob13()
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

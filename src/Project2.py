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

def plotOnePendulum(theta, n, T, name):
    print(PLOT_TEXT)
    x = np.linspace(0, T, n+1)
    y = theta[0]

    fig, ax = plt.subplots(1,1)
    ax.grid()
    ax.plot(x, y)
    ax.set_xlabel("Time")
    ax.set_ylabel("Angle")
    ax.set_title("Angle of one pendulum over time")
    fig.savefig(PLOT_PATH.format(name))
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

def prob3():
    print("\n---- Problem 3")
    prob3InitialVal = np.matrix([[np.pi/12], [0]])
    theta = euler(500, 20, prob3InitialVal, F1)
    # animateOnePendulum(theta, 500, 20, "eirikur\AnimateProb3.gif")
    plotOnePendulum(theta, 500, 20, "problem3")

def prob4():
    print("\n---- Problem 4")
    prob4InitialVal = np.matrix([[np.pi/2], [0]])
    theta = euler(500, 20, prob4InitialVal, F1)
    animateOnePendulum(theta, 500, 20, "problem4")
    plotOnePendulum(theta, 500, 20, "problem4")


def prob5():
    print("\n---- Problem 5")
    prob5InitialVal = np.matrix([[np.pi/12], [0]])
    theta = RungeKutta(500, 20, prob5InitialVal, F1)
    animateOnePendulum(theta, 500, 20, "problem5_a")

    prob5InitialVal = np.matrix([[np.pi/2], [0]])
    theta = RungeKutta(500, 20, prob5InitialVal, F1)
    animateOnePendulum(theta, 500, 20, "problem5_b")
    
def prob6():
    print("\n---- Problem 6")
    print("Not done")

def prob7():
    print("\n---- Problem 7")
    prob6InitialVal = np.matrix([[np.pi/3], [0], [np.pi/6], [0]])
    theta = euler(200, 8, prob6InitialVal, F2)
    animateTwoPendulums(theta, 200, 8, "problem7")

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
    print("Not done")

def prob10():
    print("\n---- Problem 10")
    print("Not done")
    
def prob11():
    print("\n---- Problem 11")
    print("Not done")

def prob12():
    print("\n---- Problem 12")
    print("Not done")

def prob13():
    print("\n---- Problem 13")
    print("Not done")

if __name__ == "__main__":
    prob1()
    prob2()
    prob3()
    prob4()
    prob5()
    prob6()
    prob7()
    prob8()
    prob9()
    prob10()
    prob11()
    prob12()
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

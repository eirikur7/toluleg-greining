import numpy as np
import matplotlib.pyplot as plt
from matplotlib.animation import PillowWriter
from matplotlib import animation
g = 9.81 #m/s^2 
L1 = L2 = 2 #m
m1 = m2 = 1 #kg

#Single pendulum
def prob1():
    pass

def prob2F(theta):
    L = 2
    z1 = theta[1,0]
    z2 = (-g/L)*np.sin(theta[0,0])
    return np.matrix([[z1],[z2]])

def euler(n, T, initalValues, F):
    h = T/n
    theta = np.zeros((initalValues.size, n+1))
    theta[:,0:1] = initalValues
    for i in range(n):
        theta[:,(i+1):(i+2)] = theta[:,i:(i+1)] + F(theta[:,i:(i+1)])*h

    return theta

def prob2():
    prob2InitialVal = np.matrix([[1], [0]])
    theta = euler(200, 5, prob2InitialVal, prob2F)
    print(theta)

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

def prob3():
    prob3InitialVal = np.matrix([[np.pi/12], [0]])
    theta = euler(500, 20, prob3InitialVal, prob2F)
    animatePendulum(theta, 500, 20, "prob3.gif")

def prob4():
    prob4InitialVal = np.matrix([[np.pi/2], [0]])
    theta = euler(500, 20, prob4InitialVal, prob2F)
    animatePendulum(theta, 500, 20, "prob4.gif")

def prob5():
    prob5InitialVal = np.matrix([[np.pi/12], [0]])
    theta = RungeKutta(500, 20, prob5InitialVal, prob2F)
    animatePendulum(theta, 500, 20, "prob5a.gif")

    prob5InitialVal = np.matrix([[np.pi/2], [0]])
    theta = RungeKutta(500, 20, prob5InitialVal, prob2F)
    animatePendulum(theta, 500, 20, "prob5b.gif")
    

# Double pendulum

#d^2*theta1/dt^2 = theta1'' =
#(m^2*L_1*w_1^2*sin(delta)*cos(delta) + m_2*g*sin(theta_2)*cos(delta)+m_2*L_2*w_2^2*sin(delta)-(m_1+m_2)*g*sin(theta_1))/(m_1+m_2*sin(delta)^2)

#d^2*theta2/dt^2 = theta2'' =
#(-m_2*L_2*w_2^2*sin(delta)*cos(delta) +(m_1+m_2)(g*sin(theta_1)*cos(delta)-L_1*w_1^2*sin(delta)-g*sin(theta_2)))/((m_1+m_2)L_2-m_2*L_2*cos^2(delta))
#delta = theta_1 - theta_2
#w_1 = d*theta_1/dt
#w_2 = d*theta_2/dt

def prob6():
    pass


#Supporting functions
def animatePendulum(theta, n, T, name):
    def animate(i):
        line.set_data([0,x[i]], [0,y[i]])
        ln1.set_data([x[i]], [y[i]])
    x = np.cos(theta[0]-(np.pi/2)) * L1
    y = np.sin(theta[0]-(np.pi/2)) * L1
    
    print("x:{}\ny:{}".format(x, y))
    fig,ax = plt.subplots(1,1)
    ax.grid()
    ln1, = plt.plot([], [], 'ro', markersize=10)
    line, = plt.plot([], [], 'b-', linewidth=3)
    ax.set_ylim(-2.5, 2.5)
    ax.set_xlim(-2.5, 2.5)
    ani = animation.FuncAnimation(fig, animate, frames=n, interval=T)
    #also draw trajectory of pendulum
    ani.save(name, writer='pillow', fps=25)
    print("done")



if __name__ == "__main__":
    # prob2()
    # prob3()
    prob4()
    prob5()
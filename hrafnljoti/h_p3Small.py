import numpy as np
from numpy import linalg as LA
import matplotlib.pyplot as plt


def calculateY(n, boundry, A1, An, B1, Bn):
    h = (boundry[1]-boundry[0])/(n-1)
    a = 1
    b = -(2 + (h**2))
    c = 1
    A_arr = np.zeros((n, n))
    B_arr = np.zeros((n,1))
    B_arr[0] = B1 
    B_arr[-1] = Bn 
    A_arr[0] = A1 
    A_arr[n-1] = An 
    for i in range(1, n-1):
        A_arr[i, (i-1):(i+2)] = [a, b, c]

    print(A_arr)
    print(B_arr)
    y = LA.solve(A_arr, B_arr)
    return y

def calculateActualY(n, c1, c2):
    '''Returns actual y and also x'''
    y2 = np.zeros((n, 1))
    x = np.linspace(0, 1, n)
    for i in range(n):
        y2[i,0] = c1 *(np.e**(x[i])) + c2*(np.e**(-x[i]))

    return y2, x

def prob2(n, boundry):
    A1 = np.zeros((1, n))
    A1[0, 0] = 1
    An = np.zeros((1, n))
    An[0, n-1] = 1
    B1 = np.zeros((1, 1))
    B1[0, -0] = 1
    Bn = np.zeros((1, 1))
    Bn[0, 0] = -1

    y = calculateY(n, boundry, A1, An, B1, Bn)
    y2, x = calculateActualY(n, -0.5820, 1.5820)

    fig, ax = plt.subplots(1, 1)
    ax.plot(x, y, 'b-', linewidth=2, label="Calculated")
    ax.plot(x, y2, 'r--', label="Actual Solution")
    ax.legend()
    ax.set_ylabel("y")
    ax.set_label("x")

    plt.show()

def prob3(n, boundry):
    h = (boundry[1]-boundry[0])/(n-1)
    A1 = np.zeros((1, n))
    A1[0, 0:3] = [-3, 4, -1]
    An = np.zeros((1, n))
    An[0, n-3:n] = [-1, 4, -3]
    B1 = np.zeros((1, 1))
    B1[0, 0] = 0
    Bn = np.zeros((1, 1))
    Bn[0, 0] = -(2*h)
    y = calculateY(n, boundry, A1, An, B1, Bn)
    y2, x = calculateActualY(n, 0.4255, 0.4255)

    fig, ax = plt.subplots(1, 1)
    ax.plot(x, y, 'b-', linewidth=2, label="Calculated")
    ax.plot(x, y2, 'r--', label="Actual Solution")
    ax.legend()
    ax.set_ylabel("y")
    ax.set_label("x")
    plt.show()

def prob4(n, boundry):
    h = (boundry[1]-boundry[0])/(n-1)
    A1 = np.zeros((1, n))
    A1[0, 0:3] = [-(3+2*h), 4, -1]
    An = np.zeros((1, n))
    An[0, n-3:n] = [-1, 4, -(3+2*h)]
    B1 = np.zeros((1, 1))
    B1[0, 0] = 2*h
    Bn = np.zeros((1, 1))
    Bn[0, 0] = 0
    y = calculateY(n, boundry, A1, An, B1, Bn)
    y2, x = calculateActualY(n, 0, -0.5)

    fig, ax = plt.subplots(1, 1)
    ax.plot(x, y, 'b-', linewidth=2, label="Calculated")
    ax.plot(x, y2, 'r--', label="Actual Solution")
    ax.legend()
    ax.set_ylabel("y")
    ax.set_label("x")
    plt.show()



if __name__ == "__main__":
    # prob2(100, [0,1])
    # prob3(100, [0,1])
    prob4(100, [0,1])




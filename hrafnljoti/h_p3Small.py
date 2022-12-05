import numpy as np
from numpy import linalg as LA
import matplotlib.pyplot as plt


def prob2(n, boundry):
    h = (boundry[1]-boundry[0])/(n-1)
    a = 1
    b = -(2 + (h**2))
    c = 1
    initialy1 = 1
    initialyn = -1
    row = n-2
    col = n
    x = np.linspace(0, 1, n)
    print(x.shape)
    # for i in range(n):
    #     x[i, 0] = h*(i)
    A_arr = np.zeros((n, n))
    B_arr = np.zeros((n,1))
    B_arr[0,0] = 1 # Set first B to y1
    B_arr[-1,0] = -1 #Set Last B to yn
    A_arr[0,0] = 1 #A[1,1] is just y1=y1 => A[1,1]=1
    A_arr[n-1, -1] = 1#A[n,n] is just yn=yn => A[n,n]=1
    # A_arr[1, 0:3] = [a, b, c]
    for i in range(1, n-1):
        A_arr[i, (i-1):(i+2)] = [a, b, c]

    # print(A_arr)
    # print(B_arr)
    # print(A_arr.shape)
    # print(B_arr.shape)
    y = LA.solve(A_arr, B_arr)
    # print(y)
    y2 = np.zeros((n, 1))
    for i in range(n):
        y2[i,0] = -0.5820 *(np.e**(x[i])) + 1.5820*(np.e**(-x[i]))

    fig, ax = plt.subplots(1, 1)
    ax.plot(x, y, 'b-', linewidth=2, label="Calculated")
    ax.plot(x, y2, 'r--', label="Actual Solution")
    ax.legend()
    ax.set_ylabel("y")
    ax.set_label("x")
    # plt.plot(x, y)
    # plt.plot(x, y2, 'r--')

    plt.show()

def prob3(n, boundry):
    h = (boundry[1]-boundry[0])/(n-1)
    a = 1
    b = -(2 + (h**2))
    c = 1
    x = np.linspace(0, 1, n)
    
    A_arr = np.zeros((n, n))
    B_arr = np.zeros((n,1))
    B_arr[0,0] = 0 # Set first B to y1
    B_arr[-1,0] = -2*h #Set Last B 
    A_arr[0,0:3] = [-3, 4, -1] 
    A_arr[n-1, n-3:n] = [-1, 4, -3]
    
    for i in range(1, n-1):
        A_arr[i, (i-1):(i+2)] = [a, b, c]
    print(A_arr)
    y = LA.solve(A_arr, B_arr)
    y2 = np.zeros((n, 1))
    for i in range(n):
        y2[i,0] = 0.4255 *(np.e**(x[i])) + 0.4255*(np.e**(-x[i]))

    fig, ax = plt.subplots(1, 1)
    ax.plot(x, y, 'b-', linewidth=2, label="Calculated")
    ax.plot(x, y2, 'r--', label="Actual Solution")
    ax.legend()
    ax.set_ylabel("y")
    ax.set_label("x")
    plt.show()


if __name__ == "__main__":
    # prob2(100, [0,1])
    prob3(100, [0,1])




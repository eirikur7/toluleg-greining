import numpy as np
from numpy import linalg as LA
import matplotlib.pyplot as plt


g = 9.81
L = 2

def prob1(L,y):
    return np.array([y[1],-(L/g)*np.sin(y[0])])


def Euler_step(initalval, T, n):
    h = T/n
    theta = np.zeros((n,2))
    theta[0] = initalval
    for i in range(1,n):
      theta[i] =  theta[i-1] + h*prob1(L, theta[i-1])
    return theta
    

    

if __name__ == "__main__":
    a = Euler_step([np.pi/12, 0], 20, 500)
    for i in range(len(a)):
        print(a[i])
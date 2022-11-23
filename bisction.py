from math import pow
import matplotlib as m

def f(x,y):
    return x**3 - 3*(x**2)*(y**2) + x*y + 1

def g(x):
    return f(2,x)

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

res = bisection(g, -50, 50, 0.0001)
print(res)

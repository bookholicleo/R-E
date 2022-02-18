import math
#D. E.
def f(x,y):
    return (y-1)*(y-2)*(y-6)*(y-7)*(y-4)**2 #f(x,y)=y'(x,y)
def EX(x,y):
    return 0 #write exact function if known

#Runge - Kutta 4th Order
def RK4(x,y,h):
    k1=f(x,y)
    k2=f(x+0.5*h, y+0.5*k1*h)
    k3=f(x+0.5*h, y+0.5*k2*h)   
    k4=f(x+h, y+k3*h)
    k=(k1+2*(k2+k3)+k4)/6
    return x+h, y+k*h

#input
a, b=map(float, input('initial point ').split(' '))
h=float(input('stepsize '))
N=int(input('number of steps '))

x,y=a,b
for i in range(N):
    x,y=RK4(x,y,h)
    print(x,y,EX(x,y))



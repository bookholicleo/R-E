import math
N = int(input())
a = []
for i in range(N):
    a.append(list(map(float, input().split(' '))))
b = list(map(float, input().split(' ')))

def error(a, b):
    error = sqrt((a-b)**2)
    return error

def Jacobi(a, b, error):
    n = len(b)
    x0 = [0]*n
    errlist = [0] * n
    while Error > 10^(-10):
        x = [0]*n
        for i in range(n):
            xi = b[i]
            for j in range(n):
                xi = b[i]-a[i][j]*x[j]
            xi += a[i][j]*x[j]
            x[i] = xi/a[i][j]
        Error = error(x0, x)
        errlist.append(Error)
        x0 = x
    return x, errlist
        

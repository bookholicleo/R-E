#2D FDM
#D.E. at 46
#==========solve sets of equations==========
#gaussian_elimination from CS2
def gaussian_elimination(a, b):
    assert len(a) == len(a[0]) == len(b)
    n = len(b)

    # forward elimination
    for p in range(n):
        # partial pivoting
        q = p
        for i in range(p+1, n):
            if abs(a[i][p]) > abs(a[q][p]):
                q = i
        for j in range(n): # swapping rows p/q
            a[p][j], a[q][j] = a[q][j], a[p][j]
        b[p], b[q] = b[q], b[p]

        if a[p][p] == 0:
            return None  # no solution or infinitely many solutions
        
        # zero out entries below the pivot a[p][p]
        for i in range(p+1, n):
            c = a[i][p] / a[p][p]
            for j in range(p, n):  # why range(p,n) instead of range(n)?
                a[i][j] -= c * a[p][j]
            b[i] -= c * b[p]

    # now in upper triangular form
    # back substitution
    x = [None] * n
    for i in range(n-1, -1, -1):
        sum = 0
        for j in range(i+1, n):
            sum += a[i][j] * x[j]
        x[i] = (b[i] - sum) / a[i][i]
        
    return x

#input
x0,x1,y0,y1=map(float, input('a, b, c, d: ').split(' '))
N,M=input('N,M:').split(' ')
N,M=int(N),int(M)

#==========D.E.==========
#auxx+buyy+cux+duy+eu=f
#u(a,y)=p(y)
#u(b,y)=q(y)
#u(x,c)=r(x)
#u(x,d)=s(x)
def a(x,y):
    return 0 #type a
def b(x,y):
    return 0 #type b
def c(x,y):
    return 0 #type c
def d(x,y):
    return 0 #type d
def e(x,y):
    return 0 #type e
def f(x,y):
    return 0 #type f
def p(y):
    return 0
def q(y):
    return 0
def r(x):
    return 0
def s(x):
    return 0
def ex(x,y):
    return 0 #type exact function if known
#==========compute==========


h,k=(x1-x0)/N,(y1-y0)/M
A=[[0]*((N-1)*(M-1)) for i in range((N-1)*(M-1))]
B=[0]*(N-1)*(M-1)

#setmatrix
#A
for i in range(2,N-1):
    for j in range(2,M-1):
        n=(i-1)*(M-1)+(j-1)
        x=x0+h*i
        y=y0+k*j
        A[n][n]=e(x,y)-2*(a(x,y)/(h*h))-2*(b(x,y)/(k*k))
        A[n][n-1]=(b(x,y)/(k*k))-(d(x,y)/(2*k))
        A[n][n+1]=(b(x,y)/(k*k))+(d(x,y)/(2*k))
        A[n][n-(M-1)]=(a(x,y)/(h*h))-(c(x,y)/(2*h))
        A[n][n+(M-1)]=(a(x,y)/(h*h))+(c(x,y)/(2*h))
        B[n]=f(x,y)
#boundary
i=1
for j in range(2,M-1):
    n=(i-1)*(M-1)+(j-1)
    x=x0+h*i
    y=y0+k*j
    A[n][n]=e(x,y)-2*(a(x,y)/(h*h))-2*(b(x,y)/(k*k))
    A[n][n-1]=(b(x,y)/(k*k))-(d(x,y)/(2*k))
    A[n][n+1]=(b(x,y)/(k*k))+(d(x,y)/(2*k))
    A[n][n+(M-1)]=(a(x,y)/(h*h))+(c(x,y)/(2*h))
    B[n]=f(x,y)-((a(x,y)/(h*h))-(c(x,y)/(2*h)))*p(y)
i=N-1
for j in range(2,M-1):
    n=(i-1)*(M-1)+(j-1)
    x=x0+h*i
    y=y0+k*j
    A[n][n]=e(x,y)-2*(a(x,y)/(h*h))-2*(b(x,y)/(k*k))
    A[n][n-1]=(b(x,y)/(k*k))-(d(x,y)/(2*k))
    A[n][n+1]=(b(x,y)/(k*k))+(d(x,y)/(2*k))
    A[n][n-(M-1)]=(a(x,y)/(h*h))-(c(x,y)/(2*h))
    B[n]=f(x,y)-((a(x,y)/(h*h))+(c(x,y)/(2*h)))*q(y)
    
j=1
for i in range(2,N-1):
    n=(i-1)*(M-1)+(j-1)
    x=x0+h*i
    y=y0+k*j
    A[n][n]=e(x,y)-2*(a(x,y)/(h*h))-2*(b(x,y)/(k*k))
    A[n][n+1]=(b(x,y)/(k*k))+(d(x,y)/(2*k))
    A[n][n-(M-1)]=(a(x,y)/(h*h))-(c(x,y)/(2*h))
    A[n][n+(M-1)]=(a(x,y)/(h*h))+(c(x,y)/(2*h))
    B[n]=f(x,y)-((b(x,y)/(k*k))-(d(x,y)/(2*k)))*r(x)
j=M-1
for i in range(2,N-1):
    n=(i-1)*(M-1)+(j-1)
    x=x0+h*i
    y=y0+k*j
    A[n][n]=e(x,y)-2*(a(x,y)/(h*h))-2*(b(x,y)/(k*k))
    A[n][n-1]=(b(x,y)/(k*k))-(d(x,y)/(2*k))
    A[n][n-(M-1)]=(a(x,y)/(h*h))-(c(x,y)/(2*h))
    A[n][n+(M-1)]=(a(x,y)/(h*h))+(c(x,y)/(2*h))
    B[n]=f(x,y)-((b(x,y)/(k*k))+(d(x,y)/(2*k)))*s(x)
i,j=1,1
n=(i-1)*(M-1)+(j-1)
x=x0+h*i
y=y0+k*j
A[n][n]=e(x,y)-2*(a(x,y)/(h*h))-2*(b(x,y)/(k*k))
A[n][n+1]=(b(x,y)/(k*k))+(d(x,y)/(2*k))
A[n][n+(M-1)]=(a(x,y)/(h*h))+(c(x,y)/(2*h))
B[n]=f(x,y)-((b(x,y)/(k*k))-(d(x,y)/(2*k)))*r(x)-((a(x,y)/(h*h))-(c(x,y)/(2*h)))*p(y)
i,j=1,M-1
n=(i-1)*(M-1)+(j-1)
x=x0+h*i
y=y0+k*j
A[n][n]=e(x,y)-2*(a(x,y)/(h*h))-2*(b(x,y)/(k*k))
A[n][n-1]=(b(x,y)/(k*k))-(d(x,y)/(2*k))
A[n][n+(M-1)]=(a(x,y)/(h*h))+(c(x,y)/(2*h))
B[n]=f(x,y)-((b(x,y)/(k*k))+(d(x,y)/(2*k)))*s(x)-((a(x,y)/(h*h))-c(x,y)/(2*h))*p(y)
i,j=N-1,1
n=(i-1)*(M-1)+(j-1)
x=x0+h*i
y=y0+k*j
A[n][n]=e(x,y)-2*(a(x,y)/(h*h))-2*(b(x,y)/(k*k))
A[n][n+1]=(b(x,y)/(k*k))+(d(x,y)/(2*k))
A[n][n-(M-1)]=(a(x,y)/(h*h))-(c(x,y)/(2*h))
B[n]=f(x,y)-((b(x,y)/(k*k))-(d(x,y)/(2*k)))*r(x)-((a(x,y)/(h*h))+(c(x,y)/(2*h)))*q(y)
i,j=N-1,M-1
n=(i-1)*(M-1)+(j-1)
x=x0+h*i
y=y0+k*j
A[n][n]=e(x,y)-2*(a(x,y)/(h*h))-2*(b(x,y)/(k*k))
A[n][n-1]=(b(x,y)/(k*k))-(d(x,y)/(2*k))
A[n][n-(M-1)]=(a(x,y)/(h*h))-(c(x,y)/(2*h))
B[n]=f(x,y)-((b(x,y)/(k*k))+(d(x,y)/(2*k)))*s(x)-((a(x,y)/(h*h))+(c(x,y)/(2*h)))*q(y)



for i in range((N-1)*(M-1)):
    print(A[i])
print('')
print(B)
Y=gaussian_elimination(A,B)
print(Y)

#output
Y=gaussian_elimination(A, B)
for i in range(1,N):
    for j in range(1,M):
        n=(i-1)*(M-1)+(j-1)
        x=x0+h*i
        y=y0+k*j
        print(x,y,Y[n],ex(x,y))
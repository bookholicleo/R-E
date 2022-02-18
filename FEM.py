#-(sux)x-(tuy)y+du=f
def s(x,y):
    return 1
def t(x,y):
    return 1
def d(x,y):
    return 0
def f(x,y):
    return -2

#border
def g(x,y):
    return x*x+2*y

#input
M,N=map(int, input().split(' '))
a,B=map(float, input().split(' '))
c,D=map(float, input().split(' '))
h=(B-a)/M
k=(D-c)/N
#FEM
def P(x,y,num):
    val=[(1-x)*(1-y)/4,(1+x)*(1-y)/4,(1+x)*(1+y)/4,(1-x)*(1+y)/4]
    return val[num]
def Px(x,y,num):
    val=[(y-1)/4,(1-y)/4,(1+y)/4,(-1-y)/4]
    return val[num]
def Py(x,y,num):
    val=[(x-1)/4,(-1-x)/4,(1+x)/4,(1-x)/4]
    return val[num]
#locb
def locb(i,j):
    global h
    global k
    n=N*i+j
    x=a+(i+0.5)*h                                                                   
    y=c+(j+0.5)*k
    val=(0.6)**0.5
    INT=[0]*4
    for p in range(4):
        INT[p] += 25*(P(-val,-val,p)*f(x-val*h/2,y-val*k/2)+P(val,-val,p)*f(x+val*h/2,y-val*k/2)+P(-val,val,p)*f(x-val*h/2,y+val*k/2)+P(val,val,p)*f(x+val*h/2,y+val*k/2))
        INT[p] += 40*(P(-val,0,p)*f(x-val*h/2,y)+P(val,0,p)*f(x+val*h/2,y)+P(0,-val,p)*f(x,y-val*k/2)+P(0,val,p)*f(x,y+val*k/2))
        INT[p] += 64*P(0,0,p)*f(x,y)
        INT[p]=(INT[p]/81)*h*k/4
    return INT

#locs, loct, locd
def loc(i,j):
    global h
    global k    
    n=N*i+j
    x=a+(i+0.5)*h
    y=c+(j+0.5)*k
    val=(0.6)**0.5  
    INT=[[0]*4 for x in range(4)]
    for p in range(4):
        for q in range(4):
            temp = 25*(Px(-val, -val, p)*Px(-val, -val, q)*s(x-val*h/2,y-val*k/2)+Px(val, -val, p)*Px(val, -val, q)*s(x+val*h/2,y-val*k/2)+Px(-val, val, p)*Px(-val, val, q)*s(x-val*h/2,y+val*k/2)+Px(val, val, p)*Px(val, val, q)*s(x+val*h/2,y+val*k/2))
            temp+=40*(Px(-val, 0, p)*Px(-val, 0, q)*s(x-val*h/2,y)+Px(val, 0, p)*Px(val, 0, q)*s(x+val*h/2,y)+Px(0, -val, p)*Px(0, -val, q)*s(x,y-val*k/2)+Px(0, val, p)*Px(0, val, q)*s(x,y+val*k/2))
            temp += 64*Px(0,0,p)*Px(0,0,q)*s(x,y)
            INT[p][q] += temp/81
                
            temp = 25*(Py(-val, -val, p)*Py(-val, -val, q)*t(x-val*h/2,y-val*k/2)+Py(val, -val, p)*Py(val, -val, q)*t(x+val*h/2,y-val*k/2)+Py(-val, val, p)*Py(-val, val, q)*t(x-val*h/2,y+val*k/2)+Py(val, val, p)*Py(val, val, q)*t(x+val*h/2,y+val*k/2))
            temp+=40*(Py(-val, 0, p)*Py(-val, 0, q)*t(x-val*h/2,y)+Py(val, 0, p)*Py(val, 0, q)*t(x+val*h/2,y)+Py(0, -val, p)*Py(0, -val, q)*t(x,y-val*k/2)+Py(0, val, p)*Py(0, val, q)*t(x,y+val*k/2))
            temp += 64*Py(0,0,p)*Py(0,0,q)*t(x,y)
            INT[p][q] += temp/81
                
            temp = 25*(P(-val, -val, p)*P(-val, -val, q)*d(x-val*h/2,y-val*k/2)+P(val, -val, p)*P(val, -val, q)*d(x+val*h/2,y-val*k/2)+P(-val, val, p)*P(-val, val, q)*d(x-val*h/2,y+val*k/2)+P(val, val, p)*P(val, val, q)*d(x+val*h/2,y+val*k/2))
            temp+=40*(P(-val, 0, p)*P(-val, 0, q)*d(x-val*h/2,y)+P(val, 0, p)*P(val, 0, q)*d(x+val*h/2,y)+P(0, -val, p)*P(0, -val, q)*d(x,y-val*k/2)+P(0, val, p)*P(0, val, q)*d(x,y+val*k/2))
            temp += 64*P(0,0,p)*P(0,0,q)*d(x,y)
            INT[p][q] += (temp/81)*h*k/4
    return INT

#set matrix
pos=[[-1,-1],[0,-1],[0,0],[-1,0]]
b=[0]*((M-1)*(N-1))
A=[[0]*9 for i in range((M-1)*(N-1))]
def boundary(i,j):
    return (i==-1 or i==M-1) or (j==-1 or j==N-1)

def setMatrix():    
    for i in range(M):
        for j in range(N):
            n=i*(N-1)+j
            temp=locb(i,j)
            point1=[n-N, n-1, n, n-N+1]
            for p in range(4):
                if not boundary(i+pos[p][0], j+pos[p][1]):
                    b[point1[p]] += temp[p]
                
    for i in range(M):
        for j in range(N):
            n=i*(N-1)+j
            temp=loc(i,j)
            point1=[n-N, n-1, n, n-N+1]
            point2=[[4,5,8,7],[3,4,7,6],[0,1,4,3],[1,2,5,4]]            
            for p in range(4):
                if not boundary(i+pos[p][0], j+pos[p][1]):
                    for q in range(4):
                        if boundary(i+pos[q][0], j+pos[q][1]):
                            x=a+(1+i+pos[q][0])*h
                            y=c+(1+j+pos[q][1])*k
                            b[point1[p]] -= g(x,y)*temp[p][q]
                        else:    
                            A[point1[p]][point2[p][q]] += temp[p][q]
                    
                
#solve
#multiplication with A
def mult_A(v):
    prod = [0]*(len(v))
    for i in range(M-1):
        for j in range(N-1):
            n=i*(N-1)+j
            for p in range(-1,2):
                for q in range(-1,2):
                    n2=(i+q)*(N-1)+(j+p)
                    if not boundary(i+q, j+p):
                        prod[n] += A[n][3*(p+1)+q+1]*v[n2]
    return prod

#multiplication with scalar
def mult_c(c,v):
    vect=v.copy()
    for i in range(len(v)):
        vect[i] *= c
    return vect

#subtraction
def sub(u,v):
    sub=u.copy()
    for i in range(len(u)):
        sub[i] -= v[i]
    return sub
    
#inner product
def prod_std(u,v):
    count=0
    for i in range(len(u)):
        count += v[i]*u[i]
    return count

def prod(u,v):
    return prod_std(u,mult_A(v))

#compute -gradient
def r(x):
    return sub(b,mult_A(x))

#get conjugate direction
def getp(x,prev):
    grad=r(x)
    p=grad.copy()
    for i in range(len(prev)):
        p=sub(p, mult_c(prod(prev[i],grad)/prod(prev[i],prev[i]),prev[i]))
    return p

#Conjugate Gradient method

def solve(maxloops,mindiff):
    prev=[]
    x=[0]*((M-1)*(N-1))
    for i in range(maxloops):
        grad=r(x)
        p=getp(x,prev)
        if prod_std(grad,grad) <mindiff:
            return x
        prev.append(p)
        a=prod_std(p,grad)/prod(p,p)
        x=sub(x,mult_c(-a,p))
    return x

#==========================================================================
#==========================================================================
setMatrix()
print(b)
ans=solve(100,10**-10)
for i in range(M-1):
    for j in range(N-1):
        n=i*(N-1)+j
        x,y=a+(i+1)*h, c+(j+1)*k
        print((ans[n],g(x,y)),end=' ')
    print('')
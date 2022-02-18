import time

#-(sux)x-(tuy)y+du=f
def s(x,y):
    return 1
def t(x,y):
    return 1
def d(x,y):
    return 1
def f(x,y):
    return x*x+2*y-2

#border
def g(x,y):
    return x*x+2*y


#FEM subfunction
def integral(I): #I_3
    val = (0.6)**0.5
    INTEGRAL = 25*(I(-val,-val)+I(-val,val)+I(val,-val)+I(val,val))
    INTEGRAL += 40*(I(-val,0)+I(val,0)+I(0,-val)+I(0,val))
    INTEGRAL += 64*(I(0,0))
    INTEGRAL /= 81
    return INTEGRAL

#basis function
def P(x,y,num):
    val=[(1-x)*(1-y)/4,(1+x)*(1-y)/4,(1+x)*(1+y)/4,(1-x)*(1+y)/4]
    return val[num]
def Px(x,y,num):
    val=[(y-1)/4,(1-y)/4,(1+y)/4,(-1-y)/4]
    return val[num]
def Py(x,y,num):
    val=[(x-1)/4,(-1-x)/4,(1+x)/4,(1-x)/4]
    return val[num]

def locb(M,N,i,j):
    global h, k
    n=N*i+j
    x=x1+(i+0.5)*h                                                                   
    y=y1+(j+0.5)*k
    INT=[None]*4
    for q in range(4):
        def I(u,v):
            return P(u,v,q)*f(x+u*h/2,y+v*k/2)*h*k/4
        INT[q]=integral(I)
    return INT

def loc(M,N,i,j):
    global h, k
    n=N*i+j
    x=x1+(i+0.5)*h
    y=y1+(j+0.5)*k
    INT=[[None]*4 for asdf in range(4)]
    for p in range(4):
        for q in range(4):
            def I(u,v):
                val=P(u,v,p)*P(u,v,q)*d(x+u*h/2,y+v*k/2)*h*k/4
                val += Px(u,v,p)*Px(u,v,q)*s(x+u*h/2,y+v*k/2)
                val += Py(u,v,p)*Py(u,v,q)*t(x+u*h/2,y+v*k/2)
                return val
            INT[p][q]=integral(I)
    return INT

#set matrix
pos=[[-1,-1],[0,-1],[0,0],[-1,0]]

def boundary(M,N,i,j):
    return (i==-1 or i==M-1) or (j==-1 or j==N-1)

def setMatrix(M,N):
    b=[0]*((M-1)*(N-1))
    A=[[0]*9 for i in range((M-1)*(N-1))]    
    for i in range(M):
        for j in range(N):
            n=i*(N-1)+j
            temp=locb(M,N,i,j)
            point1=[n-N, n-1, n, n-N+1]
            for p in range(4):
                if not boundary(M,N,i+pos[p][0], j+pos[p][1]):
                    b[point1[p]] += temp[p]
                
    for i in range(M):
        for j in range(N):
            n=i*(N-1)+j
            temp=loc(M,N,i,j)
            point1=[n-N, n-1, n, n-N+1]
            point2=[[4,5,8,7],[3,4,7,6],[0,1,4,3],[1,2,5,4]]            
            for p in range(4):
                if not boundary(M,N,i+pos[p][0], j+pos[p][1]):
                    for q in range(4):
                        if boundary(M,N,i+pos[q][0], j+pos[q][1]):
                            x=x1+(1+i+pos[q][0])*h
                            y=y1+(1+j+pos[q][1])*k
                            b[point1[p]] -= g(x,y)*temp[p][q]
                        else:    
                            A[point1[p]][point2[p][q]] += temp[p][q]
    return A,b
                            
#Conjugate Gradient Method
#multiplication with A
def mult_A(M,N,A,v):
    prod = [0]*(len(v))
    for i in range(M-1):
        for j in range(N-1):
            n=i*(N-1)+j
            for p in range(-1,2):
                for q in range(-1,2):
                    n2=(i+q)*(N-1)+(j+p)
                    if not boundary(M,N,i+q, j+p):
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

def prod(M,N,A,u,v):
    return prod_std(u,mult_A(M,N,A,v))

#compute -gradient
def r(M,N,A,b,x):
    return sub(b,mult_A(M,N,A,x))


#Multigrid method
def interpolate(M,N,T,x,y):
    h=(x2-x1)/M
    k=(y2-y1)/N    
    value=0
    i,j=int((x-x1)/h)-1,int((y-y1)/k)-1
    dx,dy=((x-x1)/h)%1,((y-y1)/k)%1
    n=(N-1)*i+j
    pos=[[i,j],[i+1,j],[i+1,j+1],[i,j+1]]
    idx=[n, n+N-1, n+N, n+1]
    for p in range(4):
        if boundary(M,N,pos[p][0], pos[p][1]):
            value += g(x1+(pos[p][0]+1)*h,y1+(pos[p][1]+1)*k)*P(2*dx-1,2*dy-1,p)
        else:
            value += T[idx[p]]*P(2*dx-1,2*dy-1,p)
    return value


def SOLVE(M,N,maxloops,mindiff): #conjugate gradient method
    global h,k
    h=(x2-x1)/M
    k=(y2-y1)/N    
    A,b=setMatrix(M,N)
    ans=[0]*((M-1)*(N-1))
    if len(ans) > 100: # ´ëÃæ ³·Àº ¼ö
        prev=SOLVE(M//2,N//2,maxloops,mindiff)                
        for i in range(1,M):
            for j in range(1,N):
                h=(x2-x1)/M
                k=(y2-y1)/N                
                x,y=x1+h*i, y1+k*j
                n=n=(N-1)*(i-1)+(j-1)
                ans[n]=interpolate(M//2,N//2,prev,x,y)
        h=(x2-x1)/M
        k=(y2-y1)/N                    
        
    grad,p=r(M,N,A,b,ans),r(M,N,A,b,ans)

    for i in range(maxloops):
        a=prod_std(p,grad)/prod(M,N,A,p,p)
        ans=sub(ans,mult_c(-a,p))
        b=1/prod_std(grad,grad)
        grad=sub(grad,mult_c(a,mult_A(M,N,A,p)))
        if prod_std(grad,grad) < mindiff:
            return ans        
        b *= prod_std(grad,grad)
        p = sub(grad, mult_c(-b,p))
        
        
    return ans
#==========================================================================
#==========================================================================
#input
M0,N0=map(int, input().split(' '))
x1,x2=map(float, input().split(' '))
y1,y2=map(float, input().split(' '))
t1=time.time()
ans=SOLVE(M0,N0,M0*N0, 10**-10)
t2=time.time()
for i in range(M0-1):
    for j in range(N0-1):
        n=i*(N0-1)+j
        x,y=x1+(i+1)*h, y1+(j+1)*k
        print((ans[n],g(x,y)),end=' ')
    print('')

print(t2-t1)
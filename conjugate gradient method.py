#conjugate gradient method
import numpy as np
from matplotlib import pyplot as plt
from mpl_toolkits import mplot3d
from copy import deepcopy

A = np.array([[4, 1],[1, 3]], dtype=float)
b = np.array([1, 2], dtype=float)
c = 0.
print(A.shape, b.shape, c)
n = len(A[0])

def f(x):
    """ x is a vector, of shape=[n]"""
    return (1 / 2) * np.matmul(np.matmul(x.T, A), x) - np.matmul(b.T, x) + c

def steepest_grad_decent(A, b, init, step, history=False):
    x = deepcopy(init)
    memo = [x]
    bound = 1e-7
    while step - 1:
        r = b - np.matmul(A, x)  # -graadient of f(x_i), of shape=[n]
        # find optimal learning rate for the current point.
        alpha = np.matmul(r.T, r) / np.matmul(np.matmul(r.T, A), r)  # scalar
        x = x + alpha * r
        if bound > sum(np.abs(alpha * r)): break
        if history: memo.append(x)
        step -= 1
    if not history: return x
    return x, np.array(list(zip(*memo)))

def grad_decent(f, init, step, lr=0.001, history=False):
    
    def gradient(f, x, epsilon=1e-7):
        """ numerically find gradients. 
        x: shape=[n] """
        grad = np.zeros_like(x, dtype=float)
        for i in range(len(x)):
            h = np.zeros_like(x, dtype=float)
            h[i] = epsilon
            grad[i] = (f(x + h) - f(x - h)) / (2 * h[i])
        return grad
    
    x = deepcopy(init)
    memo = [x]
    for i in range(step - 1):
        # grad = gradient(f, x)
        grad = np.matmul(A, x) - b  # -graadient of f(x_i), of shape=[n]
        x = x - lr * grad
        if history: memo.append(x)
    if not history: return x
    return x, np.array(list(zip(*memo)))

def conjugate(A, b, init, step, history=False):
    x = deepcopy(init)
    memo = [x]
    bound = 1e-7
    r = b - np.matmul(A, x)  # -graadient of f(x_i), of shape=[n]
    p = deepcopy(r) # inital search direction
    
    while step - 1:
        # line search
        pap = np.matmul(np.matmul(p.T, A), p)  # scalar
        alpha = np.matmul(r.T, r) / pap  # scalar
        
        # update estimate 
        x = x + alpha * p  # [n]
        if bound > sum(np.abs(alpha * p)): break
        
        # update residual and search direction
        r_prev = deepcopy(r)
        r = b - np.matmul(A, x)

        # r = r - np.matmul(alpha * A, p)  # [n]
        # beta = np.matmul(np.matmul(r.T, A), p) / pap  # scalar
        beta = np.matmul(r.T, r) / np.matmul(r_prev.T, r_prev)
        p = r + beta * p  # [n]
        if history: memo.append(x)
        step -= 1
    if not history: return x
    return x, np.array(list(zip(*memo)))

num_steps = 100
init = np.array([2, 1], dtype=float)
ans, history = conjugate(A, b, init, step=num_steps, history=True)
print(f(init), f(ans))
print(ans)
from numpy import *
import sympy as sy
from sympy.abc import alpha, beta, gamma, L, B, T

def A(alpha,beta,gamma,L,B,T):
    # Dimensions of cuboid should be L>B>T
    
    # Angles in radians 
    alpha = deg2rad(alpha)
    beta  = deg2rad(beta)
    gamma = deg2rad(gamma)

    # Unit vectors
    i = [1,0,0]
    j = [0,1,0]
    k = [0,0,1]
    
    # Rotational matrices
    X = sy.Matrix([i,[0,sy.cos(alpha),-sy.sin(alpha)],[0,sy.sin(alpha),sy.cos(alpha)]])
    Y = sy.Matrix([[sy.cos(beta),0,sy.sin(beta)],j,[-sy.sin(beta),0,sy.cos(beta)]])
    Z = sy.Matrix([[sy.cos(gamma),-sy.sin(gamma),0],[sy.sin(gamma),sy.cos(gamma),0],k])
    
    R = Z*Y*X
    
    # Coordinates of the vertices before and after rotation
    N1 = sy.Matrix([ L/2, B/2, T/2])
    N2 = sy.Matrix([ L/2,-B/2, T/2])
    N3 = sy.Matrix([-L/2,-B/2, T/2])
    N4 = sy.Matrix([-L/2, B/2, T/2])
    N5 = sy.Matrix([-L/2, B/2,-T/2])
    N6 = sy.Matrix([ L/2, B/2,-T/2])
    N7 = sy.Matrix([ L/2,-B/2,-T/2])
    
    R1 = R*N1
    R2 = R*N2
    R3 = R*N3
    R4 = R*N4
    R5 = R*N5
    R6 = R*N6
    R7 = R*N7
    
    l = R1-R4
    b = R1-R2
    t = R1-R6
    
    k = sy.Matrix([0,0,1])
    lp = l - l.dot(k)*k
    bp = b - b.dot(k)*k
    tp = t - t.dot(k)*k
    
    Mx = sy.Matrix([bp.transpose(),tp.transpose()])
    My = sy.Matrix([lp.transpose(),tp.transpose()])
    Mz = sy.Matrix([lp.transpose(),bp.transpose()])
    
    Mx = Mx*Mx.transpose()
    My = My*My.transpose()
    Mz = Mz*Mz.transpose()
    
    Ax = sy.sqrt(Mx.det())
    Ay = sy.sqrt(My.det())
    Az = sy.sqrt(Mz.det())

    return Ax + Ay + Az
    
n = 90/5+1
Area = zeros(shape=(n,n,n))

for a in range(0,n):
    for b in range(0,n):
        for c in range(0,n):
            Area[a,b,c] = real(complex(A(5*a,5*b,5*c,7,6,4)))

save("area.npy", Area)
        
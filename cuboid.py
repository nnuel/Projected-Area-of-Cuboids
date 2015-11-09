from __future__ import division
from area import A

from mpl_toolkits.mplot3d import Axes3D
import matplotlib.pyplot as plt

from numpy import linalg as LA
from numpy import *
from numpy import max as mx
import scipy as sp
import sympy as sy

from sympy.abc import alpha, beta, gamma, L, B, T

L = 7
B = 6
T = 4

Dxyz = sy.Matrix([ L/2, B/2, T/2])
D    = sy.Matrix([ L/2, B/2,-T/2])

N   = Dxyz.cross(D)
N   = N/N.norm()

alpha = rad2deg(math.acos(N[0]))
beta  = rad2deg(math.acos(N[1]))
gamma = rad2deg(math.acos(N[2]))

print(A(alpha,beta,gamma,L,B,T))

Area = load("area.npy")

alpha = list(xrange(0,91,5))
beta  = list(xrange(0,91,5))
gamma  = list(xrange(0,91,5))

alpha, beta = meshgrid(alpha, beta)

ax = Axes3D(plt.gcf())
ax.plot_surface(alpha, beta, Area[:,:,0])

fig = plt.figure()
ax = fig.add_subplot(111, projection='3d')
ax.plot_surface(alpha, beta, Area[:,0,:], rstride=1, cstride=1)
cset = ax.contour(alpha, beta, Area[:,0,:], zdir='x', offset=0)
cset = ax.contour(alpha, beta, Area[:,0,:], zdir='y', offset=90)
cset = ax.contour(alpha, beta, Area[:,0,:], zdir='z', offset=0)
ax.set_xlabel('alpha')
ax.set_xlim(0, 90)
ax.set_ylabel('gamma')
ax.set_ylim(0, 90)
ax.set_zlabel('area [m2]')
ax.set_zlim(0, 60)
plt.show()

plt.plot(gamma, Area[0,:,0])
plt.plot(gamma, Area[:,-1,0])
plt.plot(gamma, Area[0,0,:])
plt.plot([rad2deg(arctan(6/7)),rad2deg(arctan(6/7))],[24,38])
plt.show()

mxA = mx(Area)
A = where(Area==mxA)
print(Area[7,5,5])

# alpha = 90
# beta  = 0
# gamma = rad2deg(arctan(4/7))
# print(A(alpha,beta,gamma))

# alpha = list(xrange(0,91,5))
# beta  = list(xrange(0,91,5))
# ax = Axes3D(plt.gcf())
# ax.contourf(alpha, beta, Area[:,:,0])
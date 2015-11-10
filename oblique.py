from __future__ import division

import matplotlib.pyplot as plt
from scipy.interpolate import spline

from numpy import *
import scipy as sp
import sympy as sy

rho = 1000      # [kg/m^3]
nu  = 1.35E-6   # [m^2/s]
Re  = 4.7E4     # [-]

Ay  = 1         # [m^2]
AR  = 1         # [-] = Ax/Ay
V   = 2         # [m/s]

angle   = arange(0,90+1,5)              # [deg]
Cd_a    = array([2.2,2.1,1.8,1.3,1.9,2.1,2.2,2.3,2.4,2.4])
Cd_a    = append(Cd_a,Cd_a[:-1][::-1])

plt.figure()
plt.plot(angle,Cd_a)
plt.ylim(0, 2.5)
plt.ylabel('Drag coefficient of projected area $C_d(\\alpha) \ [-]$')
plt.xlabel('Angle of attack $\\alpha \ [\degree]$')
plt.savefig('Cd_a.eps', format='eps', dpi=1200)

angleR  = deg2rad(angle)                # [rad]
cos_a   = cos(angleR)                   # [-]
sin_a   = sin(angleR)                   # [-]
Vx      = V*cos_a                       # [m/s]
Vy      = V*sin_a                       # [m/s]
A_a     = Ay*(AR*cos_a + sin_a)         # [m^2]

Cdx     = Cd_a[0]                       # [-]
Cdy     = Cd_a[0]                       # [-]
Cdx_a   = A_a/(AR*Ay)*Cd_a              # [-]
Cdy_a   = A_a/Ay*Cd_a                   # [-]

OFdx    = 0.5*rho*Cdx*AR*Ay*Vx*abs(V)   # [N]
OFdy    = 0.5*rho*Cdy*Ay*Vy*abs(V)      # [N]
OFd     = sqrt(OFdx**2+OFdy**2)         # [N]

Fd      = 0.5*rho*Cd_a*A_a*V*abs(V)     # [N]

WFdx    = 0.5*rho*Cdx_a*AR*Ay*Vx*abs(V) # [N]
WFdy    = 0.5*rho*Cdy_a*Ay*Vy*abs(V)    # [N]
WFd     = sqrt(WFdx**2+WFdy**2)         # [N]

OFd_a   = rad2deg(arctan(OFdy/OFdx))    # [deg]
WFd_a   = rad2deg(arctan(WFdy/WFdx))    # [deg]

plt.figure()
plt.plot(angle,OFd)
plt.plot(angle,Fd)
plt.plot(angle,WFd)
plt.ylim(0, math.ceil(Fd.max()/500)*500)
plt.ylabel('Total Drag Force $F_d$ on the object $[N]$')
plt.xlabel('Angle of attack $\\alpha \ [\degree]$')
plt.savefig('Fd.eps', format='eps', dpi=1200)

plt.figure()
plt.plot(angle,OFd_a)
plt.plot(angle,WFd_a)
plt.title('Direction check of drag force $F_d \ [\degree]$')
plt.ylabel('Angle of Total Drag Force $F_d$ on the object $[\degree]$')
plt.xlabel('Angle of attack $\\alpha \ [\degree]$')
plt.savefig('Fd_a.eps', format='eps', dpi=1200)

plt.figure()
plt.plot(angle,Ay*sin_a)
plt.plot(angle,Ay*AR*cos_a)
plt.plot(angle,A_a)
plt.ylim(0, sqrt(2))
plt.title('Projected area and both components')
plt.ylabel('Projected area $[m^2]$')
plt.xlabel('Angle of attack $\\alpha \ [\degree]$')
plt.savefig('A_a.eps', format='eps', dpi=1200)

plt.show()
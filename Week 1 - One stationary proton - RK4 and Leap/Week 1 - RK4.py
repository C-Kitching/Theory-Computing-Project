#Imports
import numpy as np
import matplotlib.pyplot as plt

#Function to calculate derivative
def deriv(r,t):
    return np.array([r[2],r[3],-r[0]/((r[0]**2+r[1]**2)**(3/2)),-r[1]/((r[0]**2+r[1]**2)**(3/2))])

#Function to do single RK4 step
def rk4(r,dr,t,h):
    k1=dr(r,t)
    k2=dr(r+h/2*k1,t+h/2)
    k3=dr(r+h/2*k2,t+h/2)
    k4=dr(r+h*k3,t+h)
    r=r+h*(k1+2*k2+2*k3+k4)/6
    t=t+h
    return (t,r)

r=np.array([10,0,0,0.5])  # initial x,y displacements; initial x,y velocities
Ro=np.array([r[0],r[1],r[2],r[3]])
t=0  # time start
tmax=500  # time end
h=1  # step size
ts=np.array([t])  # store all times
rs=np.array([r])  # store all solution points

for i in range(int(tmax/h)):
    (t,r)=rk4(r,deriv,t,h)
    ts=np.append(ts,t)
    rs=np.concatenate((rs,np.array([r])))
        
[x,y,Vx,Vy]=rs.transpose()

L=((Ro[0]**2+Ro[1]**2)**(1/2))*((Ro[2]**2+Ro[3]**2)**(1/2))
initialE=(1/2*Ro[2]**2+1/2*Ro[3]**2)-1/((Ro[0]**2+Ro[1]**2)**(1/2))
print('inital Energy =')
print(initialE)
eccentricity = (1+2*(L**2)*initialE)**(1/2)
print('Eccentricity =')
print(eccentricity)

plt.figure()
plt.plot(x,y)
plt.xlabel('x-Displacement')
plt.ylabel('y-Displacement')
plt.axis('equal')
plt.title('x-y Path - RK4 - Hyperbolic')

"""
plt.figure()
plt.plot(ts,x)
plt.xlabel('time')
plt.ylabel('x displacement')
plt.title('x-Displacement vs Time - RK4 - Elliptical')

plt.figure()
plt.plot(ts,y)
plt.xlabel('time')
plt.ylabel('y displacement')
plt.title('y-Displacement vs Time - RK4 - Elliptical')
"""

plt.figure()
plt.plot(ts,y, label = 'y-displacement')
plt.plot(ts,x, label = 'x-displacement')
plt.title('Displacement vs Time - RK4 - Hyperbolic')
plt.xlabel('Time')
plt.ylabel('Displacement')
plt.legend()

plt.figure()
plt.plot(ts,Vy, label = 'y-velocity')
plt.plot(ts,Vx, label = 'x-velocity')
plt.title('Velocity vs Time - RK4 - Hyperbolic')
plt.xlabel('Time')
plt.ylabel('Velocity')
plt.legend()

V=-1/((x**2+y**2)**(1/2))
T=1/2*Vx**2+1/2*Vy**2
E=T+V

plt.figure()
plt.plot(ts,V, label = 'Potential Energy')
plt.plot(ts,T, label = 'Kinetic Energy')
plt.plot(ts,E, label = 'Total Energy')
plt.xlabel('Time')
plt.ylabel('Energy')
plt.title('Energy-Time - RK4 - Hyperbolic')
plt.legend()
plt.show()

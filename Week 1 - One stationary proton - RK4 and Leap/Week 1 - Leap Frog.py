#Imports
import numpy as np
import matplotlib.pyplot as plt

#Function to calculate derivative
def deriv(r):
    return np.array([-r[0]/((r[0]**2+r[1]**2)**(3/2)),-r[1]/((r[0]**2+r[1]**2)**(3/2))])

def LeapFrog(r,v,dt,t):
    v=v+ 1/2 *dt*deriv(r)
    r=r+dt*v
    v=v +1/2*dt*deriv(r)
    t=t+dt
    return (t,r,v)

r=np.array([10,0])  # initial x,y displacements
v=np.array([0,0.47])  #initial x,y velocities
Ro=np.array([r[0],r[1],v[0],v[1]])
t=0  # time start
tmax=500  # time end
h=1  # step size
ts=np.array([t])  # store all times
rs=np.array([r])  # store all solution points
vs=np.array([v])  # store all solution points

for i in range(int(tmax/h)):
    (t,r,v)=LeapFrog(r,v,h,t)
    ts=np.append(ts,t)
    rs=np.concatenate((rs,np.array([r])))
    vs=np.concatenate((vs,np.array([v])))

[x,y]=rs.transpose()
[Vx,Vy]=vs.transpose() 
  
L=((Ro[0]**2+Ro[1]**2)**(1/2))*((Ro[2]**2+Ro[3]**2)**(1/2))
initialE=(1/2*Ro[2]**2+1/2*Ro[3]**2)-1/((Ro[0]**2+Ro[1]**2)**(1/2))
print('inital Energy =')
print(initialE)
eccentricity = (1+2*(L**2)*initialE)**(1/2)
print('Eccentricity =')
print(eccentricity)
 
V=-1/((x**2+y**2)**(1/2))
T=1/2*Vx**2+1/2*Vy**2
E=T+V

plt.figure()
plt.plot(x,y)
plt.xlabel('x-Displacement')
plt.ylabel('y-Displacement')
plt.axis('equal')
plt.title('x-y Path - Leap Frog - Hyperbolic')

"""
plt.figure()
plt.plot(ts,y)
plt.xlabel('time')
plt.ylabel('y displacement')
plt.title('y-Displacement vs Time - Leap Frog - Parabolic')

plt.figure()
plt.plot(ts,x)
plt.xlabel('time')
plt.ylabel('x displacement')
plt.title('x-Displacement vs Time - Leap Frog - Parabolic')
"""
plt.figure()
plt.plot(ts,y, label = 'y-Displacement')
plt.plot(ts,x, label = 'x-Displacement')
plt.title('Displacement vs Time - Leap Frog - Hyperbolic')
plt.xlabel('Time')
plt.ylabel('Displacement')
plt.legend()

plt.figure()
plt.plot(ts,Vy, label = 'y-Velocity')
plt.plot(ts,Vx, label = 'x-Velocity')
plt.title('Velocity vs Time - Leap Frog - Hyperbolic')
plt.xlabel('Time')
plt.ylabel('Velocity')
plt.legend()


plt.figure()
plt.plot(ts,V, label = 'Potential Energy')
plt.plot(ts,T, label = 'Kinetic Energy')
plt.plot(ts,E, label = 'Total Energy')
plt.xlabel('Time')
plt.ylabel('Energy')
plt.title('Energy vs Time - Leap Frog - Hyperbolic')
plt.legend()
plt.show()
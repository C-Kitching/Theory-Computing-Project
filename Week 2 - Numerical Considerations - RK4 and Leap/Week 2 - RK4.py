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

r=np.array([10,0,0,0.3])  # initial x,y displacements; initial x,y velocities
Ro=r  # save initial conditions
t=0  # time start
tmax=170000  # time end
h=6  # step size
ts=np.array([t])  # store all times
rs=np.array([r])  # store all solution points

for i in range(int(tmax/h)):
    (t,r)=rk4(r,deriv,t,h)
    ts=np.append(ts,t)
    rs=np.concatenate((rs,np.array([r])))
        
[x,y,Vx,Vy]=rs.transpose()

"""
plt.figure()
plt.plot(x,y)
plt.axis('equal')
plt.xlabel('x')
plt.ylabel('y')
plt.title('Displacement')

plt.figure()
plt.plot(Vx,Vy)
plt.axis('equal')
plt.xlabel('y')
plt.ylabel('x')
plt.title('Velocity')

plt.figure()
plt.plot(x,Vx, label = 'x')
plt.plot(y,Vy, label = 'y')
plt.axis('equal')
plt.xlabel('Displacement')
plt.ylabel('Velocity')
plt.legend()
plt.title('Velocity - Displacement')

plt.figure()
plt.plot(ts,x, label='x')
plt.plot(ts,y, label = 'y')
plt.legend()
plt.xlabel('Time / s')
plt.ylabel('Displacement')
plt.title('Displacment-Time')

plt.figure()
plt.plot(ts,Vx, label='x')
plt.plot(ts,Vx, label='y')
plt.legend()
plt.xlabel('Time / s')
plt.ylabel('Velocity')
plt.title('Velocity-Time')
"""

V=-1/((x**2+y**2)**(1/2))
T=1/2*Vx**2+1/2*Vy**2
E=T+V

plt.figure()
plt.plot(ts,V, label = 'Potential Energy')
plt.plot(ts,T, label = 'Kinetic Energy')
plt.plot(ts,E, label = 'Total Energy')
plt.xlabel('Time')
plt.ylabel('Energy')
plt.title('Energy vs Time - RK4')
plt.legend()
plt.show()

"""
L=((Ro[0]**2+Ro[1]**2)**(1/2))*((Ro[2]**2+Ro[3]**2)**(1/2))
initialE=(1/2*Ro[2]**2+1/2*Ro[3]**2)-1/((Ro[0]**2+Ro[1]**2)**(1/2))
print('inital Energy =')
print(initialE)
eccentricity = (1+2*(L**2)*initialE)**(1/2)
print('Eccentricity =')
print(eccentricity)
maxPhi=np.arctan2(y[-1],x[-1])
phi=np.arange(-np.pi,np.pi,0.01)
#phi=np.arange(0,maxPhi,0.01)
R=L**2/(1+eccentricity*np.cos(phi))
xAna=R*np.cos(phi)
yAna=R*np.sin(phi)

plt.figure()
plt.xlabel('x')
plt.ylabel('y')
plt.title('Analytic vs RK4')
plt.axis('equal')
plt.plot(xAna,yAna, label = 'Analytic')
plt.plot(x,y, label = 'RK4')
plt.legend()


initialE=(1/2*Ro[2]**2+1/2*Ro[3]**2)-1/((Ro[0]**2+Ro[1]**2)**(1/2))
avE=[]
step=[]
#fig=plt.figure()

#ax = fig.add_axes([0.1, 0.1, 0.6, 0.75])
for i in range(1,2000,50):
    r=Ro  # initial x,y displacements; initial x,y velocities
    t=0  # time start
    tmax=300  # time end
    h=i/100  # step size
    ts=np.array([t])  # store all times
    rs=np.array([r])  # store all solution points

    for j in range(int(tmax/h)):
        (t,r)=rk4(r,deriv,t,h)
        ts=np.append(ts,t)
        rs=np.concatenate((rs,np.array([r])))
        
    [x,y,Vx,Vy]=rs.transpose()
    E=(1/2*Vx**2+1/2*Vy**2)-1/((x**2+y**2)**(1/2))
    avE=np.append(avE,np.mean(E))
    step=np.append(step,h)
    deltaE=initialE-E
    #ax.plot(ts,E, label = i/100)
    
#plt.xlabel('Time')
#plt.ylabel('Total Energy')
#plt.title('Timestep') 
#ax.legend(bbox_to_anchor=(1.04,1), loc="upper left", title = 'h')

plt.figure()
plt.plot(step,avE, label='RK4')
plt.plot([np.amin(step),np.amax(step)],[initialE,initialE], label = 'True value')
plt.plot([np.amin(step),np.amax(step)],[initialE*1.001,initialE*1.001], color='green', dashes=[6, 2], label = '0.1% error')
plt.plot([np.amin(step),np.amax(step)],[initialE*0.999,initialE*0.999], color='green', dashes=[6, 2], label = '0.1% error')

# Remove duplicate legends
handles, labels = plt.gca().get_legend_handles_labels()
newLabels, newHandles = [], []
for handle, label in zip(handles, labels):
  if label not in newLabels:
    newLabels.append(label)
    newHandles.append(handle)
plt.legend(newHandles, newLabels)

plt.gca().invert_xaxis()
plt.tight_layout()
plt.xlabel('Step size')
plt.ylabel('Energy Difference')
plt.title('Energy convergence with timestep - RK4')
"""







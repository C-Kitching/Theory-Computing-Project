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
v=np.array([0,0.3])  #initial x,y velocities, E<0 -> ellipse, E=0 -> parabola, E>0 -> hyperbola
Ro=np.array([r[0],r[1],v[0],v[1]])
t=0  # time start
tmax=14000  # time end
h=4  # step size
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

plt.figure()
plt.plot(x,y)
plt.axis('equal')
plt.xlabel('x')
plt.ylabel('y')
plt.title('Displacement')
"""
plt.figure()
plt.plot(Vx,Vy)
plt.axis('equal')
plt.xlabel('Vy')
plt.ylabel('Vx')
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
plt.title('Energy vs Time - Leap Frog')
plt.legend()
plt.show()

L=((Ro[0]**2+Ro[1]**2)**(1/2))*((Ro[2]**2+Ro[3]**2)**(1/2))
initialE=(1/2*Ro[2]**2+1/2*Ro[3]**2)-1/((Ro[0]**2+Ro[1]**2)**(1/2))
print('inital Energy =')
print(initialE)
eccentricity = (1+2*(L**2)*initialE)**(1/2)
print('Eccentricity =')
print(eccentricity)
#maxPhi=np.arctan2(y[-1],x[-1])
phi=np.arange(0,2*np.pi,0.01)
#phi=np.arange(0,maxPhi,0.01)
R=L**2/(1+eccentricity*np.cos(phi))
xAna=R*np.cos(phi)
yAna=R*np.sin(phi)

plt.figure()
plt.xlabel('x')
plt.ylabel('y')
plt.title('Analytic vs Leap')
plt.axis('equal')
plt.plot(xAna,yAna, label = 'Analytic')
plt.plot(x,y, label = 'RK4')
plt.legend()

initialE=(1/2*Ro[2]**2+1/2*Ro[3]**2)-1/((Ro[0]**2+Ro[1]**2)**(1/2))
avE=[]
step=[]

#fig=plt.figure()
#ax = fig.add_axes([0.1, 0.1, 0.6, 0.75])
for i in range(1,1000,10):
    r=np.array([Ro[0],Ro[1]])  # initial x,y displacements
    v=np.array([Ro[2],Ro[3]])  #initial x,y velocities
    t=0  # time start
    tmax=300  # time end
    h=i/100  # step size
    ts=np.array([t])  # store all times
    rs=np.array([r])  # store all solution points
    vs=np.array([v])  # store all solution points

    for j in range(int(tmax/h)):
        (t,r,v)=LeapFrog(r,v,h,t)
        ts=np.append(ts,t)
        rs=np.concatenate((rs,np.array([r])))
        vs=np.concatenate((vs,np.array([v])))
        
    [x,y]=rs.transpose()
    [Vx,Vy]=vs.transpose() 
    E=(1/2*Vx**2+1/2*Vy**2)-1/((x**2+y**2)**(1/2))
    avE=np.append(avE,np.mean(E))
    step=np.append(step,h)
    deltaE=initialE-E
    #ax.plot(ts,E, label = i/100)
    
#plt.xlabel('Time')
#plt.ylabel('Total Energy')
#plt.title('Timestep') 
#ax.legend(bbox_to_anchor=(1.04,1), loc="upper left", title = 'h')

fig = plt.figure()
plt.plot(step,avE, label='Leap Frog')
plt.plot([np.amin(step),np.amax(step)],[initialE,initialE], label = 'True value')
plt.plot([np.amin(step),np.amax(step)],[initialE*1.001,initialE*1.001], color='green', dashes=[6, 2], label = '0.1% error')
plt.plot([np.amin(step),np.amax(step)],[initialE*0.999,initialE*0.999], color='green', dashes=[6, 2], label = '0.1% error')
#plt.legend(bbox_to_anchor=(1.04,1), loc="upper left", title = 'h')

# Remove duplicate legends
handles, labels = plt.gca().get_legend_handles_labels()
newLabels, newHandles = [], []
for handle, label in zip(handles, labels):
  if label not in newLabels:
    newLabels.append(label)
    newHandles.append(handle)
plt.legend(newHandles, newLabels)

plt.tight_layout()
plt.gca().invert_xaxis()
plt.xlabel('Step size')
plt.ylabel('Energy Difference')
plt.title('Energy convergence with timestep - Leap Frog')



#Imports
import numpy as np
import matplotlib.pyplot as plt

#Function to calculate derivative
def derivLeap(r):
    return np.array([-r[0]/((r[0]**2+r[1]**2)**(3/2)),-r[1]/((r[0]**2+r[1]**2)**(3/2))])

def LeapFrog(r,v,dt,t):
    v=v+ 1/2 *dt*derivLeap(r)
    r=r+dt*v
    v=v +1/2*dt*derivLeap(r)
    t=t+dt
    return (t,r,v)

#Function to calculate derivative
def derivRK4(r,t):
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

r=np.array([10,0])  # initial x,y displacements
v=np.array([0,0.3])  #initial x,y velocities, E<0 -> ellipse, E=0 -> parabola, E>0 -> hyperbola
Ro=np.array([r[0],r[1],v[0],v[1]])
t=0  # time start
tmax=13000  # time end
h=10 # step size
ts=np.array([t])  # store all times
rs=np.array([r])  # store all solution points
vs=np.array([v])  # store all solution points

counterLeap=0
periodLeap=[]

for i in range(int(tmax/h)):
    (t,r,v)=LeapFrog(r,v,h,t)
    ts=np.append(ts,t)
    rs=np.concatenate((rs,np.array([r])))
    vs=np.concatenate((vs,np.array([v])))
    
    if i>1:
        if r[1]/rs[-2][1]<0:
            if rs[-1][0]>0:
                counterLeap += 1
                periodLeap=np.append(periodLeap,t)

period=periodLeap[2]-periodLeap[1]
[xLeap,yLeap]=rs.transpose()
[VxLeap,VyLeap]=vs.transpose() 


r=np.array([10,0,0,0.3])  # initial x,y displacements; initial x,y velocities
Ro=r  # save initial conditions
t=0  # time start
tmax=13000  # time end
h=10  # step size
ts=np.array([t])  # store all times
rs=np.array([r])  # store all solution points

for i in range(int(tmax/h)):
    (t,r)=rk4(r,derivRK4,t,h)
    ts=np.append(ts,t)
    rs=np.concatenate((rs,np.array([r])))
        
[xRK4,yRK4,VxRK4,VyRK4]=rs.transpose()


Vrk4=-1/((xRK4**2+yRK4**2)**(1/2))
Trk4=1/2*VxRK4**2+1/2*VyRK4**2
Erk4=Trk4+Vrk4

VLeap=-1/((xLeap**2+yLeap**2)**(1/2))
TLeap=1/2*VxLeap**2+1/2*VyLeap**2
ELeap=TLeap+VLeap



fig, (ax1, ax2) = plt.subplots(1,2,sharey=True)

ax1.plot(ts/period, VLeap)
ax1.plot(ts/period, TLeap)
ax1.plot(ts/period, ELeap)
ax1.set_title('Leapfrog', fontsize=13)

ax2.plot(ts/period, Vrk4, label='Potential energy')
ax2.plot(ts/period, Trk4, label='Kinetic energy')
ax2.plot(ts/period, Erk4, label='Total energy')
ax2.set_title('RK4', fontsize=13)

fig.text(0.5, 0.02, 'Time / orbital periods', ha='center', fontsize=12)
fig.text(0.04, 0.5, 'Energy', va='center', rotation='vertical', fontsize=12)
fig.suptitle('Energy vs time', fontsize=14)
fig.legend(bbox_to_anchor=(1,0.7))
plt.tight_layout


plt.show()













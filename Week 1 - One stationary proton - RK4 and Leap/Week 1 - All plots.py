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


#CIRCULAR
r=np.array([10,0])  # initial x,y displacements
v=np.array([0,0.32])  #initial x,y velocities
Ro=np.array([r[0],r[1],v[0],v[1]])
t=0  # time start
tmax=700  # time end
h=1  # step size
ts=np.array([t])  # store all times
rs=np.array([r])  # store all solution points
vs=np.array([v])  # store all solution points

for i in range(int(tmax/h)):
    (t,r,v)=LeapFrog(r,v,h,t)
    ts=np.append(ts,t)
    rs=np.concatenate((rs,np.array([r])))
    vs=np.concatenate((vs,np.array([v])))

[xCir,yCir]=rs.transpose()
[VxCir,VyCir]=vs.transpose() 


#ELLIPTICAL
r=np.array([10,0])  # initial x,y displacements
v=np.array([0,0.1])  #initial x,y velocities
Ro=np.array([r[0],r[1],v[0],v[1]])
t=0  # time start
tmax=300  # time end
h=0.01  # step size
ts=np.array([t])  # store all times
rs=np.array([r])  # store all solution points
vs=np.array([v])  # store all solution points

for i in range(int(tmax/h)):
    (t,r,v)=LeapFrog(r,v,h,t)
    ts=np.append(ts,t)
    rs=np.concatenate((rs,np.array([r])))
    vs=np.concatenate((vs,np.array([v])))

[xEll,yEll]=rs.transpose()
[VxEll,VyEll]=vs.transpose() 


#PARABOLIC
r=np.array([10,0])  # initial x,y displacements
v=np.array([0,0.44])  #initial x,y velocities
Ro=np.array([r[0],r[1],v[0],v[1]])
t=0  # time start
tmax=1000  # time end
h=1  # step size
ts=np.array([t])  # store all times
rs=np.array([r])  # store all solution points
vs=np.array([v])  # store all solution points

for i in range(int(tmax/h)):
    (t,r,v)=LeapFrog(r,v,h,t)
    ts=np.append(ts,t)
    rs=np.concatenate((rs,np.array([r])))
    vs=np.concatenate((vs,np.array([v])))

[xPara,yPara]=rs.transpose()
[VxPara,VyPara]=vs.transpose() 


#HYPERBOLIC
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

[xHyp,yHyp]=rs.transpose()
[VxHyp,VyHyp]=vs.transpose() 


fig=plt.figure()

ax = plt.subplot(2, 2, 1)
ax.plot(xCir, yCir)
ax.plot(0,0, 'ro')
ax.plot(10,0,'bo')
ax.axis('equal')
ax.set_title('Circular', fontsize=10)
ax.set_xlabel('x-Displacement', fontsize=8)
ax.set_ylabel('y-Displacement', fontsize=8)


ax2 = plt.subplot(2, 2, 2)
ax2.plot(xEll, yEll)
ax2.plot(0,0, 'ro')
ax2.plot(10,0,'bo')
ax2.axis('equal')
ax2.set_title('Elliptical', fontsize=10)
ax2.set_xlabel('x-Displacement', fontsize=8)
ax2.set_ylabel('y-Displacement', fontsize=8)

ax3 = plt.subplot(2, 2, 3)
ax3.plot(xPara, yPara)
ax3.plot(0,0, 'ro')
ax3.plot(10,0,'bo')
ax3.axis('equal')
ax3.set_title('Parabolic', fontsize=10)
ax3.set_xlabel('x-Displacement', fontsize=8)
ax3.set_ylabel('y-Displacement', fontsize=8)

ax4 = plt.subplot(2, 2, 4)
ax4.plot(xHyp, yHyp, label = 'Electron path')
ax4.plot(10,0,'bo', label = 'Initial electron position')
ax4.plot(0,0, 'ro', label = 'Proton position')
ax4.axis('equal')
ax4.set_title('Hyperbolic', fontsize=10)
ax4.set_xlabel('x-Displacement', fontsize=8)
ax4.set_ylabel('y-Displacement', fontsize=8)

fig.suptitle('x-y particle paths', fontsize=14)
fig.legend(bbox_to_anchor=(0.9,0.1), ncol=3)
fig.tight_layout()
fig.subplots_adjust(bottom=0.2, wspace=0.4)
fig.subplots_adjust(top=0.88)

plt.show()






















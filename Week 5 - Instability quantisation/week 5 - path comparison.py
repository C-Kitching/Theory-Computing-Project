#Imports
import numpy as np
import matplotlib.pyplot as plt

def derivStat(r,a,b):
    return np.array([-((r[0]-a[0])/(((r[0]-a[0])**2+(r[1]-a[1])**2)**(3/2))+(r[0]-b[0])/(((r[0]-b[0])**2+(r[1]-b[1])**2)**(3/2))),-((r[1]-a[1])/(((r[0]-a[0])**2+(r[1]-a[1])**2)**(3/2))+(r[1]-b[1])/(((r[0]-b[0])**2+(r[1]-b[1])**2)**(3/2)))])


def derivMov(r,ratio):
    return np.array([(-((r[0]-r[2])/(((r[0]-r[2])**2+(r[1]-r[3])**2)**(3/2))+(r[0]-r[4])/(((r[0]-r[4])**2+(r[1]-r[5])**2)**(3/2)))),\
                     (-((r[1]-r[3])/(((r[0]-r[2])**2+(r[1]-r[3])**2)**(3/2))+(r[1]-r[5])/(((r[0]-r[4])**2+(r[1]-r[5])**2)**(3/2)))),\
                     1/ratio*(-((r[2]-r[0])/(((r[2]-r[0])**2+(r[3]-r[1])**2)**(3/2))-(r[2]-r[4])/(((r[2]-r[4])**2+(r[3]-r[5])**2)**(3/2)))),\
                     1/ratio*(-((r[3]-r[1])/(((r[2]-r[0])**2+(r[3]-r[1])**2)**(3/2))-(r[3]-r[5])/(((r[2]-r[4])**2+(r[3]-r[5])**2)**(3/2)))),\
                     1/ratio*(-((r[4]-r[0])/(((r[4]-r[0])**2+(r[5]-r[1])**2)**(3/2))-(r[4]-r[2])/(((r[4]-r[2])**2+(r[5]-r[3])**2)**(3/2)))),\
                     1/ratio*(-((r[5]-r[1])/(((r[4]-r[0])**2+(r[5]-r[1])**2)**(3/2))-(r[5]-r[3])/(((r[4]-r[2])**2+(r[5]-r[3])**2)**(3/2))))])  

def LeapFrogStat(r,v,dt,t,a,b):
    v=v+ 1/2 *dt*derivStat(r,a,b)
    r=r+dt*v
    v=v +1/2*dt*derivStat(r,a,b)
    t=t+dt
    return (t,r,v)    
    
    
def LeapFrogMov(r,v,dt,t,ratio):
    v=v+ 1/2 *dt*derivMov(r,ratio)
    r=r+dt*v
    v=v +1/2*dt*derivMov(r,ratio)
    t=t+dt
    return (t,r,v)

elM=9.109e-31
proM=1.673e-27
e=1.602e-19
ratio=proM/elM

rEl=np.array([20,0])  # initial electron x,y displacements
vEl=np.array([0,0.3])  #initial electron x,y velocities, E<0 -> ellipse, E=0 -> parabola, E>0 -> hyperbola
Ro=np.array([rEl[0],rEl[1],vEl[0],vEl[1]])  # save initial electron values

a=np.array([5,0])  # 1st proton x,y displacements
aV=np.array([0,0])  # 2nd proton x,y velocities
Ao=np.array([a[0],a[1],aV[0],aV[1]])  # save initial 1st proton values
b=np.array([-5,0])  # 2nd proton x,y displacements
bV=np.array([0,0])  # 2nd proton x,y velocities
Bo=np.array([b[0],b[1],bV[0],bV[1]])  # save initial 2nd proton values

r=np.array([rEl[0],rEl[1],a[0],a[1],b[0],b[1]])
v=np.array([vEl[0],vEl[1],aV[0],aV[1],bV[0],bV[1]])

t=0  # time start
tmax=1696  # time end
h=0.1  # step size
ts=np.array([t])  # store all times
rs=np.array([r])  # store all solution points
vs=np.array([v])  # store all solution points

for i in range(int(tmax/h)):
    (t,r,v)=LeapFrogMov(r,v,h,t,ratio)
    ts=np.append(ts,t)
    rs=np.concatenate((rs,np.array([r])))
    vs=np.concatenate((vs,np.array([v])))

[xEl,yEl,xA,yA,xB,yB]=rs.transpose()
[VxEl,VyEl,VxA,VyA,VxB,VyB]=vs.transpose() 



r=np.array([20,0])  # initial electron x,y displacements (20,0)
v=np.array([0,0.3])  #initial electron x,y velocities, E<0 -> ellipse, E=0 -> parabola, E>0 -> hyperbola (0,0.15)
Ro=np.array([r[0],r[1],v[0],v[1]])  # save initial electron values

a=np.array([5,0])  # 1st proton x,y displacements
Ao=a  # save initial 1st proton values
b=np.array([-5,0])  # 2nd proton x,y displacements
Bo=b  # save initial 2ndst proton values

t=0  # time start
tmax=10000  # time end
h=0.1  # step size
ts=np.array([t])  # store all times
rs=np.array([r])  # store all solution points
vs=np.array([v])  # store all solution points

for i in range(int(tmax/h)):
    (t,r,v)=LeapFrogStat(r,v,h,t,a,b)
    ts=np.append(ts,t)
    rs=np.concatenate((rs,np.array([r])))
    vs=np.concatenate((vs,np.array([v])))

[x,y]=rs.transpose()
[Vx,Vy]=vs.transpose()




fig, (ax1, ax2) = plt.subplots(1,2)

ax1.plot(xEl,yEl, color='blue', label = 'Electron Path')
ax1.plot(Ro[0],Ro[1], 'co', label = 'Initial Electron Position')
ax1.plot(xA,yA, color='red', label = 'Proton 1 Path')
ax1.plot(Ao[0],Ao[1], 'ro', label = 'Initial Proton 1 Position')
ax1.plot(xB,yB, color='green', label = 'Proton 2 Path')
ax1.plot(Bo[0],Bo[1], 'go', label = 'Initial Proton 2 Position')
ax1.axis('equal')
ax1.set_title('3-Body Problem')

ax2.plot(x,y, color='blue')
ax2.plot(Ro[0],Ro[1], 'co')
ax2.plot(Ao[0],Ao[1], 'ro')
ax2.plot(Bo[0],Bo[1], 'go')
ax2.axis('equal')
ax2.set_title('Limiting Case of 3-Body Problem')

fig.text(0.5, 0.2, 'x-displacement', ha='center', fontsize=12)
fig.text(0.5, 0.95, 'x-y Path', ha='center', fontsize=14)
fig.text(0.02, 0.5, 'y-displacement', va='center', rotation='vertical', fontsize=12)
#fig.suptitle('Optimal Step Size', fontsize=14)
fig.legend(bbox_to_anchor=(0.9,0.15), ncol=3)
fig.tight_layout()
fig.subplots_adjust(bottom=0.3, wspace=0.4)
fig.subplots_adjust(top=0.88)


plt.show()





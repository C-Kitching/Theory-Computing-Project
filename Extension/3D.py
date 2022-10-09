#Imports
import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D

#Function to calculate derivative
def deriv(r,ratio):
    return np.array([(-((r[0]-r[3])/(((r[0]-r[3])**2+(r[1]-r[4])**2+(r[2]-r[5])**2)**(3/2))+(r[0]-r[6])/(((r[0]-r[6])**2+(r[1]-r[7])**2+(r[2]-r[8])**2)**(3/2)))),\
                     (-((r[1]-r[4])/(((r[0]-r[3])**2+(r[1]-r[4])**2+(r[2]-r[5])**2)**(3/2))+(r[1]-r[7])/(((r[0]-r[6])**2+(r[1]-r[7])**2+(r[2]-r[8])**2)**(3/2)))),\
                     (-((r[2]-r[5])/(((r[0]-r[3])**2+(r[1]-r[4])**2+(r[2]-r[5])**2)**(3/2))+(r[2]-r[8])/(((r[0]-r[6])**2+(r[1]-r[7])**2+(r[2]-r[8])**2)**(3/2)))),\
                     1/ratio*(-((r[3]-r[0])/(((r[0]-r[3])**2+(r[1]-r[4])**2+(r[2]-r[5])**2)**(3/2))-(r[3]-r[6])/(((r[3]-r[6])**2+(r[4]-r[7])**2+(r[5]-r[8])**2)**(3/2)))),\
                     1/ratio*(-((r[4]-r[1])/(((r[0]-r[3])**2+(r[1]-r[4])**2+(r[2]-r[5])**2)**(3/2))-(r[4]-r[7])/(((r[3]-r[6])**2+(r[4]-r[7])**2+(r[5]-r[8])**2)**(3/2)))),\
                     1/ratio*(-((r[5]-r[2])/(((r[0]-r[3])**2+(r[1]-r[4])**2+(r[2]-r[5])**2)**(3/2))-(r[5]-r[8])/(((r[3]-r[6])**2+(r[4]-r[7])**2+(r[5]-r[8])**2)**(3/2)))),\
                     1/ratio*(-((r[6]-r[0])/(((r[0]-r[6])**2+(r[1]-r[7])**2+(r[2]-r[8])**2)**(3/2))-(r[6]-r[3])/(((r[3]-r[6])**2+(r[4]-r[7])**2+(r[5]-r[8])**2)**(3/2)))),\
                     1/ratio*(-((r[7]-r[1])/(((r[0]-r[6])**2+(r[1]-r[7])**2+(r[2]-r[8])**2)**(3/2))-(r[7]-r[4])/(((r[3]-r[6])**2+(r[4]-r[7])**2+(r[5]-r[8])**2)**(3/2)))),\
                     1/ratio*(-((r[8]-r[2])/(((r[0]-r[6])**2+(r[1]-r[7])**2+(r[2]-r[8])**2)**(3/2))-(r[8]-r[5])/(((r[3]-r[6])**2+(r[4]-r[7])**2+(r[5]-r[8])**2)**(3/2))))])  

def LeapFrog(r,v,dt,t,ratio):
    v=v+ 1/2 *dt*deriv(r,ratio)
    r=r+dt*v
    v=v +1/2*dt*deriv(r,ratio)
    t=t+dt
    return (t,r,v)

elM=9.109e-31
proM=1.673e-27
e=1.602e-19
ratio=proM/elM

k=(1.602e-19)**2/(4*np.pi*8.854e-12)

rEl=np.array([20,0,0])  # initial electron x,y displacements
vEl=np.array([0,0,0.3])  #initial electron x,y velocities, E<0 -> ellipse, E=0 -> parabola, E>0 -> hyperbola
Ro=np.array([rEl[0],rEl[1],rEl[2],vEl[0],vEl[1],vEl[2]])  # save initial electron values

a=np.array([5,0,0])  # 1st proton x,y displacements
aV=np.array([0,0.01,0])  # 2nd proton x,y velocities
Ao=np.array([a[0],a[1],a[2],aV[0],aV[1],aV[2]])  # save initial 1st proton values
b=np.array([-5,0,0])  # 2nd proton x,y displacements
bV=np.array([0,-0.01,0])  # 2nd proton x,y velocities
Bo=np.array([b[0],b[1],b[2],bV[0],bV[1],bV[2]])  # save initial 2nd proton values

r=np.array([rEl[0],rEl[1],rEl[2],a[0],a[1],a[2],b[0],b[1],b[2]])
v=np.array([vEl[0],vEl[1],vEl[2],aV[0],aV[1],aV[1],bV[0],bV[1],bV[2]])

t=0  # time start
tmax=20000  # time end
h=1  # step size
ts=np.array([t])  # store all times
rs=np.array([r])  # store all solution points
vs=np.array([v])  # store all solution points

Counter=0
Period=[]

for i in range(int(tmax/h)):
    (t,r,v)=LeapFrog(r,v,h,t,ratio)
    ts=np.append(ts,t)
    rs=np.concatenate((rs,np.array([r])))
    vs=np.concatenate((vs,np.array([v])))
    
"""    
    if i>1:
        if r[1]/rs[-2][1]<0:
            if rs[-1][0]>a[0]:
                Counter += 1
                Period=np.append(Period,t)
            
period=Period[2]-Period[1]
"""

[xEl,yEl,zEl,xA,yA,zA,xB,yB,zB]=rs.transpose()
[VxEl,VyEl,VzEl,VxA,VyA,VzA,VxB,VyB,VzB]=vs.transpose() 

fig=plt.figure(figsize=(10,10))  # Plotting a figure with size (10,10)
ax=fig.gca(projection='3d')  # Define a new environment to use a 3D axis

ax.plot(xEl,yEl,zEl, label='Electron path')
ax.scatter(Ro[0],Ro[1],Ro[2], s=20, label='Initial electron position')
ax.plot(xA,yA,zA, label='Proton 1 path')
ax.scatter(Ao[0],Ao[1],Ao[2], s=20, label='Initial proton 1 position')
ax.plot(xB,yB,zB, label='Proton 2 path')
ax.scatter(Bo[0],Bo[1],Bo[2], s=20, label='Initial proton 2 position')
ax.legend()

ax.set_xlabel('x-displacement')
ax.set_ylabel('y-displacement')
ax.set_zlabel('z-displacement')
plt.title('3-D particle paths')







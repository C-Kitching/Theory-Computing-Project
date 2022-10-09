#Imports
import numpy as np
import matplotlib.pyplot as plt

#Function to calculate derivative
def deriv(r,ratio):
    return np.array([(-((r[0]-r[2])/(((r[0]-r[2])**2+(r[1]-r[3])**2)**(3/2))+(r[0]-r[4])/(((r[0]-r[4])**2+(r[1]-r[5])**2)**(3/2)))),\
                     (-((r[1]-r[3])/(((r[0]-r[2])**2+(r[1]-r[3])**2)**(3/2))+(r[1]-r[5])/(((r[0]-r[4])**2+(r[1]-r[5])**2)**(3/2)))),\
                     1/ratio*(-((r[2]-r[0])/(((r[2]-r[0])**2+(r[3]-r[1])**2)**(3/2))-(r[2]-r[4])/(((r[2]-r[4])**2+(r[3]-r[5])**2)**(3/2)))),\
                     1/ratio*(-((r[3]-r[1])/(((r[2]-r[0])**2+(r[3]-r[1])**2)**(3/2))-(r[3]-r[5])/(((r[2]-r[4])**2+(r[3]-r[5])**2)**(3/2)))),\
                     1/ratio*(-((r[4]-r[0])/(((r[4]-r[0])**2+(r[5]-r[1])**2)**(3/2))-(r[4]-r[2])/(((r[4]-r[2])**2+(r[5]-r[3])**2)**(3/2)))),\
                     1/ratio*(-((r[5]-r[1])/(((r[4]-r[0])**2+(r[5]-r[1])**2)**(3/2))-(r[5]-r[3])/(((r[4]-r[2])**2+(r[5]-r[3])**2)**(3/2))))])  

def LeapFrog(r,v,dt,t,ratio):
    v=v+ 1/2*dt*deriv(r,ratio)
    r=r+dt*v
    v=v +1/2*dt*deriv(r,ratio)
    t=t+dt
    return (t,r,v)

def LeapFrogVar(r,v,dt,t,ratio):
    v=v+ 1/2*(dt/(((deriv(r,ratio)[0])**2+(deriv(r,ratio)[1])**2)**(1/2))**(1/2))*deriv(r,ratio)
    r=r+(dt/(((deriv(r,ratio)[0])**2+(deriv(r,ratio)[1])**2)**(1/2))**(1/2))*v
    v=v +1/2*(dt/(((deriv(r,ratio)[0])**2+(deriv(r,ratio)[1])**2)**(1/2))**(1/2))*deriv(r,ratio)
    t=t+dt
    return (t,r,v)

elM=9.109e-31
proM=1.673e-27
e=1.602e-19
ratio=proM/elM

k=(1.602e-19)**2/(4*np.pi*8.854e-12)

rEl=np.array([20,0])  # initial electron x,y displacements
vEl=np.array([0,0.15])  #initial electron x,y velocities, E<0 -> ellipse, E=0 -> parabola, E>0 -> hyperbola
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
tmax=10000  # time end
h=0.1  # step size
ts=np.array([t])  # store all times
rs=np.array([r])  # store all solution points
vs=np.array([v])  # store all solution points

for i in range(int(tmax/h)):
    
    if ((r[0]-r[2])**2+(r[1]-r[3])**2)**(1/2) < 1 or ((r[0]-r[4])**2+(r[1]-r[5])**2)**(1/2) < 1 :
        h=0.001
    else:
        h=0.1
        
    (t,r,v)=LeapFrog(r,v,h,t,ratio)
    ts=np.append(ts,t)
    rs=np.concatenate((rs,np.array([r])))
    vs=np.concatenate((vs,np.array([v])))
    
    """
    if i==0:
        (t,r,v)=LeapFrog(r,v,h,t,ratio)
        ts=np.append(ts,t)
        rs=np.concatenate((rs,np.array([r])))
        vs=np.concatenate((vs,np.array([v])))
    else:
        (t,r,v)=LeapFrogVar(r,v,h,t,ratio)
        ts=np.append(ts,t)
        rs=np.concatenate((rs,np.array([r])))
        vs=np.concatenate((vs,np.array([v])))
    """
    
    
[xEl,yEl,xA,yA,xB,yB]=rs.transpose()
[VxEl,VyEl,VxA,VyA,VxB,VyB]=vs.transpose() 

plt.figure()
plt.plot(xEl,yEl,color='dodgerblue', label = 'Electron path')
plt.plot(Ro[0],Ro[1], 'bo', label= 'Initial electron position')
plt.plot(xA,yA, color='red', label = 'Proton 1 path')
plt.plot(Ao[0],Ao[1], 'ro', label= 'Initial proton 1 position')
plt.plot(xB,yB, color='green', label = 'Proton 2 path')
plt.plot(Bo[0],Bo[1], 'go', label= 'Initial proton 2 position')
plt.axis('equal')
plt.xlabel('x-Displacement')
plt.ylabel('y-Displacement')
plt.title('x-y particle paths')

plt.legend(bbox_to_anchor=(1.1,-0.2), ncol=3)
plt.subplots_adjust(bottom=0.3, wspace=0.4)
plt.subplots_adjust(top=0.88)









V=-1/(((xEl-xA)**2+(yEl-yA)**2)**(1/2))-1/(((xEl-xB)**2+(yEl-yB)**2)**(1/2))+1/(((xB-xA)**2+(yB-yA)**2)**(1/2))
T=1/2*(VxEl**2+VyEl**2)+1/2*ratio*(VxA**2+VyA**2)+1/2*ratio*(VxB**2+VyB**2)
E=T+V

plt.figure()
plt.plot(ts,V, label = 'Potential Energy')
plt.plot(ts,T, label = 'Kinetic Energy')
plt.plot(ts,E, label = 'Total Energy')
plt.xlabel('Time')
plt.ylabel('Energy')
plt.title('Energy vs Time')
plt.legend()
plt.show()
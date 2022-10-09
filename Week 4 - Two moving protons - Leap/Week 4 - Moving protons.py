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
tmax=2000  # time end
h=0.1  # step size
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
    
    if i>1:
        if r[1]/rs[-2][1]<0:
            if rs[-1][0]>a[0]:
                Counter += 1
                Period=np.append(Period,t)
            
period=Period[2]-Period[1]

[xEl,yEl,xA,yA,xB,yB]=rs.transpose()
[VxEl,VyEl,VxA,VyA,VxB,VyB]=vs.transpose() 

plt.figure()
plt.plot(xEl,yEl, color='dodgerblue', label = 'Electron path')
plt.plot(Ro[0],Ro[1], 'bo', label = 'Initial electron position')
plt.plot(xA,yA, color='red', label = 'Proton 1 path')
plt.plot(Ao[0],Ao[1], 'ro', label = 'Initial proton 1 position')
plt.plot(xB,yB, color='green', label = 'Proton 2 path')
plt.plot(Bo[0],Bo[1], 'go', label = 'Initial proton 2 position')
plt.axis('equal')
plt.xlabel('x-Displacement')
plt.ylabel('y-Displacement')
plt.legend(bbox_to_anchor=(1.1,-0.2), ncol=3)
plt.subplots_adjust(bottom=0.3, wspace=0.4)
plt.subplots_adjust(top=0.88)
plt.title('x-y particle paths')

"""
plt.figure()
plt.plot(ts,(xEl**2+yEl**2)**(1/2)-(xA**2+yA**2)**(1/2))
"""

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
plt.plot(ts,xEl, label='Electron x-Displacement')
plt.plot(ts,yEl, label = 'Electron y-Displacement')
plt.plot(ts,xA, label='Proton 1 x-Displacement')
plt.plot(ts,yA, label = 'Proton 1 y-Displacement')
plt.plot(ts,xB, label='Proton 2 x-Displacement')
plt.plot(ts,yB, label = 'Proton 2 y-Displacement')
plt.legend(bbox_to_anchor=(1.04,1), loc="upper left")  # shift legend position
#plt.legend()
plt.tight_layout()
plt.xlabel('Time')
plt.ylabel('Displacement')
plt.title('Displacment vs Time')

plt.figure()
plt.plot(ts,VxEl, label='Electron x-Velocity')
plt.plot(ts,VyEl, label='Electron y-Velocity')
plt.plot(ts,VxA, label='Proton 1 x-Velocity')
plt.plot(ts,VyA, label='Proton 1 y-Velocity')
plt.plot(ts,VxB, label='Proton 2 x-Velocity')
plt.plot(ts,VyB, label='Proton 2 y-Velocity')
plt.legend(bbox_to_anchor=(1.04,1), loc="upper left")  # shift legend position
#plt.legend()
plt.tight_layout()
plt.xlabel('Time')
plt.ylabel('Velocity')
plt.title('Velocity vs Time')

"""
V=-1/(((xEl-xA)**2+(yEl-yA)**2)**(1/2))-1/(((xEl-xB)**2+(yEl-yB)**2)**(1/2))+1/(((xB-xA)**2+(yB-yA)**2)**(1/2))
T=1/2*(VxEl**2+VyEl**2)+1/2*ratio*(VxA**2+VyA**2)+1/2*ratio*(VxB**2+VyB**2)
E=T+V

plt.figure()
plt.plot(ts/period,V, label = 'Potential energy')
plt.plot(ts/period,T, label = 'Kinetic energy')
plt.plot(ts/period,E, label = 'Total energy')
plt.xlabel('Time / Orbital periods')
plt.ylabel('Energy')
plt.title('Energy vs time')
plt.legend()
plt.show()
"""

L=((Ro[0]**2+Ro[1]**2)**(1/2))*((Ro[2]**2+Ro[3]**2)**(1/2))
print('Angular Momentum:')
print(L)
initialE=(1/2*Ro[2]**2+1/2*Ro[3]**2)-1/(((Ro[0]-Ao[0])**2+(Ro[1]-Ao[1])**2)**(1/2))-1/(((Ro[0]-Bo[0])**2+(Ro[1]-Bo[1])**2)**(1/2))
print('Initial Energy:')
print(initialE)
eccentricity = (1+2*(L**2)*initialE)**(1/2)
print('Eccentricity:')
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
plt.title('Analytic vs Leap Frog')
plt.axis('equal')
plt.plot(xAna,yAna, label = 'Analytic')
plt.plot(x,y, label = 'Leap Frog')
plt.legend()
"""
"""
initialE=1/2*(Ro[2]**2+Ro[3]**2)+1/2*ratio*(Ao[2]**2+Ao[3]**2)+1/2*ratio*(Bo[2]**2+Bo[3]**2)\
    +(-1/(((Ro[0]-Ao[0])**2+(Ro[1]-Ao[1])**2)**(1/2))-1/(((Ro[0]-Bo[0])**2+(Ro[1]-Bo[1])**2)**(1/2))+1/(((Bo[0]-Ao[0])**2+(Bo[1]-Ao[1])**2)**(1/2)))
print('Initial Energy:')
print(initialE)
avE=[]
step=[]

#fig=plt.figure()
#ax = fig.add_axes([0.1, 0.1, 0.6, 0.75])
for i in range(10,100,10):
    
    r=np.array([Ro[0],Ro[1],Ao[0],Ao[1],Bo[0],Bo[1]])
    v=np.array([Ro[2],Ro[3],Ao[2],Ao[3],Bo[2],Bo[3]])

    t=0  # time start
    tmax=300  # time end
    h=i/100  # step size
    ts=np.array([t])  # store all times
    rs=np.array([r])  # store all solution points
    vs=np.array([v])  # store all solution points

    for i in range(int(tmax/h)):
        (t,r,v)=LeapFrog(r,v,h,t,ratio)
        ts=np.append(ts,t)
        rs=np.concatenate((rs,np.array([r])))
        vs=np.concatenate((vs,np.array([v])))

    [xEl,yEl,xA,yA,xB,yB]=rs.transpose()
    [VxEl,VyEl,VxA,VyA,VxB,VyB]=vs.transpose() 
        
    E=1/2*(VxEl**2+VyEl**2)+1/2*ratio*(VxA**2+VyA**2)+1/2*ratio*(VxB**2+VyB**2)\
        +(-1/(((xEl-xA)**2+(yEl-yA)**2)**(1/2))-1/(((xEl-xB)**2+(yEl-yB)**2)**(1/2))+1/(((xB-xA)**2+(yB-yA)**2)**(1/2)))
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

#plt.legend(bbox_to_anchor=(1.04,1), loc="upper left", title = 'h')  # shift legend position
#plt.legend()

# Remove duplicate legends
handles, labels = plt.gca().get_legend_handles_labels()
newLabels, newHandles = [], []
for handle, label in zip(handles, labels):
  if label not in newLabels:
    newLabels.append(label)
    newHandles.append(handle)
plt.legend(newHandles, newLabels)

plt.gca().invert_xaxis()  # invert x-axis

plt.xlabel('Step size')
plt.ylabel('Energy Difference')
plt.title('Optimal Step Size')
"""


"""
OneBodyEnergy=1/2*elM*(Ro[2]**2+Ro[3]**2)**(1/2)-2*k/((Ro[0]**2+Ro[1]**2)**(1/2))
TwoBodyEnergy=1/2*elM*(Ro[2]**2+Ro[3]**2)**(1/2)-k/(((Ro[0]-Ao[0])**2+(Ro[1]-Ao[1])**2)**(1/2))-k/(((Ro[0]-Bo[0])**2+(Ro[1]-Bo[1])**2)**(1/2))
print('Initial OneBodyEnergy:')
print(OneBodyEnergy)
print('Initial TwoBodyEnergy:')
print(TwoBodyEnergy)

E1body=[]
E2body=[]
E3body=[]
EnergyDif=[]
radius=[]

#fig=plt.figure()
#ax = fig.add_axes([0.1, 0.1, 0.6, 0.75])
for i in range(50,500,50):
    for j in range(50,500,50):
        r=np.array([i,j])  # initial x,y displacements
        v=np.array([Ro[2],Ro[3]])  #initial x,y velocities
        a=Ao  # set 1st proton position
        b=Bo  # set 2nd proton position
        E3=1/2*elM*(v[0]**2+v[1]**2)+1/2*proM*(a[2]**2+a[3]**2)+1/2*proM*(b[2]**2+b[3]**2)\
            +k*(-1/(((r[0]-a[0])**2+(r[1]-a[1])**2)**(1/2))-1/(((r[0]-b[0])**2+(r[1]-b[1])**2)**(1/2))+1/(((b[0]-a[0])**2+(b[1]-a[1])**2)**(1/2)))
        E2=1/2*elM*(v[0]**2+v[1]**2)**(1/2)-k/(((r[0]-a[0])**2+(r[1]-a[1])**2)**(1/2))-k/(((r[0]-b[0])**2+(r[1]-b[1])**2)**(1/2))
        E1=1/2*elM*(v[0]**2+v[1]**2)**(1/2)-2*k/((r[0]**2+r[1]**2)**(1/2))
        E1body=np.append(E1body,E1)
        E2body=np.append(E2body,E2)
        E3body=np.append(E3body,E3)
        EnergyDif=np.append(EnergyDif,E2-E3)
        radius=np.append(radius,(r[0]**2+r[1]**2)**(1/2))
    
fig = plt.figure()
#EnergyDif=np.sort(EnergyDif)
#radius=np.sort(radius)
plt.plot(radius,EnergyDif, 'ro')
plt.plot([np.amin(radius),np.amax(radius)],[0,0], color='green', dashes=[6, 2], label = 'Zero Energy Difference')
#plt.legend(bbox_to_anchor=(1.04,1), loc="upper left", title = 'h')  # shift legend position
#plt.legend()

# Remove duplicate legends
handles, labels = plt.gca().get_legend_handles_labels()
newLabels, newHandles = [], []
for handle, label in zip(handles, labels):
  if label not in newLabels:
    newLabels.append(label)
    newHandles.append(handle)
plt.legend(newHandles, newLabels)

plt.xlabel('Initial Displacement / m')
plt.ylabel('Energy Difference / J')
plt.title('Two body convergence')

plt.figure()
plt.plot(radius,E1body,'ro', label='One Body')
plt.plot(radius,E2body,'go', label='Two Body')
plt.plot(radius,E3body,'bo', label='Three Body')

plt.legend()
plt.xlabel('Initial Displacement / m')
plt.ylabel('Energy / J')
plt.title('Multiple Body Energies')
"""
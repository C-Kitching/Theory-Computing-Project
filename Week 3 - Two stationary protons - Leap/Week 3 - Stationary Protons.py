#Imports
import numpy as np
import matplotlib.pyplot as plt

#Function to calculate derivative
def deriv(r,a,b):
    return np.array([-((r[0]-a[0])/(((r[0]-a[0])**2+(r[1]-a[1])**2)**(3/2))+(r[0]-b[0])/(((r[0]-b[0])**2+(r[1]-b[1])**2)**(3/2))),-((r[1]-a[1])/(((r[0]-a[0])**2+(r[1]-a[1])**2)**(3/2))+(r[1]-b[1])/(((r[0]-b[0])**2+(r[1]-b[1])**2)**(3/2)))])

def LeapFrog(r,v,dt,t,a,b):
    v=v+ 1/2 *dt*deriv(r,a,b)
    r=r+dt*v
    v=v +1/2*dt*deriv(r,a,b)
    t=t+dt
    return (t,r,v)

r=np.array([20,0])  # initial electron x,y displacements (20,0)
v=np.array([0,0.15])  #initial electron x,y velocities, E<0 -> ellipse, E=0 -> parabola, E>0 -> hyperbola (0,0.15)
Ro=np.array([r[0],r[1],v[0],v[1]])  # save initial electron values

a=np.array([5,0])  # 1st proton x,y displacements
Ao=a  # save initial 1st proton values
b=np.array([-5,0])  # 2nd proton x,y displacements
Bo=b  # save initial 2ndst proton values

t=0  # time start
tmax=20000  # time end
h=1  # step size
ts=np.array([t])  # store all times
rs=np.array([r])  # store all solution points
vs=np.array([v])  # store all solution points

for i in range(int(tmax/h)):
    (t,r,v)=LeapFrog(r,v,h,t,a,b)
    ts=np.append(ts,t)
    rs=np.concatenate((rs,np.array([r])))
    vs=np.concatenate((vs,np.array([v])))

[x,y]=rs.transpose()
[Vx,Vy]=vs.transpose() 

plt.figure()
plt.plot(x,y, color='dodgerblue',label= 'Electron path')
plt.plot(Ro[0],Ro[1], 'bo', label = 'Initial electron position')
plt.plot(Ao[0],Ao[1], 'ro', label='Proton 1 position')
plt.plot(Bo[0],Bo[1], 'go', label = 'Proton 2 position ')
plt.legend(bbox_to_anchor=(0.9,-0.2), ncol=2)
plt.subplots_adjust(bottom=0.3, wspace=0.4)
plt.subplots_adjust(top=0.88)
plt.axis('equal')
plt.xlabel('x-Displacement')
plt.ylabel('y-Displacement')
plt.title('x-y particle paths')

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
plt.plot(ts,x, label='x-Displacement')
plt.plot(ts,y, label = 'y-Displacement')
plt.legend()
plt.xlabel('Time')
plt.ylabel('Displacement')
plt.title('Displacment vs Time')

plt.figure()
plt.plot(ts,Vx, label='x-Velocity')
plt.plot(ts,Vy, label='y-Velocity')
plt.legend()
plt.xlabel('Time')
plt.ylabel('Velocity')
plt.title('Velocity vs Time')

V=-1/(((x-a[0])**2+(y-a[1])**2)**(1/2))-1/(((x-b[0])**2+(y-b[1])**2)**(1/2))
T=1/2*Vx**2+1/2*Vy**2
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
initialE=(1/2*Ro[2]**2+1/2*Ro[3]**2)-1/(((Ro[0]-Ao[0])**2+(Ro[1]-Ao[1])**2)**(1/2))-1/(((Ro[0]-Bo[0])**2+(Ro[1]-Bo[1])**2)**(1/2))
print('Initial Energy:')
print(initialE)
avE=[]
step=[]

#fig=plt.figure()
#ax = fig.add_axes([0.1, 0.1, 0.6, 0.75])
for i in range(10,30,1):
    r=np.array([Ro[0],Ro[1]])  # initial x,y displacements
    v=np.array([Ro[2],Ro[3]])  #initial x,y velocities
    a=Ao  # set 1st proton position
    b=Bo  # set 2nd proton position
    t=0  # time start
    tmax=3000  # time end
    h=i/100  # step size
    ts=np.array([t])  # store all times
    rs=np.array([r])  # store all solution points
    vs=np.array([v])  # store all solution points

    for j in range(int(tmax/h)):
        (t,r,v)=LeapFrog(r,v,h,t,a,b)
        ts=np.append(ts,t)
        rs=np.concatenate((rs,np.array([r])))
        vs=np.concatenate((vs,np.array([v])))
        
    [x,y]=rs.transpose()
    [Vx,Vy]=vs.transpose() 
    E=(1/2*Vx**2+1/2*Vy**2)-1/(((x-a[0])**2+(y-a[1])**2)**(1/2))-1/(((x-b[0])**2+(y-b[1])**2)**(1/2))
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
plt.tight_layout()
plt.xlabel('Step size')
plt.ylabel('Energy Difference')
plt.title('Energy convergence with stepsize')
"""
"""
fig=plt.figure()
ax = fig.add_axes([0.1, 0.1, 0.6, 0.75])
 
Quad=[]
R=[]

for i in range(30,100,10):
    #for j in range(40,70,10):
    
        r=np.array([i,0])  # initial electron x,y displacements
        v=np.array([0,0.15])  #initial electron x,y velocities, E<0 -> ellipse, E=0 -> parabola, E>0 -> hyperbola
        Ro=np.array([r[0],r[1],v[0],v[1]])  # save initial electron values

        a=np.array([5,0])  # 1st proton x,y displacements
        Ao=a  # save initial 1st proton values
        b=np.array([-5,0])  # 2nd proton x,y displacements
        Bo=b  # save initial 2ndst proton values

        t=0  # time start
        tmax=3000  # time end
        h=1  # step size
        ts=np.array([t])  # store all times
        rs=np.array([r])  # store all solution points
        vs=np.array([v])  # store all solution points

        for k in range(int(tmax/h)):
            (t,r,v)=LeapFrog(r,v,h,t,a,b)
            ts=np.append(ts,t)
            rs=np.concatenate((rs,np.array([r])))
            vs=np.concatenate((vs,np.array([v])))

        [x,y]=rs.transpose()
        [Vx,Vy]=vs.transpose() 
        
        Radius = (Ro[0]**2+Ro[1]**2)**(1/2)
        ax.plot(x,y,label= Radius)
        
plt.xlabel('x')
plt.ylabel('y')
plt.title('Displacement')
ax.legend(bbox_to_anchor=(1.04,1), loc="upper left", title = 'Inital Displacement')  # shift legend position

"""
OneBodyEnergy=(1/2*Ro[2]**2+1/2*Ro[3]**2)-2/((Ro[0]**2+Ro[1]**2)**(1/2))
print('OneBodyEnergy:')
print(OneBodyEnergy)

E1body=[]
E2body=[]
EnergyDif=[]
radius=[]

protonSep=((Ao[0]-Bo[0])**2+(Ao[1]-Bo[1])**2)**(1/2)

#fig=plt.figure()
#ax = fig.add_axes([0.1, 0.1, 0.6, 0.75])
for i in range(50,500,50):
    for j in range(50,500,50):
        r=np.array([i,j])  # initial x,y displacements
        v=np.array([Ro[2],Ro[3]])  #initial x,y velocities
        a=Ao  # set 1st proton position
        b=Bo  # set 2nd proton position
        E2=(1/2*v[0]**2+1/2*v[1]**2)-1/(((r[0]-a[0])**2+(r[1]-a[1])**2)**(1/2))-1/(((r[0]-b[0])**2+(r[1]-b[1])**2)**(1/2))
        E1=(1/2*v[0]**2+1/2*v[1]**2)-2/((r[0]**2+r[1]**2)**(1/2))
        E1body=np.append(E1body,E1)
        E2body=np.append(E1body,E2)
        EnergyDif=np.append(EnergyDif,E1-E2)
        radius=np.append(radius,(r[0]**2+r[1]**2)**(1/2))
    
fig = plt.figure()
#EnergyDif=np.sort(EnergyDif)
#radius=np.sort(radius)
plt.plot(radius/protonSep,EnergyDif, 'ro')
plt.plot([np.amin(radius)/protonSep,np.amax(radius)/protonSep],[0,0], color='green', dashes=[6, 2], label = 'Zero energy difference')
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

plt.tight_layout()
plt.xlabel('Ratio between inital electron displacement and proton separation')
plt.ylabel('Total energy difference')
plt.title('Convergence of the two-body problem to a one-body problem')































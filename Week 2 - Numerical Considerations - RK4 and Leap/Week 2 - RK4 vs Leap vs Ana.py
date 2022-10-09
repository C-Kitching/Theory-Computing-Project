#Imports
import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D  # Import 3D axes support module

#Function to calculate derivative
def derivRK4(r,t):
    return np.array([r[2],r[3],-r[0]/((r[0]**2+r[1]**2)**(3/2)),-r[1]/((r[0]**2+r[1]**2)**(3/2))])

def derivLeap(r):
    return np.array([-r[0]/((r[0]**2+r[1]**2)**(3/2)),-r[1]/((r[0]**2+r[1]**2)**(3/2))])

#Function to do single RK4 step
def rk4(r,dr,t,h):
    k1=dr(r,t)
    k2=dr(r+h/2*k1,t+h/2)
    k3=dr(r+h/2*k2,t+h/2)
    k4=dr(r+h*k3,t+h)
    r=r+h*(k1+2*k2+2*k3+k4)/6
    t=t+h
    return (t,r)

def LeapFrog(r,v,dt,t):
    v=v+ 1/2 *dt*derivLeap(r)
    r=r+dt*v
    v=v +1/2 *dt*derivLeap(r)
    t=t+dt
    return (t,r,v)


#Master Conditions
rM=np.array([10,0,0,0.3]) # 10,0,0,0.3
tmax=20000  # time end
h=9  # step size

#Counters
counterRK4=0  # RK4 counter
counterLeap=0  # Leap frog counter


r=np.array([rM[0],rM[1]])  # initial x,y displacements
v=np.array([rM[2],rM[3]])  #initial x,y velocities, E<0 -> ellipse, E=0 -> parabola, E>0 -> hyperbola
Ro=np.array([r[0],r[1],v[0],v[1]])
t=0  # time start
tsRK4=np.array([t])  # store all times
rs=np.array([r])  # store all solution points
vs=np.array([v])  # store all solution points

periodRK4=[]

for i in range(int(tmax/h)):
    (t,r,v)=LeapFrog(r,v,h,t)
    tsRK4=np.append(tsRK4,t)
    rs=np.concatenate((rs,np.array([r])))
    vs=np.concatenate((vs,np.array([v])))
    
    if i>1:
        if r[1]/rs[-2][1]<0:
            if rs[-1][0]>0:
                counterRK4 += 1
                periodRK4=np.append(periodRK4,t)
                
    
[xLeap,yLeap]=rs.transpose()
[VxLeap,VyLeap]=vs.transpose()

r=np.array([rM[0],rM[1],rM[2],rM[3]])  # initial x,y displacements; initial x,y velocities
Ro=r  # save initial conditions
t=0  # time start
tsLeap=np.array([t])  # store all times
rs=np.array([r])  # store all solution points

periodLeap=[]

for i in range(int(tmax/h)):
    (t,r)=rk4(r,derivRK4,t,h)
    tsLeap=np.append(tsLeap,t)
    rs=np.concatenate((rs,np.array([r])))
    
    if i>1:
        if r[1]/rs[-2][1]<0:
            if rs[-1][0]>0:
                counterLeap += 1
                periodLeap=np.append(periodLeap,t)
        
[xRK4,yRK4,VxRK4,VyRK4]=rs.transpose()

LeapPeriod=periodLeap[2]-periodLeap[1]
RK4Period=periodRK4[2]-periodRK4[1]


L=((Ro[0]**2+Ro[1]**2)**(1/2))*((Ro[2]**2+Ro[3]**2)**(1/2))
initialE=(1/2*Ro[2]**2+1/2*Ro[3]**2)-1/((Ro[0]**2+Ro[1]**2)**(1/2))
print('Initial Energy:')
print(initialE)
eccentricity = (1+2*(L**2)*initialE)**(1/2)
print('Eccentricity:')
print(eccentricity)
lastPeriodPhi=np.arctan2(yLeap[-1],xLeap[-1])
firstPeriodPhi=np.arctan2(yLeap[0],xLeap[0])
phi=np.arange(0,2*np.pi*counterLeap-1+firstPeriodPhi+lastPeriodPhi,(2*np.pi*counterLeap-1+firstPeriodPhi+lastPeriodPhi)/((tmax/h)+1))
#phi=np.arange(0,maxPhi,maxPhi*h/(tmax+1))
R=L**2/(1+eccentricity*np.cos(phi))
xAna=R*np.cos(phi)
yAna=R*np.sin(phi)
xAna=xAna[:-1]
yAna=yAna[:-1]

plt.figure()
plt.xlabel('x-Displacement')
plt.ylabel('y-Displacement')
plt.title('Analytic vs RK4 vs Leap Frog')
plt.axis('equal')
plt.plot(xAna,yAna, label = 'Analytic')
plt.plot(xRK4,yRK4, label = 'RK4')
plt.plot(xLeap,yLeap, label = 'Leap Frog')
plt.legend()

DifAnaRK4=(xAna**2+yAna**2)**(1/2)-(xRK4**2+yRK4**2)**(1/2)
DifAnaLeap=(xAna**2+yAna**2)**(1/2)-(xLeap**2+yLeap**2)**(1/2)
plt.figure()
plt.plot(tsRK4/RK4Period,DifAnaRK4, label = 'RK4', color='orange')
plt.plot(tsLeap/LeapPeriod,DifAnaLeap, label = 'Leapfrog', color='green')
plt.xlabel('Time / Orbital periods')
plt.ylabel('DIfference in electron displacement with respect to the \n origin between the analytic and numerical solutions')
plt.title('Analytic solution vs numerical approximations')
plt.legend()




"""
AnaLeap=[]
AnaRK4=[]
Vel=[]
Dis=[]


for i in range(20,100,20):
    for j in range(20,100,20):
        for k in range(20,100,20):
            for l in range(20,100,20):
                
                r=np.array([i,j,k/100,l/100])
                Ro=r  # save initial conditions
                t=0  # time start
                tmax=500  # time end
                h=1  # step size
                tsRK4=np.array([t])  # store all times
                rsRK4=np.array([r])  # store all solution points
                
                for m in range(int(tmax/h)):
                    (t,r)=rk4(r,derivRK4,t,h)
                    tsRK4=np.append(tsRK4,t)
                    rsRK4=np.concatenate((rsRK4,np.array([r])))

                [xRK4,yRK4,VxRK4,VyRK4]=rsRK4.transpose()
                
                r=np.array([i,j])  # initial x,y displacements
                v=np.array([k/100,l/100])  #initial x,y velocities, E<0 -> ellipse, E=0 -> parabola, E>0 -> hyperbola
                Ro=np.array([r[0],r[1],v[0],v[1]])
                t=0  # time start
                tmax=500  # time end
                h=1  # step size
                tsLeap=np.array([t])  # store all times
                rsLeap=np.array([r])  # store all solution points
                vsLeap=np.array([v])  # store all solution points
                
                for n in range(int(tmax/h)):
                    (t,r,v)=LeapFrog(r,v,h,t)
                    tsLeap=np.append(tsLeap,t)
                    rsLeap=np.concatenate((rsLeap,np.array([r])))
                    vsLeap=np.concatenate((vsLeap,np.array([v])))
    
                [xLeap,yLeap]=rsLeap.transpose()
                [VxLeap,VyLeap]=vsLeap.transpose()

                L=((Ro[0]**2+Ro[1]**2)**(1/2))*((Ro[2]**2+Ro[3]**2)**(1/2))
                initialE=(1/2*Ro[2]**2+1/2*Ro[3]**2)-1/((Ro[0]**2+Ro[1]**2)**(1/2))
                eccentricity = (1+2*(L**2)*initialE)**(1/2)
                
                if (eccentricity >= 1):
                    maxPhi=np.arctan2(yLeap[-1],xLeap[-1])
                    phi=np.arange(0,maxPhi,0.01)
                    R=L**2/(1+eccentricity*np.cos(phi))
                    xAna=R*np.cos(phi)
                    yAna=R*np.sin(phi)
                
                else:
                    phi=np.arange(0,2*np.pi,0.01)
                    R=L**2/(1+eccentricity*np.cos(phi))
                    xAna=R*np.cos(phi)
                    yAna=R*np.sin(phi)
                
                AnaRK4dif=np.mean((xAna**2+yAna**2)**(1/2))-np.mean((xRK4**2+yRK4**2)**(1/2))
                AnaRK4=np.append(AnaRK4,AnaRK4dif)
                
                AnaLeapdif=np.mean((xAna**2+yAna**2)**(1/2))-np.mean((xLeap**2+yLeap**2)**(1/2))
                AnaLeap=np.append(AnaLeap,AnaLeapdif)
                
                Vel=np.append(Vel,((k/100)**2+(l/100)**2)**(1/2))
                Dis=np.append(Dis,(i**2+j**2)**(1/2))

fig=plt.figure(figsize=(10,10))  # Plotting a figure with size (10,10)
ax=fig.gca(projection='3d')  # Define a new environment to use a 3D axis

ax.scatter(Vel,Dis,abs(AnaLeap),zdir='z', s=20, c='c', marker='o',label='Leap Frog', depthshade=True) # Plot all Leap points as cyan circles, opacity is a 3D effect
ax.scatter(Vel,Dis,abs(AnaRK4),zdir='z', s=20, c='r', marker='^',label='RK4', depthshade=True) # Plot all RK4 points as red triangles, opacity is a 3D effect    

ax.set_xlabel('Initial Velocity Magnitude')  #
ax.set_ylabel('Initial Displacement Magnitude')  
ax.set_zlabel('Displacement Difference')  
plt.title('Analytic vs Numerical Approximation')  

ax.legend()  # Create legend
plt.legend(bbox_to_anchor=(0.3, 0.8),bbox_transform=plt.gcf().transFigure)  # Shift legend position to top left

plt.show()  # Display the plot
"""

















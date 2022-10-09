#Imports
import numpy as np
import matplotlib.pyplot as plt

#Function to calculate derivative
def derivMoving(r,elM,proM):
    return np.array([(-((r[0]-r[2])/(((r[0]-r[2])**2+(r[1]-r[3])**2)**(3/2))+(r[0]-r[4])/(((r[0]-r[4])**2+(r[1]-r[5])**2)**(3/2)))),\
                     (-((r[1]-r[3])/(((r[0]-r[2])**2+(r[1]-r[3])**2)**(3/2))+(r[1]-r[5])/(((r[0]-r[4])**2+(r[1]-r[5])**2)**(3/2)))),\
                     elM/proM*(-((r[2]-r[0])/(((r[2]-r[0])**2+(r[3]-r[1])**2)**(3/2))-(r[2]-r[4])/(((r[2]-r[4])**2+(r[3]-r[5])**2)**(3/2)))),\
                     elM/proM*(-((r[3]-r[1])/(((r[2]-r[0])**2+(r[3]-r[1])**2)**(3/2))-(r[3]-r[5])/(((r[2]-r[4])**2+(r[3]-r[5])**2)**(3/2)))),\
                     elM/proM*(-((r[4]-r[0])/(((r[4]-r[0])**2+(r[5]-r[1])**2)**(3/2))-(r[4]-r[2])/(((r[4]-r[2])**2+(r[5]-r[3])**2)**(3/2)))),\
                     elM/proM*(-((r[5]-r[1])/(((r[4]-r[0])**2+(r[5]-r[1])**2)**(3/2))-(r[5]-r[3])/(((r[4]-r[2])**2+(r[5]-r[3])**2)**(3/2))))])  

def LeapFrogMoving(r,v,dt,t,elM,proM):
    v=v+ 1/2 *dt*derivMoving(r,elM,proM)
    r=r+dt*v
    v=v +1/2*dt*derivMoving(r,elM,proM)
    t=t+dt
    return (t,r,v) 

def deriv(r,a,b):
    return np.array([-((r[0]-a[0])/(((r[0]-a[0])**2+(r[1]-a[1])**2)**(3/2))+(r[0]-b[0])/(((r[0]-b[0])**2+(r[1]-b[1])**2)**(3/2))),-((r[1]-a[1])/(((r[0]-a[0])**2+(r[1]-a[1])**2)**(3/2))+(r[1]-b[1])/(((r[0]-b[0])**2+(r[1]-b[1])**2)**(3/2)))])

def LeapFrog(r,v,dt,t,a,b):
    v=v+ 1/2 *dt*deriv(r,a,b)
    r=r+dt*v
    v=v +1/2*dt*deriv(r,a,b)
    t=t+dt
    return (t,r,v)


#Constants
elM=9.109e-31
proM=1.673e-27
e=1.602e-19
ratio=proM/elM
k=(1.602e-19)**2/(4*np.pi*8.854e-12)


fig, ax= plt.subplots(2)

for k in range(1,10,1):
    
    #Master Conditions
    rM=np.array([20,0,0,0.3])
    aM=np.array([5,0,0,0])
    bM=np.array([-5,0,0,0])
    tmax=2000  # time end
    h=k/100  # step size

    #Counters
    statCounter=0
    movCounter=0
    

    r=np.array([rM[0],rM[1]])  # initial electron x,y displacements (20,0)
    v=np.array([rM[2],rM[3]])  #initial electron x,y velocities, E<0 -> ellipse, E=0 -> parabola, E>0 -> hyperbola (0,0.15)
    Ro=np.array([r[0],r[1],v[0],v[1]])  # save initial electron values

    a=np.array([aM[0],aM[1]])  # 1st proton x,y displacements
    Ao=a  # save initial 1st proton values
    b=np.array([bM[0],bM[1]])  # 2nd proton x,y displacements
    Bo=b  # save initial 2ndst proton values

    t=0  # time start
    ts=np.array([t])  # store all times
    rs=np.array([r])  # store all solution points
    vs=np.array([v])  # store all solution points

    statProtonPeriod=[]

    for i in range(int(tmax/h)):
        (t,r,v)=LeapFrog(r,v,h,t,a,b)
        ts=np.append(ts,t)
        rs=np.concatenate((rs,np.array([r])))
        vs=np.concatenate((vs,np.array([v])))
    
        if i>1:
            if r[1]/rs[-2][1]<0:
                if rs[-1][0]>a[0]:
                    statCounter += 1
                    statProtonPeriod=np.append(statProtonPeriod,t)

    [x,y]=rs.transpose()
    [Vx,Vy]=vs.transpose()

    period=statProtonPeriod[1]-statProtonPeriod[0]



    #MOVING PROTONS
    rEl=np.array([rM[0],rM[1]])  # initial electron x,y displacements
    vEl=np.array([rM[2],rM[3]])  #initial electron x,y velocities, E<0 -> ellipse, E=0 -> parabola, E>0 -> hyperbola
    Ro=np.array([rEl[0],rEl[1],vEl[0],vEl[1]])  # save initial electron values

    a=np.array([aM[0],aM[1]])  # 1st proton x,y displacements
    aV=np.array([aM[2],aM[3]])  # 2nd proton x,y velocities
    Ao=np.array([a[0],a[1],aV[0],aV[1]])  # save initial 1st proton values
    b=np.array([bM[0],bM[1]])  # 2nd proton x,y displacements
    bV=np.array([bM[2],bM[3]])  # 2nd proton x,y velocities
    Bo=np.array([b[0],b[1],bV[0],bV[1]])  # save initial 2nd proton values

    r=np.array([rEl[0],rEl[1],a[0],a[1],b[0],b[1]])
    v=np.array([vEl[0],vEl[1],aV[0],aV[1],bV[0],bV[1]])

    t=0  # time start
    ts=np.array([t])  # store all times
    rs=np.array([r])  # store all solution points
    vs=np.array([v])  # store all solution points

    movingProtonPeriod=[]

    for j in range(int(tmax/h)):
        (t,r,v)=LeapFrogMoving(r,v,h,t,elM,proM)
        ts=np.append(ts,t)
        rs=np.concatenate((rs,np.array([r])))
        vs=np.concatenate((vs,np.array([v])))
    
        if j>1:
            if r[1]/rs[-2][1]<0:
                if rs[-1][0]>rs[-1][2]:
                    movCounter += 1
                    movingProtonPeriod=np.append(movingProtonPeriod,t)
                

    [xEl,yEl,xA,yA,xB,yB]=rs.transpose()
    [VxEl,VyEl,VxA,VyA,VxB,VyB]=vs.transpose() 


    DifEl=(xEl**2+yEl**2)**(1/2)-(x**2+y**2)**(1/2)
    DifPro=((xA-xB)**2+(yA-yB)**2)**(1/2)-(Ao[0]-Bo[0])
    ax[0].plot(ts/period,DifEl, label = k/100)
    ax[1].plot(ts/period,DifPro)


ax[0].set_xlabel('Time / Orbital Periods')
ax[0].set_ylabel('Electron displacment')
ax[1].set_xlabel('Time / Orbital Periods')
ax[1].set_ylabel('Proton separation')
ax[0].legend(bbox_to_anchor=(1.05, 1), loc=2, borderaxespad=0.,title = 'Step Size')
fig.suptitle('Instability with Step Size')
fig.tight_layout()
fig.subplots_adjust(top=0.9)

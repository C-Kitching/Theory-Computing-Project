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

initialE=(1/2*Ro[2]**2+1/2*Ro[3]**2)-1/((Ro[0]**2+Ro[1]**2)**(1/2))


avELeap=[]
stepLeap=[]

for i in range(1,900,10):
    r=np.array([Ro[0],Ro[1]])  # initial x,y displacements
    v=np.array([Ro[2],Ro[3]])  #initial x,y velocities
    t=0  # time start
    tmax=1000  # time end
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
    avELeap=np.append(avELeap,np.mean(E))
    stepLeap=np.append(stepLeap,h)



avERK4=[]
stepRK4=[]

for i in range(1,900,10):
    r=Ro  # initial x,y displacements; initial x,y velocities
    t=0  # time start
    tmax=1000  # time end
    h=i/100  # step size
    ts=np.array([t])  # store all times
    rs=np.array([r])  # store all solution points

    for j in range(int(tmax/h)):
        (t,r)=rk4(r,derivRK4,t,h)
        ts=np.append(ts,t)
        rs=np.concatenate((rs,np.array([r])))
        
    [x,y,Vx,Vy]=rs.transpose()
    E=(1/2*Vx**2+1/2*Vy**2)-1/((x**2+y**2)**(1/2))
    avERK4=np.append(avERK4,np.mean(E))
    stepRK4=np.append(stepRK4,h)


"""
fig, (ax1, ax2) = plt.subplots(1,2)

ax1.plot(stepLeap,avELeap, label='Leap Frog')
ax1.plot([np.amin(stepLeap),np.amax(stepLeap)],[initialE,initialE], label = 'True value')
ax1.plot([np.amin(stepLeap),np.amax(stepLeap)],[initialE*1.001,initialE*1.001], color='green', dashes=[6, 2], label = '0.1% error')
ax1.plot([np.amin(stepLeap),np.amax(stepLeap)],[initialE*0.999,initialE*0.999], color='green', dashes=[6, 2])
ax1.set_title('Leap Frog', fontsize=13)
ax1.invert_xaxis()

ax2.plot(stepRK4,avERK4, color='red', label='RK4')
ax2.plot([np.amin(stepRK4),np.amax(stepRK4)],[initialE,initialE], color='orange')
ax2.plot([np.amin(stepRK4),np.amax(stepRK4)],[initialE*1.001,initialE*1.001], color='green', dashes=[6, 2])
ax2.plot([np.amin(stepRK4),np.amax(stepRK4)],[initialE*0.999,initialE*0.999], color='green', dashes=[6, 2])
ax2.set_title('RK4', fontsize=13)
ax2.invert_xaxis()

fig.text(0.5, 0.1, 'Step Size', ha='center', fontsize=12)
fig.text(0.5, 0.95, 'Optimal Step Size', ha='center', fontsize=14)
fig.text(0.04, 0.5, 'Energy', va='center', rotation='vertical', fontsize=12)
#fig.suptitle('Optimal Step Size', fontsize=14)
fig.legend(bbox_to_anchor=(0.75,0.08), ncol=4)
fig.tight_layout()
fig.subplots_adjust(bottom=0.2, wspace=0.4)
fig.subplots_adjust(top=0.88)
"""

fig=plt.figure()
plt.plot(stepRK4,avERK4, color='red', label='RK4')
plt.plot(stepLeap,avELeap, color='blue', label='Leapfrog')
plt.plot([np.amin(stepRK4),np.amax(stepRK4)],[initialE,initialE], color='orange', label = 'True value')
plt.plot([np.amin(stepRK4),np.amax(stepRK4)],[initialE*1.001,initialE*1.001], color='green', dashes=[6, 2], label = '0.1% error')
plt.plot([np.amin(stepRK4),np.amax(stepRK4)],[initialE*0.999,initialE*0.999], color='green', dashes=[6, 2])
plt.title('Optimum time step')
plt.xlabel('Step-size')
plt.ylabel('Energy')
plt.gca().invert_xaxis()
plt.legend()
fig.tight_layout()









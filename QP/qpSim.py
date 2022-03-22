#QP Simulation
#Eric Fackelman
#March 2022

#20.3 update: having difficulty bc the number of steps changes the trajectories.
#Need to make trajectories independent of steps!

import numpy as np
import matplotlib.pyplot as plt
import constants_qp as c
import sys

#sys.argv[0] = python code name
#sys.argv[1] must be # of particles

#in order to do several particles (z), do a range and save values to a csv
w=1 #counter which isn't used elsewhere
z=int(sys.argv[1]) #number of particles
columns = [] #columns array for stacking the values for saving
header = [] #header row
filename = "../CSVs/FODOattempt_{:.0f}".format(c.B0*1e6)+"uT,m-2_"+str(z)+".csv"
with open(filename,"w") as f: #must make new/overwrite past file
  f.write("FODO arrays:\n")
print(filename)

#make constants and arrays outside loop
#print(c.ess_i/1e6*c.ess_t/c.e)
p0=np.zeros(3)
f = [0,0,0] #force initial vector
s=10
N=500 #don't change this!
n=N*s #number of steps
tinQP=c.l/c.v
dt = tinQP/N #step length based on time in QP
#print(dt)
dz = c.v*dt
Bf=[0,0,0] #[T/m]
Ef = [0,0,0] #[V/m]

#especially for the 3D arrays
pos = np.zeros([z,3,n]) #position, velocity, time, energy arrays for graphing
vels = np.zeros([z,3,n])
ts = np.zeros([z,n])
eng = np.zeros([z,n])
Bs = np.zeros([z,3,n])
xprod = np.zeros(3) #needed for cross product calculation
#for storing values in graphing arrays!

for w in range(z):
  #add the particle number to the header row, and this keeps columns lined up
  header += ["ts_"+str(w)+"      ",'Bx_'+str(w)+"      ",'By_'+str(w)+"      ",'Bz_'+str(w)+"      ",
        'x_'+str(w)+"       ",'y_'+str(w)+"       ",'z_'+str(w)+"       "] 

  class Particle:
    def __init__(self,pos,q,v,m):
      self.pos = pos #[m]
      self.q = q #[C]
      self.v = v #[m/s]
      self.m = m #[kg]
  p1 = Particle([(w+1)*1e-6,(w+1)*-1e-6,0],c.e,[0,0,c.v],c.mp) #actual particle 
  for j in range(3):
    #Bf[0] = - c.g147 * p1.pos[1]
    #Bf[1] = - c.g147 * p1.pos[0]
    Bf[j] = Bf[j] + c.dBydx[j]*(p1.pos[j]-p0[j]) #from L4 p13
  xprod = np.cross(p1.v,Bf)
  #print('v x B =',xprod)

  #print('Time in QP is',tinQP)
  t=0
  i=0
  while t <= s*tinQP:
    if p1.pos[2] <= c.l: #for 1st QP
      for j in range(3):
        Bf[j] = Bf[j] - c.dBydx[j]*(p1.pos[j]-p0[j]) #from L4 p13
      xprod = np.cross(p1.v,Bf)
    elif p1.pos[2] >= c.l and p1.pos[2] <= c.qp1start: #for drift between QPs
      for j in range(3):
        Bf[j] = 0 #from L4 p13
      xprod = np.cross(p1.v,Bf)
    elif p1.pos[2] >= c.qp2start and p1.pos[2] <= (c.qp2start+c.l): #for 2nd QP
      for j in range(3):
        Bf[j] = Bf[j] - c.dBydx2[j]*(p1.pos[j]-p0[j]) #from L4 p13
      xprod = np.cross(p1.v,Bf)
    elif p1.pos[2] > (c.qp2start + c.l): #for drift after QPs
      for j in range(3):
        Bf[j] = 0 #from L4 p13
      xprod = np.cross(p1.v,Bf)
  
    for j in range(3): #condense the lines of code
      f[j] = p1.q/p1.m * xprod[j] + p1.q/p1.m * Ef[j] #calculate force at this velocity
      p1.v[j] = p1.v[j] + f[j] * dt #calculate the velocity change from the force
      p1.pos[j] = p1.pos[j] + p1.v[j]*dt #move the particle

      #save values before calculations in order to save first ones
      vels[w,j,i] = p1.v[j] #save velocities for Larmor radius calc
      pos[w,j,i] = p1.pos[j] #save positions for graphing
      Bs[w,j,i] = Bf[j] #save Bfield for graphing
    en = 0.5 * np.sqrt(p1.v[0]**2 + p1.v[1]**2 + p1.v[2]**2)
    ts[w,i] = t #save time before changing
    eng[w,i] = en
    t = t + dt
    i = i + 1
  columns += [ts[w],Bs[w,0],Bs[w,1],Bs[w,2],pos[w,0],pos[w,1],pos[w,2]] #stack the values to save
  print(w,t,i)

#outside w for loop
#print(np.shape(Bs),np.shape(pos))
with open(filename,"a") as f: #append header
  np.savetxt(f,np.column_stack(header),fmt='%2s') 
with open(filename,"a") as f: #append values
  np.savetxt(f,np.column_stack(columns),delimiter=',',fmt='%+.3e') #force +/- and 3 decimals

plottitle = 'Particle through Quadrupole, B0={:.1f}'.format(c.B0*1e6)+r'$/mu$m'
#define fig outside loop
from mpl_toolkits.mplot3d import Axes3D
fig = plt.figure(tight_layout=True,figsize=(12,12))
s1 = fig.add_subplot(111, projection='3d')
#fig = plt.figure(figsize=(16,12))
#plt.subplots_adjust(hspace=0.25,wspace=0.35)
##s1 = fig.add_subplot(2,2,(1,2)) #4x4 plots
#s1 = fig.add_subplot(2,3,(1,3)) #6x6 plots
##s2 = fig.add_subplot(2,2,3) #4x4 plots
#s2 = fig.add_subplot(2,3,4) #6x6 plots
##s3 = fig.add_subplot(2,2,4) #4x4 plots
#s3 = fig.add_subplot(2,3,5) #6x6 plots
#s4 = fig.add_subplot(2,3,6)

for w in range(z):
  #plot inside loops to go over the particles

  #3D plot:
  s1.plot3D(pos[w,2],pos[w,0]*1e6,pos[w,1]*1e6)

  #2Dplots
  #s1.plot(ts[w]*1e9,pos[w,0]*1e6,label='X{:d}'.format(w))
  #s1.plot(ts[w]*1e9,pos[w,1]*1e6,label='Y{:d}'.format(w))
  
  #s2.plot(pos[w,2],pos[w,0]*1e6)
  
  #s3.plot(pos[w,2],pos[w,1]*1e6)

  #s1.plot(ts*1e6,pos[2]*1e6,'r.',label='Z')
  #print(ts[w,4],pos[w,0,4])
  #s2.plot(pos[0]*1e6,pos[1]*1e6)
  #s4.plot(ts*1e9,eng)
  #s4.grid(visible=True)
  #s4.set_xlabel('Time [ns]',fontsize=20)
  #s4.set_ylabel('Energy [eV]',fontsize=20)

  #s4 = fig.add_subplot(2,3,6) #attempt at plotting B, WIP
  #s4.plot(pos[2],Bf[0])
  #s4.plot(pos[2],Bf[1])
  #s4.grid(visible=True)
  #s4.set_xlabel('Z Position [m]',fontsize=20)
  #s4.set_ylabel('Magnetic Field [T/m]',fontsize=20)

#add options AFTER plotting
s1.grid(visible=True)
s1.set_title(r'Particle through Quadrupole, $\frac{d B}{d x,y}$='+'{:.1f}'.format(c.B0*1e6)+r'$\frac{\mu T}{m^2}$',fontsize=24,pad=20)
#s1.set_xlabel('Time [ns]',fontsize=20)
s1.set_ylabel('Particle X Position [$\mu$m]',fontsize=20,labelpad=10)
s1.set_zlabel('Particle Y Position [$\mu$m]',fontsize=20,labelpad=10)
s1.set_xlabel('Particle Z Position [m]',fontsize=20,labelpad=10)
#s1.legend(loc='lower left',bbox_to_anchor=(-0.01,0.1),frameon=False)
#s1.text(-0.5,440,,fontsize=20)

#s2.set_xlabel('Particle X Position [$\mu$m]',fontsize=20)
#s2.set_xlabel('Particle Z Position [m]',fontsize=20)
#s2.set_ylabel('Particle X Position [$\mu$m]',fontsize=20)
#s2.grid(visible=True)

#s3.set_xlabel('Particle Z Position [m]',fontsize=20)
#s3.set_ylabel('Particle Y Position [$\mu$m]',fontsize=20)
#s3.grid(visible=True)

#print(np.mean(eng))
savename='../pictures/FODO_{:.0f}ns_{:.0f}e_P0x{:.1f},y{:.1f}um_{:.1f}mTm-2__N{:d}_p{:d}.pdf'.format(t*1e9,p1.q*1/c.e,pos[0,0,0]*1e6,pos[0,1,0]*1e6,c.B0*1e3,N,z)
print(savename)
#plt.savefig(savename)
plt.show()
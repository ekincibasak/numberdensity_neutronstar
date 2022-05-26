#calculate number density of a neutron star
import matplotlib.pyplot as plt
import numpy as np
import h5py
import scidata.carpet.hdf5 as hdf5
# Import math Library
import math 
from pylab import *

import pandas as pd
import timeit
import sys
reload(sys)
sys.setdefaultencoding('utf-8')

rc('font', weight='bold')
cons=6.06E-6
#1.234E-6
#row:y:axis1
#column:x:axis0
#plane:z:axis2

#im column index





######


#This function uses scidata-routines in order to
#read the data. It takes as parameters the file, iteration,
#refinement-level (which grid) and the component.
def get_data(df,it,tl,rl):
    grid=df.get_reflevel(iteration=it,timelevel=tl,reflevel=rl)
    xxx,yxx,zxx=grid.mesh()
#    xxy,yxy=grid.mesh()
    datxx=df.get_reflevel_data(grid,iteration=it)
#    datxy=df.get_reflevel_data(grid,iteration=it)
    return xxx,yxx,zxx,datxx
######################################################################




cs=(2.99*(10**23))*((1/197)**4)
gammac=4.68*(10**(-19))
delta=6
x_1=189/(367*((math.pi)**2))
x_2=21/((367*math.pi)**4)
x_3=3/((1835*math.pi)**6)
############################################################
ktt = []
####################################################################
it=1024
old=0
tl=0
rl=4
#while it <292864 :



datafiletemp = hdf5.dataset("/home/basak/basak/data0/temperature.xyz.all.h5")
datafilerho = hdf5.dataset("/home/basak/basak/data0/rho.xyz.all.h5")
datafileye = hdf5.dataset("/home/basak/basak/data0/Y_e.xyz.all.h5")
datafilevel1= hdf5.dataset("/home/basak/basak/data0/vel1.xyz.all.h5")
datafilevel2 = hdf5.dataset("/home/basak/basak/data0/vel2.xyz.all.h5")
datafilevel0 = hdf5.dataset("/home/basak/basak/data0/vel0.xyz.all.h5")


initial=timeit.default_timer()   



kt=((it*0.028235)*4.926)*0.001
time = float("{0:.2f}".format(kt))
f = open('eos{}.txt'.format(time), 'wb')
fn = open('vect{}.txt'.format(time), 'wb')
fall = open('all{}.txt'.format(time), 'wb')
ktt.append(kt)




xxx,yxx,zxx,rho = get_data(datafilerho,it,tl,rl)#6.column
  
xxx,yxx,zxx,Y_e = get_data(datafileye,it,tl,rl) 

xxx,yxx,zxx,velx = get_data(datafilevel0,it,tl,rl) #0.column
xxx,yxx,zxx,vely = get_data(datafilevel1,it,tl,rl) #1.column
xxx,yxx,zxx,velz = get_data(datafilevel2,it,tl,rl) #2.column   
xxx,yxx,zxx,Y_e = get_data(datafileye,it,tl,rl) 
xxx,yxx,zxx,Tempp = get_data(datafiletemp,it,tl,rl) 
getd=timeit.default_timer()   


denc=369.52 
aa=0.376
im=rho.shape[0] #column
jm=rho.shape[1] #row
km=rho.shape[2] #plane
dt=1024



ndi= np.empty((im, jm,km)) 
ndf = ndi.copy()
     
#dff = np.loadtxt('eos.table')
#dff = pd.read_csv (eos.table, sep='\s+')
#bden =dff[:,1]/369.52
#bden =dff[:,1]/369.52 
#Temp=dff[:,0]
#dmu =dff[:,4]
#arr_2d = np.reshape(arr, (2, 5), order='F')
 


#denpp = np.reshape(bden, (im, jm, km), order='F')
#tempp = np.reshape(Temp, (im, jm, km), order='F')
#yerp = np.reshape(dmu , (im, jm, km), order='F')

#To convert the mass density in geometrized units to a baryon density in units of nuclear saturation density


    
     
#print(dff.shape)
print(im*jm*km)  



# initial conditions for the initial grid


#initial condition
ndi[:,:,:]=(1-Y_e[:,:,:])*(rho[:,:,:])*denc  
initial=timeit.default_timer()   
   
#im=im-1
#jm=jm-1
#km=km-1           
            
 
#boundray conditions 
ndi[0,:,:]=(1-Y_e[0,:,:])*(rho[0,:,:])*denc   
ndi[jm-1,:,:]=(1-Y_e[jm-1,:,:])*(rho[jm-1,:,:])*denc    
 
 
ndi[:,:,0]=(1-Y_e[:,:,0])*(rho[:,:,0])*denc   
ndi[:,:,km-1]=(1-Y_e[:,:,km-1])*(rho[:,:,km-1])*denc  


ndi[:,0,:]=(1-Y_e[:,0,:])*(rho[:,0,:])*denc   
ndi[:,im-1,:]=(1-Y_e[:,im-1,:])*(rho[:,im-1,:])*denc  

bound=timeit.default_timer()   

# evolution time step
ndf[1:-1,1:-1,1:-1]=ndi[1:-1,1:-1,1:-1]+\
                -(dt/(2*aa))*((ndi[2:,1:-1,1:-1]-ndi[:-2,1:-1,1:-1])*velx[1:-1,1:-1,1:-1]\
                -(ndi[1:-1,2:,1:-1]-ndi[1:-1,:-2,1:-1])*vely[1:-1,1:-1,1:-1]\
                -(ndi[1:-1,1:-1,2:]-ndi[1:-1,1:-1,:-2])*velz[1:-1,1:-1,1:-1])\
                -(dt/(2*aa))*((velx[2:,1:-1,1:-1]-velx[:-2,1:-1,1:-1])*ndi[1:-1,1:-1,1:-1]\
                -(vely[1:-1,2:,1:-1]-vely[1:-1,:-2,1:-1])*ndi[1:-1,1:-1,1:-1]\
                -(velz[1:-1,1:-1,2:]-velz[1:-1,1:-1,:-2])*ndi[1:-1,1:-1,1:-1])
evola=timeit.default_timer()   
ndf=ndi.copy()              
#             f.write("%e %e %e \n" % (Tempw,  Tempp[j][i][k], pdn ))
denpp = np.reshape(rho, ((im)*(jm)*(km)), order='F')
for z in range((im-2)*(jm-2)*(km-2)):
#for row in denpp:

    fn.write("%e  \n" % (denpp[z] ))

fn.close() 

#             fall.write("%e %e %e %e %e %e %e %e \n" % (Tempw, Tempp[j][i][k], (rho[j][i][k])*denc,denpp[j][i][k]*369.52 , pdn, 




f.close()              
alld=timeit.default_timer()   
print('alld='+str((alld))+'s')
print('init='+str((initial-getd))+'s')
print('boundray'+str((bound-initial))+'s')
print('evulation'+str((evola-bound))+'s')



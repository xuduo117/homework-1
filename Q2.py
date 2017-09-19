#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Tue Sep 12 13:35:42 2017

@author: xuduo
"""



import numpy as np
import matplotlib.pyplot as plt
from astropy.time import Time

times = ['2017-08-01T00:00:00', '2017-12-31T23:59:59']
t = Time(times, format='isot', scale='utc')
start_time,end_time=t.jd1-2454424.857-111.43637*31

###Newton-Raphson

def derivative(f, x, h):
      return (f(x+h) - f(x-h)) / (2.0*h)  

def quadratic(x):
    return x-eccentricity*np.sin(x)-mean_anomaly     

def solve(f, x0, h):
    lastX = x0
    nextX = lastX + 10* h 
    while (abs(lastX - nextX) > h):  
        newY = f(nextX)                     
#        print "f(", nextX, ") = ", newY     
        lastX = nextX
        nextX = lastX - newY / derivative(f, lastX, h) 
    return nextX

##### angle calculation
  
def angle_xy(x_pos,y_poz):
    angle=np.mod(np.arctan(y_poz/x_pos)*180/np.pi,360)
    if x_pos <0:
        angle=np.mod(np.arctan(y_poz/x_pos)*180/np.pi+180,360)
    return angle


######constant cgs

c=29979245800.0
G=6.674079999999999e-08
R_sun=69570000000.0
R_jup=7149200000.0
M_sun=1.9884754153381438e+33
M_earth=5.972364730419773e+27
M_jup=1.8981871658715508e+30
AU=14959787070000.0
yr=3.15576e7
day=86400.0
#######input paprameter: time, the orbital elements, and the object masses 

M_1=M_sun*0.98
M_2=M_jup*3.94
#M_2=M_sun*1.0

#t_all=np.linspace(0,3e7,num=100000)
t_all=np.linspace(start_time*day,end_time*day,num=10000)

RV_z1=[]
RV_z2=[]
seperation_dis=[]
position_angle=[]
energy_all=[]
angular_momentum=[]


for ctt_t in range(len(t_all)):
    t=t_all[ctt_t]
    
    a=AU*0.449
    eccentricity=0.93366
    I=89.26/180.0*np.pi     ###Inclination	
    Omega=160.98/180.0*np.pi       ###Longitude of ascending node -19.02 or +160.98
    omega=300.651/180.0*np.pi       ###Argument of perihelion
    t_0=0     ### time of periastron passage
    
    ####  calculat eccentric anomaly 
    n=(G*(M_1+M_2)/a**3)**0.5
    mean_anomaly=n*(t-t_0)
    eccentric_anomaly = solve(quadratic, np.pi/2, 1e-8)    # call the solver
    #print "solution: x = ", eccentric_anomaly        # print the result
    E=eccentric_anomaly*1.0
    e=eccentricity*1.0
    ####  rotational array
    
    P=np.array([[np.cos(Omega)*np.cos(omega)-np.sin(Omega)*np.sin(omega)*np.cos(I)],
                [np.sin(Omega)*np.cos(omega)+np.cos(Omega)*np.sin(omega)*np.cos(I)],
                [np.sin(omega)*np.sin(I)]])
    
    Q=np.array([[-np.cos(Omega)*np.sin(omega)-np.sin(Omega)*np.cos(omega)*np.cos(I)],
                [-np.sin(Omega)*np.sin(omega)+np.cos(Omega)*np.cos(omega)*np.cos(I)],
                [np.cos(omega)*np.sin(I)]])

#    P=np.array([[np.cos(Omega)*np.cos(omega)-np.sin(Omega)*np.sin(omega)*np.cos(I)],
#                [-np.sin(Omega)*np.cos(omega)-np.cos(Omega)*np.sin(omega)*np.cos(I)],
#                [np.sin(omega)*np.sin(I)]])
#    
#    Q=np.array([[np.cos(Omega)*np.sin(omega)+np.sin(Omega)*np.cos(omega)*np.cos(I)],
#                [-np.sin(Omega)*np.sin(omega)+np.cos(Omega)*np.cos(omega)*np.cos(I)],
#                [-np.cos(omega)*np.sin(I)]])
        
    
    r_vector=np.dot(a*(np.cos(E)-e),P)+np.dot(a*(1-e**2)**0.5*np.sin(E),Q)
    r_scale=np.linalg.norm(r_vector)
    v_vector=np.dot(-a**2*n/r_scale*np.sin(E),P)+np.dot(a**2*n/r_scale*(1-e**2)**0.5*np.cos(E),Q)
    
    
    #print r_vector/AU, r_scale/AU
    #print v_vector/1e5  
    
    r_vector=np.dot(a*(np.cos(E)-e),P)+np.dot(a*(1-e**2)**0.5*np.sin(E),Q)
    r_scale=np.linalg.norm(r_vector)
    v_vector=np.dot(-a**2*n/r_scale*np.sin(E),P)+np.dot(a**2*n/r_scale*(1-e**2)**0.5*np.cos(E),Q)
     
    
    #print r_vector/AU, r_scale/AU
    #print v_vector/1e5  
    
    r_COM=r_vector*M_2/(M_1+M_2)
    v_COM=v_vector*M_2/(M_1+M_2)
    
    r_M1_COM=0-r_COM
    r_M2_COM=r_vector-r_COM
    v_M1_COM=0-v_COM
    v_M2_COM=v_vector-v_COM
    RV_z2.append(-v_M2_COM[2][0])
    RV_z1.append(-v_M1_COM[2][0])    
    seperation_dis.append(np.linalg.norm(r_vector[0:2]))
    energy_all.append(0.5*np.linalg.norm(v_vector)**2-G*(M_1+M_2)/r_scale)
    angular_momentum.append(np.linalg.norm(np.cross(r_vector[:,0],v_vector[:,0])))
    
#print r_M1_COM, r_M2_COM,v_M1_COM,v_M2_COM

#print 	np.linalg.norm(r_vector[1:])/AU
plt.figure(1)
plt.plot(t_all/day-start_time,np.asarray(RV_z1)/1e5)
inter_month=1.0 ##653.59477124183
labels=['8/1','9/1', '10/1', '11/1', '12/1','1/1']
x_lable=[inter_month*0,inter_month*31,inter_month*61,inter_month*92,inter_month*122,inter_month*153]
plt.xticks(x_lable, labels)
plt.xlabel('Date')
plt.ylabel('V (km/s)')
plt.title('RV-star')
plt.savefig('latex/Q2_rv_star'+'.pdf',bbox_inches='tight')
plt.clf()

print np.argsort(RV_z1)[0:10]
print np.argsort(RV_z1)[-10:]

time_rv_max=[t_all[1576]/day+2457879.38447,t_all[8864]/day+2457879.38447]
time_rv_min=[t_all[1465]/day+2457879.38447,t_all[8754]/day+2457879.38447]

t_max=Time(time_rv_max, format='jd')
t_min=Time(time_rv_min, format='jd')

print t_max.datetime
print t_min.datetime

#plt.xlabel([])

#plt.figure(2)
#plt.plot(t_all,seperation_dis)
#
#plt.figure(3)
#plt.plot(t_all,position_angle)

#print angle_yz(r_vector[1],r_vector[2])
#print r_vector
#t_all=t_all/yr
    
#plt.plot(t_all,energy_all)
#plt.ylim([np.nanmin(energy_all)*1.5,0])
#plt.show()

#plt.plot(t_all,angular_momentum)
#plt.ylim([np.nanmin(energy_all)*1.5,0])
#plt.show()





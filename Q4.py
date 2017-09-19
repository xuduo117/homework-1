#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Sun Sep 17 20:11:25 2017

@author: xuduo
"""




import numpy as np
import matplotlib.pyplot as plt
from astropy.time import Time

times = ['2017-01-01T00:00:00', '2022-01-01T00:00:00']
t = Time(times, format='isot', scale='utc')
start_time,end_time=t.jd1-2454424.857-111.43637*29
duration_time=end_time-start_time

def time_utc2jd(days):
    t=Time([2454424.857+111.43637*31+days], format='jd')
    return t.datetime


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
pc=206265*AU
deg2mas=180/np.pi*3600*1000
#######input paprameter: time, the orbital elements, and the object masses 

M_1=M_sun*0.98
M_2=M_jup*3.94
#M_2=M_sun*1.0

#t_all=np.linspace(0,3e7,num=100000)
#t_all=np.linspace(start_time*day,end_time*day,num=100000)
t_all=np.random.uniform(low=start_time*day,high=end_time*day,size=100)
t_all=np.sort(t_all)

RV_z1=[]
RV_z2=[]
seperation_dis=[]
seperation_dis_x=[]
seperation_dis_y=[]
position_angle=[]
energy_all=[]
angular_momentum=[]
position_z=[]
t_transite=[]
seperation_dis_transite=[]
x_ra_earth=[]
y_dec_earth=[]
PM_ra=[]
PM_dec=[]

for ctt_t in range(len(t_all)):
    t=t_all[ctt_t]
    
    radius_M1=0.98*R_sun
    radius_M2=R_jup*0.921
    a=AU*0.449
    eccentricity=0.93366
    I=89.285/180.0*np.pi     ###Inclination	
#    I=90.0/180.0*np.pi     ###Inclination	
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
    
    r_COM=r_vector*M_2/(M_1+M_2)
    v_COM=v_vector*M_2/(M_1+M_2)
    
    r_M1_COM=0-r_COM
    r_M2_COM=r_vector-r_COM
    v_M1_COM=0-v_COM
    v_M2_COM=v_vector-v_COM
    RV_z2.append(-v_M2_COM[2][0])
    RV_z1.append(-v_M1_COM[2][0])    
    seperation_dis.append(np.linalg.norm(r_vector[0:2]))
    seperation_dis_x.append(r_M1_COM[0][0])
    seperation_dis_y.append(r_M1_COM[1][0])
    position_angle.append(angle_xy(r_vector[0],r_vector[1]))
    energy_all.append(0.5*np.linalg.norm(v_vector)**2-G*(M_1+M_2)/r_scale)
    angular_momentum.append(np.linalg.norm(np.cross(r_vector[:,0],v_vector[:,0])))
    position_z.append(r_vector[2,0])
    if (r_vector[2,0] >0) & (np.linalg.norm(r_vector[0:2])<radius_M1):
        t_transite.append(t)
        seperation_dis_transite.append(np.linalg.norm(r_vector[0:2]))
    x_ra_earth.append(17.13*np.cos(np.mod(t/day-318.9882700001567-91.25,365)/365.0*2*np.pi))
    y_dec_earth.append(17.13*np.sin(np.mod(t/day-318.9882700001567-91.25,365)/365.0*2*np.pi))
    PM_ra.append((t/day-start_time)/365.0*45.76)
    PM_dec.append((t/day-start_time)/365.0*16.56)

    
ra_gaia_err=np.random.normal(loc=0.0, scale=0.133,size=100)
dec_gaia_err=np.random.normal(loc=0.0, scale=0.133,size=100)

plt.figure(2)
plt.plot(ra_gaia_err,dec_gaia_err,'*')
plt.xlabel(r'$\Delta \alpha$ (mas)')
plt.ylabel(r'$\Delta \delta$ (mas)')
plt.title('Gaia uncertainty')
#plt.legend()
plt.savefig('latex/Q4_gaia_radec'+'.pdf',bbox_inches='tight')
plt.clf()



year_day=duration_time/5.0

plt.figure(2)
#plt.xlim([start_time,end_time])
plt.plot(t_all/day,np.asarray(seperation_dis_y)/58./pc*deg2mas,'*')
labels=['2017','2018', '2019', '2020', '2021','2022']
x_lable=start_time+np.arange(0,5,1)*year_day
plt.xticks(x_lable, labels)
plt.xlabel(r'Year')
plt.ylabel(r'$\Delta \alpha$ (mas)')
plt.title('RA offset')
#plt.legend()
plt.savefig('latex/Q4_1_ra'+'.pdf',bbox_inches='tight')
plt.clf()


plt.figure(2)
#plt.xlim([start_time,end_time])
plt.plot(t_all/day,np.asarray(seperation_dis_x)/58./pc*deg2mas,'*')
labels=['2017','2018', '2019', '2020', '2021','2022']
x_lable=start_time+np.arange(0,5,1)*year_day
plt.xticks(x_lable, labels)
plt.xlabel(r'Year')
plt.ylabel(r'$\Delta \delta$ (mas)')
plt.title('Dec offset')
#plt.legend()
plt.savefig('latex/Q4_1_dec'+'.pdf',bbox_inches='tight')
plt.clf()


plt.figure(2)
#plt.xlim([start_time,end_time])
plt.plot(np.asarray(seperation_dis_y)/58./pc*deg2mas,np.asarray(seperation_dis_x)/58./pc*deg2mas,'*')
plt.xlabel(r'$\Delta \alpha$ (mas)')
plt.ylabel(r'$\Delta \delta$ (mas)')
plt.title('planet influence')
plt.savefig('latex/Q4_1_radec'+'.pdf',bbox_inches='tight')
plt.clf()

plt.figure(2)
#plt.xlim([start_time,end_time])
plt.plot(t_all/day,np.asarray(seperation_dis_y)/58./pc*deg2mas+ra_gaia_err,'*')
labels=['2017','2018', '2019', '2020', '2021','2022']
x_lable=start_time+np.arange(0,5,1)*year_day
plt.xticks(x_lable, labels)
plt.xlabel(r'Year')
plt.ylabel(r'$\Delta \alpha$ (mas)')
plt.title('RA offset')
#plt.legend()
plt.savefig('latex/Q4_1_ra_err'+'.pdf',bbox_inches='tight')
plt.clf()


plt.figure(2)
#plt.xlim([start_time,end_time])
plt.plot(t_all/day,np.asarray(seperation_dis_x)/58./pc*deg2mas+dec_gaia_err,'*')
labels=['2017','2018', '2019', '2020', '2021','2022']
x_lable=start_time+np.arange(0,5,1)*year_day
plt.xticks(x_lable, labels)
plt.xlabel(r'Year')
plt.ylabel(r'$\Delta \delta$ (mas)')
plt.title('Dec offset')
#plt.legend()
plt.savefig('latex/Q4_1_dec_err'+'.pdf',bbox_inches='tight')
plt.clf()


plt.figure(2)
#plt.xlim([start_time,end_time])
plt.plot(np.asarray(seperation_dis_y)/58./pc*deg2mas+ra_gaia_err,np.asarray(seperation_dis_x)/58./pc*deg2mas+dec_gaia_err,'*')
plt.xlabel(r'$\Delta \alpha$ (mas)')
plt.ylabel(r'$\Delta \delta$ (mas)')
plt.title('planet influence+Gaia uncertainty')
plt.savefig('latex/Q4_1_radec_err'+'.pdf',bbox_inches='tight')
plt.clf()



#########  2 

#t = Time('2017-08-10T23:15:00', format='isot', scale='utc')
#start_time,end_time=t.jd1-2454424.857-111.43637*29


plt.figure(2)
#plt.xlim([start_time,end_time])
plt.plot(t_all/day,np.asarray(x_ra_earth)+np.asarray(seperation_dis_y)/58./pc*deg2mas+ra_gaia_err,'*')
labels=['2017','2018', '2019', '2020', '2021','2022']
x_lable=start_time+np.arange(0,5,1)*year_day
plt.xticks(x_lable, labels)
plt.xlabel(r'Year')
plt.ylabel(r'$\Delta \alpha$ (mas)')
plt.title('RA offset')
#plt.legend()
plt.savefig('latex/Q4_2_ra_err'+'.pdf',bbox_inches='tight')
plt.clf()


plt.figure(2)
#plt.xlim([start_time,end_time])
plt.plot(t_all/day,np.asarray(y_dec_earth)+np.asarray(seperation_dis_x)/58./pc*deg2mas+dec_gaia_err,'*')
labels=['2017','2018', '2019', '2020', '2021','2022']
x_lable=start_time+np.arange(0,5,1)*year_day
plt.xticks(x_lable, labels)
plt.xlabel(r'Year')
plt.ylabel(r'$\Delta \delta$ (mas)')
plt.title('Dec offset')
#plt.legend()
plt.savefig('latex/Q4_2_dec_err'+'.pdf',bbox_inches='tight')
plt.clf()


plt.figure(2)
#plt.xlim([start_time,end_time])
plt.plot(np.asarray(x_ra_earth)+np.asarray(seperation_dis_y)/58./pc*deg2mas+ra_gaia_err,np.asarray(y_dec_earth)+np.asarray(seperation_dis_x)/58./pc*deg2mas+dec_gaia_err,'*')
plt.xlabel(r'$\Delta \alpha$ (mas)')
plt.ylabel(r'$\Delta \delta$ (mas)')
plt.title('above+parallax motion')
plt.savefig('latex/Q4_2_radec_err'+'.pdf',bbox_inches='tight')
plt.clf()


#########  3


plt.figure(2)
#plt.xlim([start_time,end_time])
plt.plot(t_all/day,np.asarray(PM_ra)+np.asarray(x_ra_earth)+np.asarray(seperation_dis_y)/58./pc*deg2mas+ra_gaia_err,'*')
labels=['2017','2018', '2019', '2020', '2021','2022']
x_lable=start_time+np.arange(0,5,1)*year_day
plt.xticks(x_lable, labels)
plt.xlabel(r'Year')
plt.ylabel(r'$\Delta \alpha$ (mas)')
plt.title('RA offset')
#plt.legend()
plt.savefig('latex/Q4_3_ra_err'+'.pdf',bbox_inches='tight')
plt.clf()


plt.figure(2)
#plt.xlim([start_time,end_time])
plt.plot(t_all/day,np.asarray(PM_dec)+np.asarray(y_dec_earth)+np.asarray(seperation_dis_x)/58./pc*deg2mas+dec_gaia_err,'*')
labels=['2017','2018', '2019', '2020', '2021','2022']
x_lable=start_time+np.arange(0,5,1)*year_day
plt.xticks(x_lable, labels)
plt.xlabel(r'Year')
plt.ylabel(r'$\Delta \delta$ (mas)')
plt.title('Dec offset')
#plt.legend()
plt.savefig('latex/Q4_3_dec_err'+'.pdf',bbox_inches='tight')
plt.clf()


plt.figure(2)
#plt.xlim([start_time,end_time])
plt.plot(np.asarray(PM_ra)+np.asarray(x_ra_earth)+np.asarray(seperation_dis_y)/58./pc*deg2mas+ra_gaia_err,np.asarray(PM_dec)+np.asarray(y_dec_earth)+np.asarray(seperation_dis_x)/58./pc*deg2mas+dec_gaia_err,'*')
plt.xlabel(r'$\Delta \alpha$ (mas)')
plt.ylabel(r'$\Delta \delta$ (mas)')
plt.title('above all+proper motion')
plt.savefig('latex/Q4_3_radec_err'+'.pdf',bbox_inches='tight')
plt.clf()


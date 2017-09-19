#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Sat Sep  9 16:17:09 2017

@author: xuduo
"""



import numpy as np
import matplotlib.pyplot as plt


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

M_1=M_sun*1.0
M_2=M_jup*1.0
#M_2=M_sun*1.0

t=3.58e9

a=AU*5.2026
eccentricity=0.048498
I=1.303/180.0*np.pi     ###Inclination	
Omega=100.464/180.0*np.pi       ###Longitude of ascending node
omega=14.75385/180.0*np.pi       ###Argument of perihelion
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
                [-np.sin(Omega)*np.cos(omega)-np.cos(Omega)*np.sin(omega)*np.cos(I)],
                [np.sin(omega)*np.sin(I)]])
    
Q=np.array([[np.cos(Omega)*np.sin(omega)+np.sin(Omega)*np.cos(omega)*np.cos(I)],
                [-np.sin(Omega)*np.sin(omega)+np.cos(Omega)*np.cos(omega)*np.cos(I)],
                [-np.cos(omega)*np.sin(I)]])
    

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


#print r_M1_COM, r_M2_COM,v_M1_COM,v_M2_COM

#print 	np.linalg.norm(r_vector[1:])/AU



print angle_xy(r_vector[0],r_vector[1])
print r_vector









# -*- coding: utf-8 -*-
"""
Created on Tue Jun 16 15:17:07 2015
@author: Jerry Hong and Sathvik Nair
"""

import math
import numpy
#import sys
import matplotlib
matplotlib.use('Agg')
from matplotlib import pyplot as plt
from scipy.interpolate import interp1d as t1d
from scipy.interpolate import spline

plt.rc('text', usetex=True)
font = {'family' : 'serif',
        'weight' : 'bold',
        'size'   : 16}

plt.rc('font', **font)

#Constants
H = 6.626068 * 10**-27 # erg s
K = 1.38066 * 10**-16 # erg/k
c = 3*10**10 #speed of light
b = 2.8977721*10**(-3) #wien's law constant (m*K)

#initial parameters and arrays, step sizes
#Radius is in centimeters
MIN_LOG_NU = 12.5
MAX_LOG_NU = 16
MIN_LOG_R = 12 #Moves SED to right/left (higher frequency)
MAX_LOG_R = 18.94

log_lum_norm = 43.69-14.77 #point on observed SED used to normalize
log_nu_norm = 14.77
log_nu_peak = 15.45 #used to find temperature at which SED peaks at this log_nu, assume temperature is constant with radius.
flux_lum_const = 43.65

dlog_nu = 0.01
dlog_r= 0.01

log_r_arr = numpy.arange(MIN_LOG_R, MAX_LOG_R, dlog_r)
log_nu_arr = numpy.arange(MIN_LOG_NU, MAX_LOG_NU, dlog_nu)
log_lum_arr = []
log_lum_nu_arr = []
lum_arr = []
log_flux_arr = []

emm_arr = []
max_log_r_arr = []

log_r_arr_trad2 = []
log_t_arr_trad2 =[]

#observed SEDs arrays
log_nu_matt = [14.5, 15, 15.8, 16.2, 16.95, 18.3]
log_lum_matt = [45.3, 45.3, 45.75, 45.75, 44.15, 44.6]
log_nu_zrq = [14.3, 14.5, 15, 16.8, 17.7, 18.5]
log_lum_zrq = [45.15, 45.05, 45.4, 44.8, 44.1, 44.1]
log_nu_zrl = [14.3, 14.5, 15, 16.7, 17.7, 18.3]
log_lum_zrl = [45.15, 45.05, 45.4,44.4, 44.2, 44.4]
log_nu_stev = [15.78, 15.46, 14.74]
log_flux_stev = []
log_fluxnu_wink = [-0.15, -0.1, 0.1, 0.3, 0.6, 0.65]
log_nu_wink = [14.6, 14.7, 14.8, 14.9, 15.2, 15.4]
log_flux_wink = map(lambda x,y: x-y, log_fluxnu_wink, log_nu_wink)
log_lum_wink = map(lambda x: x+flux_lum_const, log_flux_wink)


#Tau versus beta versus alpha
beta_arr = [-1.9, -1.8, -1.7, -1.6, -1.5, -1.4, -1.33, -1.3 ,-1.2 ,-1.1,-1,-0.9,-0.8,-0.75,-0.7,-0.65,-0.6,-0.55,-0.5,-0.45,-0.4]
beta_arr = map(lambda x:-x, beta_arr)
alpha_arr = [1.92 ,1.87,1.82,1.75,1.66,1.571,1.5,1.46,1.33333,1.1818,1,0.7777,0.5,0.33,0.14,-0.08,-0.33,-0.64,-1,-1.44,-2]
tau_arr = [0.03,0.06,0.08,0.11,0.14,0.17,0.19,0.2,0.25,0.29,0.35,0.42,0.52,0.58,0.64,0.74,0.86,0.99,1.22,1.5,1.95]
log_tau_arr = map(lambda x:math.log10(x/0.62), tau_arr)

plt.plot(beta_arr, log_tau_arr)
plt.plot((3,-3),(0,0),'k-')
plt.xlabel(r"$\beta$")
plt.ylabel(r"$\log(\tau_{relative})$")
plt.xlim(0.4, 1.9)
#plt.savefig('beta_tau.eps', bbox_inches='tight')
plt.show()

plt.plot(alpha_arr, log_tau_arr)
plt.plot((3,-3),(0,0),'k-')
plt.xlabel(r"$\alpha$")
plt.ylabel(r"$\log(\tau_{relative})$")
plt.xlim(-2,2)
#plt.savefig('alpha_tau.eps', bbox_inches='tight')
plt.show()

#normalization for ngc5548
flux_4392_edel = 1.5 
flux_5468_edel = 1.2
flux_5100_edel = 1.5 + (5100-4392)*(1.2-1.5)/(5468-4392)
lambda_flux_5100_edel = flux_5100_edel*5100*10**-14
edel_const = 4*math.pi*(2.44*10**26)**2
lum_nu_5100_edel = lambda_flux_5100_edel*edel_const
print math.log10(lum_nu_5100_edel)

"""
#Inflection point versus Tau(5100) plot
log_r_inflect_arr = [17.2,16.9,16.6,16.3,16.,15.73,15.48,15.35,15.05,14.98,14.89,14.63,13.96,13.29]
log_tau_arr = [1.32,1.32,1.32,1.32,1.35,1.42,1.17,0.88,0.67,0.66,0.65,0.64,0.63,0.62]
plt.plot(log_r_inflect_arr, log_tau_arr)
plt.show()
"""


#Defining parameters for T_rad5 (temp structure parameters ex. slope)
t_inflect = 4.655 #increasing moves SED to right
inflect = 14.8 #Up and Down
beta_in = 0.55
beta_out = 0.51

"""
#Trad6
alpha_1 = 0.75
r_inflect_1 = 13.8
t_inflect_1 = 7
alpha_2 = 0.75
r_inflect_2 = 55.8
alpha_3 = alpha_2
"""

"""
#use wien's law to find temperature at peak
wavelength_peak = c/10**(log_nu_peak)*10**-2
t_inflect = math.log(b/wavelength_peak,10)
"""

"""
#read t for log_R, log_T for T_rad2
txt = open('temp_structure.txt', 'r')
for columns in ( raw.strip().split() for raw in txt ):  
    log_r_arr_trad2.append(float(columns[0]))
    log_t_arr_trad2.append(float(columns[1]))
txt.close()
"""
#power function
def T_rad(log_r):
    T_grad = -0.58
    t_scale = 10**log_t_scale
    r = 10**log_r
    t = t_scale*(r**(T_grad))
    log_t = math.log(t, 10)
    return log_t
    
#linear interpolation
def T_rad2(log_r):
#   log_T = [6.1, 3.3, 2.1, 5.61, 1.73, 1.86, 1.21, 2.8, 5.2, 1.1, .5]
   i_lo = 0
   i_high = 1
   for i in xrange(len(log_r_arr_trad2)):
       if log_r_arr_trad2[i] <= log_r <= log_r_arr_trad2[i + 1]:
           i_lo = i
           i_high = i_lo+1
   log_t = log_t_arr_trad2[i_lo] + (log_t_arr_trad2[i_high]-log_t_arr_trad2[i_lo])*(log_r-log_t_arr_trad2[i_lo])/(log_r_arr_trad2[i_high]-log_r_arr_trad2[i_lo])
   return log_t

#3 point linear interpolation   
def T_rad3(log_r):
    log_R1= [3, 5, 15, 17]
    log_T1 = [11, 9, 4.4, 3.44]
    #log_T1 = [8.5, 4.67, 4.2, 2.5]
    #log_T1 = [9.2, 5.27, 4.8, 3.1]
    if log_r <= log_R1[1]:
        log_t = log_T1[0] + (log_T1[1]-log_T1[0])*(log_r-log_R1[0])/(log_R1[1]-log_R1[0])
    elif log_R1[1] < log_r < log_R1[2]:
        log_t = log_T1[1] + (log_T1[2]-log_T1[1])*(log_r-log_R1[1])/(log_R1[2]-log_R1[1])
    elif log_r > log_R1[2]:
        log_t = log_T1[2] + (log_T1[3]-log_T1[2])*(log_r-log_R1[2])/(log_R1[3]-log_R1[2])
    return log_t
    
"""
#plot T_rad3
log_t_arr_graph = []    
for log_r in log_r_arr:
    log_t = T_rad3(log_r)
    log_t_arr.append(log_t)
plt.scatter (log_R1, log_T1, color = "red")
plt.scatter(log_r_arr, log_t_arr, color = "blue")
plt.xlabel("log(r)")
plt.ylabel("log(t)")
"""

#Quadratic Interpolation
def T_rad4(log_r):
   log_R = [13.0, 13.5, 14.0, 14.5, 15.0, 15.5, 16.0, 16.5, 17.0]
   log_T = [8, 6.32, 5.95, 3.57, 2.20, 1.82, 0.45, -.07, -.70]
   if log_r == log_R[0]:
       return log_T[0]
   elif log_r == log_R[-1]:
       return log_T[-1]
   else:
       f2 = t1d(log_R, log_T, kind='quadratic')
       log_t = f2(log_r)
       return log_t
       
#Piecewise function with lines of two slopes 
def T_rad5(log_r, inflect, beta_in, beta_out, t_inflect):     
    if log_r > inflect:
        log_T = t_inflect-beta_out*(log_r-inflect)
        return log_T
    elif log_r <= inflect:
        log_T=t_inflect-beta_in*(log_r-inflect)
        return log_T
        
def T_rad6(log_r, alpha_1, r_inflect_1, t_inflect_1, r_inflect_2, alpha_2, alpha_3):
    t_inflect_2 = t_inflect_1-alpha_2*(r_inflect_2-r_inflect_1)
    if log_r > r_inflect_2:
        log_T = t_inflect_2-alpha_3*(log_r-r_inflect_2)
        return log_T
    elif r_inflect_1 <= log_r <= r_inflect_2:
        log_T = t_inflect_1-alpha_2*(log_r-r_inflect_1)
        return log_T
    elif log_r < r_inflect_1:
        log_T=t_inflect_1-alpha_1*(log_r-r_inflect_1)
        return log_T


#Planck Function
def planck(log_r, nu):
    log_temp = T_rad5(log_r,inflect,beta_in,beta_out,t_inflect) 
#    log_temp = T_rad3(log_r)   
    temp = 10**(log_temp)
    return ((4*math.pi*H*nu**3/c**2)/((math.e)**(H*nu/(K*temp))-1)) #Correction for 2pi, Stefan-Boltzmann in sterradians

#Area of Annulus
def area(r):
    return math.pi*float(r**2-(r/math.pow(10, dlog_r))**2)
        
#Find luminosity at a certain log_nu
def lum(log_nu_in):
    sum = 0
    for log_r in log_r_arr:
        r = (10.0)**(log_r)
        nu = 10**log_nu_in        
        rar = area(r)
        pl = planck(log_r, nu)
        sum += rar*pl 
    return sum
    
#Normalization for SED
norm_calc = math.log(lum(log_nu_norm), 10)
r_scale = 10**((log_lum_norm-norm_calc)/2)
log_r_scale = ((log_lum_norm-norm_calc)/2)
LOG_R_MIN = MIN_LOG_R + log_r_scale
LOG_R_MAX = MAX_LOG_R + log_r_scale


def areaf(r):
    return math.pi*float (((r*r_scale)**2)-((r*r_scale/math.pow(10, dlog_r))**2))

def lumf(log_nu_in):
    sum = 0
    for log_r in log_r_arr:
        r = (10.0)**(log_r)
        nu = 10**log_nu_in        
        ar = areaf(r)
        pl = planck(log_r, nu)
        sum += ar*pl 
    return sum      
    
"""
#plot T_rad5
log_T_arr = []        
for log_r in log_r_arr:
    log_T_arr.append(T_rad5(log_r, inflect, beta_in, beta_out, t_inflect))
plt.plot(log_r_arr, log_T_arr, color = "blue")
plt.xlabel("log(r)")
plt.ylabel("log(t)")
plt.show()
"""

#generate log_nu_stev
for log_nu in log_nu_stev:
    if log_nu > 15.46:
         log_flux_stev.append(-14.8 +(log_nu-15.46)*(-1.32))
    elif log_nu < 15.46:
        log_flux_stev.append(-14.8 + (log_nu-15.46) * (-0.74))
    else: 
        log_flux_stev.append(-14.8)
log_lum_stev = map(lambda x: x+flux_lum_const, log_flux_stev)

#Iterate along log_nu, generate SED
print "Normalized at log_nu =", str(log_nu_norm), "=", str("{:6.0f}".format(3*10**10/10**log_nu_norm*10**8)), "angstroms", "log_nu_lum =", str(log_nu_norm + log_lum_norm)
#print "temp = ", str(10**t_inflect_1), "log_temp =", str(t_inflect_1)
print "log_r_scale = ", str("{:6.2}".format(log_r_scale)), "r_scale = ", str("{:6.2}".format(10**log_r_scale))
print "Maximum radius =", str("{:6.2f}".format(LOG_R_MAX)), "=", str("{:6.2f}".format(10**LOG_R_MAX/(2.59*10**15))), "ld"
for log_nu in log_nu_arr:
    #sum is luminosity at that log_nu
    sum = lumf(log_nu)
#    total = flux(log_nu)
    nu = 10**log_nu       
    log_lum_arr.append(math.log(sum,10))
#    log_lum_nu_arr.append(math.log((sum*nu),10))
#    lum_arr.append(nu*sum)
#    log_flux_arr.append(math.log(total, 10))
    #Determine calculad luminosity  
    if abs(log_nu-log_nu_norm) < 0.01:
        print "{:6.2f}".format(log_nu)
        log_nu += dlog_nu
        print "Calculated log_nu_lum: ", str("{:6.2f}".format(log_nu + math.log(sum, 10))), "Observed log_nu_lum:" + str(log_lum_norm+log_nu)
"""
    diff = numpy.diff(log_lum_nu_arr)
    log_lum_max = 0
    log_nu_max = 0       
    for v in diff:       
        if v <= 0.0:
            for x in xrange(len(diff)):
                if v == diff[x]:
                    log_nu_max = log_nu_arr[x]
                    log_lum_max = log_lum_nu_arr[x]
                    print str("{:6.2f}".format(log_nu_max))
            break        
        else:
            continue
#All these pitiful lines for determining t_scale, given t_grad and luminosity at a particular nu
log_t_scale = input("Enter a guess for the log_t_scale:")
log_t_scale_arr = [0.0, 0.0]
delta_arr = [0.0, 0.0]
#Generate two deltas and two log_t_scales to interpolate
for v in range(len(delta_arr)):   
    sum = lum(log_nu_lum)
#    log_lum_arr.append(math.log(sum,10))
    #Determine calculated luminosity
#    print "{:6.2f}".format(log_nu_lum)
#        print "{:6.2f}".format(math.log(sum, 10))
    nulum = log_nu_lum + math.log(sum, 10)
#    nulum_arr.append(nulum)
    print "Calculated Luminosity of random: ", str("{:6.2f}".format(nulum)), "(observed =" + str("{:6.2f}".format(obslum)) +  ")" 
    log_t_scale_arr[v] = log_t_scale
    delta_arr[v] = nulum - obslum
    if delta_arr[v] <= 0:
        log_t_scale += 0.1
    elif delta_arr[v] >= 0:
        log_t_scale -= 0.1  
#Generate first weighted average, which is not quite right...            
log_t_scale = (log_t_scale_arr[0]*(delta_arr[1])+log_t_scale_arr[1]*(-delta_arr[0]))/(-delta_arr[0]+delta_arr[1])
print "Random log_t_scale_arr picked:" + str(log_t_scale_arr)
print "Calculated luminosity minus observed:" + str(delta_arr)
#print log_t_scale
print "Rough scale factor:" + str(10**log_t_scale)
 #Do it all again, but replace the element of the array log_t_scale_arr that is farther from log_t_scale
if abs(log_t_scale_arr[0]-log_t_scale) < abs(log_t_scale_arr[1]-log_t_scale):
    index_replace = 1    
else:
    index_replace = 0  
log_t_scale_arr[index_replace] = log_t_scale
sum = lum(log_nu_lum)
#    log_lum_arr.append(math.log(sum,10))
    #Determine calculated luminosity
#    print "{:6.2f}".format(log_nu_lum)
nulum = log_nu_lum + math.log(sum, 10)
#    nulum_arr.append(nulum)
print "Calculated Luminosity of rough: ", str("{:6.2f}".format(nulum)), "(observed =" + str("{:6.2f}".format(obslum)) +  ")" 
delta_arr[index_replace] = nulum - obslum
#Determine final scale factor:
log_t_scale = (log_t_scale_arr[0]*(delta_arr[1])+log_t_scale_arr[1]*(-delta_arr[0]))/(-delta_arr[0]+delta_arr[1])
print "Final scale factor: " + str(10**log_t_scale)
print "Final scale factor (logged): " + str(log_t_scale)
#plotting (log L) versus (log V), i.e. SED
plt.plot(log_nu_arr, log_lum_arr)
plt.xlabel('log(nu)')
plt.ylabel('log(lum)')
plt.show
#Plotting Log vL versus (log V)
fig1,ax1=plt.subplots(1,1)
plt.plot(log_nu_arr, log_lum_nu_arr)
plt.xlabel('log(nu)')
plt.ylim(44,46)
plt.ylabel('log(nu*lum)')
plt.plot(log_nu_matt, log_lum_matt)
plt.show()
#plt.savefig('SED_temp5_mathews.pdf', bbox_inches='tight')
 
#Plotting vL versus (log V)
fig1,ax1=plt.subplots(1,1)
plt.plot(log_nu_arr, lum_arr)
plt.plot(log_nu_stev, lum_stev)
plt.xlabel('log(nu)')
plt.ylabel('Nu*L')
plt.show


#Plotting lum versus log_nu
fig1, ax1 = plt.subplots (1,1)
plt.plot (log_nu_arr, log_lum_arr)
plt.plot(log_nu_stev, log_lum_stev)
plt.plot(log_nu_wink, log_lum_wink)
plt.xlabel("Log(Nu)")
plt.ylabel("Log(Lum)")
#plt.xlim(14.5,15.8)
#plt.ylim(28.6,29.0)
plt.show



"""
#Emmissivity Graphs
log_r_femm = numpy.arange(LOG_R_MIN,LOG_R_MAX, dlog_r)
log_r_emm = map(lambda x: x-log_r_scale, log_r_femm)
"""
plt.plot(log_r_femm, log_T_arr, color = "red")
plt.xlabel("log(r)")
plt.ylabel("log(t)")
plt.ylim(0,20)
plt.show()
"""    

wavelength = [800, 1367, 1928, 2246, 2600, 3000, 3465, 4000, 4392,5100, 5468, 6000]
colors = ['r', 'g', 'c', 'm', 'y', 'b', 'p', 'b', 'y', 'r', 'g', 'm']
#fig1,ax2=plt.subplots(1,1)
fig1, ax1 = plt.subplots()
for num in xrange(len(wavelength)): 
    emm_arr = []
    log_emm_arr = []
    nu = c/(wavelength[num]*10**-8)
    for log_r in log_r_emm:
        r = (10.0)**(log_r)
        ar = areaf(r)
        pl = planck(log_r, nu)
        # construct planck function at that nu
        emm_arr.append(ar*pl)
#        log_emm_arr.append(math.log10(ar*pl))
#finding the where the functions reaches maximum
    diff = numpy.diff(emm_arr)
    ar_pl_max = 0
    log_r_max = 0       
    for v in diff:       
        if v <= 0.0:
            for x in xrange(len(diff)):
                if v == diff[x]:
                    log_r_max = log_r_emm[x]
                    ar_pl_max = emm_arr[x]
                    max_log_r_LD = (10**(log_r_max+log_r_scale))/(2.59*10**15)
                    max_log_r_arr.append(max_log_r_LD)
                    print "Log R(" + str(wavelength[num]) + ")" "= " + str("{:6.2f}".format(log_r_max+log_r_scale)) + " =" + str("{:6.2f}".format(max_log_r_LD)) + " ld" + str("{:6.2f}".format((max_log_r_LD)-(max_log_r_arr[0])))
            break        
        else:
            continue
    #Plotting E(log R)
    emm_arr = map(lambda x: x*0.88/10**27, emm_arr) #small scale factor adjustment
    ax1.plot(log_r_femm, emm_arr, colors[num]+'-')
    ax1.set_xlabel(r"$\log(R)$")
    ax1.set_ylabel(r"Emmissivity")

log_T_arr_575 = []
for log_r in log_r_femm:
    log_T_arr_575.append(T_rad5(log_r-log_r_scale,inflect,beta_in,beta_out,t_inflect))
ax2 = ax1.twinx()
ax2.plot(log_r_femm, log_T_arr_575, "k-")

plt.xlim(12,17)
plt.show()

#data from edelson et al. 2015
lambda_5548 = [4.4, 25.3, 1367, 1928, 2246, 2600, 3465, 4392, 5468]
tau_5548 = [-0.66, 0.08, 0, 0.40, 0.35, 0.61, 1.35, 1.23, 1.56]
tau_err_5548 = [0.46, 0.52, 0.25, 0.17, 0.16, 0.20, 0.24, 0.29, 0.50]
plt.errorbar(lambda_5548, tau_5548, yerr=tau_err_5548, fmt = 'o', color = 'k')
#plt.plot(wavelength, max_log_r_arr, 'k')
plt.ylabel(r"Lag (days)")
plt.xlabel(r"Wavelength (\r{A})")
x_smooth = numpy.linspace(wavelength[0],wavelength[-1], 200)
y_smooth = spline(wavelength, max_log_r_arr, x_smooth)
plt.plot(x_smooth, y_smooth, 'k')
#plt.show()

beta_out = 0.75
beta_in = 0.75
norm_calc = math.log(lum(log_nu_norm), 10)
r_scale = 10**((log_lum_norm-norm_calc)/2)
log_r_scale = ((log_lum_norm-norm_calc)/2)
LOG_R_MIN = MIN_LOG_R + log_r_scale
LOG_R_MAX = MAX_LOG_R + log_r_scale


max_log_r_arr = []
for num in xrange(len(wavelength)): 
    emm_arr = []
    log_emm_arr = []
    nu = c/(wavelength[num]*10**-8)
    for log_r in log_r_emm:
        r = (10.0)**(log_r)
        ar = areaf(r)
        pl = planck(log_r, nu)
        # construct planck function at that nu
        emm_arr.append(ar*pl)
#        log_emm_arr.append(math.log10(ar*pl))
#finding the where the functions reaches maximum
    diff = numpy.diff(emm_arr)
    ar_pl_max = 0
    log_r_max = 0       
    for v in diff:       
        if v <= 0.0:
            for x in xrange(len(diff)):
                if v == diff[x]:
                    log_r_max = log_r_emm[x]
                    ar_pl_max = emm_arr[x]
                    max_log_r_LD = (10**(log_r_max+log_r_scale))/(2.59*10**15)
                    max_log_r_arr.append(max_log_r_LD)
                    print "Log R(" + str(wavelength[num]) + ")" "= " + str("{:6.2f}".format(log_r_max+log_r_scale)) + " =" + str("{:6.2f}".format(max_log_r_LD)) + " ld" + str("{:6.2f}".format((max_log_r_LD)-(max_log_r_arr[0])))
            break        
        else:
            continue
    #Plotting E(log R)
    emm_arr = map(lambda x: x*0.88/10**27, emm_arr) #small scale factor adjustment
    ax1.plot(log_r_femm, emm_arr, colors[num]+'-')
    ax1.set_xlabel(r"$\log(R)$")
    ax1.set_ylabel(r"Emmissivity")

x_smooth = numpy.linspace(wavelength[0],wavelength[-1], 200)
y_smooth = spline(wavelength, max_log_r_arr, x_smooth)
plt.plot(x_smooth, y_smooth, 'k--')
#plt.show()










#Augment luminostiy and see if there's a better result

log_lum_norm += 0.25

beta_in = 0.55
beta_out = 0.55

norm_calc = math.log(lum(log_nu_norm), 10)
r_scale = 10**((log_lum_norm-norm_calc)/2)
log_r_scale = ((log_lum_norm-norm_calc)/2)
LOG_R_MIN = MIN_LOG_R + log_r_scale
LOG_R_MAX = MAX_LOG_R + log_r_scale

log_r_femm = numpy.arange(LOG_R_MIN,LOG_R_MAX, dlog_r)
log_r_emm = map(lambda x: x-log_r_scale, log_r_femm)

max_log_r_arr = []
for num in xrange(len(wavelength)): 
    emm_arr = []
    log_emm_arr = []
    nu = c/(wavelength[num]*10**-8)
    for log_r in log_r_emm:
        r = (10.0)**(log_r)
        ar = areaf(r)
        pl = planck(log_r, nu)
        # construct planck function at that nu
        emm_arr.append(ar*pl)
#        log_emm_arr.append(math.log10(ar*pl))
#finding the where the functions reaches maximum
    diff = numpy.diff(emm_arr)
    ar_pl_max = 0
    log_r_max = 0     
    for v in diff:       
        if v <= 0.0:
            for x in xrange(len(diff)):
                if v == diff[x]:
                    log_r_max = log_r_emm[x]
                    ar_pl_max = emm_arr[x]
                    max_log_r_LD = (10**(log_r_max+log_r_scale))/(2.59*10**15)
                    max_log_r_arr.append(max_log_r_LD)
                    print "Log R(" + str(wavelength[num]) + ")" "= " + str("{:6.2f}".format(log_r_max+log_r_scale)) + " =" + str("{:6.2f}".format(max_log_r_LD)) + " ld" + str("{:6.2f}".format((max_log_r_LD)-(max_log_r_arr[0])))
            break        
        else:
            continue
    #Plotting E(log R)
#    emm_arr = map(lambda x: x*0.88/10**27, emm_arr) #small scale factor adjustment
#    ax1.plot(log_r_femm, emm_arr, colors[num]+'-')
#    ax1.set_xlabel(r"$\log(R)$")
#    ax1.set_ylabel(r"Emmissivity")

log_T_arr_575 = []
for log_r in log_r_femm:
    log_T_arr_575.append(T_rad5(log_r-log_r_scale,inflect,beta_in,beta_out,t_inflect))
ax2 = ax1.twinx()
ax2.plot(log_r_femm, log_T_arr_575, "k-")

#plt.xlim(12,17)
#plt.show()

x_smooth = numpy.linspace(wavelength[0],wavelength[-1], 200)
y_smooth = spline(wavelength, max_log_r_arr, x_smooth)
plt.plot(x_smooth, y_smooth, 'r')


beta_in = 0.75
beta_out = 0.75
norm_calc = math.log(lum(log_nu_norm), 10)
r_scale = 10**((log_lum_norm-norm_calc)/2)
log_r_scale = ((log_lum_norm-norm_calc)/2)
LOG_R_MIN = MIN_LOG_R + log_r_scale
LOG_R_MAX = MAX_LOG_R + log_r_scale

log_r_femm = numpy.arange(LOG_R_MIN,LOG_R_MAX, dlog_r)
log_r_emm = map(lambda x: x-log_r_scale, log_r_femm)


max_log_r_arr = []
for num in xrange(len(wavelength)): 
    emm_arr = []
    log_emm_arr = []
    nu = c/(wavelength[num]*10**-8)
    for log_r in log_r_emm:
        r = (10.0)**(log_r)
        ar = areaf(r)
        pl = planck(log_r, nu)
        # construct planck function at that nu
        emm_arr.append(ar*pl)
#        log_emm_arr.append(math.log10(ar*pl))
#finding the where the functions reaches maximum
    diff = numpy.diff(emm_arr)
    ar_pl_max = 0
    log_r_max = 0       
    for v in diff:       
        if v <= 0.0:
            for x in xrange(len(diff)):
                if v == diff[x]:
                    log_r_max = log_r_emm[x]
                    ar_pl_max = emm_arr[x]
                    max_log_r_LD = (10**(log_r_max+log_r_scale))/(2.59*10**15)
                    max_log_r_arr.append(max_log_r_LD)
                    print "Log R(" + str(wavelength[num]) + ")" "= " + str("{:6.2f}".format(log_r_max+log_r_scale)) + " =" + str("{:6.2f}".format(max_log_r_LD)) + " ld" + str("{:6.2f}".format((max_log_r_LD)-(max_log_r_arr[0])))
            break        
        else:
            continue
    #Plotting E(log R)
#    emm_arr = map(lambda x: x*0.88/10**27, emm_arr) #small scale factor adjustment
#    ax1.plot(log_r_femm, emm_arr, colors[num]+'-')
#    ax1.set_xlabel(r"$\log(R)$")
#    ax1.set_ylabel(r"Emmissivity")

x_smooth = numpy.linspace(wavelength[0],wavelength[-1], 200)
y_smooth = spline(wavelength, max_log_r_arr, x_smooth)
plt.plot(x_smooth, y_smooth, 'r--')

plt.savefig('timedelay_5548_fig.eps', bbox_inches='tight')

log_nu_slope = []
log_lum_slope = []
for x in range(len(log_nu_arr)):
#    print "(",str(log_nu_arr[x]), ",", str(log_lum_arr[x]), ")"
    if abs(log_nu_arr[x]-14.77) <= 0.001 or abs(log_nu_arr[x]-14.78) <= 0.001:
        log_nu_slope.append(log_nu_arr[x])
        log_lum_slope.append(log_lum_arr[x])
print "Slope of SED at 5100 = ", str((log_lum_slope[1]-log_lum_slope[0])/(log_nu_slope[1]-log_nu_slope[0]))

log_nu_slope = []
log_lum_slope = []
for x in range(len(log_nu_arr)):
#    print "(",str(log_nu_arr[x]), ",", str(log_lum_arr[x]), ")"
    if abs(log_nu_arr[x]-15.36) <= 0.001 or abs(log_nu_arr[x]-15.37) <= 0.001:
        log_nu_slope.append(log_nu_arr[x])
        log_lum_slope.append(log_lum_arr[x])
print "Slope of SED at 1310 = ", str((log_lum_slope[1]-log_lum_slope[0])/(log_nu_slope[1]-log_nu_slope[0]))


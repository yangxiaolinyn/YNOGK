import os
import matplotlib.pyplot as plt
from matplotlib import font_manager as fm, rcParams 
import numpy as np 
#import myplot 
import openpyxl
from openpyxl.styles import Font  
import math

font2 = {'family': 'STIXGeneral',
         'style': 'normal',
         'weight': 'normal', 
        'size': 16,
        }


#~~~~~ Set the observer's viewing angle ~~~~~~~~~~~~~~~~~~~~~~~~

psi = 45.0
the = 20.0
phi = 110.0

lws = 0.15


file1 = 'sphmotion.txt'
  


fig, ax = plt.subplots(figsize=(8, 8))
  
#~~~~~~~~~~~~~~~~~Readin figure data~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ 
data = np.loadtxt(file1, usecols = (0, 1, 2))
xx = data[::, 0]
yy = data[::, 1]
zz = data[::, 2]

Pi = 3.141592653589793

psi1 = psi / 180.*math.pi
theta1 = the / 180.*math.pi
phi1 = phi / 180.*math.pi
 
#~~~~~~~~~~~~~~~~~Rotation the coordinate and Draw the curve ~~~~~~~~~~~~~~~~ 

x0 = xx * math.cos(psi1) - yy * math.sin(psi1)
y0 = xx * math.sin(psi1) + yy * math.cos(psi1)
z0 = zz

z1 = z0 * math.cos(theta1) - x0 * math.sin(theta1)
x1 = z0 * math.sin(theta1) + x0 * math.cos(theta1)
y1 = y0

x2 = x1 * math.cos(phi1) - y1 * math.sin(phi1)
y2 = x1 * math.sin(phi1) + y1 * math.cos(phi1)
z2 = z1

xp = np.array(x2)
y3 = np.array(y2)
zp = np.array(z2)
 

posi = np.where( y3 >= 0.0)
neg = np.where( y3 < 0.0)


xn = np.array(x2)
zn = np.array(z2)
 

for i in posi:
    xn[i] = float('nan')
    zn[i] = float('nan')
#xn[posi] = 0#float('nan')
#zn[posi] = 0#float('nan')
 

ax.plot(xn, zn, ':', c = 'black', linewidth = 0.5 )

#print('xn = ', xn)

#print('zn = ', zn)

 
 
for i in neg:
    xp[i] = float('nan')
    zp[i] = float('nan')
 
ax.plot(xp, zp, '-', c = 'black', linewidth = 1 )

#print('xp = ', xp)

#print('zp = ', zp)

#plt.style.use('seaborn')
#plt.style.use('fivethirtyeight')
#plt.style.use('seaborn-ticks')
#plt.style.use('seaborn-darkgrid')
#plt.style.use('seaborn-dark')


#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~


#~~~~~~~~~~~~~~~~~~Set the axises ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#ax.set_title("figures", fontsize = 24)
ax.set_xlabel("y", fontdict=font2 )
ax.set_ylabel("x", fontdict=font2 )


ax.tick_params(axis = 'both', labelsize = 10)
 
#--snip--

sz1 = 2.0
ax.axis([-sz1, sz1, -sz1, sz1] )
 

rg = 1.403
a = 0.980
rh = math.sqrt(rg*rg + a*a)

num = 1000
xsp = np.arange(num, dtype = np.float32)
ysp = np.arange(num, dtype = np.float32)
for i in range(1000):
   xsp[i] = rh * math.cos(i * 2.*math.pi/1000.0)
   ysp[i] = rh * math.sin(i * 2.*math.pi/1000.0)

ax.plot(xsp, ysp, '-', c = 'black', linewidth = 0.5 )
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~


#~~~~~~~~~~~~draw the equatorial circle~~~~~~~~~~~~~~~~~~~~~~~~~~~
num = 500
xep = np.arange(num, dtype = np.float32)
yep = np.arange(num, dtype = np.float32)
zep = np.arange(num, dtype = np.float32)
for i in range(500):
   xep[i] = rh * math.cos(i * 2.*math.pi/500.0)
   yep[i] = rh * math.sin(i * 2.*math.pi/500.0)
   zep[i] = 0.0

#ax.plot(xep, yep, '-', c = 'red', linewidth = 0.5 )

x0 = xep * math.cos(psi1) - yep * math.sin(psi1)
y0 = xep * math.sin(psi1) + yep * math.cos(psi1)
z0 = zep

z1 = z0 * math.cos(theta1) - x0 * math.sin(theta1)
x1 = z0 * math.sin(theta1) + x0 * math.cos(theta1)
y1 = y0

x2 = x1 * math.cos(phi1) - y1 * math.sin(phi1)
y2 = x1 * math.sin(phi1) + y1 * math.cos(phi1)
z2 = z1

xp = np.array(x2)
y3 = np.array(y2)
zp = np.array(z2)
 

posi = np.where( y3 >= 0.0)
neg = np.where( y3 < 0.0)
 
 
for i in neg:
    xp[i] = float('nan')
    zp[i] = float('nan')
 
ax.plot(xp, zp, '-', c = 'blue', linewidth = 0.5 )

#print('xp = ', xp)

#print('zp = ', zp)

xn = np.array(x2)
zn = np.array(z2)
 

for i in posi:
    xn[i] = float('nan')
    zn[i] = float('nan')
#xn[posi] = 0#float('nan')
#zn[posi] = 0#float('nan')
 

ax.plot(xn, zn, ':', c = 'blue', linewidth = 0.5 )
#ax.text(xp[300], zp[300], r"equatorial plane", fontsize = 10, fontdict=font2)
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

#plt.show()

fig.savefig('fig5.eps', dpi=120, bbox_inches = 'tight')




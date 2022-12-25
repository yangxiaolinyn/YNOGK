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
the = 40.0
phi = 110.0


file1 = 'rayxx.txt'
file2 = 'rayyy.txt'
  


fig, ax = plt.subplots(figsize=(8, 8))
  
#~~~~~~~~~~~~~~~~~Readin figure data~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ 
xx = np.loadtxt(file1, usecols = (0, 1)) 

Num = 41

tmpx = np.array(xx[::, 0])
N_data = int(len(tmpx)/(Num))
Xdata = tmpx.reshape(Num, N_data )


tmpy = np.array(xx[::, 1])
Ydata = tmpy.reshape(Num, N_data )

#~~~~~~~~~~~~~~~~~Readin figure data~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ 
yy = np.loadtxt(file2, usecols = (0, 1)) 
 

tmpx = np.array(yy[::, 0])
N_data = int(len(tmpx)/(Num))
Xdata2 = tmpx.reshape(Num, N_data )


tmpy = np.array(yy[::, 1])
Ydata2 = tmpy.reshape(Num, N_data )




Pi = 3.141592653589793

psi1 = psi / 180.*math.pi
theta1 = the / 180.*math.pi
phi1 = phi / 180.*math.pi
 
#~~~~~~~~~~~~~~~~~Rotation the coordinate and Draw the curve ~~~~~~~~~~~~~~~~ 
  
N2 = Num - 1
for i in range(Num): 
    ax.plot( -Ydata[i, ::], -Xdata[i, ::], '-', c = 'black', linewidth = 0.6 )
    ax.plot( -Ydata2[i, ::], -Xdata2[i, ::], '-', c = 'black', linewidth = 0.6 )
    #ax.plot( -Ydata[::, i], -Xdata[::, i], '-', c = 'black', linewidth = 0.6 )
 
 


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

ssz = 10
ax.axis([-ssz, ssz, -ssz, ssz] )
 
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
 

#plt.show() 

fig.savefig('fig11.eps', dpi=120, bbox_inches = 'tight')




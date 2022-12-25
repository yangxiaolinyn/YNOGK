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


file1 = 'rays.txt'
  


fig, ax = plt.subplots(figsize=(8, 8))
  
#~~~~~~~~~~~~~~~~~Readin figure data~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ 
xx = np.loadtxt(file1, usecols = (0, 1, 2)) 

Num = 21*21

tmpx = np.array(xx[::, 0])
N_data = int(len(tmpx)/(Num))
Xdata = tmpx.reshape(Num, N_data )


tmpy = np.array(xx[::, 1])
Ydata = tmpy.reshape(Num, N_data )

tmpz= np.array(xx[::, 2])
Zdata = tmpz.reshape(Num, N_data )
 
 


Pi = 3.141592653589793

psi1 = psi / 180.*math.pi
theta1 = the / 180.*math.pi
phi1 = phi / 180.*math.pi
 
#~~~~~~~~~~~~~~~~~Rotation the coordinate and Draw the curve ~~~~~~~~~~~~~~~~ 
  
N2 = Num - 1
for i in range(Num): 
    ax.plot( Xdata[i, ::], Ydata[i, ::], '-', c = 'black', linewidth = 0.6 )
    #ax.plot( -Ydata2[i, ::], -Xdata2[i, ::], '-', c = 'black', linewidth = 0.6 )
    #ax.plot( -Ydata[::, i], -Xdata[::, i], '-', c = 'black', linewidth = 0.6 )
 
 


#plt.style.use('seaborn')
#plt.style.use('fivethirtyeight')
#plt.style.use('seaborn-ticks')
#plt.style.use('seaborn-darkgrid')
#plt.style.use('seaborn-dark')


#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~


#~~~~~~~~~~~~~~~~~~Set the axises ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#ax.set_title("figures", fontsize = 24)
ax.set_xlabel("X [$r_g$]", fontdict=font2 )
ax.set_ylabel("Y [$r_g$]", fontdict=font2 )


ax.tick_params(axis = 'both', labelsize = 10)
 
#--snip--

ssz = 5
ax.axis([-ssz, ssz, -ssz, ssz] )
 
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
 

#plt.show() 

fig.savefig('fig14.eps', dpi=120, bbox_inches = 'tight')




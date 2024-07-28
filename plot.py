import pandas as pd
import numpy as np
#Data is loaded
rfilename = 'Rhistory.dat'
xfilename = 'Xhistory.dat'
yfilename = 'Yhistory.dat'
ffilename = 'Fhistory.dat'

rdata = open(rfilename)
xdata = open(xfilename)
ydata = open(yfilename)
fdata = open(ffilename)
rdata_raw = rdata.readline()
xdata_raw = xdata.readlines()
ydata_raw = ydata.readlines()
fdata_raw = fdata.readlines()
iterations = len(xdata_raw)
rdata_raw = list(map(float,(rdata_raw.split(','))[:-1]))
for i in range(iterations):
    xdata_raw[i] = list(map(float,(xdata_raw[i].split(','))[:-1]))
    ydata_raw[i] = list(map(float,(ydata_raw[i].split(','))[:-1]))
    fdata_raw[i] = list(map(float,(fdata_raw[i].split(','))[:-1]))
rdata.close()
xdata.close()
ydata.close()
fdata.close()

#Animation in paralel of scalar fields
import matplotlib.animation as animation
import matplotlib.pyplot as plt

R = rdata_raw
IPF = 5 #Iterations per frame
FRAMES = int(iterations/IPF)

fig, ax = plt.subplots(3,1,figsize=(12,12),dpi=100)

Fplot = ax[0].plot(R,fdata_raw[0])[0]
Xplot = ax[1].plot(R,xdata_raw[0])[0]
Yplot = ax[2].plot(R,ydata_raw[0])[0]
ax[1].set_ylim(-ax[0].get_ylim()[1],ax[0].get_ylim()[1])
ax[0].set_ylim(ax[0].get_ylim())
ax[2].set_ylim(ax[0].get_ylim())
ax[0].set_title('$\phi$')
ax[1].set_title('$\Phi$')
ax[2].set_title('$\Pi$')

#Create a function that changes the data every frame
def animate(i):
    it = i*IPF
    Fplot.set_data(R,fdata_raw[it])
    Xplot.set_data(R,xdata_raw[it])
    Yplot.set_data(R,ydata_raw[it])
    return Xplot,Yplot

#Make the animation
anim = animation.FuncAnimation(fig, animate, frames=FRAMES,
                                interval=100, repeat_delay=3000)

anim.save('newtest.gif',writer='pillow')


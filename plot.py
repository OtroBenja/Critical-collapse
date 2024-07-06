import pandas as pd
import numpy as np
#Se importan los datos
rfilename = 'Rhistory.dat'
xfilename = 'Xhistory.dat'
yfilename = 'Yhistory.dat'

rdata = open(rfilename)
xdata = open(xfilename)
ydata = open(yfilename)
rdata_raw = rdata.readline()
xdata_raw = xdata.readlines()
ydata_raw = ydata.readlines()
iterations = len(xdata_raw)
rdata_raw = list(map(float,(rdata_raw.split(','))[:-1]))
for i in range(iterations):
    xdata_raw[i] = list(map(float,(xdata_raw[i].split(','))[:-1]))
    ydata_raw[i] = list(map(float,(ydata_raw[i].split(','))[:-1]))
rdata.close()
xdata.close()
ydata.close()

#Animación de X e Y en paralelo
import matplotlib.animation as animation
import matplotlib.pyplot as plt

R = rdata_raw
IPF = 5 #Iterations per frame
FRAMES = int(iterations/IPF)

fig, ax = plt.subplots(1,2,figsize=(12,6),dpi=100)

Xplot = ax[0].plot(R,xdata_raw[0])[0]
Yplot = ax[1].plot(R,ydata_raw[0])[0]
ax[0].set_ylim(-ax[0].get_ylim()[1],ax[0].get_ylim()[1])
ax[1].set_ylim(ax[0].get_ylim())

#Se crea una función que cambia los datos necesarios en cada frame
def animate(i):
    it = i*IPF
    Xplot.set_data(R,xdata_raw[it])
    Yplot.set_data(R,ydata_raw[it])
    return Xplot,Yplot

  #Se crea la animación
anim = animation.FuncAnimation(fig, animate, frames=FRAMES,
                                interval=100, repeat_delay=3000)

anim.save('newtest.gif',writer='pillow')


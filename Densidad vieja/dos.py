#instalar paquetes
#para instalar paquete nuevo abrir "cmd" o "simbolo del sistema"
#pip3 install nombre del paquete
from matplotlib import pyplot as plt
#from scipy.signal import savgol_filter

import math
import numpy as np
#import scipy.signal




#---------------------
fig,ax=plt.subplots() #esto hace algo no se que
#fig = plt.figure()
#---------------------
infile=open("DOSc0g0001.dat","r")
x=[]
y=[]
for line in infile:
	words=line.split()
	x.append(float(words[0])) #columna eje x (python empieza por 0)
	y.append(float(words[1])) #columna eje y (python empieza por 0)
infile.close()
x=np.array(x)
y=np.array(y)
#yhat = scipy.signal.savgol_filter(y, 9, 2) # window size 51, polynomial order 3

plt.plot(x,y,color="k",markersize=1,linestyle="solid",linewidth=1,label="$\Gamma$=0.001")
plt.plot(x*-1,y1,color="k",markersize=1,linestyle="solid",linewidth=2,label="")
#plt.fill_between(x,y1, color= "grey")

#---------------------
infile=open("DOSc0g0005.dat","r") #abrir tu archivo de datos "nombre"
x=[]
y=[]
for line in infile:
	words=line.split()
	x.append(float(words[0])) #columna eje x (python empieza por 0)
	y.append(float(words[1])) #columna eje y (python empieza por 0)
infile.close()
x=np.array(x)
y=np.array(y)
#yhat = scipy.signal.savgol_filter(y, 9, 2) # window size 51, polynomial order 3

plt.plot(x,y,color="m",markersize=1,linestyle="dotted",linewidth=1,label="$\Gamma$=0.005")
plt.plot(x*-1,y2,color="k",markersize=1,linestyle="dotted",linewidth=2,label="")
#---------------------	
#---------------------
infile=open("DOSc0g0010.dat","r") #abrir tu archivo de datos "nombre"
x=[]
y=[]
for line in infile:
	words=line.split()
	x.append(float(words[0])) #columna eje x (python empieza por 0)
	y.append(float(words[1])) #columna eje y (python empieza por 0)
infile.close()
x=np.array(x)
y=np.array(y)
#yhat = scipy.signal.savgol_filter(y, 9, 2) # window size 51, polynomial order 3
plt.plot(x,y,color="g",markersize=1,linestyle="dashdot",linewidth=1,label="$\Gamma$=0.010")
plt.plot(x*-1,y1,color="k",markersize=1,linestyle="solid",linewidth=2,label="")
#plt.fill_between(x,y1, color= "grey")

#---------------------
infile=open("DOSc0g0015.dat","r") #abrir tu archivo de datos "nombre"
x=[]
y=[]
for line in infile:
	words=line.split()
	x.append(float(words[0])) #columna eje x (python empieza por 0)
	y.append(float(words[1])) #columna eje y (python empieza por 0)
infile.close()
x=np.array(x)
y=np.array(y)
#yhat = scipy.signal.savgol_filter(y, 9, 2) # window size 51, polynomial order 3
plt.plot(x,y,color="b",markersize=1,linestyle="dashed",linewidth=1,label="$\Gamma$=0.015")
plt.plot(x*-1,y2,color="k",markersize=1,linestyle="dotted",linewidth=2,label="")
#---------------------	

infile=open("DOSc0g0020.dat","r") #abrir tu archivo de datos "nombre"
x=[]
y=[]
for line in infile:
	words=line.split()
	x.append(float(words[0])) #columna eje x (python empieza por 0)
	y.append(float(words[1])) #columna eje y (python empieza por 0)
infile.close()
x=np.array(x)
y=np.array(y)
#yhat = scipy.signal.savgol_filter(y, 9, 3) # window size 51, polynomial order 3
plt.plot(x,y,color="r",markersize=1,linestyle="dotted",linewidth=2,label="$\Gamma$=0.020")
plt.plot(x*-1,y5,color="k",markersize=1,linestyle="solid",linewidth=1,label="")
#plt.fill_between(x*-1,y5, color= "gainsboro")"""
#---------------------	

#---------------------
plt.tight_layout()
plt.legend(loc="upper right", prop={'size': 10})
#ax.text(31.28, 2.196, "$\Gamma^{+}$=0,15 meV",fontsize=10)
plt.style.use("classic")
plt.ylabel("Im[$\widetilde{\omega}$($\omega$+ $0^{+}$)](meV)",fontsize=12)
plt.xlabel("$\omega$(meV)",fontsize=12)
ax.tick_params(direction='in', length=6, width=1.0, grid_alpha=0.5)
plt.xlim(0,2)
plt.ylim(0,3)
plt.tight_layout()
#---------------------

plt.savefig("dos.pdf")
plt.show()
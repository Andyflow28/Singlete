import numpy as np
import matplotlib.pyplot as plt

file1 = np.loadtxt('DOSc0g0001.dat')
file2 = np.loadtxt('DOSc0g0005.dat')
file3 = np.loadtxt('DOSc0g0010.dat')
file4 = np.loadtxt('DOSc0g0015.dat')
file5 = np.loadtxt('DOSc0g0020.dat')

def plot_everything():
    fig,ax=plt.subplots()
    
    x,y = file1.T
    plt.plot(x,y,color="k",markersize=1,linestyle="solid",linewidth=1,label="$\Gamma$=0.001")
    plt.plot(x*(-1),y,color="k",markersize=1,linestyle="solid",linewidth=1,label="$\Gamma$=0.001")
#     plt.plot(x*-1,y1,color="k",markersize=1,linestyle="solid",linewidth=2,label="")
    
    x,y = file2.T
    plt.plot(x,y,color="m",markersize=1,linestyle="dotted",linewidth=1,label="$\Gamma$=0.005")
#     plt.plot(x*-1,y2,color="k",markersize=1,linestyle="dotted",linewidth=2,label="")
    
    x,y = file3.T
    plt.plot(x,y,color="g",markersize=1,linestyle="dashdot",linewidth=1,label="$\Gamma$=0.010")
#     plt.plot(x*-1,y1,color="k",markersize=1,linestyle="solid",linewidth=2,label="")
    
    x,y = file4.T
    plt.plot(x,y,color="b",markersize=1,linestyle="dashed",linewidth=1,label="$\Gamma$=0.015")
#     plt.plot(x*-1,y2,color="k",markersize=1,linestyle="dotted",linewidth=2,label="")
    
    x,y = file5.T
    plt.plot(x,y,color="r",markersize=1,linestyle="dotted",linewidth=2,label="$\Gamma$=0.020")
#     plt.plot(x*-1,y5,color="k",markersize=1,linestyle="solid",linewidth=1,label="")
    
    plt.legend(loc="upper right", prop={'size': 10})
    #ax.text(31.28, 2.196, "$\Gamma^{+}$=0,15 meV",fontsize=10)
    plt.style.use("classic")
    plt.ylabel("Im[$\widetilde{\omega}$($\omega$+ $0^{+}$)](meV)",fontsize=12)
    plt.xlabel("$\omega$(meV)",fontsize=12)
    ax.tick_params(direction='in', length=6, width=1.0, grid_alpha=0.5)
    plt.xlim(0,2)
    plt.ylim(0,3)
    plt.tight_layout()
    plt.savefig("dos.pdf")
    plt.show()


plot_everything()


file1[:,0] = file1[:,0]*(-1)
file2[:,0] = file2[:,0]*(-1)
file3[:,0] = file3[:,0]*(-1)
file4[:,0] = file4[:,0]*(-1)
file5[:,0] = file5[:,0]*(-1)

np.savetxt('DOSc0g0001_x_negativo.dat',file1)
np.savetxt('DOSc0g0005_x_negativo.dat',file2)
np.savetxt('DOSc0g0010_x_negativo.dat',file3)
np.savetxt('DOSc0g0015_x_negativo.dat',file4)
np.savetxt('DOSc0g0020_x_negativo.dat',file5)
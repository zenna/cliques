import matplotlib.pyplot as plt
import numpy as np
import mpl_toolkits.mplot3d.axes3d as p3
import matplotlib.colors
from ipdb import set_trace

space = np.loadtxt('coords.mat')
colours = np.loadtxt('stabs.mat')
normalise = matplotlib.colors.Normalize()
normalise.autoscale(colours)
cmap = matplotlib.colors.Colormap('Spectral')


if space.shape[1] == 2:
    plt.scatter(space[:,0], space[:,1], cmap=matplotlib.cm.jet(colours), c=matplotlib.cm.spectral(colours))
else:
    fig=plt.figure()
    ax = p3.Axes3D(fig)
    ax.scatter3D(space[:,0],space[:,2],space[:,4], c=matplotlib.cm.spectral(colours))
    ax.set_xlabel('X')
    ax.set_ylabel('Y')
    ax.set_zlabel('Z')

plt.show()
        
        
#plt.show()

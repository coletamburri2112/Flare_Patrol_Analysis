
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.animation as animation
import sunpy.visualization.colormaps as cm
import matplotlib

def storeSequence(data, movieName, dpi=300, write=True):
  fig =plt.figure(dpi=300)
  max1=np.max(data[0])
  min1=np.min(data[0])
  #norm = (data[0]-min1)/(max1-min1)
  norm=data[0]/max1
  im = plt.imshow(norm,cmap=matplotlib.colormaps['afmhot'],vmin=0,vmax=.98)
  plt.tick_params(
    axis='x',          # changes apply to the x-axis
    which='both',      # both major and minor ticks are affected
    bottom=False,      # ticks along the bottom edge are off
    top=False,         # ticks along the top edge are off
    labelbottom=False) # labels along the bottom edge are off
  
  plt.tick_params(
    axis='y',          # changes apply to the x-axis
    which='both',      # both major and minor ticks are affected
    left=False,      # ticks along the bottom edge are off
    right=False,         # ticks along the top edge are off
    labelleft=False)
  fig.tight_layout()

  def animate(n):
    print(n)
    norm=data[n]/max1
    #norm=(data[n]-min1)/(max1-min1)
    im.set_data(norm)
    return im


  ani = animation.FuncAnimation(fig, animate, frames=data.shape[0], interval=100)
  
  if write:
    writer = animation.writers['ffmpeg'](fps=10)
    ani.save(movieName, writer=writer, dpi=dpi)
  plt.close(fig)  
  
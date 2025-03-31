
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.animation as animation
import sunpy.visualization.colormaps as cm
import matplotlib

def storeSequence(data, movieName, dpi=300, write=True):
  fig =plt.figure()
  im = plt.imshow(data[0,:,:], cmap=matplotlib.colormaps['sdoaia1600'], interpolation='none')
  
  fig.tight_layout()

  def animate(n):
    im.set_data(data[n,:,:])
    return im


  ani = animation.FuncAnimation(fig, animate, frames=data.shape[0], interval=100)
  
  if write:
    writer = animation.writers['ffmpeg'](fps=10)
    ani.save(movieName, writer=writer, dpi=dpi)
  plt.close(fig)  
  
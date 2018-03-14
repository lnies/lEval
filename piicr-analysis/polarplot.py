import matplotlib
import numpy as np
from matplotlib.pyplot import figure, show, rc, grid
import matplotlib.pyplot as plt

# radar green, solid grid lines
rc('grid', color='black', linewidth=1, linestyle='-')
rc('xtick', labelsize=15)
rc('ytick', labelsize=15)

# force square figure and square axes looks better for polar, IMO
width, height = matplotlib.rcParams['figure.figsize']
size = min(width, height)
# make a square figure
fig = figure(figsize=(size, size))
ax = fig.add_axes([0.1, 0.1, 0.8, 0.8], polar=True, alpha=0)
ax.set_theta_zero_location('W', offset=10)

#r = np.arange(0, 10.0, 0.01)
#theta = 2*np.pi*r
#ax.plot(theta, r, color='#ee8d18', lw=3)
ax.set_rmax(5.0)
#grid(True)

ax.set_title("", fontsize=20)
show()

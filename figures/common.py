import matplotlib as mpl
import matplotlib.pyplot as plt

from matplotlib.patches import Patch
from matplotlib.lines import Line2D

import numpy as np

prop_cycle = plt.rcParams['axes.prop_cycle']
colors = prop_cycle.by_key()['color']

titlesize = 18
mpl.rcParams['xtick.labelsize'] = 14
mpl.rcParams['ytick.labelsize'] = 14

mpl.rcParams['axes.labelsize'] = titlesize
mpl.rcParams['legend.fontsize'] = 11

mpl.rcParams['figure.titlesize'] = titlesize

import sys, os
sys.path.append(os.path.abspath('../'))
from common import *

import pandas as pd
import seaborn as sns

grid = 21
L = 5
n = 100

fig, axs = plt.subplots(3, 2, figsize=(10, 12))

L = 10

def get_minmax(label):
    data = []
    for model in ['null', 'full']:
        data.append(pd.read_csv(f'data/L{L}_{label}_{model}.dat', sep='\t')[label])
    return (np.amin(data), np.amax(data))


for i, (label, stat_name) in enumerate([('mean', 'mean'), ('var', 'variance'), ('gamma', 'gamma epistasis')]):
    min, max = get_minmax(label)
    for j, (model, model_name) in enumerate([('full', 'ecological'), ('null', 'non-ecological')]):
        title = fr'{stat_name}, {model_name}'
        data = pd.read_csv(f'data/L{L}_{label}_{model}.dat', sep='\t')

        d = np.flip(data[label].to_numpy().reshape((grid, grid)).T)
        plot = sns.heatmap(d, ax=axs[i, j], square=True, cbar_kws={"shrink": 0.8}, vmin=0, vmax=max, rasterized=True)

        plot.set_xticks(np.arange(0, grid+1, 5.25))
        plot.set_xticklabels([0, 0.25, 0.5, 0.75, 1])

        plot.set_yticks(np.arange(0, grid+1, 5.25))
        plot.set_yticklabels([1, 0.5, 0, -0.5, -1])

        axs[i, j].set_xlabel(r'epistasis $e$')
        axs[i, j].set_ylabel(r'trait correlation $\rho$')

        axs[i, j].set_title(title)

fig.subplots_adjust(wspace=0.45, hspace=0.45)

plt.savefig(f'../figS1.png', dpi=300)
plt.savefig(f'../figS1.pdf', dpi=300)

import sys, os
sys.path.append(os.path.abspath('../'))
from common import *

import pandas as pd
import seaborn as sns

grid = 21

fig, axs = plt.subplots(2, 2, figsize=(10, 8))

label = 'mean_phenotypic_distance'

for i, L in enumerate([5, 10]):
    max_value = 0
    for j, (model, model_name) in enumerate([('full', 'ecological'), ('null', 'non-ecological')]):
        title = fr'L = {L}, {model_name}'
        data = pd.read_csv(f'data/L{L}_{label}_{model}.dat', sep='\t')

        d = np.flip(data[label].to_numpy().reshape((grid, grid)).T)

        if model == 'full':
            max_value = np.max(d)

        plot = sns.heatmap(d, ax=axs[i, j], square=True, cbar_kws={"shrink": 0.8}, vmin=0, vmax=max_value, rasterized=True)

        plot.set_xticks(np.arange(0, grid+1, 5.25))
        plot.set_xticklabels([0, 0.25, 0.5, 0.75, 1])

        plot.set_yticks(np.arange(0, grid+1, 5.25))
        plot.set_yticklabels([1, 0.5, 0, -0.5, -1])

        axs[i, j].set_xlabel(r'epistasis $e$')
        axs[i, j].set_ylabel(r'trait correlation $\rho$')

        axs[i, j].set_title(title)

fig.suptitle('mean phenotypic distance')

fig.subplots_adjust(wspace=0.45, hspace=0.45)

plt.savefig(f'../fig5.png', dpi=300)
plt.savefig(f'../fig5.pdf', dpi=300)

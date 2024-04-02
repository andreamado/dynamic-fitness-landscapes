import pandas as pd
import seaborn as sns
from glob import glob
import re

import sys, os
sys.path.append(os.path.abspath('../'))
from common import *

mpl.rcParams['xtick.labelsize'] = 12
mpl.rcParams['ytick.labelsize'] = 12

n_pop = 1000

cs = [[0.00100, 0.07070], [0.01000, 0.07035], [0.02500, 0.06846], [0.05000, 0.06123], [0.07500, 0.04677], [0.09900, 0.00997], [0.09990, 0.00316], [0.09999, 0.00099]]

correlations = np.array([-0.9, 0, 0.9])
sigma = 0.1

fig, axs = plt.subplots(2, 3, figsize=(10, 8))

for (i, L) in enumerate([5, 10]):
    for (j, correlation) in enumerate(correlations):
        cax = axs[i, j]
        data = {
            'data': [],
            'null': [],
            'sa':   [],
            'sb':   [],
            'correlation': []
        }

        for file in glob(f'processed_data/c_const/haplotype_diversity_L{L}_RMF_S2_mu0_ca*_cb*_c{correlation:.5f}_N{n_pop}_null.dat'):
            file = file.replace('\\', '/')
            for null in [True, False]:
                d = None
                if null:
                    d = np.loadtxt(file)[:, 3]
                else:
                    file = file.replace('null', 'full')
                    d = np.loadtxt(file)[:, 3]

                p = re.compile(f'processed_data/c_const/haplotype_diversity_L{L}_RMF_S2_mu0_ca([0-9.]+)_cb([0-9.]+)_c{correlation:.5f}_N{n_pop}_{"n" if null else "f"}ull.dat')

                m = p.match(file)
                sa = float(m.group(1)) / sigma
                sb = float(m.group(2)) * np.sqrt(2) / sigma

                data['null'].extend(np.ones(len(d), dtype=bool) if null else np.zeros(len(d), dtype=bool))
                data['data'].extend(d)

                data['sa'].extend(np.ones(len(d)) * round(sa, 3))
                data['sb'].extend(np.ones(len(d)) * round(sb, 3))

                data['correlation'].extend(np.ones(len(d)) * correlation)

        data = pd.DataFrame(data)

        sns.boxplot(
            x="sb",
            y='data',
            data=data,
            hue='null',
            fliersize=2,
            ax=cax
        )

        cax.set_xlabel(r'epistasis $e$')
        cax.set_ylabel(r'hap. diversity $h$')
        cax.set_ylim(0, 0.7)

        cax.tick_params(axis='x', rotation=45)

        cax.get_legend().remove()

        cax.set_title(fr'L={L}, correlation $\rho = {correlation}$')


handles = [
    Patch(color=colors[0], label='ecological model'),
    Patch(color=colors[1], label='non-ecological model'),
]
axs[1, 0].legend(
    ncol=2,
    handles=handles,
    frameon=False,
    bbox_to_anchor=(-0.05, -0.7),
    columnspacing=2,
    loc='lower left'
)


plt.tight_layout()
fig.subplots_adjust(bottom=0.2, wspace=0.45, hspace=0.65)

plt.savefig(f'../fig4.png', dpi=300)
plt.savefig(f'../fig4.pdf')
plt.close('all')

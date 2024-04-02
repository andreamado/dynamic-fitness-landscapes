import numpy as np
import matplotlib as mpl
import matplotlib.pyplot as plt

from glob import glob

t_max = 250

L = 5
n_genotypes = 2**L
n_pop = 1000

l = 0

min_alpha = 0.2

min_radius = 0.2
max_radius = 4.

titlesize = 18
labelsize = 18
mpl.rcParams['xtick.labelsize'] = 14
mpl.rcParams['ytick.labelsize'] = 14

def transform_frequency(frequency):
    return (min_radius + frequency) / (1. + min_radius) * max_radius

id_null = '9985'
id_full = '2564'

filename = glob(f'data/data_full/landscape_data_l{l}_{id_full}_{0:06d}.dat')[0].replace('\\', '/')

genotypes = np.loadtxt(filename)[:, :L]
for i in range(L-1, -1, -1):
    genotypes = genotypes[genotypes[:, i].argsort()]

cm = plt.get_cmap('viridis')
colors = [cm(x) for x in np.linspace(0, 1, n_genotypes)]


fig, axs = plt.subplots(2, 2, figsize=(10, 7))

for i, model in enumerate(['null', 'full']):
    id = id_full if model == 'full' else id_null
    title = 'ecological' if model == 'full' else 'non-ecological'

    filename = glob(f'data/fitness_L{L}_{model}_l{l}_{id}.dat')
    fitness = np.loadtxt(filename[0])

    filename = glob(f'data/frequency_L{L}_{model}_l{l}_{id}.dat')
    frequency = np.loadtxt(filename[0])

    for t in range(t_max):
        for g in range(n_genotypes):
            axs[0, i].scatter(t, np.log(fitness[t, g]), s=transform_frequency(frequency[t, g])**2, color=colors[g], alpha=(min_alpha + (1-min_alpha)*frequency[t, g]), edgecolors='none', rasterized=True)

    axs[0, i].set_title(title, fontsize=titlesize)
    axs[0, i].set_xlabel('generations', fontsize=labelsize)
    axs[0, i].set_ylabel('rel. fitness', fontsize=labelsize)

    axs[0, i].set_ylim(-0.5, 0.5)

    del fitness
    del frequency

################################################################################
# Haplotype diversity
max = 0.
for model in ['null', 'full']:
    id = id_full if model == 'full' else id_null
    filename = glob(f'data/data_{model}/L{L}*_{id}.dat')
    data = np.loadtxt(filename[0])

    max_local = np.amax(data[:t_max, 5])
    if max < max_local:
        max = max_local

    label = 'ecological' if model == 'full' else 'non-ecological'
    axs[1, 0].plot(data[:t_max, 3], data[:t_max, 5], label=label)

    del data

axs[1, 0].set_xlabel('generations', fontsize=labelsize)
axs[1, 0].set_ylabel('haplotype diversity $h$', fontsize=labelsize)

axs[1, 0].set_ylim(0, max + 0.1)

axs[1, 0].legend(
    loc='lower left',
    bbox_to_anchor=(-0.05, -0.75),
    frameon=False,
    fontsize=titlesize,
    ncol=2
)

################################################################################
# Mean phenotypic distance
max = 0.
for model in ['null', 'full']:
    id = id_full if model == 'full' else id_null
    filename = glob(f'data/data_{model}/L{L}*_{id}.dat')
    data = np.loadtxt(filename[0])

    max_local = np.amax(data[:t_max, 16])
    if max < max_local:
        max = max_local

    axs[1, 1].plot(data[:t_max, 3], data[:t_max, 16])

    del data

axs[1, 1].set_xlabel('generations', fontsize=labelsize)
axs[1, 1].set_ylabel('phenotypic distance $D$', fontsize=labelsize)

axs[1, 1].set_ylim(0, max + 0.1)

################################################################################
# Save figure

plt.tight_layout()
fig.subplots_adjust(bottom=0.2, wspace=0.45, hspace=0.45)

plt.savefig(f'../fig2.pdf', dpi=600)
plt.savefig(f'../fig2.png', dpi=300)
plt.close('all')

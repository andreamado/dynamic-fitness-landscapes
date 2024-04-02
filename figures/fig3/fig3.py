import sys, os
sys.path.append(os.path.abspath('../'))
from common import *

L = 10
norm = 100

fig, axs = plt.subplots(1, 2, figsize=(10, 4))

a = np.array([0., 0.7071, 1.])
e = np.flip(np.around(np.sqrt(1 - a**2), 4))
cs = [-0.999, -0.8, -0.6, -0.4, -0.2, 0, 0.2, 0.4, 0.6, 0.8, 0.999]
for fi, file in enumerate([f'a_const_L{L}_a09999', f'a_const_L{L}_asqrt05', f'a_const_L{L}_a00001']):
    f = np.loadtxt(f'data/{file}_full.dat')
    n = np.loadtxt(f'data/{file}_null.dat')

    N = len(f[:, 1])

    accf = np.sum((f[:, 2:]/norm)**2, axis=1)
    accn = np.sum((n[:, 2:]/norm)**2, axis=1)

    axs[0].plot(cs, accf, '-',  color=colors[fi])
    axs[0].plot(cs, accn, '--', color=colors[fi])

axs[0].set_xlabel(r'trait correlation $\rho$')

axs[0].set_ylim((0, 1))
axs[0].set_ylabel('repeatability')

handles = [
    Patch(label='Model', color='none'),
    Patch(label=r'Epistasis $e$', color='none'),
    Line2D([0], [0], color='grey', lw=1, label='eco', linestyle='-'),
    Patch(color=colors[0], label=e[0]),
    Line2D([0], [0], color='grey', lw=1, label='non-eco', linestyle='--'),
    Patch(color=colors[1], label=e[1]),
    Patch(label=r'', color='none'),
    Patch(color=colors[2], label=e[2]),
]

axs[0].legend(
    ncol=4, handles=handles,
    frameon=False,
    bbox_to_anchor=(-0.25, -0.55),
    columnspacing=1,
    loc='lower left'
)



cs = [-0.9, 0., 0.9]
full = np.loadtxt(f'data/c_const_L{L}_full.dat')
null = np.loadtxt(f'data/c_const_L{L}_null.dat')

a_vals = np.array([0.01, 0.1, 0.25, 0.5, 0.75, 0.99, 0.999, 0.9999])
e_vals = np.sqrt(1 - a_vals**2)
for fi, c in enumerate(cs):
    f = full[full[:, 1] == c]
    n = null[null[:, 1] == c]

    N = len(f[:, 0])

    accf = np.sum((f[:, 2:]/norm)**2, axis=1)
    accn = np.sum((n[:, 2:]/norm)**2, axis=1)

    axs[1].plot(e_vals, accf, '-',  color=colors[fi], label=fr'$\rho$ = {c} (+ eco)')
    axs[1].plot(e_vals, accn, '--', color=colors[fi], label=fr'$\rho$ = {c} (- eco)')

axs[1].set_xlabel('epistasis $e$')

axs[1].set_ylim((0, 1))
axs[1].set_ylabel('repeatability')


handles = [
    Patch(label='Model', color='none'),
    Patch(label=r'Correlation $\rho$', color='none'),
    Line2D([0], [0], color='grey', lw=1, label='eco', linestyle='-'),
    Patch(color=colors[0], label=cs[0]),
    Line2D([0], [0], color='grey', lw=1, label='non-eco', linestyle='--'),
    Patch(color=colors[1], label=cs[1]),
    Patch(label=r'', color='none'),
    Patch(color=colors[2], label=cs[2]),
]

axs[1].legend(
    ncol=4, handles=handles,
    frameon=False,
    bbox_to_anchor=(-0.25, -0.55),
    columnspacing=1,
    loc='lower left'
)

fig.subplots_adjust(bottom=0.3, wspace=0.45)

plt.savefig(f'../fig3.pdf')
plt.savefig(f'../fig3.png', dpi=300)
plt.close('all')

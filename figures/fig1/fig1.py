import pickle
import numpy as np
import matplotlib.pyplot as plt
import os

class PhenotypicLandscape():
    """docstring for PhenotypicLandscape."""

    def __init__(
        self, n_resources, L, cov_matrix_additive, cov_matrix_epistatic, mean,
        multiplicative_phenotype=True
    ):
        self.L = L
        self.size = 2**L
        self.n_resources = n_resources

        self.cov_matrix_additive  = cov_matrix_additive
        self.cov_matrix_epistatic = cov_matrix_epistatic
        self.mean = mean

        self.additive_effect = np.random.multivariate_normal(
            self.mean, self.cov_matrix_additive, size=L)
        self.epistatic_component = np.random.multivariate_normal(
            np.zeros(n_resources), self.cov_matrix_epistatic, size=self.size)

        self.seqs = np.flip(np.unpackbits(np.array([[i] for i in range(
            self.size)], dtype=np.uint8), axis=1, count=L, bitorder='little'), axis=1)
        self.idx = np.flip(2**np.arange(L))

        self.phenotype = \
            np.dot(self.seqs, self.additive_effect) + self.epistatic_component
        if multiplicative_phenotype:
            self.phenotype = np.exp(self.phenotype)

    def save(self, filename):
        with open(filename, 'wb') as f:
            pickle.dump(self, f)

    def load(filename):
        return pickle.load(filename)

    def plot(self, ax, title=None):
        genotype_pairs = []
        for i in range(2**self.L):
            gen0 = self.seqs[i]
            for j in range(self.L):
                gen1 = gen0.copy()
                gen1[j] = 1 - gen1[j]

                genotype_pairs.append([gen0, gen1])

        for pair in genotype_pairs:
            idx0 = pair[0].dot(self.idx)
            idx1 = pair[1].dot(self.idx)

            (x0, y0) = self.phenotype[idx0]
            (x1, y1) = self.phenotype[idx1]

            ax.plot([x0, x1], [y0, y1], color='gray', linewidth=0.5)

        ax.plot(self.phenotype.T[0], self.phenotype.T[1], '.', color='black')

        if title:
            ax.set_title(title, weight='bold', fontsize=16)

        ax.set_xlabel('trait 1', fontsize=16)
        ax.set_ylabel('trait 2', fontsize=16)

        ax.set_aspect('equal')


# parameters
L = 4
n_resources = 2

sigma = 0.1
a = 1./np.sqrt(2.)
mean = np.zeros(n_resources)
rs = [-0.75, 0, 0.75]

# create landscapes
landscapes = []
for r in rs:
    covariance_matrix = np.array([[1,   r],
                                  [r,   1]])

    landscape = None
    filename = f'fig_traits_r{r}'.replace('-', 'm') + '.pickle'
    if os.path.exists(filename):
        with open(filename, 'rb') as f:
            landscape = PhenotypicLandscape.load(f)
    else:
        landscape = PhenotypicLandscape(
            n_resources,
            L,
            sigma**2 * a**2         * covariance_matrix,
            sigma**2 * (1 - a**2)/2 * covariance_matrix,
            mean,
            multiplicative_phenotype=False
        )
        landscape.save(filename)

    landscapes.append(landscape)


fig, axs = plt.subplots(1, 3, figsize=(10, 4))

min = 1.1 * np.min([landscape.phenotype for landscape in landscapes])
max = 1.1 * np.max([landscape.phenotype for landscape in landscapes])

for i in range(3):
    r = rs[i]
    title = f'$(\\rho = {r})$'
    if r < 0.:
        title = 'negative correlation\n' + title
    elif r > 0.:
        title = 'positive correlation\n' + title
    else:
        title = 'no correlation\n' + title

    landscapes[i].plot(axs[i], title)
    axs[i].set_xlim(min, max)
    axs[i].set_ylim(min, max)

    axs[i].set_xticks([-2*sigma, 0, 2*sigma])
    axs[i].set_yticks([-2*sigma, 0, 2*sigma])

    axs[i].tick_params(labelsize=14)

plt.tight_layout()

plt.savefig('../fig1.pdf')
plt.savefig('../fig1.png', dpi=300)
plt.close('all')

import os
import datetime
from classes.snp_database import SnpDatabase
from classes.tfam import Tfam
from classes.intensities import Intensities
import numpy as np
import matplotlib.pyplot as plt
from sklearn.neighbors import KernelDensity

MISSING_POSITIVE = [
    "SNP_A-4249117",
    "SNP_A-8311817",
    "SNP_A-4284077",
    "SNP_A-2199193",
    "SNP_A-1933659",
    "SNP_A-2015495",
    "SNP_A-8331382",
    "SNP_A-8565877",
    "SNP_A-2305037"
]

MISSING_NEGATIVE = [
    "SNP_A-4256494",
    "SNP_A-2194348",
    "SNP_A-8657565",
    "SNP_A-1851845",
    "SNP_A-1962720",
    "SNP_A-1850450",
    "SNP_A-8552096",
    "SNP_A-2301699",
    "SNP_A-8311383"
]

y_offset = 1e-5
kde_bandwith = 75
min_x = 0
max_x = 2250
parent_mark = 'ro'
children_mark = 'b+'


def plot_a_b_intensities(probe_id, intensities, parents_index, offspring_index):
    A = intensities.get(probe_id + '-A')
    B = intensities.get(probe_id + '-B')

    num_parents = len(parents_index)
    num_children = len(offspring_index)
    total = len(intensities.get(probe_id + '-A'))

    assert total == num_children + num_parents

    fig, axes = plt.subplots(nrows=3, ncols=2, sharex='all', sharey='all')

    y_parents = np.zeros(num_parents) + y_offset
    y_children = np.zeros(num_children) + y_offset

    for col, intens in enumerate([A, B]):
        axes[0, col].plot(intens[parents_index], y_parents, parent_mark)
        fit_gaussian_kde(intens[parents_index], axes[0, col])

        axes[1, col].plot(intens[offspring_index], y_children, children_mark)
        fit_gaussian_kde(intens[offspring_index], axes[1, col])

        axes[2, col].plot(intens[parents_index], y_parents, parent_mark, intens[offspring_index], y_children, children_mark)
        fit_gaussian_kde(intens, axes[2, col])

    plt.show()


def fit_gaussian_kde(intens, ax):
    X = intens[:, np.newaxis]
    X_plot = np.linspace(min_x, max_x, 1000)[:, np.newaxis]
    kde = KernelDensity(kernel='gaussian', bandwidth=kde_bandwith).fit(X)
    log_dens = kde.score_samples(X_plot)
    ax.plot(X_plot[:, 0], np.exp(log_dens), 'k-')


if __name__ == "__main__":
    print 'Started:', datetime.datetime.now().isoformat()

    birdseed_base = '/home/victor/Escritorio/matesanz2015'
    birdseed2_base = '/home/victor/Escritorio/Genotipado_Alternativo/data/birdseed_out'

    tfam = Tfam(os.path.join(birdseed_base, 'matesanz2015.tfam'))
    db = SnpDatabase()
    db.load_from_birdseed(birdseed_base, '8090939')
    intensities = Intensities.load_birdseed_summary_intensities(
        os.path.join(birdseed2_base, 'birdseed-dev.summary.txt.gz'), limit=500)

    parents_index = np.array(tfam.get_parents_index(intensities['subjects']))
    offspring_index = np.array(tfam.get_offspring_index(intensities['subjects']))

    probe_id = 'SNP_A-1975121'  # MISSING_POSITIVE[0]
    plot_a_b_intensities(probe_id, intensities['intensities'], parents_index, offspring_index)

    print 'Finished:', datetime.datetime.now().isoformat()

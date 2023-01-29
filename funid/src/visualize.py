from Bio import SeqIO
from Bio import pairwise2
from Bio.pairwise2 import format_alignment
from sklearn.decomposition import PCA
from functools import partial
import matplotlib.pyplot as plt
from matplotlib.patches import Rectangle
import pandas as pd
import matplotlib
from matplotlib.pyplot import cm
import multiprocessing as mp
import numpy as np
import os

# for visualizing statistical analysis on groupal clustering


def plot_pca(name, group_range, seq_matrix, path):

    matplotlib.use("cairo")

    plt.clf()
    plt.cla()
    X = seq_matrix
    pca = PCA(n_components=3)
    X3d = pca.fit_transform(X)

    fig, ax = plt.subplots()

    len_f = 0

    for n, group in enumerate(group_range.keys()):
        i = group_range[group][0]
        j = group_range[group][1]
        ax.scatter(X3d[i:j], X3d[i:j], cmap=cm.get_cmap("rainbow", n))

    plt.legend()

    plt.savefig(
        f"{path.result}/pca_{name}.svg", dpi=150, bbox_inches="tight", pad_inches=0
    )
    plt.close()


def plot_heatmap(name, group_range, seq_matrix, path, vmax=100):

    # matplotlib.use("cairo")

    plt.clf()
    plt.cla()
    fig = plt.figure()
    ax = fig.add_subplot(1, 1, 1)
    plt.pcolor(seq_matrix, vmax=vmax)
    plt.colorbar()
    for key in group_range.keys():
        ax.add_patch(
            Rectangle(
                (group_range[key][0], group_range[key][0]),
                group_range[key][1] - group_range[key][0],
                group_range[key][1] - group_range[key][0],
                fill=False,
                edgecolor="white",
                lw=0.5,
            )
        )
    ax.set_xticks(
        [(group_range[x][0] + group_range[x][1]) / 2 for x in group_range.keys()]
    )
    ax.set_xticklabels([x for x in group_range.keys()], fontsize=2)
    # plt.savefig(f"{path.result}/heatmap_{name}.svg", dpi=1200)
    plt.savefig(f"{path.result}/heatmap_{name}.png", dpi=1200)

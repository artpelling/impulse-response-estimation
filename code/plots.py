#!/usr/bin/env python3

import matplotlib.pyplot as plt
import numpy as np
from pathlib import Path


single_column_width = 5.90666
double_column_width = 3.4402
ratio = (1 + np.sqrt(5)) / 2 # golden ratio

w = double_column_width
h = w/ratio

plt.rcParams.update({
    'figure.figsize': (w, h),
    'figure.dpi': 200,
    'figure.constrained_layout.use': True,
    'lines.linewidth': 0.8,
    'pgf.texsystem': 'pdflatex',
    'font.family': 'serif',
    'text.usetex': True,
    'text.latex.preamble': r'\usepackage{amsfonts}',
    'font.size': 10.0,
    'pgf.rcfonts': False,
    'savefig.bbox': 'tight',
    'savefig.format': 'pdf',
    'contour.linewidth': 0.1,
})

import_dir = Path(__file__).parents[1] / 'data'
export_dir = Path(__file__).parents[1] / 'latex' / 'figures'


def plot():
    name = 'plot'
    fig, ax = plt.subplots()
    plt.savefig(export_dir / name)


if __name__ == '__main__':
    plt.close('all')
    plot()

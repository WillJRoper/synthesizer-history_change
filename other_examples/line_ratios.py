"""
This example generates line-ratios vs. Z for all available line ratios
"""

import matplotlib.pyplot as plt
import numpy as np
from synthesizer.grid import Grid, get_available_lines
from synthesizer.line import LineRatios

if __name__ == '__main__':

    grid_dir = '/Users/stephenwilkins/Dropbox/Research/data/synthesizer/grids'

    grid_names = [
        'bpass-2.2.1-bin_chabrier03-0.1,300.0_cloudy-geometryspherical',
        'bpass-2.2.1-bin_chabrier03-0.1,300.0_cloudy-grainsFalse',
        'bpass-2.2.1-bin_chabrier03-0.1,300.0_cloudy-d2m0.0-grainsFalse',
        'bpass-2.2.1-bin_chabrier03-0.1,300.0_cloudy',
    ]

    fig, axes = plt.subplots(len(LineRatios.available_ratios), 1, sharex=True, figsize=(4, 10))
    plt.subplots_adjust(bottom=0.1, wspace=0.0, hspace=0.0, top=0.95, left=0.2, right=0.95)

    for grid_name in grid_names:

        grid = Grid(grid_name, grid_dir=grid_dir, read_spectra=False, read_lines=True)

        Zsun = grid.metallicities/0.0124

        ia = 0  # 1 Myr

        lines = grid.get_lines_info(grid.line_list, (ia, 0))

        limits = {}

        # or loop over availalable ratios
        for ratio_id, ax in zip(lines.available_ratios, axes):

            ratios = []
            for iZ, Z in enumerate(grid.metallicities):
                grid_point = (ia, iZ)
                lines = grid.get_lines_info(grid.line_list, grid_point)
                ratios.append(lines.get_ratio(ratio_id))

            limits[ratio_id] = (np.min(ratios)*0.5, np.max(ratios)*2)

            ax.plot(Zsun, ratios, label=grid_name)

    for ratio_id, ax in zip(lines.available_ratios, axes):

        ax.set_xlim([0.01, 1])
        ax.set_ylim(limits[ratio_id])
        ax.set_xscale('log')
        ax.set_yscale('log')
        ax.set_ylabel(ratio_id)

    axes[0].legend(fontsize=7)

    fig.savefig('figs/line_ratios.pdf')

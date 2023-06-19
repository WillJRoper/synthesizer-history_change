



import numpy as np
import matplotlib.pyplot as plt
import h5py
import cmasher as cmr

from synthesizer.utils import flatten
from synthesizer.line import LineRatios
from synthesizer.cloudy import Ions
from synthesizer.agn import Feltre16
from synthesizer.ggrid import Grid, flatten_linelist
from synthesizer.plt import single
from unyt import Angstrom


if __name__ == '__main__':

    """
    Plot BPT diagram like Feltre16.

    """

    grid_dir = '/Users/sw376/Dropbox/Research/data/synthesizer/grids'
    grid_name = 'agn_feltre16_alpha_Z_U_nH'
    grid_name += '_spherical'

    # load grid
    grid = Grid(grid_name, grid_dir=grid_dir, read_lines=True)

    # print grid summary
    print(grid)


    # this should *roughly* correspond to the second row, middle panel, though the metallicity is slightly different (0.01 vs. 0.008 in Feltre 16)
   
    params_ = {
        'log10Z': -2.1,
        'log10n_H': 3.,
    }
    

    # initialise figure
    fig, ax = single((4., 4.))


    # add extreme starburst line from Kewley+01

    x = np.arange(-3, 2, 0.01) # log10(NII/Ha)
    y = 0.61/(x-0.47)+1.19# log10(OIII5007/Hb)
    ax.plot(10**x,10**y,lw=2,c='k',alpha=0.1)

    # generate BPT-NII diagram
    diagram_id = 'BPT-NII'

    # get flattened line list
    line_ids = list(flatten(LineRatios().diagrams[diagram_id]))
    print(line_ids)

    # fixed alpha
    for j, alpha in enumerate(grid.bin_centres['alpha']):

        x = []
        y = []

        for i, log10U in enumerate(grid.bin_centres['log10U']):


            params = params_ | {'log10U': log10U, 'alpha': alpha} 

            grid_point_indices = grid.get_grid_point([params[k] for k in grid.axes])

            lines = grid.get_lines_info(line_ids, grid_point_indices)
            x_, y_ = lines.get_diagram(diagram_id)

            print(params, grid_point_indices, x_, y_)

            x.append(x_)
            y.append(y_)

        ax.plot(x, y, label = f'alpha={alpha}', ls = '--')

    
    
    # fixed log10U

    for i, log10U in enumerate(grid.bin_centres['log10U']):
    
        x = []
        y = []

        for j, alpha in enumerate(grid.bin_centres['alpha']):

            params = params_ | {'log10U': log10U, 'alpha': alpha} 

            grid_point_indices = grid.get_grid_point([params[k] for k in grid.axes])

            lines = grid.get_lines_info(line_ids, grid_point_indices)
            x_, y_ = lines.get_diagram(diagram_id)

            print(params, grid_point_indices, x_, y_)

            x.append(x_)
            y.append(y_)

        ax.plot(x, y, label = f'log10U={log10U}', ls='-')



    
    ax.set_xlim([0.003, 5])
    ax.set_ylim([0.05, 50])
    ax.set_xscale('log')
    ax.set_yscale('log')
    ax.legend()

    # grab x and y labels, this time use "fancy" label ids
    xlabel, ylabel = lines.get_diagram_label(diagram_id, fancy=True)

    ax.set_xlabel(xlabel)
    ax.set_ylabel(ylabel)
    plt.show()



    
    
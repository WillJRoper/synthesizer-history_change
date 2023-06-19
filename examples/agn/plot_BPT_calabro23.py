



import numpy as np
import matplotlib.pyplot as plt
import h5py
import cmasher as cmr

from synthesizer.utils import flatten
from synthesizer.lline import LineRatios
from synthesizer.cloudy import Ions
from synthesizer.agn import Feltre16
from synthesizer.ggrid import Grid, flatten_linelist
from synthesizer.plt import single
from unyt import Angstrom


if __name__ == '__main__':

    """
    Plot BPT diagram like Calabro 2023.

    """

    grid_dir = '/Users/sw376/Dropbox/Research/data/synthesizer/grids'
    grid_name = 'agn_cloudy_calabro23-llines'

    # load grid
    grid = Grid(grid_name, grid_dir=grid_dir, read_lines=True, read_spectra=False)

    # print grid summary
    print(grid)

    # initialise figure
    fig, ax = single((4., 4.))

    # add extreme starburst line from Kewley+01
    x = np.arange(-3, 2, 0.01) # log10(NII/Ha)
    y = 0.61/(x-0.47)+1.19# log10(OIII5007/Hb)
    ax.plot(x,y,lw=2,c='k',alpha=0.1)

    # generate BPT-NII diagram
    diagram_id = 'BPT-NII'

    # get flattened line list
    line_ids = list(flatten(LineRatios().diagrams[diagram_id]))
    print(line_ids)

    # calabro colours by log10U and then plots every other model

    colours = cmr.take_cmap_colors('cmr.cosmic', len(grid.bin_centres['log10U']))
     
    for i, (log10U, c) in enumerate(zip(grid.bin_centres['log10U'], colours)):

        for log10Z in  grid.bin_centres['log10Z']:
       
            for aox in grid.bin_centres['aox']:

                for log10n_H in grid.bin_centres['log10n_H']:

                    params = {'log10U': log10U, 'aox': aox, 'log10Z': log10Z, 'log10n_H': log10n_H} 

                    grid_point_indices = grid.get_grid_point([params[k] for k in grid.axes])

                    lines = grid.get_lines_info(line_ids, grid_point_indices)
                    x_, y_ = lines.get_diagram(diagram_id)

                    ax.scatter(np.log10(x_), np.log10(y_), c=c)

    
   

    
    ax.set_xlim([-2.7, 1.3])
    ax.set_ylim([-3.0, 2.0])

    ax.legend()

    # grab x and y labels, this time use "fancy" label ids
    xlabel, ylabel = lines.get_diagram_label(diagram_id, fancy=True)

    ax.set_xlabel(xlabel)
    ax.set_ylabel(ylabel)
    plt.show()



    
    
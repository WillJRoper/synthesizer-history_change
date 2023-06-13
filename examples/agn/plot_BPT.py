



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


    grid_dir = '/Users/sw376/Dropbox/Research/data/synthesizer/grids'
    # grid_name = 'agn_cloudy'
    grid_name = 'agn_feltre16'

    # load grid
    grid = Grid(grid_name, grid_dir=grid_dir, read_lines=True)

    # print grid summary
    print(grid)

    # define parameters of the grid point we're interested in
    if grid_name == 'agn_cloudy':
        params = {
            'log10T': 7.,
            'log10Z': -2.,
            'log10U': -2.,
            'log10n_H': 4.,
        }

    if grid_name == 'agn_feltre16':
        params = {
            'alpha': -1.5,
            'log10Z': -2.,
            'log10U': -2.,
            'log10n_H': 3.,
        }
    

    # get the grid point indices
    grid_point_indices = grid.get_grid_point([params[k] for k in grid.axes])
    print(grid_point_indices)



    # generate BPT-NII diagram
    diagram_id = 'BPT-NII'

    # get flattened line list
    line_ids = list(flatten(LineRatios().diagrams[diagram_id]))
    print(line_ids)

    x = []
    y = []
    for i, p in enumerate(grid.bin_centres['log10U']):

        # this is not the most efficient way to do this
        # params_ = params | {'log10Z': log10Z} 
        params_ = params | {'log10U': p} 
        grid_point_indices = grid.get_grid_point([params_[k] for k in grid.axes])

        lines = grid.get_lines_info(line_ids, grid_point_indices)
        x_, y_ = lines.get_diagram(diagram_id)

        print(p, grid_point_indices, x_, y_)

        x.append(x_)
        y.append(y_)

    plt.plot(x, y)
    plt.xlim([0.01, 10])
    plt.ylim([0.05, 20])
    plt.xscale('log')
    plt.yscale('log')

    # grab x and y labels, this time use "fancy" label ids
    xlabel, ylabel = lines.get_diagram_label(diagram_id, fancy=True)

    plt.xlabel(xlabel)
    plt.ylabel(ylabel)
    plt.show()



    
    # plot lines as a function of metallicity for different temperatures in BPT space

    # fig, ax = single((3.5, 3.5))

    # spec_name = 'incident'

    # colours = cmr.take_cmap_colors('cmr.ember', len(grid.bin_centres['log10T'][:]))

    # for i, (log10T, c) in enumerate(zip(grid.bin_centres['log10T'][:], colours)):

    #     for iZ


    #     # update grid point indices to select the correct spectra
    #    
    #     # get the specific spectra at the desired grid point
    #     lnu = grid.spectra[spec_name][tuple(gpi)] # this is not an Sed object

    #     ax.plot(np.log10(grid.lam), np.log10(lnu), label = log10T, c=c)

    
    
    
    # ax.set_xlim([2., 4.5])
    # ax.set_ylim([-20, -14.5])
    # ax.legend(fontsize = 8, labelspacing = 0.0, title=r'$\rm log_{10}(T/K)=$')
    # ax.set_xlabel(r'$\rm log_{10}(\lambda/\AA)$')
    # ax.set_ylabel(r'$\rm log_{10}(L_{\nu})$')

    # plt.show()



import numpy as np
import matplotlib.pyplot as plt
import h5py
import cmasher as cmr

from synthesizer.cloudy import Ions
from synthesizer.agn import Feltre16
from synthesizer.ggrid import Grid
from synthesizer.plt import single
from unyt import Angstrom


if __name__ == '__main__':


    grid_dir = '/Users/sw376/Dropbox/Research/data/synthesizer/grids'
    grid_name = 'agn_cloudy'

    # load grid
    grid = Grid(grid_name, grid_dir=grid_dir)

    # print grid summary
    print(grid)

    # define parameters of the grid point we're interested in
    params = {
        'log10T': 5.,
        'log10Z': -2.,
        'log10U': -2.,
        'log10n_H': 2.,
    }

    # get the grid point indices
    grid_point_indices = grid.get_grid_point([params[k] for k in grid.axes])
    print(grid_point_indices)


    # plot the spectra at this grid point

    fig, ax = single((6., 3.5))

    for spec_name in grid.spec_names:

        # get the specific spectra at the desired grid point
        lnu = grid.spectra[spec_name][tuple(grid_point_indices)] # this is not an Sed object

        ax.plot(np.log10(grid.lam), np.log10(lnu), label = spec_name)

    ax.set_xlim([1., 4.5])
    ax.set_ylim([-20, -14])
    ax.legend(fontsize = 8, labelspacing = 0.0)
    ax.set_xlabel(r'$\rm log_{10}(\lambda/\AA)$')
    ax.set_ylabel(r'$\rm log_{10}(L_{\nu})$')

    plt.show()

    
    # plot the spectra for the different temperatures at this grid point

    fig, ax = single((6., 3.5))

    spec_name = 'incident'

    colours = cmr.take_cmap_colors('cmr.ember', len(grid.bin_centres['log10T'][:]))


    # plot ionisation energies

    # for ion, wavelength in Ions.wavelength.items(): # all ions
    for ion in ['HI', 'HeII', 'OII']:
        wavelength = Ions.wavelength[ion]
        l = wavelength.to('angstrom')
        ax.axvline(np.log10(l), c='k', lw=1, alpha=0.1)
        ax.text(np.log10(l), -14.5, ion, fontsize=7, rotation=90)

    for i, (log10T, c) in enumerate(zip(grid.bin_centres['log10T'][:], colours)):

        # update grid point indices to select the correct spectra
        gpi = list(grid_point_indices)
        gpi[0] = i
        gpi = tuple(gpi)

        # get the specific spectra at the desired grid point
        lnu = grid.spectra[spec_name][tuple(gpi)] # this is not an Sed object

        ax.plot(np.log10(grid.lam), np.log10(lnu), label = log10T, c=c)

    ax.set_xlim([2., 4.5])
    ax.set_ylim([-20, -14.5])
    ax.legend(fontsize = 8, labelspacing = 0.0, title=r'$\rm log_{10}(T/K)=$')
    ax.set_xlabel(r'$\rm log_{10}(\lambda/\AA)$')
    ax.set_ylabel(r'$\rm log_{10}(L_{\nu})$')



    plt.show()
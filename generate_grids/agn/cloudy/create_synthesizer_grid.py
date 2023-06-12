
"""
This reads in a cloudy grid of models

"""


from scipy import integrate
import os
import shutil
from synthesizer.utils import read_params
from synthesizer.cloudy import read_wavelength, read_continuum, read_lines
from synthesizer.sed import calculate_Q
from unyt import eV
import argparse
import numpy as np
import h5py
import yaml





def get_grid_properties(hf, verbose = True):


    """ 
    Get the properties of the grid including the dimensions etc.
    
    """

    # the grid axes
    axes = hf.attrs['grid_axes']
    if verbose: print(f'axes: {axes}')

    # number of axes
    n_axes = len(axes)
    if verbose: print(f'number of axes: {n_axes}')

    # the shape of the grid (useful for creating outputs)
    shape = list([len(hf[axis][:]) for axis in axes])
    if verbose: print(f'shape: {shape}')

    # determine number of models
    n_models = np.prod(shape)
    if verbose: print(f'number of models to run: {n_models}')

    # create the mesh of the grid
    mesh = np.array(np.meshgrid(*[np.array(hf[axis][:]) for axis in axes]))

    # create the list of the models 
    model_list = mesh.T.reshape(n_models, n_axes)
    if verbose: 
        print('model list:')
        print(model_list)

    # create a list of the indices

    index_mesh = np.array(np.meshgrid(*[range(n) for n in shape]))

    index_list =  index_mesh.T.reshape(n_models, n_axes)
    if verbose: 
        print('index list:')
        print(index_list)


    return axes, n_axes, shape, n_models, mesh, model_list, index_list





def check_cloudy_runs(grid_name, synthesizer_data_dir, replace=False):
    """
    Check that all the cloudy runs have run properly

    Parameters
    ----------
    grid_name : str
        Name of the grid
    synthesizer_data_dir : str
        Directory where synthesizer data is kept.
    replace : boolean
        If a run has failed simply replace the model with the previous one
    """

    # open the new grid
    with h5py.File(f'{synthesizer_data_dir}/grids/{grid_name}.hdf5', 'r') as hf:

        # Get the properties of the grid including the dimensions etc.
        axes, n_axes, shape, n_models, mesh, model_list, index_list = get_grid_properties(hf)

        # list of failed models
        failed_list = []

        for i, grid_params_ in enumerate(model_list):

            infile = f"{synthesizer_data_dir}/{grid_name.replace('_','/')}/{i}"

            failed = False

            # check if files exist
            if not os.path.isfile(infile+'.cont'):  # attempt to open run.
                failed = True   
            if not os.path.isfile(infile+'.lines'):  # attempt to open run.
                failed = True
            
            # if they exist also check they have size >0
            if not failed:
                if os.path.getsize(infile+'.cont') == 0:
                    failed = True
                if os.path.getsize(infile+'.lines') == 0:
                    failed = True

            if failed:
                print(i)    
                failed_list.append(i)

        return failed_list


def add_spectra(grid_name, synthesizer_data_dir):
    """
    Open cloudy spectra and add them to the grid

    Parameters
    ----------
    grid_name : str
        Name of the grid
    synthesizer_data_dir : str
        Directory where synthesizer data is kept.
    """

    # spec_names = ['incident','transmitted','nebular','nebular_continuum','total','linecont']
    #  the cloudy spectra to save (others can be generated later)
    spec_names = ['incident', 'transmitted', 'nebular', 'linecont']

    # open the new grid
    with h5py.File(f'{synthesizer_data_dir}/grids/{grid_name}.hdf5', 'a') as hf:

        # Get the properties of the grid including the dimensions etc.
        axes, n_axes, shape, n_models, mesh, model_list, index_list = get_grid_properties(hf)

        # read first spectra from the first grid point to get length and wavelength grid
        lam = read_wavelength(f"{synthesizer_data_dir}/{grid_name.replace('_','/')}/0")

        if 'spectra' in hf:
            del hf['spectra']

        spectra = hf.create_group('spectra')  # create a group holding the spectra in the grid file
        spectra.attrs['spec_names'] = spec_names  # save list of spectra as attribute

        spectra['wavelength'] = lam  # save the wavelength
        nu = 3E8 / (lam*1E-10)

        nlam = len(lam)  # number of wavelength points

        # make spectral grids and set them to zero
        for spec_name in spec_names:
            spectra[spec_name] = np.zeros((*shape, nlam))

        # array for holding the normalisation which is calculated below and used by lines
        normalisation = np.zeros(shape)

        for i, indices in enumerate(index_list):

            indices = tuple(indices)

            # define the infile
            infile = f"{synthesizer_data_dir}/{grid_name.replace('_','/')}/{i}"

            # read the continuum file containing the spectra 
            spec_dict = read_continuum(infile, return_dict=True)

            # for an arbitrary grid, we should normalise by the bolometric luminosity of the incident spectra
            norm = np.trapz(spec_dict['incident'][::-1], x=nu[::-1])

            # save normalisation for later use (rescaling lines)
            normalisation[indices] = norm

            # save the normalised spectrum to the correct grid point 
            for spec_name in spec_names:
                spectra[spec_name][indices] = spec_dict[spec_name] / norm

    return normalisation


def get_default_line_list(interesting=True):

    with open('default_lines.dat') as f:
        line_list = f.read().splitlines()

    if interesting:

        with open('interesting_lines.dat') as f:
            line_list += f.read().splitlines()

    return line_list





def add_lines(grid_name, synthesizer_data_dir, normalisation, lines_to_include):
    """
    Open cloudy lines and add them to the HDF5 grid

    Parameters
    ----------
    grid_name: str
        Name of the grid.
    synthesizer_data_dir: str
        Directory where synthesizer data is kept.
    dlog10Q
        The difference between the original and cloudy log10Q used for rescaling the cloudy spectra
    """

    # open the new grid
    with h5py.File(f'{synthesizer_data_dir}/grids/{grid_name}.hdf5', 'a') as hf:

        # Get the properties of the grid including the dimensions etc.
        axes, n_axes, shape, n_models, mesh, model_list, index_list = get_grid_properties(hf)

        # delete lines group if it already exists
        if 'lines' in hf:
            del hf['lines']

        # create group for holding lines
        lines = hf.create_group('lines')
        # lines.attrs['lines'] = list(lines_to_include)  # save list of spectra as attribute

        # set up output arrays
        for line_id in lines_to_include:
            lines[f'{line_id}/luminosity'] = np.zeros(shape)
            lines[f'{line_id}/stellar_continuum'] = np.zeros(shape)
            lines[f'{line_id}/nebular_continuum'] = np.zeros(shape)
            lines[f'{line_id}/continuum'] = np.zeros(shape)

        for i, indices in enumerate(index_list):

            # convert indices array to tuple
            indices = tuple(indices)

            # define the infile
            infile = f"{synthesizer_data_dir}/{grid_name.replace('_','/')}/{i}"

            # get TOTAL continuum spectra
            nebular_continuum = spectra['nebular'][indices] - spectra['linecont'][indices]
            continuum = spectra['transmitted'][indices] + nebular_continuum

            # get line quantities
            id, blend, wavelength, intrinsic, emergent = read_lines(infile)

            # identify lines we want to keep
            s = np.nonzero(np.in1d(id, np.array(lines_to_include)))[0]

            for id_, wavelength_, emergent_ in zip(id[s], wavelength[s], emergent[s]):

                line = lines[id_]

                # save line wavelength
                line.attrs['wavelength'] = wavelength_

                # calculate line luminosity and save it. Uses normalisation from spectra.
                line['luminosity'][indices] = 10**(emergent_)/normalisation[indices]  # erg s^-1
                
                # calculate stellar continuum at the line wavelength and save it. 
                line['stellar_continuum'][indices] = np.interp(
                    wavelength_, lam, spectra['transmitted'][indices])  # erg s^-1 Hz^-1
                
                # calculate nebular continuum at the line wavelength and save it. 
                line['nebular_continuum'][indices] = np.interp(
                    wavelength_, lam, nebular_continuum)  # erg s^-1 Hz^-1
                
                # calculate total continuum at the line wavelength and save it. 
                line['continuum'][indices] = np.interp(
                    wavelength_, lam, continuum)  # erg s^-1 Hz^-1


if __name__ == "__main__":

    parser = argparse.ArgumentParser(description=('Create synthesizer HDF5 grid '
                                                  'for a given grid.'))

    parser.add_argument("-synthesizer_data_dir", type=str, required=True) # path to synthesizer_data_dir
    parser.add_argument("-grid", "--grid", type=str, required=True)

    args = parser.parse_args()

    synthesizer_data_dir = args.synthesizer_data_dir
    grid_name = args.grid

    # check cloudy runs
    failed_list = check_cloudy_runs(grid_name, synthesizer_data_dir)

    print(failed_list)

    # if failed prompt to re-run
    if len(failed_list)>0:

        print(f'ERROR: {len(failed_list)} cloudy runs have failed. You should re-run these with command:')
        print(f'  qsub -t 1:{len(failed_list)}  run_grid.job')

        # replace input_names with list of failed runs
        with open(f"{synthesizer_data_dir}/{grid_name.replace('_','/')}/input_names.txt", "w") as myfile:
            myfile.write('\n'.join(map(str, failed_list)))

    #if not failed, go ahead and add spectra and lines
    else:
        
        # add spectra
        normalisation = add_spectra(grid_name, synthesizer_data_dir)

        # get list of lines
        lines_to_include = get_default_line_list()

        # add lines
        add_lines(grid_name, synthesizer_data_dir, normalisation, lines_to_include)

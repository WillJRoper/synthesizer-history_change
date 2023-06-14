
"""
This reads in a cloudy grid of models

"""

from scipy import integrate
import os
import shutil
from unyt import eV
import argparse
import numpy as np
import h5py
import yaml

from synthesizer.utils import read_params
from synthesizer.cloudy import read_wavelength, read_continuum, read_lines
from synthesizer.sed import calculate_Q

from utils import get_grid_properties





def get_grid_properties_hf(hf, verbose=True):

    axes = hf.attrs['grid_axes'] # list of axes
    axes_values = {axis: hf[axis][:] for axis in axes} # dictionary of axis grid points
    
    # Get the properties of the grid including the dimensions etc.
    return get_grid_properties(axes, axes_values, verbose=verbose)



def check_cloudy_runs(grid_name, synthesizer_data_dir, replace=False, include_spectra = True):
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
        axes, n_axes, shape, n_models, mesh, model_list, index_list = get_grid_properties_hf(hf)

        # list of failed models
        failed_list = []

        for i, grid_params_ in enumerate(model_list):

            infile = f"{synthesizer_data_dir}/cloudy/{grid_name}/{i}"

            failed = False

            # check if files exist
            if include_spectra:
                if not os.path.isfile(infile+'.cont'):  # attempt to open run.
                    print('failed cont')
                    failed = True   
            if not os.path.isfile(infile+'.lines'):  # attempt to open run.
                failed = True
                print('failed lines')
            
            # if they exist also check they have size >0
            if not failed:
                if include_spectra:
                    if os.path.getsize(infile+'.cont') < 1000:
                        print('failed cont size')
                        failed = True
                if os.path.getsize(infile+'.lines') < 1000:
                    failed = True
                    print('failed lines size')

            if failed:
                print(i, model_list[i])    
                failed_list.append(i)

                # if replace is specified, instead replace the grid point
                if replace:
                    shutil.copyfile(f"{synthesizer_data_dir}/cloudy/{grid_name}/{i-1}.lines", infile+'.lines')
                    if include_spectra:
                        shutil.copyfile(f"{synthesizer_data_dir}/cloudy/{grid_name}/{i-1}.cont", infile+'.cont')
                    
        if replace:
            failed_list = []


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
        axes, n_axes, shape, n_models, mesh, model_list, index_list = get_grid_properties_hf(hf)

        # read first spectra from the first grid point to get length and wavelength grid
        lam = read_wavelength(f"{synthesizer_data_dir}/cloudy/{grid_name}/0")

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
        spectra['normalisation'] = np.zeros(shape)

        for i, indices in enumerate(index_list):

            indices = tuple(indices)

            # define the infile
            infile = f"{synthesizer_data_dir}/cloudy/{grid_name}/{i}"

            # read the continuum file containing the spectra 
            spec_dict = read_continuum(infile, return_dict=True)

            # for an arbitrary grid, we should normalise by the bolometric luminosity of the incident spectra
            norm = np.trapz(spec_dict['incident'][::-1], x=nu[::-1])

            # save normalisation for later use (rescaling lines)
            spectra['normalisation'][indices] = norm

            # save the normalised spectrum to the correct grid point 
            for spec_name in spec_names:
                spectra[spec_name][indices] = spec_dict[spec_name] / norm


def get_default_line_list(interesting=True):

    with open('default_lines.dat') as f:
        line_list = f.read().splitlines()

    if interesting:

        with open('interesting_lines.dat') as f:
            line_list += f.read().splitlines()

    return line_list





def add_lines(grid_name, synthesizer_data_dir, lines_to_include, include_spectra = True):
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
        axes, n_axes, shape, n_models, mesh, model_list, index_list = get_grid_properties_hf(hf)

        # delete lines group if it already exists
        if 'lines' in hf:
            del hf['lines']

        # define spectra
        spectra = hf['spectra']
        normalisation = spectra['normalisation'][:]
        lam = spectra['wavelength'][:]

        # create group for holding lines
        lines = hf.create_group('lines')
        # lines.attrs['lines'] = list(lines_to_include)  # save list of spectra as attribute

        # set up output arrays
        for line_id in lines_to_include:
            lines[f'{line_id}/luminosity'] = np.zeros(shape)
            lines[f'{line_id}/intrinsic_luminosity'] = np.zeros(shape)
            if include_spectra:
                lines[f'{line_id}/stellar_continuum'] = np.zeros(shape)
                lines[f'{line_id}/nebular_continuum'] = np.zeros(shape)
                lines[f'{line_id}/continuum'] = np.zeros(shape)

        for i, indices in enumerate(index_list):

            # convert indices array to tuple
            indices = tuple(indices)

            # define the infile
            infile = f"{synthesizer_data_dir}/cloudy/{grid_name}/{i}"

            # get TOTAL continuum spectra
            nebular_continuum = spectra['nebular'][indices] - spectra['linecont'][indices]
            continuum = spectra['transmitted'][indices] + nebular_continuum

            # get line quantities
            id, blend, wavelength, intrinsic, emergent = read_lines(infile)

            # identify lines we want to keep
            s = np.nonzero(np.in1d(id, np.array(lines_to_include)))[0]

            for id_, wavelength_, emergent_, intrinsic_ in zip(id[s], wavelength[s], emergent[s], intrinsic[s]):

                line = lines[id_]

                # save line wavelength
                line.attrs['wavelength'] = wavelength_

                # calculate line luminosity and save it. Uses normalisation from spectra.
                line['luminosity'][indices] = 10**(emergent_)/normalisation[indices]  # erg s^-1
                line['intrinsic_luminosity'][indices] = 10**(intrinsic_)/normalisation[indices]  # erg s^-1
                
                if include_spectra:

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
    parser.add_argument("-grid_name", "--grid_name", type=str, required=True)
    parser.add_argument("-replace", "--replace", type=bool, default=False, required=False)
    parser.add_argument("-include_spectra", "--include_spectra", type=bool, default=True)

    args = parser.parse_args()

    include_spectra = args.include_spectra

    synthesizer_data_dir = args.synthesizer_data_dir
    grid_name = args.grid_name

    # check cloudy runs
    failed_list = check_cloudy_runs(grid_name, synthesizer_data_dir, replace = args.replace, include_spectra = include_spectra)

    print(failed_list)

    # if failed prompt to re-run
    if len(failed_list)>0:

        print(f'ERROR: {len(failed_list)} cloudy runs have failed. You should re-run these with command:')
        print(f'  qsub -t 1:{len(failed_list)}  run_grid.job')

        # replace input_names with list of failed runs
        with open(f"{synthesizer_data_dir}/cloudy/{grid_name}/input_names.txt", "w") as myfile:
            myfile.write('\n'.join(map(str, failed_list)))

    #if not failed, go ahead and add spectra and lines
    else:
        
        print('- passed checks')

        # add spectra
        if include_spectra:
            add_spectra(grid_name, synthesizer_data_dir)
            print('- spectra added')

        # get list of lines
        lines_to_include = get_default_line_list()

        # add lines
        add_lines(grid_name, synthesizer_data_dir, lines_to_include)
        print('- lines added')
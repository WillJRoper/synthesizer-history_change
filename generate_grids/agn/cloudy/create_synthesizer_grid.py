
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

        grid_axes = hf.attrs['grid_axes']
        print(grid_axes)

        n_axes = len(grid_axes)

        grid = np.array(np.meshgrid(*[np.array(hf[grid_axis][:]) for grid_axis in grid_axes]))

        # determine number of models
        N = 1
        for dim in np.shape(grid): N *= dim
        N /= len(grid)
        N = int(N)

        print(f'number of models to run: {N}')

        grid_list = grid.T.reshape(N, n_axes)

        print(grid_list)

        failed_list = []

        for i, grid_params_ in enumerate(grid_list):

            infile = f'{synthesizer_data_dir}/cloudy/{grid_name}/{i}'

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

        # --- short hand for later
        metallicities = hf['metallicities']
        log10ages = hf['log10ages']

        nZ = len(metallicities)  # number of metallicity grid points
        na = len(log10ages)  # number of age grid points

        # --- read first spectra from the first grid point to get length and wavelength grid
        lam = read_wavelength(f'{synthesizer_data_dir}/cloudy/{grid_name}/0_0')

        spectra = hf.create_group('spectra')  # create a group holding the spectra in the grid file
        spectra.attrs['spec_names'] = spec_names  # save list of spectra as attribute

        spectra['wavelength'] = lam  # save the wavelength

        nlam = len(lam)  # number of wavelength points

        # --- make spectral grids and set them to zero
        for spec_name in spec_names:
            spectra[spec_name] = np.zeros((na, nZ, nlam))

        # --- now loop over meallicity and ages

        dlog10Q = np.zeros((na, nZ))

        for iZ, Z in enumerate(metallicities):

            # print(f'{iZ+1}/{len(metallicities)}')

            for ia, log10age in enumerate(log10ages):

                infile = f'{synthesizer_data_dir}/cloudy/{grid_name}/{ia}_{iZ}'

                spec_dict = read_continuum(infile, return_dict=True)
                for spec_name in spec_names:
                    spectra[spec_name][ia, iZ] = spec_dict[spec_name]

                # --- we need to rescale the cloudy spectra to the original SPS spectra. This is done here using the ionising photon luminosity, though could in principle by done at a fixed wavelength.

                # calculate log10Q_HI

                Q = calculate_Q(lam, spectra['incident'][ia, iZ, :],
                                ionisation_energy=13.6 * eV)

                dlog10Q[ia, iZ] = hf['log10Q/HI'][ia, iZ] - np.log10(Q)

                # renormalise each spectra
                for spec_name in spec_names:
                    spectra[spec_name][ia, iZ] *= 10**dlog10Q[ia, iZ]

        return dlog10Q


def get_default_line_list(interesting=True):

    with open('default_lines.dat') as f:
        line_list = f.read().splitlines()

    if interesting:

        with open('interesting_lines.dat') as f:
            line_list += f.read().splitlines()

    return line_list


def get_line_list(grid_name, synthesizer_data_dir, threshold_line='H 1 4862.69A', relative_threshold=2.0):
    """
    Get a list of lines meeting some threshold at the reference age and metallicity

    NOTE: the updated base grid HDF5 file must have been created first.
    NOTE: changing the threshold to 2.5 doubles the number of lines and produces repeats that will need to be merged.

    Parameters
    ----------
    grid_name : str
        Name of the grid.
    synthesizer_data_dir : str
        Directory where synthesizer data is kept.
    threshold : float
        The log threshold relative to Hbeta for which lines should be kept.
        Default = 2.0 which implies L > 0.01 * L_Hbeta

    Returns
    ----------
    list
        list of the lines meeting the threshold criteria
    """

    # open the new grid
    with h5py.File(f'{synthesizer_data_dir}/grids/{grid_name}.hdf5', 'a') as hf:

        # get the reference metallicity and age grid point
        ia = hf.attrs['ia_ref']
        iZ = hf.attrs['iZ_ref']

    reference_filename = f'{synthesizer_data_dir}/cloudy/{grid_name}/{ia}_{iZ}'

    cloudy_ids, blends, wavelengths, intrinsic, emergent = read_lines(
        reference_filename)

    threshold = emergent[cloudy_ids == threshold_line] - relative_threshold

    s = (emergent > threshold) & (blends == False) & (wavelengths < 50000)

    line_list = cloudy_ids[s]

    # print(f'number of lines: {np.sum(s)}')
    # print(line_list)
    # print(len(list(set(line_list))))

    return line_list


def add_lines(grid_name, synthesizer_data_dir, dlog10Q, lines_to_include):
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

        # --- short hand for later
        metallicities = hf['metallicities']
        log10ages = hf['log10ages']
        spectra = hf['spectra']
        lam = spectra['wavelength']

        nZ = len(metallicities)  # number of metallicity grid points
        na = len(log10ages)  # number of age grid points

        # -- get list of lines

        lines = hf.create_group('lines')
        # lines.attrs['lines'] = list(lines_to_include)  # save list of spectra as attribute

        for line_id in lines_to_include:
            lines[f'{line_id}/luminosity'] = np.zeros((na, nZ))
            lines[f'{line_id}/stellar_continuum'] = np.zeros((na, nZ))
            lines[f'{line_id}/nebular_continuum'] = np.zeros((na, nZ))
            lines[f'{line_id}/continuum'] = np.zeros((na, nZ))

        for iZ, Z in enumerate(metallicities):

            for ia, log10age in enumerate(log10ages):

                infile = f'{synthesizer_data_dir}/cloudy/{grid_name}/{ia}_{iZ}'

                # --- get TOTAL continuum spectra
                nebular_continuum = spectra['nebular'][ia, iZ] - spectra['linecont'][ia, iZ]
                continuum = spectra['transmitted'][ia, iZ] + nebular_continuum

                # --- get line quantities
                id, blend, wavelength, intrinsic, emergent = read_lines(infile)

                s = np.nonzero(np.in1d(id, np.array(lines_to_include)))[0]

                for id_, wavelength_, emergent_ in zip(id[s], wavelength[s], emergent[s]):

                    line = lines[id_]

                    line.attrs['wavelength'] = wavelength_

                    line['luminosity'][ia, iZ] = 10**(emergent_ + dlog10Q[ia, iZ])  # erg s^-1
                    line['stellar_continuum'][ia, iZ] = np.interp(
                        wavelength_, lam, spectra['transmitted'][ia, iZ])  # erg s^-1 Hz^-1
                    line['nebular_continuum'][ia, iZ] = np.interp(
                        wavelength_, lam, nebular_continuum)  # erg s^-1 Hz^-1
                    line['continuum'][ia, iZ] = np.interp(
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

    # if failed prompt to re-run
    if len(failed_list)>0:

        print(f'ERROR: {len(failed_list)} cloudy runs have failed. You should re-run these with command:')
        print(f'  qsub -t 1:{len(failed_list)}  run_grid.job')

        # replace input_names with list of failed runs
        with open(f"{output_dir}/input_names.txt", "a") as myfile:
            myfile.write('\n'.join(failed_list))

    # if not failed, go ahead and add spectra and lines
    # else:
        
    #     # add spectra
    #     dlog10Q = add_spectra(grid_name, synthesizer_data_dir)

    #     # add lines
    #     lines_to_include = get_default_line_list()
    #     add_lines(grid_name, synthesizer_data_dir, dlog10Q, lines_to_include)

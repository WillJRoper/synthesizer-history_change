"""
Given an SED (from an SPS model, for example), generate a cloudy atmosphere
grid. Can optionally generate an array of input files for selected parameters.
"""

import yaml
import numpy as np
from scipy import integrate
from pathlib import Path
import argparse
from unyt import c, h, angstrom, eV, erg, s, Hz, unyt_array
import h5py

from synthesizer.agn import Feltre16
from synthesizer.cloudy import create_cloudy_input, ShapeCommands
from synthesizer.abundances import Abundances
from write_submission_script import (apollo_submission_script,
                                     cosma7_submission_script)






def load_cloudy_parameters(param_file='default.yaml', default_param_file='default.yaml'):
    """
    Load CLOUDY parameters from a YAML file

    Parameters
    ----------
    param_file : str
        Location of YAML file.
    default_param_file : str
        Location of default_parameter_file

    Returns
    -------
    dict
        Dictionary of cloudy parameters
    """

    cloudy_params = {}

    # open paramter file
    with open(param_file, "r") as stream:
        try:
            cloudy_params_ = yaml.safe_load(stream)
        except yaml.YAMLError as exc:
            print(exc)

    # open default parameter file
    with open(default_param_file, "r") as stream:
        try:
            default_cloudy_params = yaml.safe_load(stream)
        except yaml.YAMLError as exc:
            print(exc)

    # combine default parameters with 
    cloudy_params = default_cloudy_params | cloudy_params_

    return cloudy_params


if __name__ == "__main__":

    parser = argparse.ArgumentParser(description='Run a grid of Cloudy AGN models')
    parser.add_argument("-machine", type=str, required=True) # machine (for submission script generation)
    parser.add_argument("-synthesizer_data_dir", type=str, required=True) # path to synthesizer_data_dir
    parser.add_argument("-grid_name", type=str, required=True) # grid_name
    parser.add_argument("-grid_params", type=str, required=False, default='grid4D.yaml') # cloudy parameter file
    parser.add_argument("-cloudy_params", type=str, required=False, default='default.yaml') # cloudy parameter file
    parser.add_argument("-cloudy", type=str, required=True) # path to cloudy executable
    parser.add_argument("-dry_run", type=bool, required=False, default=False) # boolean for dry run
    args = parser.parse_args()

    grid_name = args.grid_name
    machine = args.machine
    output_dir = f"{args.synthesizer_data_dir}/agn/{grid_name}"
    cloudy = args.cloudy

    # load cloudy parameters
    cloudy_params = load_cloudy_parameters(param_file = args.cloudy_params)

    # read in the grid parameters
    with open(args.grid_params, "r") as stream:
        try:
            grid_params = yaml.safe_load(stream)
        except yaml.YAMLError as exc:
            print(exc)

    print(machine)
    print(grid_name)
    print(output_dir)
    print(cloudy)
    print(cloudy_params)
    print(grid_params)
    
    
    grid_axes = list(grid_params.keys())
    print(grid_axes)

    n_axes = len(grid_axes)

    grid = np.array(np.meshgrid(*[np.array(v) for k,v in grid_params.items()]))

    # determine number of models
    N = 1
    for dim in np.shape(grid): N *= dim
    N /= len(grid)
    N = int(N)

    print(f'number of models to run: {N}')

    grid_list = grid.T.reshape(N, n_axes)

    print(grid_list)



    if not args.dry_run:

        # create path for cloudy runs
        Path(output_dir).mkdir(parents=True, exist_ok=True)

        # for submission system output files
        Path(f'{output_dir}/output').mkdir(parents=True, exist_ok=True)

        
        # open the new grid
        with h5py.File(f'{args.synthesizer_data_dir}/grids/agn_{grid_name}.hdf5', 'w') as hf:

            # add attribute with the grid axes
            hf.attrs['grid_axes'] = grid_axes

            # add the grid axis values
            for grid_axis in grid_axes:
                hf[grid_axis] = grid_params[grid_axis]


    for i, grid_params_ in enumerate(grid_list):
    
        grid_params = dict(zip(grid_axes, grid_params_))

        params = cloudy_params | grid_params

        print(i, params)

        if not args.dry_run:

            abundances = Abundances(10**params['log10Z'])

            # generate incident spectra
            lam = np.arange(1, 20000, 1) * angstrom

            lnu = Feltre16.intrinsic(lam, alpha=params['alpha'])

            # create shape commands
            shape_commands = ShapeCommands.table_sed(str(i), lam, lnu,  output_dir=output_dir)
            
            # create input file
            create_cloudy_input(str(i), shape_commands, abundances, output_dir=output_dir, **params)

            # write out input file
            with open(f"{output_dir}/input_names.txt", "a") as myfile:
                myfile.write(f'{i}\n')

    if machine == 'apollo':
        apollo_submission_script(N, output_dir, cloudy)
    elif machine == 'cosma7':
        cosma7_submission_script(N, output_dir, cloudy,
                                cosma_project='cosma7',
                                cosma_account='dp004')

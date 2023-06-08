"""
Given an SED (from an SPS model, for example), generate a cloudy atmosphere
grid. Can optionally generate an array of input files for selected parameters.
"""

import numpy as np
from scipy import integrate
from pathlib import Path
import argparse
from unyt import c, h, angstrom, eV, erg, s, Hz, unyt_array
import h5py

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
    with open(grid_params, "r") as stream:
        try:
            grid_params = yaml.safe_load(stream)
        except yaml.YAMLError as exc:
            print(exc)



    if args.dry_run:
        print(machine)
        print(grid_name)
        print(output_dir)
        print(cloudy)
        print(cloudy_params)

    if not args.dry_run:

        # create path for cloudy runs
        Path(output_dir).mkdir(parents=True, exist_ok=True)

        # for submission system output files
        Path(f'{output_dir}/output').mkdir(parents=True, exist_ok=True)

        # log10TBB grid
        log10TBBs = np.arange(4., 7., 0.2)

        # metallicity grid
        log10Zs = np.arange(-5., -1., 0.5)

        # log10U
        log10Us = np.array([-4., -3, -2, -1., 0., 1.])

        # log10nH
        log10n_Hs = np.array([0.,1.,2.,3.,4.])

        # total number of models
        N = len(log10Ts)*len(log10Zs)*len(log10Us)*len(log10n_Hs)

        # open the new grid
        with h5py.File(f'{synthesizer_data_dir}/grids/{grid_name}.hdf5', 'w') as hf:

            # add attribute with the grid axes
            hf.attrs['grid_axes'] = ['log10TBB', 'log10Z', 'log10U', 'log10n_H']

            
            hf['log10TBB'] = log10Ts
            hf['log10Z'] = log10Zs
            hf['log10U'] = log10Us
            hf['log10n_H'] = log10n_Hs


        for iZ, log10Z in enumerate(log10Zs):
            abundances = Abundances(10**log10Z)
            for iT, log10T in enumerate(log10Ts):
                for iU, log10U in enumerate(log10Us):
                    for inH, log10n_H in enumerate(log10n_Hs):
        
                        # cloudy_params

                        cloudy_params['log10U'] = log10U
                        cloudy_params['log10n_H'] = log10n_H

                        # specify model_name
                        model_name = f"{iT}_{iZ}_{iU}_{inH}"
            
                        # create shape commands
                        TBB = 10**log10T
                        cinputs = ShapeCommands.cloudy_agn(TBB)
                        
                        # create input file
                        create_cloudy_input(model_name, shape_commands, abundances, output_dir=output_dir, **cloudy_params)

                        # write out input file
                        with open(f"{output_dir}/input_names.txt", "a") as myfile:
                            myfile.write(f'{model_name}\n')
        
        if machine == 'apollo':
            apollo_submission_script(N, output_dir, cloudy)
        elif machine == 'cosma7':
            cosma7_submission_script(N, output_dir, cloudy,
                                    cosma_project='cosma7',
                                    cosma_account='dp004')

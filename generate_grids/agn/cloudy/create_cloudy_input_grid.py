"""
Given an SED (from an SPS model, for example), generate a cloudy atmosphere
grid. Can optionally generate an array of input files for selected parameters.
"""

import numpy as np
from scipy import integrate
from pathlib import Path
from unyt import c, h, angstrom, eV, erg, s, Hz, unyt_array
import h5py

from synthesizer.cloudy import calculate_Q_from_U
from synthesizer.abundances import Abundances
from write_submission_script import (apollo_submission_script,
                                     cosma7_submission_script)




def make_cloudy_input_agn_cloudy(TBB, aox=-1.4, auv=-0.5, ax=-1.):

    """
    Define the command intitialising the cloudy AGN model. See 6.2 Hazy1.pdf.
    """

    cinput = []
    cinput.append(f'AGN T ={TBB} k, a(ox) = {aox}, a(uv)= {auv} a(x)={ax} \n')

    return 



if __name__ == "__main__":

    grid_name = 'agn_cloudy'

    # replace with arguments
    machine = 'apollo'
    synthesizer_data_dir = "/research/astrodata/highz/synthesizer/"
    output_dir = f"{synthesizer_data_dir}/cloudy/agn_cloudy"
    cloudy = '/its/home/sw376/flare/software/cloudy/c17.03/source/cloudy.exe'

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
    log10nHs = np.array([0.,1.,2.,3.,4.])

    # total number of models
    N = len(log10Ts)*len(log10Zs)*len(log10Us)

    # open the new grid
    with h5py.File(f'{synthesizer_data_dir}/grids/{grid_name}.hdf5', 'w') as hf:

        # add attribute with the grid axes for future when using >2D grid or AGN grids
        hf.attrs['grid_axes'] = ['log10TBB', 'log10Z', 'log10U']

        
        hf['log10TBB'] = log10Ts
        hf['log10Z'] = log10Zs
        hf['log10U'] = log10Us
        hf['log10nH'] = log10Us


    for iZ, log10Z in enumerate(log10Zs):
        abundances = Abundances(10**log10Z)
        for iT, log10T in enumerate(log10Ts):
            for iU, log10U in enumerate(log10Us):
    
                model_name = f"{iT}_{iZ}_{iU}"
    
                # this will need changing
                
    
        
    
                with open(f"{output_dir}/input_names.txt", "a") as myfile:
                    myfile.write(f'{model_name}\n')
    
    if machine == 'apollo':
        apollo_submission_script(N, output_dir, cloudy)
    elif machine == 'cosma7':
        cosma7_submission_script(N, output_dir, cloudy,
                                 cosma_project='cosma7',
                                 cosma_account='dp004')

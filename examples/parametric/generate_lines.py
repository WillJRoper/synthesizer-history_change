"""
Example for generating a emission lines for a parametric galaxy. This example will:
- show the available lines to a grid
- build a parametric galaxy (see make_sfzh and make_sed)
- calculate intrinsic line properties
- calculate dust-attenuated line properties
"""
import os

from synthesizer.units import Units
from synthesizer.grid import get_available_lines, Grid
from synthesizer.parametric.sfzh import SFH, ZH, generate_sfzh
from synthesizer.parametric.galaxy import Galaxy
from unyt import yr, Myr


if __name__ == '__main__':

    # Get the location of this script, __file__ is the absolute path of this
    # script, however we just want to directory
    script_path = os.path.abspath(os.path.dirname(__file__))

    # Define the grid
    grid_name = "test_grid"
    grid_dir = script_path + "/../../tests/test_grid/"

    # to see what lines are available in a grid we can use this helper function
    available_lines = get_available_lines(grid_name, grid_dir)
    print(available_lines)

    # list of lines. Lines in nested lists (or tuples) denote doublets for which the combined line properties are calculated
    # should result in the same behaviour as above
    line_ids = ['H 1 4862.69A', 'O 3 4960.29A', 'O 3 5008.24A', ['O 3 4960.29A', 'O 3 5008.24A']]
    grid = Grid(grid_name, grid_dir=grid_dir, read_spectra=False, read_lines=line_ids)

    # --- define the parameters of the star formation and metal enrichment histories
    sfh_p = {'duration': 100 * Myr}
    Z_p = {'log10Z': -2.0}  # can also use linear metallicity e.g. {'Z': 0.01}

    # --- define the functional form of the star formation and metal enrichment histories
    sfh = SFH.Constant(sfh_p)  # constant star formation
    Zh = ZH.deltaConstant(Z_p)  # constant metallicity

    # --- get the 2D star formation and metal enrichment history for the given SPS grid. This is (age, Z).
    sfzh = generate_sfzh(grid.log10ages, grid.metallicities, sfh, Zh)

    # --- create the Galaxy object
    galaxy = Galaxy(sfzh)

    # --- print a summary of the Galaxy object
    print(galaxy)

    # --- create the Lines dictionary which contains line objects
    lines = galaxy.get_line_intrinsic(grid, line_ids)
    print('-'*50)
    print('INTRINSIC')
    for line_id, line in lines.items():
        print(line)

    # --- calculate attenuated line properties assuming uniform dust (should leave EW unchanged)
    lines = galaxy.get_line_screen(grid, line_ids, tauV=0.5)
    print('-'*50)
    print('SCREEN')
    for line_id, line in lines.items():
        print(line)

    # --- calculate attenuated line properties assuming different dust affecting stellar and nebular components
    lines = galaxy.get_line_attenuated(grid, line_ids, tauV_stellar=0.1, tauV_nebular=0.5)
    print('-'*50)
    print('ATTENUATED')
    for line_id, line in lines.items():
        print(line)

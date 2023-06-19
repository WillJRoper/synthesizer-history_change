

"""
Run a single cloudy model based on an SPS grid
"""

import os
import sys
import argparse
import numpy as np

from synthesizer.abundances import Abundances
from synthesizer.grid import Grid
from synthesizer.cloudy import create_cloudy_input, ShapeCommands




# ---- initialise abundances object
abundances = Abundances(Z=0.3*0.02) # abundances object

print(abundances)

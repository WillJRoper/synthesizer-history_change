import matplotlib.pyplot as plt
import numpy as np
from unyt import matplotlib_support
from unyt import Angstrom

from synthesizer import dust

import cmasher as cmr

models = ['power_law', 'Calzetti2000', 'GrainsWD01', 'GrainsWD01', 'GrainsWD01']
params = [{'slope': -1.}, {'slope': 0., 'x0': 0.2175, 'ampl': 1.},
          {'model': 'MW'}, {'model': 'SMC'}, {'model': 'LMC'}]

colors = cmr.take_cmap_colors('cmr.guppy', len(models))

lam = np.arange(1000, 10000, 10)*Angstrom

for ii, (model, param) in enumerate(zip(models, params)):
    emodel = getattr(dust, model)(params=param)

    plt.plot(lam, emodel.tau(lam), color=colors[ii], label=F'{model}, {param}')

plt.xlabel(r'$\lambda/(\AA)$', fontsize=12)
plt.ylabel(r'A$_{\lambda}/$A$_{V}$', fontsize=12)
plt.yticks(np.arange(0,10))
plt.xlim(np.min(lam), np.max(lam))

plt.legend(frameon=False, fontsize=10)
plt.grid()

plt.show()

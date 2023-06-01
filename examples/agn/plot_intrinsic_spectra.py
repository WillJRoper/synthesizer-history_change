


import numpy as np
import matplotlib.pyplot as plt

from synthesizer.agn import Feltre16
from synthesizer.plt import single
from unyt import Angstrom


if __name__ == '__main__':


    lam = 10**(np.arange(1., 6., 0.1)) * Angstrom

    fig, ax = single((6., 3.5))

    lnu = Feltre16.intrinsic(lam, -1.8)

    ax.plot(np.log10(lam), np.log10(lnu))

    # ax.set_xlim([3., 4.5])
    # ax.set_ylim([15, 19])
    ax.legend(fontsize = 8, labelspacing = 0.0)
    ax.set_xlabel(r'$\rm log_{10}(\lambda/\AA)$')
    ax.set_ylabel(r'$\rm log_{10}(L_{\nu}/L_{\nu}()}^1)$')

    plt.show()

    
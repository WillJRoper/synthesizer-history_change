

# --- general
import h5py
import copy
import numpy as np
from scipy import integrate
from scipy.interpolate import InterpolatedUnivariateSpline
from unyt import yr
import time


import matplotlib.pyplot as plt
from matplotlib.colors import LogNorm
import cmasher as cmr

from .. import exceptions
from ..stats import weighted_median, weighted_mean
from ..plt import single_histxy, mlabel


class BinnedSFZH:

    """
    This is a simple class for holding a binned star formation and metal
    enrichment history. This can be extended with other methods.
    """

    def __init__(self, log10ages, metallicities, sfzh, sfh_f=None, Zh_f=None,
                 parent_sfzh=None):

        # Compute the individual star formation and metallicity histories
        self.sfzh = sfzh  # 2D star formation and metal enrichment history
        self.sfh = np.sum(self.sfzh, axis=1)  # 1D star formation history
        self.Z = np.sum(self.sfzh, axis=0)  # metallicity distribution

        # Assign ancillary data for this SFZH
        if parent_sfzh is None:
            self.log10ages = log10ages
            self.ages = 10**log10ages
            self.log10ages_lims = [self.log10ages[0], self.log10ages[-1]]
            self.metallicities = metallicities
            self.metallicities_lims = [
                self.metallicities[0], self.metallicities[-1]]
            self.log10metallicities = np.log10(metallicities)
            self.log10metallicities_lims = [
                self.log10metallicities[0], self.log10metallicities[-1]]

            # Function used to generate the star formation history if given
            self.sfh_f = sfh_f

            # Function used to generate the metallicity history/distribution
            # if given
            self.Zh_f = Zh_f

        else:

            self._inherit_props(parent_sfzh)

        # Check if metallicities on regular grid in log10metallicity or
        # metallicity or not at all (e.g. BPASS)
        if len(set(self.metallicities[:-1]-self.metallicities[1:])) == 1:
            self.metallicity_grid = 'Z'
        elif len(set(self.log10metallicities[:-1]-self.log10metallicities[1:])) == 1:
            self.metallicity_grid = 'log10Z'
        else:
            self.metallicity_grid = None

    def _inherit_props(self, parent):
        """
        A helper method to inherit properties from a parent SFZH and thus
        avoid double up computations during addition.

        Parameters
        ----------
        parent: obj (BinnedSFZH)
            The parent SFZH object from which to inherit properties.
        """

        # Assign the inherited properties
        self.log10ages = parent.log10ages
        self.ages = parent.ages
        self.log10ages_lims = parent.log10ages_lims
        self.metallicities = parent.metallicities
        self.metallicities_lims = parent.metallicities_lims
        self.log10metallicities = parent.log10metallicities
        self.log10metallicities_lims = parent.log10metallicities_lims
        self.sfh_f = parent.sfh_f
        self.Zh_f = parent.Zh_f

    def calculate_median_age(self):
        """ calculate the median age """

        return weighted_median(self.ages, self.sfh) * yr

    def calculate_mean_age(self):
        """ calculate the mean age """

        return weighted_mean(self.ages, self.sfh) * yr

    def calculate_mean_metallicity(self):
        """ calculate the mean metallicity """

        return weighted_mean(self.metallicities, self.Z)

    def __str__(self):
        """ print basic summary of the binned star formation and metal enrichment history """

        pstr = ''
        pstr += '-'*10 + "\n"
        pstr += 'SUMMARY OF BINNED SFZH' + "\n"
        pstr += f'median age: {self.calculate_median_age().to("Myr"):.2f}' + "\n"
        pstr += f'mean age: {self.calculate_mean_age().to("Myr"):.2f}' + "\n"
        pstr += f'mean metallicity: {self.calculate_mean_metallicity():.4f}' + \
            "\n"
        pstr += '-'*10 + "\n"
        return pstr

    def __add__(self, second_sfzh):
        """ Add two SFZH histories together """

        if second_sfzh.sfzh.shape == self.sfzh.shape:

            new_sfzh = self.sfzh + second_sfzh.sfzh

            return BinnedSFZH(self.log10ages, self.metallicities, new_sfzh,
                              sfh_f=self.sfh_f, Zh_f=self.Zh_f,
                              parent_sfzh=self)

        else:

            exceptions.InconsistentAddition('SFZH must be the same shape')

    def plot(self, show=True, plot_sfh=False):
        """ Make a nice plots of the binned SZFH """

        # Define the figure grid
        fig = plt.figure()
        gs = fig.add_gridspec(nrows=2, ncols=2, wspace=0.0, hspace=0.0,
                              height_ratios=[1, 4], width_ratios=[4, 1])

        # Define the axes
        ax = fig.add_subplot(gs[1, 0])
        haxy = fig.add_subplot(gs[1, 1])
        haxx = fig.add_subplot(gs[0, 0])
        haxx.axis('off')
        haxy.axis('off')

        # This is technically incorrect because metallicity is not on a an
        # actual grid.
        ax.imshow(self.sfzh.T, origin='lower',
                  extent=[
                      *self.log10ages_lims,
                      self.log10metallicities[0],
                      self.log10metallicities[-1]
                  ], cmap=cmr.sunburst, aspect='auto')

        # Add binned Z to right of the plot
        haxy.fill_betweenx(self.log10metallicities, self.Z/np.max(self.Z),
                           step='mid', color='k', alpha=0.3)

        # Add binned SFH to top of the plot
        haxx.fill_between(self.log10ages, self.sfh/np.max(self.sfh),
                          step='mid', color='k', alpha=0.3)

        # Add SFR to top of the plot
        if plot_sfh:
            x = np.linspace(*self.log10ages_lims, 1000)
            y = self.sfh_f.sfr(10**x)
            haxx.plot(x, y/np.max(y))

        # Set axes limits
        haxy.set_xlim([0., 1.2])
        haxy.set_ylim(self.log10metallicities[0], self.log10metallicities[-1])
        haxx.set_ylim([0., 1.2])
        haxx.set_xlim(self.log10ages_lims)

        # Set axis labels
        ax.set_xlabel(mlabel('log_{10}(age/yr)'))
        ax.set_ylabel(mlabel('log_{10}Z'))

        if show:
            plt.show()

        return fig, ax


def generate_sfh(ages, sfh_, log10=False):

    if log10:
        ages = 10**ages

    SFH = np.zeros(len(ages))

    min_age = 0
    for ia, age in enumerate(ages[:-1]):
        max_age = int(np.mean([ages[ia+1], ages[ia]]))  #  years
        sf = integrate.quad(sfh_.sfr, min_age, max_age)[0]
        SFH[ia] = sf
        min_age = max_age

    # --- normalise
    SFH /= np.sum(SFH)

    return SFH


def generate_instant_sfzh(log10ages, metallicities, log10age, metallicity, stellar_mass=1):
    """ simply returns the SFZH where only bin is populated corresponding to the age and metallicity """

    sfzh = np.zeros((len(log10ages), len(metallicities)))
    ia = (np.abs(log10ages - log10age)).argmin()
    iZ = (np.abs(metallicities - metallicity)).argmin()
    sfzh[ia, iZ] = stellar_mass

    return BinnedSFZH(log10ages, metallicities, sfzh)


def generate_sfzh(log10ages, metallicities, sfh, Zh, stellar_mass=1.):
    """ return an instance of the BinnedSFZH class """

    ages = 10**log10ages

    sfzh = np.zeros((len(log10ages), len(metallicities)))

    # Evaluate sfrs
    sfrs = sfh.sfr_arr(ages)

    if Zh.dist == 'delta':
        f = InterpolatedUnivariateSpline(ages, sfrs, k=3)
        min_age = 0
        for ia, age in enumerate(ages[:-1]):
            max_age = (ages[ia+1] + ages[ia]) / 2  # years
            sf = f.integral(min_age, max_age)
            iZ = (np.abs(metallicities - Zh.Z(age))).argmin()
            sfzh[ia, iZ] = sf
            min_age = max_age

    if Zh.dist == 'dist':
        print('WARNING: NOT YET IMPLEMENTED')

    # --- normalise
    sfzh /= np.sum(sfzh)
    sfzh *= stellar_mass

    return BinnedSFZH(log10ages, metallicities, sfzh, sfh_f=sfh, Zh_f=Zh)


def generate_sfzh_from_array(log10ages, metallicities, sfh, Zh, stellar_mass=1.):
    """
    Generated a BinnedSFZH from an array instead of function
    """

    if not isinstance(Zh, np.ndarray):
        iZ = np.abs(metallicities - Zh).argmin()
        Zh = np.zeros(len(metallicities))
        Zh[iZ] = 1.0

    sfzh = sfh[:, np.newaxis] * Zh

    # --- normalise
    sfzh /= np.sum(sfzh)
    sfzh *= stellar_mass

    return BinnedSFZH(log10ages, metallicities, sfzh)


class ZH:

    """ A collection of classes describing the metallicity history (and distribution) """

    class Common:

        def __str__(self):
            """ print basic summary of the parameterised star formation history """

            pstr = ''
            pstr += '-'*10 + "\n"
            pstr += 'SUMMARY OF PARAMETERISED METAL ENRICHMENT HISTORY' + "\n"
            pstr += str(self.__class__) + "\n"
            for parameter_name, parameter_value in self.parameters.items():
                pstr += f'{parameter_name}: {parameter_value}' + "\n"
            pstr += '-'*10 + "\n"
            return pstr

    class deltaConstant(Common):

        """ return a single metallicity as a function of age. """

        def __init__(self, parameters):

            self.name = 'Constant'
            self.dist = 'delta'  # set distribution type
            self.parameters = parameters
            if 'Z' in parameters.keys():
                self.Z_ = parameters['Z']
                self.log10Z_ = np.log10(self.Z_)
            elif 'log10Z' in parameters.keys():
                self.log10Z_ = parameters['log10Z']
                self.Z_ = 10**self.log10Z_

        def Z(self, age):
            return self.Z_

        def log10Z(self, age):
            return self.log10Z_


class SFH:

    """ A collection of classes describing parametric star formation histories """

    class Common:

        def sfr(self, age):

            if isinstance(age, np.ndarray) | isinstance(age, list):
                return np.array([self.sfr_(a) for a in age])
            else:
                return self.sfr_(age)

        def calculate_sfh(self, t_range=[0, 1E10], dt=1E6):
            """ calcualte the age of a given star formation history """

            t = np.arange(*t_range, dt)
            sfh = self.sfr(t)
            return t, sfh

        def calculate_median_age(self, t_range=[0, 1E10], dt=1E6):
            """ calcualte the median age of a given star formation history """

            t, sfh = self.calculate_sfh(t_range=t_range, dt=dt)

            return weighted_median(t, sfh) * yr

        def calculate_mean_age(self, t_range=[0, 1E10], dt=1E6):
            """ calcualte the median age of a given star formation history """

            t, sfh = self.calculate_sfh(t_range=t_range, dt=dt)

            return weighted_mean(t, sfh) * yr

        def calculate_moment(self, n):
            """ calculate the n-th moment of the star formation history """

            print('WARNING: not yet implemnted')
            return

        def __str__(self):
            """ print basic summary of the parameterised star formation history """

            pstr = ''
            pstr += '-'*10 + "\n"
            pstr += 'SUMMARY OF PARAMETERISED STAR FORMATION HISTORY' + "\n"
            pstr += str(self.__class__) + "\n"
            for parameter_name, parameter_value in self.parameters.items():
                pstr += f'{parameter_name}: {parameter_value}' + "\n"
            pstr += f'median age: {self.calculate_median_age().to("Myr"):.2f}' + "\n"
            pstr += f'mean age: {self.calculate_mean_age().to("Myr"):.2f}' + "\n"
            pstr += '-'*10 + "\n"
            return pstr

    class Constant(Common):

        """
        A constant star formation history
            sfr = 1; t<=duration
            sfr = 0; t>duration
        """

        def __init__(self, parameters):
            self.name = 'Constant'
            self.parameters = parameters
            self.duration = self.parameters['duration'].to('yr').value

        def sfr_(self, age):
            if age <= self.duration:
                return 1.0
            else:
                return 0.0

        def sfr_arr(self, age):
            """ age is lookback time """

            okinds = age < self.duration

            sfr = np.zeros(age.shape)
            sfr[okinds] = 1.0

            return sfr

    class Exponential(Common):

        """
        An exponential star formation history
        """

        def __init__(self, parameters):
            self.name = 'Exponential'
            self.parameters = parameters
            self.tau = self.parameters['tau'].to('yr').value

        def sfr_(self, age):

            return np.exp(-age/self.tau)

    class TruncatedExponential(Common):

        """
        A truncated exponential star formation history
        """

        def __init__(self, parameters):
            self.name = 'Truncated Exponential'
            self.parameters = parameters
            self.tau = self.parameters['tau'].to('yr').value

            if 'max_age' in self.parameters:
                self.max_age = self.parameters['max_age'].to('yr').value
            else:
                self.max_age = self.parameters['duration'].to('yr').value

        def sfr_(self, age):

            if age < self.max_age:
                return np.exp(-age/self.tau)
            else:
                return 0.0

    class SteveLogNormal(Common):
        """
        A log-normal star formation history
        """

        def __init__(self, parameters):
            self.name = 'Log Normal'
            self.parameters = parameters
            self.peak_age = self.parameters['peak_age'].to('yr').value
            self.tau = self.parameters['tau']
            self.max_age = self.parameters['max_age'].to('yr').value

            self.tpeak = self.max_age - self.peak_age
            self.T0 = np.log(self.tpeak) + self.tau ** 2

        def sfr_(self, age):
            """ age is lookback time """

            if age < self.max_age:
                return (1./(self.max_age-age))*np.exp(-(np.log(self.max_age-age)-self.T0)**2/(2*self.tau**2))
            else:
                return 0.0

        def sfr_arr(self, age):
            """ age is lookback time """

            okinds = age < self.max_age

            sfr = np.zeros(age.shape)
            sfr[okinds] = (1./(self.max_age-age[okinds])) * \
                np.exp(-(np.log(self.max_age -
                                age[okinds])-self.T0)**2/(2*self.tau**2))

            return sfr

    class LogNormal(Common):
        """
        A log-normal star formation history
        """

        def __init__(self, parameters):
            self.name = 'Log Normal'
            self.parameters = parameters
            self.peak_age = self.parameters['peak_age'].to('yr').value
            self.tau = self.parameters['tau']
            self.max_age = self.parameters['max_age'].to('yr').value

            self.tpeak = self.max_age - self.peak_age
            self.T0 = np.log(self.tpeak) + self.tau ** 2

        def sfr_(self, age):
            """ age is lookback time """

            if age < self.max_age:
                return (1./(self.max_age-age))*np.exp(-(np.log(self.max_age-age)-self.T0)**2/(2*self.tau**2))
            else:
                return 0.0

        def sfr_arr(self, age):
            """ age is lookback time """

            okinds = age < self.max_age

            sfr = np.zeros(age.shape)
            sfr[okinds] = (1./(self.max_age-age[okinds])) * \
                np.exp(-(np.log(self.max_age -
                                age[okinds])-self.T0)**2/(2*self.tau**2))

            return sfr

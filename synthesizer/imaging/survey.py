"""
Survey functionality 
"""
import numpy as np
import synthesizer.exceptions as exceptions
from synthesizer.imaging import images, spectral_cubes
from synthesizer.galaxy.particle import ParticleGalaxy
from synthesizer.galaxy.parametric import ParametricGalaxy
from synthesizer.utils import m_to_flux, flux_to_luminosity


class Instrument:
    """
    This class describes an instrument used to make a set of observations.

    Attributes
    ----------

    Methods
    -------

    """

    def __init__(self, resolution, filters, psfs=None, depths=None,
                 aperture=None, snrs=None, noises=None,
                 resolving_power=None, lam=None):
        """
        Initialise the Observatory.

        Parameters
        ----------

        """

        # Basic metadata
        self.instrument = filters.filter_codes[0].split(".")[0]

        # Store some basic instrument properties
        self.resolution = resolution
        self.filters = filters
        self.psfs = psfs

        # Store some basic spectral information for this observation.
        self.resolving_power = resolving_power
        self.lams = None

        # Intilaise noise properties which can be populated by the outside.
        self.aperture = aperture
        self.depths = depths
        self.noises = noises
        self.snrs = snrs

        # Unit information
        if resolution is not None:
            self.spatial_unit = resolution.units
        else:
            self.spatial_unit = None

    def _check_obs_args(self):
        """
        Ensures we have valid inputs.

        Parameters
        ----------

        Raises
        ------

        """
        pass

    def get_lam_from_R(self):
        """
        Calculates the wavelengths of a spectrum based on this observations
        resolving power.

        Parameters
        ----------

        Raises
        ------

        """
        pass


class Survey:
    """
    A Survey helper object which defines all the properties of a survey. This enables the defintion of a survey approach and then the production of photometry and/or images based on that instruments/filters making up the survey. This can handle PSFs and depths which vary across bands and intruments.


    Attributes
    ----------

    Methods
    -------

    """

    def __init__(self, galaxies=(), fov=None, super_resolution_factor=2):
        """
        Initialise the Survey.

        Parameters
        ----------

        """

        # Basic information
        self.ninstruments = 0
        self.nfilters = 0
        self.ngalaxies = 0

        # Information about the field/collection of images being observered
        self.fov = fov
        self.super_resolution_factor = super_resolution_factor

        # Observation configurations are held in a dict, initialise it.
        self.instruments = {}

        # Store the galaxies we are making images of
        self.galaxies = galaxies

        # Intialise somewhere to keep survey images, this is populated later
        self.imgs = None

        # Initialise somewhere to keep galaxy photometry. This is the
        # integrated flux/luminosity in each band in Survey.filters.
        self.photometry = {}

    def _check_survey_args(self):
        """
        Ensures we have valid inputs.

        Parameters
        ----------

        Raises
        ------

        """
        pass

    def add_photometric_instrument(self, filters, label, resolution=None,
                                   psfs=None, depths=None, apertures=None,
                                   snrs=None, noises=None):
        """
        Adds an instrument and all it's filters to the Survey.

        Parameters
        ----------

        Raises
        ------
        InconsistentArguments
            If the arguments do not constitute a valid combination for an
            instrument an error is thrown.
        """

        # How many filters do we have?
        nfilters = len(filters)

        # Check our inputs match the number of filters we have
        if isinstance(psfs, dict):
            if nfilters != len(psfs):
                raise exceptions.InconsistentArguments(
                    "Inconsistent number of entries in instrument dictionaries"
                    " len(filters)=%d, len(psfs)=%d)" % (nfilters, len(psfs))
                )
        if isinstance(depths, dict):
            if nfilters != len(depths):
                raise exceptions.InconsistentArguments(
                    "Inconsistent number of entries in instrument dictionaries"
                    " len(filters)=%d, len(depths)=%d)" % (nfilters,
                                                           len(depths))
                )
        if isinstance(apertures, dict):
            if nfilters != len(apertures):
                raise exceptions.InconsistentArguments(
                    "Inconsistent number of entries in instrument dictionaries"
                    " len(filters)=%d, len(apertures)=%d)" % (nfilters,
                                                              len(apertures))
                )
        if isinstance(snrs, dict):
            if nfilters != len(snrs):
                raise exceptions.InconsistentArguments(
                    "Inconsistent number of entries in instrument dictionaries"
                    " len(filters)=%d, len(snrs)=%d)" % (nfilters,
                                                         len(snrs))
                )
        if isinstance(noises, dict):
            if nfilters != len(noises):
                raise exceptions.InconsistentArguments(
                    "Inconsistent number of entries in instrument dictionaries"
                    " len(filters)=%d, len(noises)=%d)" % (nfilters,
                                                           len(noises))
                )

        # Create this observation configurations
        self.instruments[label] = Instrument(
            resolution=resolution, filters=filters, psfs=psfs,
            depths=depths, aperture=apertures, snrs=snrs, noises=noises
        )

        # Record that we included another insturment and count the filters
        self.ninstruments += 1
        self.nfilters += len(filters)

    def add_spectral_instrument(self, resolution, resolving_power,
                                psf=None, depth=None, aperture=None):
        pass

    def add_galaxies(self, galaxies):
        """
        Adds galaxies to this survey

        Parameters
        ----------
        galaxies : list
            The galaxies to include in this Survey.

        """

        # If we have no galaxies just add them
        if len(self.galaxies) == 0:
            self.galaxies = galaxies

        # Otherwise, we have to add them on to what we have, handling whether
        # we are adding 1 galaxy...
        elif (len(self.galaxies) > 0 and
              (isinstance(galaxies, ParticleGalaxy) or
               isinstance(galaxies, ParametricGalaxy))):

            # Double check galaxies is a list
            self.galaxies = list(self.galaxies)

            # Include the new galaxies
            self.galaxies.append(galaxies)

        # ... or multiple galaxies
        else:

            # Double check galaxies is a list
            self.galaxies = list(self.galaxies)

            # Include the new galaxies
            self.galaxies.extend(galaxies)

        # Count how many galaxies we have
        self.ngalaxies = len(self.galaxies)

    def convert_mag_depth_to_fnu0(self, redshift, cosmo):
        """
        Converts depths defined in absolute magnitude to the units of
        luminosity (erg / s /Hz).

        This is a helpful wrapper to handle the use of different terminology in
        the SED object.

        Parameters
        ----------
        redshift : float
            The redshift of the observation.
        cosmo : obj (astropy.cosmology)
            The cosmology object used to do the cosmological calculations
            involved in the conversion.
        """
        self.convert_mag_depth_to_lnu(redshift, cosmo)

    def convert_mag_depth_to_lnu(self, redshift, cosmo):
        """
        Converts depths defined in apparent magnitude to the units of
        luminosity (erg / s /Hz).

        Parameters
        ----------
        redshift : float
            The redshift of the observation.
        cosmo : obj (astropy.cosmology)
            The cosmology object used to do the cosmological calculations
            involved in the conversion.
        """

        # Convert the depths, looping over them if we have to.
        for inst in self.instruments:
            if isinstance(self.instruments[inst].depths, dict):
                for key in self.instruments[inst].depths:
                    flux = m_to_flux(self.instruments[inst].depths[key])
                    self.instruments[inst].depths[key] = flux_to_luminosity(
                        flux, cosmo, redshift
                    )
            else:
                flux = m_to_flux(self.instruments[inst].depths)
                self.instruments[inst].depths = flux_to_luminosity(
                    flux, cosmo, redshift
                )

    def convert_mag_depth_to_fnu(self):
        """
        Converts depths defined in apparent magnitude to the units of
        flux (nJy).
        """

        # Convert the depths, looping over them if we have to.
        for inst in self.instruments:
            if isinstance(self.instruments[inst].depths, dict):
                for key in self.instruments[inst].depths:
                    self.instruments[inst].depths[key] = m_to_flux(
                        self.instruments[inst].depths[key]
                    )
            else:
                self.instruments[inst].depths = m_to_flux(
                    self.instruments[inst].depths
                )

    def get_photometry(self, spectra_type, cosmo=None, redshift=None, igm=None):
        """

        Parameters
        ----------

        """

        # We need to handle whether different types of spectra exist.
        if spectra_type == "intrinsic":
            pass
        elif spectra_type == "stellar":
            pass
        elif spectra_type == "attenuated":
            raise exceptions.UnimplementedFunctionality(
                "Attenuated spectra coming soon!"
            )
        else:
           # TODO: make a UnknownSpectralType error
            raise exceptions.InconsistentArguments(
                "Unrecognised spectra_type!")

        # Loop over each instrument
        for key in self.instruments:

            # Loop over filters in this instrument
            for f in self.instruments[key].filters:

                # Make an entry in the photometry dictionary for this filter
                self.photometry[f.filter_code] = np.zeros(self.ngalaxies)

                # Loop over each galaxy
                for igal in range(self.ngalaxies):

                    if spectra_type not in self.galaxies[igal].spectra:
                        pass

                    # Are we getting flux or luminosity?
                    sed = self.galaxies[igal].spectra[spectra_type]
                    if cosmo is not None:
                        if sed.fnu is None:
                            sed.get_fnu(cosmo, redshift, igm)
                        spectra = sed.fnu
                        lams = sed.lamz
                    else:
                        spectra = sed.lnu
                        lams = sed.lam

                    # Calculate the photometry in this band
                    phot = f.apply_filter(spectra, lams)

                    # Store this photometry
                    self.photometry[f.filter_code][igal] = phot

        return self.photometry

    def make_field_image(self, centre):
        """

        Parameters
        ----------

        """
        pass

    def make_images(self, img_type, spectra_type, kernel_func=None,
                    rest_frame=False, cosmo=None, igm=None):
        """

        Parameters
        ----------

        """

        # Make a dictionary in which to store our image objects, within
        # this dictionary imgs are stored in list ordered by galaxy.
        self.imgs = {}

        # Loop over instruments and make images for each galaxy using each
        # instrument
        for key in self.instruments:

            # Extract the instrument
            inst = self.instruments[key]

            # Create entry in images dictionary
            self.imgs[inst] = []

            # Loop over galaxies
            for gal in self.galaxies:

                # Get images of this galaxy with this instrument
                img = gal.make_image(
                    inst.resolution, fov=self.fov, img_type=img_type,
                    sed=gal.spectra_array[spectra_type], filters=inst.filters,
                    psfs=inst.psfs, depths=inst.depths,  aperture=inst.aperture,
                    snrs=inst.snrs, kernel_func=kernel_func,
                    rest_frame=rest_frame, cosmo=cosmo, igm=igm,
                    super_resolution_factor=self.super_resolution_factor
                )

                # Store this result
                self.imgs[inst].append(img)

        return self.imgs

    def make_field_ifu(self, centre):
        """

        Parameters
        ----------

        """
        pass

    def make_ifus(self):
        """

        Parameters
        ----------

        """
        pass

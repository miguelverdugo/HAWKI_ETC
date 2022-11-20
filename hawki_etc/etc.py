
import numpy as np
import astropy.units as u

from spextra import Spextrum, Passband

import scopesim as sim

from .utils import *

#from scopesim_templates.basic.basic import empty_sky
#from scopesim_templates.b
"""
Expected usage:

initialize

etc = ETC()



"""

SED_dict = dict(template=template,
                MARCS=MARCS,
                emission_line=emission_line,
                black_body=black_body,
                powerlaw=powerlaw,
                uniform=flat_spec)

FLUX_DISTRO_dict = dict(sersic=sersic,
                        point_source=point_source,
                        uniform=uniform_flux)






#class ETC:
#    common class for setting parameters
#class ETC_HAWKI(ETC)
#class ETC_MICADO(ETC)
#class ETC_METIS(ETC)


class ETC_base:
    """
    Base class for the common interface across instruments

    it simply stores values to be used in calculations later

    """


    def __init__(self, mode, pixel_size, ao_mode=None, target_separation=0, dx=0, dy=0):

        self.mode = mode
        self.pixel_size = pixel_size
        self.ao_mode = ao_mode
        self.target_separation = target_separation
        self.dx = dx
        self.dy = dy

        self.sed_type = None
        self.sed_params = None
        self.source = None
        self.source_params = None
        self.ao_params = None
        self.sky_params = None
        self.obs_params = None

    def set_sed(self, spectrum_type, **kwargs):
        """

        Parameters
        ----------
        spectrum_type:
        kwargs

        Returns
        -------

        """

        if spectrum_type not in SED_dict.keys():
            raise ValueError("SED not available")
        else:
            func = SED_dict[spectrum_type]
            func_params = dict(kwargs)

        self.sed_type = spectrum_type
        self.sed_params = check_func_params(func, func_params)
        print("SED: '%s' with parameters %s" % (self.sed_type, self.sed_params))

    def set_source(self, distribution, **kwargs):
        """
        Same here
        Parameters
        ----------
        distribution
        kwargs

        Returns
        -------

        """
        if distribution not in FLUX_DISTRO_dict.keys():
            raise ValueError("FLUX distribution not available")
        else:
            func = FLUX_DISTRO_dict[distribution]
            func_params = dict(kwargs)

        self.source = distribution
        self.source_params = check_func_params(func, func_params)
        self.source_params.update(dict(filter_name=self.sed_params["filter_curve"],
                                       magnitude=self.sed_params["magnitude"],
                                       redshift=self.sed_params["redshift"]))
        print("Source '%s' with parameters %s" % (self.source, self.source_params))

    def set_sky_conditions(self, airmass=1.5, moon_phase=0.5, pwv=10, turbulence=100, iq=0.8):
        """
        TODO: Check how ScopeSim can read this
        TODO: turbulence and iq affect the PSF!
        """

        self.sky_params = dict(airmass=airmass,
                               moon_phase=moon_phase,
                               pwv=pwv,
                               turbulence=turbulence,
                               iq=iq)
        print("Sky parameters:", self.sky_params)

    def set_ao(self, profile_name="EsoQ4", zenDist=0, seeing=0.8):
        """
        TODO: Conditions here are redundant with self.set_sky_conditions()
        """

        self.ao_params = dict(profile_name=profile_name,
                              zenDist=zenDist,
                              seeing=seeing)

        print("AO parameters:", self.ao_params)

    def set_setup_obs(self, dit=60, ndit=1, filter_name="Ks"):

        self.obs_params = dict(dit=dit, ndit=ndit, filter_name=filter_name)
        print("Observing parameters:", self.obs_params)

    def _get_sed(self):
        try:
            sed_func = SED_dict[self.sed_type]
        except KeyError as error:
            print("sed type %s unknown" % error)
            raise

        params = self.sed_params
        sp = sed_func(**params)

        return sp

    def _get_source(self):
        src_func = FLUX_DISTRO_dict[self.source]
        params = self.source_params

        params.update({"sed": self._get_sed()
                       })

        src = src_func(**params)
        src.shift(dx=self.dx,
                  dy=self.dy)
        return src

    def _get_sky_conditions(self):

        skycalc_dict  = {"!ATMO.pwv": self.sky_params["pwv"],
                         "!ATMO.airmass": self.sky_params["airmass"],
                         "!ATMO.moon_sun_sep": self.sky_params["moon_phase"],
                         "!ATMO.wmax": 300000
                         }
        return skycalc_dict


class HAWKI_ETC(ETC_base):

   # sim.server.database.download_package(["locations/Paranal.zip", # TODO: Yaml files have some problems, update irdb
   #                                       "telescopes/VLT.zip",
   #                                       "instruments/HAWKI.zip"])

    def __init__(self, mode="imaging", pixel_size=0.106, ao_mode="no_ao"):

        super(HAWKI_ETC, self).__init__(mode=mode, pixel_size=pixel_size, ao_mode=ao_mode)

        self.set_instrument()

    def set_instrument(self):   # This is very instrument specific
        if self.ao_mode.lower() not in ["no_ao", "ao"]:
            raise ValueError

    def _get_psf(self):
        psf_effect = sim.effects.psfs.SeeingPSF(fwhm=self.sky_params["iq"])

        return psf_effect

    def _define_effects(self):


        psf_effect = self._get_psf()
        cmd = sim.UserCommands(use_instrument="HAWKI")


        cmd["!OBS.filter_name"] = self.obs_params["filter_name"]  # observing filter
        cmd["!INST.pixel_scale"] = self.pixel_size
        cmd["!OBS.modes"] = [self.ao_mode.upper(), self.mode]

        cmd["!DET.y"] = self.dy / self.pixel_size
        cmd["!DET.x"] = self.dx / self.pixel_size
        cmd["!ATMO.pwv"] = self.sky_params["pwv"]
        cmd["!ATMO.moon_sun_sep"] = self.sky_params["moon_phase"]
        cmd["!ATMO.airmass"] = self.sky_params["airmass"]
       # cmd[!SIM.spectral.spectral_bin]
        print(cmd)

        hawki = sim.OpticalTrain(cmd)

        hawki["paranal_atmo_skycalc_ter_curve"].include = True
        hawki["paranal_atmo_default_ter_curve"].include = False
        hawki['detector_linearity'].include = False
        hawki.effects["vlt_generic_psf"] = psf_effect

   #     hawki["relay_psf"].include = False
        print(hawki.effects)
   #     hawki.optics_manager["default_ro"].add_effect(psf_effect)

        return hawki

    def run(self, filename=None):
        """
        This returns the results


        Returns
        -------

        """
        src = self._get_source()

        hawki = self._define_effects()

        hawki.observe(src)
        noiseless_obj = hawki.image_planes[0].data
        signal_to_noise = np.sqrt(noiseless_obj)

        sky_src = empty_sky()
        hawki.observe(sky_src)

        #    signal_to_noise = np.sqrt(noiseless_image)
        hdu_obj = hawki.readout(filename=filename)
        observed_image = hdu_obj[0][1]


    def report(self):

        pass




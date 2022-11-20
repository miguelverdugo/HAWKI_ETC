import inspect
import warnings

import numpy as np
from astropy.table import Table

from spextra import Spextrum


#from scopesim_templates.basic.basic import empty_sky
#from scopesim_templates.basic.stars import star
from scopesim import Source
from anisocado import AnalyticalScaoPsf
try:  # temporal holder because conflict versions.
    from scopesim_templates.extragalactic.galaxies import galaxy
except ImportError:
    from scopesim_templates.basic.galaxy import galaxy


def check_func_params(func, params):
    """
    Small function to check for the parameters before evaluating.
    """
    args = inspect.getfullargspec(func).args
    defaults = inspect.getfullargspec(func).defaults
    func_params = dict(zip(args, defaults))

    for k in params.keys():
        if k not in func_params.keys():
            warnings.warn("Parameter %s not part of function %s parameters" % (k, func.__name__),
                          stacklevel=2)
        else:
            func_params[k] = params[k]

    return func_params



# SEDs

def template(template_name='pickles/a0v', magnitude=20, filter_curve="V", redshift=0):
    sp = Spextrum(template_name=template_name).redshift(redshift)
    return sp.scale_to_magnitude(amplitude=magnitude, filter_curve=filter_curve)


def from_file(filename, magnitude, filter_name, redshift):

    sp = Spextrum.from_file(filename=filename)
    sp = sp.redshift(z=redshift).scale_to_magnitude(amplitude=magnitude, filter_curve=filter_name)
    return sp


def MARCS():
    """
    TODO: Implement in spextra
    Parameters
    ----------
    self

    Returns
    -------

    """

    return NotImplementedError


def emission_line(center=9000, fwhm=2, flux=1e-16):
    sp = Spextrum.emission_line_spectra(center=center, fwhm=fwhm, flux=flux)
    return sp


def black_body(temperature=5000, magnitude=15, filter_curve="Ks"):
    sp = Spextrum.black_body_spectrum(temperature=temperature, amplitude=magnitude, filter_curve=filter_curve)
    return sp


def powerlaw(alpha=1, magnitude=15, filter_curve="Ks"):
    sp = Spextrum.powerlaw(alpha=alpha, amplitude=magnitude, filter_curve=filter_curve)
    return sp


def flat_spec(magnitude=15):
    sp = Spextrum.flat_spectrum(amplitude=magnitude)
    return sp


# Sources


def point_source(sed="pickles/a0v", magnitude=15, filter_name="V", redshift=0):
    if isinstance(sed, Spextrum):
        sp1 = sed.redshift(z=redshift)
        sp = sp1.scale_to_magnitude(amplitude=magnitude, filter_curve=filter_name)

    else:
        sp = Spextrum(sed).redshift(z=redshift).scale_to_magnitude(amplitude=magnitude, filter_curve=filter_name)

    tbl = Table(names=["x", "y", "ref", "weight", "spec_types"],
                data=[[0], [0], [0], [1], [sed]])

    src = Source(spectra=[sp], table=tbl)
    # src.meta.update()

    return src


def sersic(sed, magnitude=15, filter_name="V",  redshift=0, r_eff=5, n=4, ellip=0.1, theta=0, extend=3, x=0, y=0):
    src = galaxy(sed=sed,
                 z=redshift,
                 amplitude=magnitude,
                 filter_curve=filter_name,
                 pixel_scale=0.0015,  # !!!!
                 r_eff=r_eff,
                 n=n,
                 ellip=ellip,
                 theta=theta,
                 extend=extend)
    return src


def empty_sky():
    """
    Returns an empty source so that instrumental fluxes can be simulated

    Returns
    -------
    sky : Source

    """
   # params = {"function_call": function_call_str(empty_sky, {}),
   #           "object": "empty sky"}

    sky = Source(lam=np.array([0.7, 2.5]),
                 spectra=np.array([0, 0]), x=[0], y=[0], ref=[0], weight=[0])
  #  sky.meta.update(params)

    return sky


def uniform_flux():
    pass





def scao_psf(wavelength=2.15, profile_name="EsoQ4", zenDist=0, seeing=1, x=0, y=0):
    psf = AnalyticalScaoPsf(N=512, wavelength=wavelength, profile_name=profile_name,
                            zenDist=zenDist, seeing=seeing)
    kernel = psf.shift_off_axis(x, y)

    return kernel




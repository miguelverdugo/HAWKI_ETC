import pytest
import inspect
import os

import numpy as np
import astropy.units as u
from spextra import Spextrum

from ..etc import ETC_base
from ..utils import check_func_params


def test_eval_func_params():
    def little_function(a=1, b=2, c=3):
        return a + b + c

    params = dict(a=10, b=20, d=60)
    result = check_func_params(little_function, params)

    assert result == 60




class TestETCBasic:

    def test_initialization(self):
        etc = ETC_base(mode="imaging", pixel_size=0.1)
        assert isinstance(etc, ETC_base)

    def test_initialization_ii(self):
        etc = ETC_base(mode="imaging", pixel_size=0.0015, ao_mode="MCAO")
        assert etc.pixel_size == 0.0015

    def test_set_instrument(self):
        etc = ETC_base(mode="spectroscopy", pixel_size=0.004, ao_mode="SCAO")
        assert etc.mode =="spectroscopy"

    def test_set_sed(self):
        etc = ETC_base(mode="imaging", pixel_size=0.1)
        params = dict(template_name="kc96/s0", magnitude=18, filter_curve="R", redshift=1)

        etc.set_sed("template", **params)

        assert etc.sed_type == "template"
        assert etc.sed_params.keys() == params.keys()

    def test_get_sed(self):
        etc = ETC_base(mode="imaging", pixel_size=0.1)
        etc.set_sed("template", template_name="kc96/s0", magnitude=15, filter_curve="V", redshift=0)
        print(etc.sed_type, etc.sed_params)
        assert isinstance(etc._get_sed(), Spextrum)


    def test_get_source(self):
        etc = ETC_base(mode="imaging", pixel_size=0.1)
        etc.set_sed("template", template_name="kc96/s0")
        etc.set_source("point_source")



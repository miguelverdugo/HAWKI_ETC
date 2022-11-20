import os
import inspect

from .etc import HAWKI_ETC, ETC_base

__pkg_dir__ = os.path.dirname(inspect.getfile(inspect.currentframe()))
__data_dir__ = os.path.join(__pkg_dir__, "data")

# MICADO ETC 

This is a meta package that will help to generate the validation datasets for the MICADO ETC. 
It is based on the `ScopeSim` framework, but it will use a subset of its capabilities. 

The package is modeled after current ESO ETCs. It will generate simulated data and from that the 
proper calculations of S/N, etc.

## Usage

The proposed usage is the following

```python
from etc import MICADO_ETC

etc = ETC("instrument configuration")
etc.set_sed("sed selection and settings")
etc.set_source("source properties")
etc.set_atmo("atmospheric conditions")
etc.obs_setup("exposure time, etc")
etc.run()  # runs the simulation and produces a report
```



## Requirements

* numpy
* astropy
* ScopeSim
* irdb
* ScopeSim_Templates
* speXtra
* AnisoCADO
* skycalc_ipy

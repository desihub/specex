#!/usr/bin/env python

import desispec.bootcalib

nist_list=desispec.bootcalib.load_arcline_list(camera="all", vacuum="True", lamps=None)
for ion,wave,intensity in zip(nist_list["Ion"],nist_list["wave"],nist_list["RelInt"]) :
    print ion,wave,1,intensity


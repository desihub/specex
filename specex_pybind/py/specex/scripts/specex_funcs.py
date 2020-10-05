# Licensed under a 3-clause BSD style license - see LICENSE.rst
# -*- coding: utf-8 -*-
"""
specex.scripts.specex_funcs
==============================

High-level functions for running assignment.

"""
from __future__ import absolute_import, division, print_function

import os
import sys
import argparse
import re

from ..myargs import (MYARGS_INT, MyArgs, test_myargs)

from ..mymath import (MYMATH_INT, MyMath, test_mymath)

from specex import (specex_desi_psf_fit_main)
def test_specex(args):
    """test specex interface to c++

    """

    i=specex_desi_psf_fit_main(args)
    
    return 

def test_testclasses(args):
    """test test class interface to c++

    """
    test_myargs(args)

    a=3; b=31
    c=test_mymath(a,b)

    print('sum of ',a,' and ',b,' = ',c)

    return 

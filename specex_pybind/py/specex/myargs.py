# Licensed under a 3-clause BSD style license - see LICENSE.rst
# -*- coding: utf-8 -*-
"""
specex.myargs
=====================

Functions for testing argument passing to c++ functions

"""
from __future__ import absolute_import, division, print_function

import numpy as np

import fitsio

from ._internal import (MYARGS_INT, MyArgs)

def str_to_myargs_type(input):
    if input == "domyargs":
        return MYARGS_INT
    else:
        raise ValueError("unknown target type '{}'".format(input))
    return None

def test_myargs(list_of_args=[]):
    ma = MyArgs(MYARGS_INT)

    ma.get_args(list_of_args)
    ma.print_args()
    
    return None

def test_mymath(a,b):
    ma = MyMath(MYMATH_INT)

    ma.get_args(list_of_args)
    ma.print_args()
    
    return None

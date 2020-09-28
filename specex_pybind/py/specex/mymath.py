# Licensed under a 3-clause BSD style license - see LICENSE.rst
# -*- coding: utf-8 -*-
"""
specex.mymath
=====================

Functions for testing trivial math in c++

"""
from __future__ import absolute_import, division, print_function

import numpy as np

import fitsio

from ._internal import (MYMATH_INT, MyMath)

def str_to_mymath_type(input):
    if input == "domymath":
        return MYMATH_INT
    else:
        raise ValueError("unknown target type '{}'".format(input))
    return None

def test_mymath(a,b):
    mm = MyMath(10)
    return mm.add(a,b)


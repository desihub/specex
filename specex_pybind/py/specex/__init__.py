#
# See top-level LICENSE.rst file for Copyright information
#
# -*- coding: utf-8 -*-
"""
specex
==============

Tools for psf fitting with specex

.. _specex: https://github.com/desihub/specex
.. _DESI: http://desi.lbl.gov
.. _Python: http://python.org
"""

from __future__ import absolute_import, division, print_function

from ._version import __version__

from ._internal import specex_desi_psf_fit_main
from ._internal import Options

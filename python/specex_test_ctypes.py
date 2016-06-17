#!/usr/bin/env python

import ctypes as ct
from ctypes.util import find_library

import numpy as np


libspecex = None
try:
    libspecex = ct.CDLL('libspecex.so')
except:
    path = find_library('specex')
    if path is not None:
        libspecex = ct.CDLL(path)

if libspecex is not None:
    libspecex.cspecex_desi_psf_fit.restype = ct.c_int
    libspecex.cspecex_desi_psf_fit.argtypes = [
        ct.c_int,
        ct.POINTER(ct.POINTER(ct.c_char))
    ]


def main():
    com = ['specex_desi_psf_fit']
    com.extend(['-a', 'blah.fits'])

    argc = len(com)

    arg_buffers = [ct.create_string_buffer(com[i]) for i in range(argc)]
    addrlist = [ ct.cast(x, ct.POINTER(ct.c_char)) for x in map(ct.addressof, arg_buffers) ]
    arg_pointers = (ct.POINTER(ct.c_char) * argc)(*addrlist)

    retval = libspecex.cspecex_desi_psf_fit(argc, arg_pointers)

    return
        


if __name__ == '__main__':
    main()
    

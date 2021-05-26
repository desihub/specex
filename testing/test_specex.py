
from specex.specex import run_specex
import sys

psffile=sys.argv[1]

first_fiber='100'
last_fiber='104'

com=['desi_psf_fit', '-a', '/global/cfs/cdirs/desi/spectro/redux/blanc/preproc/20201216/00068217/preproc-b1-00068217.fits', '--in-psf', '/global/cfs/cdirs/desi/spectro/redux/blanc/exposures/20201216/00068217/shifted-input-psf-b1-00068217.fits', '--out-psf',psffile,'--lamp-lines', '/global/common/software/desi/cori/desiconda/20200801-1.4.0-spec/code/specex/0.6.7/data/specex_linelist_desi.txt', '--first-bundle', '4', '--last-bundle', '4', '--first-fiber', first_fiber, '--last-fiber', last_fiber, '--legendre-deg-wave', '3','--fit-continuum']#,'--debug']

run_specex(com)

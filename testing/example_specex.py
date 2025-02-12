
from specex.specex import run_specex
import sys

psffile=sys.argv[1]
#inputroot=sys.argv[2].rstrip('\n')
inputroot_pre='/global/cfs/cdirs/desi/spectro/redux/blanc/preproc/20201216/00068217/'
inputroot_exp='/global/cfs/cdirs/desi/spectro/redux/blanc/exposures/20201216/00068217/'
first_fiber='100'
last_fiber='104'

com=['desi_psf_fit', '-a', inputroot_pre+'/preproc-b1-00068217.fits', '--in-psf', inputroot_exp+'/shifted-input-psf-b1-00068217.fits', '--out-psf',psffile,'--lamp-lines', 'specex-dev/code/lib/python3.8/site-packages/specex-0.7.0.dev637-py3.8-linux-x86_64.egg/specex/data/specex_linelist_desi.txt', '--first-bundle', '4', '--last-bundle', '4', '--first-fiber', first_fiber, '--last-fiber', last_fiber, '--legendre-deg-wave', '3','--fit-continuum','--debug']

run_specex(com)
#run_specex(['desi_psf_fit','-a','foo'])

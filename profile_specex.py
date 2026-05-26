
import os
import sys
import cProfile
import pstats
import io

# Add build and py directories to PYTHONPATH
current_dir = os.getcwd()
sys.path.insert(0, os.path.join(current_dir, 'py'))
sys.path.insert(0, os.path.join(current_dir, 'build'))

from specex.specex import run_specex

def main():
    lamp_lines = os.path.join(current_dir, 'py/specex/data/specex_linelist_desi.txt')
    out_psf = os.path.join(os.environ.get('SCRATCH', '/tmp'), 'fit-psf-z8-00344649_05.fits')
    
    com = [
        'desi_psf_fit',
        '-a', '/dvs_ro/cfs/cdirs/desi/spectro/redux/matterhorn/preproc/20260401/00344649/preproc-z8-00344649.fits.gz',
        '--in-psf', '/dvs_ro/cfs/cdirs/desi/spectro/redux/matterhorn/exposures/20260401/00344649/shifted-input-psf-z8-00344649.fits',
        '--lamp-lines', lamp_lines,
        '--out-psf', out_psf,
        '--first-bundle', '5',
        '--last-bundle', '5',
        '--first-fiber', '125',
        '--last-fiber', '149',
        '--legendre-deg-wave', '3',
        '--fit-continuum',
        '--broken-fibers', '473,474'
    ]
    
    print(f"Running command: {' '.join(com)}")
    
    pr = cProfile.Profile()
    pr.enable()
    
    retval = run_specex(com)
    
    pr.disable()
    s = io.StringIO()
    sortby = pstats.SortKey.CUMULATIVE
    ps = pstats.Stats(pr, stream=s).sort_stats(sortby)
    ps.print_stats(50)
    print(s.getvalue())
    
    return retval

if __name__ == "__main__":
    sys.exit(main())

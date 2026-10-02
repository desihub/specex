from desispec.io.xytraceset import read_xytraceset
from desiutil.log import get_logger
import numpy as np

def trace_psf_qa(psf_filename, broken_fiber_list):
    """Check an output PSF FITS file for neighboring fiber pairs whose fitted X-position traces cross (fiber i+1's trace not everywhere above fiber i's, over the shared wavelength range), logging a warning for broken-fiber-involving crossings and an error for unexpected ones.

    Args:
        psf_filename (str): path to the output PSF FITS file to check (read
            via desispec's read_xytraceset).
        broken_fiber_list (str): comma-separated broken fiber IDs; crossings
            involving one of these are logged as warnings, not counted as
            failures.

    Returns:
        int: failcount, the number of crossing pairs where neither fiber is in
        broken_fiber_list.

    Status: DEAD -- adapted from upstream specex#91; imported by
    specex.run_specex() (as specex_psf_qa) but the one call site there is
    commented out. io.py's write_python_psf has its own independently-adapted
    inline reimplementation of this same trace-crossing check (STATUS=4
    flagging), which is what actually runs in the production pipeline.
    """

    log = get_logger()

    if len(broken_fiber_list) > 0:
        brokenfibers = list(map(int,broken_fiber_list.split(",")))
    else:
        brokenfibers = []

    fibertraces = read_xytraceset(psf_filename)
    ww = np.arange(fibertraces.wavemin, fibertraces.wavemax)
    
    failcount=0
    for fiber in range(0,fibertraces.nspec-1):
        correct = np.all(fibertraces.x_vs_wave(fiber+1, ww) > fibertraces.x_vs_wave(fiber, ww))
        if not correct:
            if fiber in brokenfibers:
                log.warning("broken fiber {} overlaps {} in {}".format(fiber, fiber+1, psf_filename))
            if fiber+1 in brokenfibers:
                log.warning("broken fiber {} overlaps {} in {}".format(fiber+1, fiber, psf_filename))
            if fiber not in brokenfibers and fiber+1 not in brokenfibers:
                log.error("overlapping traces for fibers {} and {} in {}".format(fiber, fiber+1, psf_filename))
                failcount += 1

    return failcount

def specex_psf_qa(opts):

    # trace QA
    """Run trace_psf_qa on the PSF file/broken-fiber list named in a C++-style options object.

    Args:
        opts: an options object with `.output_fits_filename` (str, PSF file
            path) and `.broken_fibers_string` (str, comma-separated fiber IDs)
            attributes -- matches the C++ pybind11 PyOptions interface used by
            specex.run_specex().

    Returns:
        int: failcount from trace_psf_qa.

    Status: DEAD -- see trace_psf_qa; the one call site (specex.run_specex())
    has this call commented out.
    """
    psf_filename = opts.output_fits_filename
    broken_fiber_list = opts.broken_fibers_string

    failcount = 0

    failcount += trace_psf_qa(psf_filename, broken_fiber_list)

    return failcount

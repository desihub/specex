from desispec.io.xytraceset import read_xytraceset
from desiutil.log import get_logger
import numpy as np

def trace_psf_qa(psf_filename, broken_fiber_list):
    """
    QA PSF by reporting overlapping fiber traces.

    Args:
        psf_filename: string, input PSF file
        broken_fiber_list: string, comma separated list of broken fibers

    Returns:
        failcount: int, number of neighboring fibers with overlapping traces
                   where neither is on in broken_fiber_list
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
    psf_filename = opts.output_fits_filename
    broken_fiber_list = opts.broken_fibers_string

    failcount = 0

    failcount += trace_psf_qa(psf_filename, broken_fiber_list)

    return failcount

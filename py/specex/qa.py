from desispec.io.xytraceset import read_xytraceset
from desiutil.log import get_logger
import numpy as np

def psf_qa(opts):
    log = get_logger()

    brokenfibers = list(map(int,opts.broken_fibers_string.split(",")))
    fibertraces = read_xytraceset(opts.output_fits_filename)
    ww = np.arange(fibertraces.wavemin, fibertraces.wavemax)
    
    failcount=0
    for fiber in range(0,fibertraces.nspec-1):
        correct = np.all(fibertraces.x_vs_wave(fiber+1, ww) > fibertraces.x_vs_wave(fiber, ww))
        if not correct:
            if fiber in brokenfibers:
                log.warning("broken fiber {} overlaps {} in {}".format(fiber, fiber+1, opts.output_fits_filename))
            if fiber+1 in brokenfibers:
                log.warning("broken fiber {} overlaps {} in {}".format(fiber+1, fiber, opts.output_fits_filename))
            if fiber not in brokenfibers and fiber+1 not in brokenfibers:
                log.error("overlapping traces for fibers {} and {} in {}".format(fiber, fiber+1, opts.output_fits_filename))
                failcount += 1

    return failcount

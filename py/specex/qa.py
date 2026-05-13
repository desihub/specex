from desispec.io.xytraceset import read_xytraceset
from desiutil.log import get_logger
import numpy as np
from astropy.io import fits

def trace_psf_qa(psf_filename, broken_fiber_list):
    """
    QA PSF by reporting overlapping fiber traces.

    Args:
        psf_filename: string, input PSF file
        broken_fiber_list: string, comma separated list of broken fibers

    Returns:
        failcount: int, number of neighboring fiber pairs with overlapping traces
                   where neither is in broken_fiber_list
        bad_fibers: set of fiber indices (0-based relative to FIBERMIN) involved
                    in overlapping traces where neither fiber is in broken_fiber_list
    """

    log = get_logger()

    if len(broken_fiber_list) > 0:
        brokenfibers = list(map(int,broken_fiber_list.split(",")))
    else:
        brokenfibers = []

    fibertraces = read_xytraceset(psf_filename)
    ww = np.arange(fibertraces.wavemin, fibertraces.wavemax)

    failcount = 0
    bad_fibers = set()
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
                bad_fibers.add(fiber)
                bad_fibers.add(fiber+1)

    return failcount, bad_fibers

def specex_psf_qa(opts):

    # trace QA
    psf_filename = opts.output_fits_filename
    broken_fiber_list = opts.broken_fibers_string

    failcount, bad_fibers = trace_psf_qa(psf_filename, broken_fiber_list)

    if bad_fibers:
        file = fits.open(psf_filename)
        params = file['PSF']['PARAM'][:]
        coeff_all = file['PSF']['COEFF'][:]
        for i, param in enumerate(params):
            if param.strip() == 'STATUS':
                for fiber in bad_fibers:
                    index = np.where(file['PSF'].data["PARAM"]=="STATUS")[0][0]
                    file['PSF'].data["COEFF"][index][fiber]=4
                    # coeff_all[i, fiber, 0] = 4
                # file['PSF'].write_column('COEFF', coeff_all)
                break
        file.close()

    return failcount

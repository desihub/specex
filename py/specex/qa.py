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
    bad_fibers = []
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
                log.info(f'Failcount is now {failcount}')
                bad_fibers.append(fiber)
                log.info(f'Bad fibers are now {bad_fibers}')
                bad_fibers.append(fiber+1)
                log.info(f'Bad fibers are now {bad_fibers}')
                log.info(f'Fibers {fiber} and {fiber+1} are flagged')
    log.info(f'Failcount is {failcount}')

    return failcount, bad_fibers

def specex_psf_qa(opts):
    log=get_logger()
    # trace QA
    psf_filename = opts.output_fits_filename
    broken_fiber_list = opts.broken_fibers_string

    failcount, bad_fibers = trace_psf_qa(psf_filename, broken_fiber_list)
    if len(bad_fibers)>0:
        log.info(f'Setting Status of Bad Fibers {bad_fibers} to 4 in PSF file {psf_filename}')
        other_psf_hdulist=fits.open(psf_filename)
        i=np.where(other_psf_hdulist["PSF"].data["PARAM"]=="STATUS")[0][0]
        status_of_fibers = \
            other_psf_hdulist["PSF"].data["COEFF"][i][:,0].astype(int)
        log.info(f'Status of fibers before: {status_of_fibers}')
        with fits.open(psf_filename, mode='update',memmap=False) as file:
            index = np.where(file['PSF'].data["PARAM"]=="STATUS")[0][0]
            # log.info(f'Index of STATUS parameter is {index}')
            for fiber in bad_fibers:
                file['PSF'].data["COEFF"][index][fiber, 0] = 4
            # status_of_fibers = \
                # file["PSF"].data["COEFF"][i][:,0].astype(int)
            # log.info(f'Status of fibers: {status_of_fibers}')
            log.info(f'Bad fibers {bad_fibers} have been set to 4 in PSF file {psf_filename}')
            file.flush()
        log.info(f'File now closed')

    log.info(f'QA complete for PSF file {psf_filename} with failcount {failcount} and bad fibers {bad_fibers}')
    other_psf_hdulist=fits.open(psf_filename,memmap=False)

    # look at what fibers where actually fit
    i=np.where(other_psf_hdulist["PSF"].data["PARAM"]=="STATUS")[0][0]
    status_of_fibers = \
        other_psf_hdulist["PSF"].data["COEFF"][i][:,0].astype(int)
    log.info(f'Status of fibers after: {status_of_fibers}')
    other_psf_hdulist.close()
    return failcount

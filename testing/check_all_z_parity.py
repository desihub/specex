import numpy as np
import os
import re
from specex.io import read_lamp_lines, load_python_psf, read_image
from specex.fitter import get_bundle_spots

def check_all_z():
    cams = [f'z{i}' for i in range(10)]
    night = '20260401'
    expid = '00344649'
    arc_b = '/dvs_ro/cfs/cdirs/desi/spectro/redux/matterhorn/preproc'
    psf_b = '/dvs_ro/cfs/cdirs/desi/spectro/redux/matterhorn/exposures'
    lines = read_lamp_lines('py/specex/data/specex_linelist_desi.txt')
    
    print(f"Cam | C++ Spots | Py (3.15) | Parity %")
    print("-" * 40)
    
    for cam in cams:
        # Python Fit
        arc = f'{arc_b}/{night}/{expid}/preproc-{cam}-{expid}.fits.gz'
        psf_f = f'{psf_b}/{night}/{expid}/shifted-input-psf-{cam}-{expid}.fits'
        
        try:
            opts = type('Opts',(object,),{'arc_image_filename':arc,'input_psf_filename':psf_f})()
            d = read_image(arc)
            psf = load_python_psf(psf_f, opts)
            broken = '367' if cam == 'z0' else '473,474'
            spots = get_bundle_spots(psf, 125, 149, lines, image=d['image'].T, weight=d['ivar'].T, sn_threshold=3.15, broken_fibers=broken)
            py_count = len(spots)
        except Exception as e:
            py_count = -1
            
        # Scrape C++ Spot Count from log
        cpp_count = -1
        log_file = f'debug_{cam}.txt'
        if os.path.exists(log_file):
            with open(log_file, 'r') as f:
                content = f.read()
                matches = re.findall(r"selected (\d+) spots out of 1750 with S/N>3", content)
                if matches:
                    cpp_count = int(matches[-1])
        
        parity = "Match" if py_count == cpp_count else f"{(py_count/cpp_count-1)*100:+.2f}%" if cpp_count > 0 else "N/A"
        print(f"{cam:<3} | {cpp_count:>9} | {py_count:>9} | {parity}")

if __name__ == "__main__":
    check_all_z()

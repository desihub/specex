import subprocess
import os
import re
import time

def get_cpp_counts():
    cams = [f'z{i}' for i in range(10)]
    results = {}
    
    for cam in cams:
        print(f"Testing {cam}...")
        log_file = f"debug_{cam}_selected.txt"
        cmd = f"module load libfabric && source env_setup.sh && desi_psf_fit -a /dvs_ro/cfs/cdirs/desi/spectro/redux/matterhorn/preproc/20260401/00344649/preproc-{cam}-00344649.fits.gz --in-psf /dvs_ro/cfs/cdirs/desi/spectro/redux/matterhorn/exposures/20260401/00344649/shifted-input-psf-{cam}-00344649.fits --lamp-lines py/specex/data/specex_linelist_desi.txt --out-psf test_cpp_{cam}.fits --first-bundle 5 --last-bundle 5 --legendre-deg-wave 3 --fit-continuum --broken-fibers 473,474"
        
        # Run until we see "selected ... with S/N>3"
        proc = subprocess.Popen(cmd, shell=True, stdout=subprocess.PIPE, stderr=subprocess.STDOUT, text=True, executable='/bin/bash')
        
        count = -1
        try:
            # Poll stdout
            start_t = time.time()
            while True:
                line = proc.stdout.readline()
                if not line: break
                # print(line.strip())
                if "selected" in line and "S/N>3" in line:
                    match = re.search(r"selected (\d+) spots", line)
                    if match:
                        count = int(match.group(1))
                        break
                if time.time() - start_t > 120: # 2 minute timeout per camera
                    print(f"Timeout for {cam}")
                    break
        finally:
            proc.terminate()
            proc.wait()
            
        results[cam] = count
        print(f"  {cam}: {count}")
    
    print("\nFinal C++ Counts:")
    for cam, count in results.items():
        print(f"  {cam}: {count}")

if __name__ == "__main__":
    get_cpp_counts()

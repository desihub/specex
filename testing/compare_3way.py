import sys
import numpy as np
import fitsio

def analyze_bundle(filename, bid=5):
    f = fitsio.FITS(filename)
    xt = f['XTRACE'].read()[bid*25:(bid+1)*25, 0] # Constant term of trace
    yt = f['YTRACE'].read()[bid*25:(bid+1)*25, 0]
    psf = f['PSF'].read()
    names = [n.strip() for n in psf['PARAM']]
    idx_sigx = names.index('GHSIGX')
    idx_gh00 = names.index('GH-0-0')
    
    sigx = psf['COEFF'][idx_sigx, bid*25:(bid+1)*25, 0]
    gh00 = psf['COEFF'][idx_gh00, bid*25:(bid+1)*25, 0]
    
    return {
        'xt': xt, 'yt': yt, 'sigx': sigx, 'gh00': gh00
    }

def compare_3way(bid=5):
    file_cpp = '/dvs_ro/cfs/cdirs/desi/spectro/redux/matterhorn/exposures/20260401/00344649/fit-psf-z8-00344649.fits'
    file_gpu = 'python-gpu-fit-z8-00344649.fits'
    file_cpu = f'python-cpu-fit-bundle-{bid}.fits'
    
    print(f"--- 3-Way Comparison for Bundle {bid} ---")
    
    data = {}
    data['CPP'] = analyze_bundle(file_cpp, bid)
    try: data['GPU'] = analyze_bundle(file_gpu, bid)
    except: print("GPU file not ready.")
    try: data['CPU'] = analyze_bundle(file_cpu, bid)
    except: print("CPU file not ready.")
    
    for key in ['sigx', 'gh00', 'xt']:
        print(f"\nParameter: {key.upper()}")
        print(f"{'Method':<8} | {'Mean':<12} | {'RMS (vs CPP)':<12}")
        print("-" * 35)
        v_cpp = data['CPP'][key]
        print(f"{'C++':<8} | {np.mean(v_cpp):>12.6f} | {'Baseline':>12}")
        
        for m in ['CPU', 'GPU']:
            if m in data:
                v = data[m][key]
                rms = np.sqrt(np.mean((v - v_cpp)**2))
                print(f"{m:<8} | {np.mean(v):>12.6f} | {rms:>12.6e}")

if __name__ == "__main__":
    compare_3way(5)

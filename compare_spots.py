import numpy as np

def load_spots(path):
    data = []
    with open(path, 'r') as f:
        for line in f:
            parts = line.strip().split(',')
            if len(parts) < 6: continue
            # fiber, wave, xc, yc, flux, eflux/snr
            data.append([float(p) for p in parts])
    return np.array(data)

cpp_path = '/pscratch/sd/c/cdwarner/specex/testing/fit-psf-z8-00344649_05_cpp_debug_v3.rawspots.txt'
py_path = '/pscratch/sd/c/cdwarner/specex/testing/pyfit-psf-z8-00344649_05_cppstyle_v6.pyrawspots.txt'

cpp_data = load_spots(cpp_path)
py_data = load_spots(py_path)

print(f"CPP spots: {len(cpp_data)}")
print(f"PY spots: {len(py_data)}")

if len(cpp_data) != len(py_data):
    print("Count mismatch!")
else:
    # Compare flux and eflux/snr
    # C++: flux=4, eflux=5
    # PY: flux=4, snr=5
    diff_flux = np.abs(cpp_data[:, 4] - py_data[:, 4])
    diff_eflux = np.abs(cpp_data[:, 5] - py_data[:, 5])
    
    print(f"Flux mean diff: {np.mean(diff_flux):.6f}, max: {np.max(diff_flux):.6f}")
    print(f"eFlux/SNR mean diff: {np.mean(diff_eflux):.6f}, max: {np.max(diff_eflux):.6f}")

    # Selection check
    cpp_sel = cpp_data[:, 5] >= 5.0 # This is wrong if index 5 is eflux, but we are testing
    # Actually C++ check is flux/eflux >= 5
    cpp_snr = cpp_data[:, 4] / cpp_data[:, 5]
    py_snr = py_data[:, 5] # In PY, index 5 is already SNR
    
    cpp_pass = cpp_snr >= 5.0
    py_pass = py_snr >= 5.0
    
    print(f"CPP pass count: {np.sum(cpp_pass)}")
    print(f"PY pass count: {np.sum(py_pass)}")
    
    mismatch = np.where(cpp_pass != py_pass)[0]
    print(f"Mismatches: {len(mismatch)}")
    if len(mismatch) > 0:
        print("\nFirst 10 mismatches:")
        for i in mismatch[:10]:
            print(f"Idx {i}: CPP(f={cpp_data[i,4]:.4f}, e={cpp_data[i,5]:.4f}, snr={cpp_snr[i]:.4f}) vs PY(f={py_data[i,4]:.4f}, snr={py_snr[i]:.4f})")


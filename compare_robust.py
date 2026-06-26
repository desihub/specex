import numpy as np

def load_spots(path):
    data = []
    with open(path, 'r') as f:
        for line in f:
            parts = line.strip().split(',')
            if len(parts) == 4:
                data.append([float(x) for x in parts])
    return np.array(data)

def find_matches(set_a, set_b, tol=1e-7):
    # set_a: (N, 4), set_b: (M, 4)
    # Matches based on fiber (index 0) and wavelength (index 1)
    matches = []
    # Optimization: only compare within same fiber
    fibers_a = np.unique(set_a[:, 0])
    for f in fibers_a:
        idx_a = np.where(set_a[:, 0] == f)[0]
        idx_b = np.where(set_b[:, 0] == f)[0]
        if len(idx_b) == 0: continue
        
        # Compare wavelengths for these fibers
        waves_a = set_a[idx_a, 1]
        waves_b = set_b[idx_b, 1]
        
        for i, w_a in enumerate(waves_a):
            # Find first match in waves_b
            diffs = np.abs(waves_b - w_a)
            match_idx = np.where(diffs < tol)[0]
            if len(match_idx) > 0:
                # Store the indices relative to the original arrays
                matches.append((idx_a[i], idx_b[match_idx[0]]))
                
    return np.array(matches)

def main():
    cpp_raw_path = '/pscratch/sd/c/cdwarner/specex/cpp-fit-psf-z8-00344649_05.rawspots.txt'
    py_raw_path = '/pscratch/sd/c/cdwarner/specex/pyfit-psf-z8-00344649_05.pyrawspots.txt'
    cpp_sel_path = '/pscratch/sd/c/cdwarner/specex/cpp-fit-psf-z8-00344649_05.cpp_cp2_pass4.txt'
    py_sel_path = '/pscratch/sd/c/cdwarner/specex/pyfit-psf-z8-00344649_05.pyspots.txt'

    cpp_raw = load_spots(cpp_raw_path)
    py_raw = load_spots(py_raw_path)
    cpp_sel = load_spots(cpp_sel_path)
    py_sel = load_spots(py_sel_path)

    print(f"Counts: Raw C++={len(cpp_raw)}, Raw Py={len(py_raw)}, Sel C++={len(cpp_sel)}, Sel Py={len(py_sel)}")

    # 1. Verify raw spots match
    raw_matches = find_matches(cpp_raw, py_raw)
    print(f"Raw matches (tol=1e-7): {len(raw_matches)} / 1700")

    # 2. Verify selected are subsets of raw
    cpp_sub_matches = find_matches(cpp_sel, cpp_raw)
    py_sub_matches = find_matches(py_sel, py_raw)
    print(f"C++ Selected is subset of C++ Raw: {len(cpp_sub_matches) == len(cpp_sel)}")
    print(f"Py Selected is subset of Py Raw: {len(py_sub_matches) == len(py_sel)}")

    # 3. Compare selected spots
    sel_matches = find_matches(cpp_sel, py_sel)
    print(f"Selected matches (tol=1e-7): {len(sel_matches)}")

    if len(sel_matches) > 0:
        # Calculate Centroid RMS
        # Indices of matches in cpp_sel and py_sel
        idx_cpp = sel_matches[:, 0]
        idx_py = sel_matches[:, 1]
        
        coord_cpp = cpp_sel[idx_cpp, 2:]
        coord_py = py_sel[idx_py, 2:]
        
        dist_sq = np.sum((coord_cpp - coord_py)**2, axis=1)
        rms = np.sqrt(np.mean(dist_sq))
        print(f"Relative Centroid RMS: {rms:.6f} pixels")

if __name__ == '__main__':
    main()

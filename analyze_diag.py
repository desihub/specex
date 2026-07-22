import numpy as np

def load_diag(path):
    data = []
    with open(path, 'r') as f:
        next(f) # skip header
        for line in f:
            parts = line.strip().split(',')
            if len(parts) == 7:
                data.append([float(x) for x in parts])
    return np.array(data)

def load_spots(path):
    data = []
    with open(path, 'r') as f:
        for line in f:
            parts = line.strip().split(',')
            if len(parts) >= 4:
                data.append([float(x) for x in parts[:4]])
    return np.array(data)

def find_matches(set_a, set_b, tol=1e-7):
    matches = []
    fibers_a = np.unique(set_a[:, 0])
    for f in fibers_a:
        idx_a = np.where(set_a[:, 0] == f)[0]
        idx_b = np.where(set_b[:, 0] == f)[0]
        if len(idx_b) == 0: continue
        waves_a = set_a[idx_a, 1]
        waves_b = set_b[idx_b, 1]
        for i, w_a in enumerate(waves_a):
            diffs = np.abs(waves_b - w_a)
            match_idx = np.where(diffs < tol)[0]
            if len(match_idx) > 0:
                matches.append((idx_a[i], idx_b[match_idx[0]]))
    return np.array(matches)

def main():
    diag_path = '/pscratch/sd/c/cdwarner/specex/pyfit-psf-z8-00344649_05.py_selection_diag.txt'
    py_sel_path = '/pscratch/sd/c/cdwarner/specex/pyfit-psf-z8-00344649_05.pyspots.txt'
    cpp_sel_path = '/global/cfs/cdirs/desicollab/users/cdwarner/code/specex/fit_cpp_z8_b5.cppspots_pass4.txt'
    
    diag = load_diag(diag_path)
    py_sel = load_spots(py_sel_path)
    cpp_sel = load_spots(cpp_sel_path)
    
    # Matches
    matches = find_matches(py_sel, cpp_sel)
    matched_py_indices = matches[:, 0]
    matched_cpp_indices = matches[:, 1]
    
    # Ghost spots (Py Only)
    ghost_py_indices = np.setdiff1d(np.arange(len(py_sel)), matched_py_indices)
    
    # Missing spots (C++ Only)
    missing_cpp_indices = np.setdiff1d(np.arange(len(cpp_sel)), matched_cpp_indices)
    
    # To link py_sel and diag: py_sel is a subset of diag.
    # We need to find where each py_sel spot is in the original diag list.
    # py_sel indices in diag:
    py_sel_in_diag = []
    for s in py_sel:
        # Search diag for Fiber and Wave
        match = np.where((diag[:, 0] == s[0]) & (np.abs(diag[:, 1] - s[1]) < 1e-7))[0]
        if len(match) > 0:
            py_sel_in_diag.append(match[0])
        else:
            py_sel_in_diag.append(-1)
    
    py_sel_in_diag = np.array(py_sel_in_diag)
    
    # Ghosts' diagnostics
    ghost_diag_indices = py_sel_in_diag[ghost_py_indices]
    ghost_stats = diag[ghost_diag_indices]
    
    # Missing's diagnostics (find them in diag)
    missing_diag_indices = []
    for s in cpp_sel[missing_cpp_indices]:
        match = np.where((diag[:, 0] == s[0]) & (np.abs(diag[:, 1] - s[1]) < 1e-7))[0]
        if len(match) > 0:
            missing_diag_indices.append(match[0])
        else:
            missing_diag_indices.append(-1)
    
    missing_stats = diag[np.array(missing_diag_indices)]
    
    print("--- GHOST STATS ---")
    print("Fiber,Wave,Flux,SNR,Chi2,Conv,Shift")
    for row in ghost_stats:
        print(",".join(map(lambda x: f"{x:.6f}", row)))
    
    print("\n--- MISSING STATS ---")
    print("Fiber,Wave,Flux,SNR,Chi2,Conv,Shift")
    for row in missing_stats:
        print(",".join(map(lambda x: f"{x:.6f}", row)))
    
    # Analysis
    print("\n--- SUMMARY ---")
    print(f"Ghost Chi2 Mean: {np.mean(ghost_stats[:, 4]):.4f}, Min: {np.min(ghost_stats[:, 4]):.4f}, Max: {np.max(ghost_stats[:, 4]):.4f}")
    print(f"Ghost SNR Mean: {np.mean(ghost_stats[:, 3]):.4f}, Min: {np.min(ghost_stats[:, 3]):.4f}, Max: {np.max(ghost_stats[:, 3]):.4f}")
    print(f"Missing Chi2 Mean: {np.mean(missing_stats[:, 4]):.4f}, Min: {np.min(missing_stats[:, 4]):.4f}, Max: {np.max(missing_stats[:, 4]):.4f}")
    print(f"Missing SNR Mean: {np.mean(missing_stats[:, 3]):.4f}, Min: {np.min(missing_stats[:, 3]):.4f}, Max: {np.max(missing_stats[:, 3]):.4f}")
    print(f"Missing Converged Rate: {np.mean(missing_stats[:, 5]):.2%}")

if __name__ == '__main__':
    main()

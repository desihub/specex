import numpy as np

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
    # Update paths to the actual files generated in the previous run
    py_sel_path = '/pscratch/sd/c/cdwarner/specex/pyfit-psf-z8-00344649_05.pyspots.txt'
    cpp_sel_path = '/global/cfs/cdirs/desicollab/users/cdwarner/code/specex/fit_cpp_z8_b5.cppspots_pass4.txt'
    
    py_sel = load_spots(py_sel_path)
    cpp_sel = load_spots(cpp_sel_path)
    
    print(f"Counts: Py={len(py_sel)}, C++={len(cpp_sel)}")
    
    # Find matches to isolate "ghost" spots in Python
    matches = find_matches(py_sel, cpp_sel)
    matched_py_indices = matches[:, 0]
    
    # Ghost spots are those in Py but not in C++
    ghost_indices = np.setdiff1d(np.arange(len(py_sel)), matched_py_indices)
    ghost_spots = py_sel[ghost_indices]
    
    print(f"Matches: {len(matches)}")
    print(f"Ghost spots: {len(ghost_indices)}")
    
    # Log ghost spots for analysis
    with open('ghost_spots_analysis.txt', 'w') as f:
        f.write("Fiber,Wavelength,X,Y\n")
        for spot in ghost_spots:
            f.write(f"{spot[0]},{spot[1]},{spot[2]},{spot[3]}\n")
    
    print("Ghost spots written to ghost_spots_analysis.txt")

if __name__ == '__main__':
    main()

import numpy as np
import argparse
from astropy.io import fits
from specex.io import read_lamp_lines
from specex.math import Legendre1DPol

def calculate_rms(diffs):
    if len(diffs) == 0: return np.nan
    return np.sqrt(np.mean(np.square(diffs)))

def load_spots(path):
    """
    Loads spots from a CSV-like file: fiber,wave,x,y
    """
    try:
        data = np.loadtxt(path, delimiter=',')
        return {
            'fiber': data[:, 0].astype(int),
            'wave': data[:, 1],
            'x': data[:, 2],
            'y': data[:, 3]
        }
    except Exception as e:
        print(f"Error loading spots from {path}: {e}")
        return None

def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--py-psf", required=True)
    parser.add_argument("--cpp-psf", required=True)
    parser.add_argument("--py-spots", required=True)
    parser.add_argument("--cpp-spots", required=True)
    parser.add_argument("--lamp-lines", required=True)
    parser.add_argument("--bundle", type=int, required=True)
    args = parser.parse_args()

    # Load Truth
    truth_lines = read_lamp_lines(args.lamp_lines)
    
    # Extract Centroids from sidecar files
    py_data = load_spots(args.py_spots)
    cpp_data = load_spots(args.cpp_spots)

    if py_data is None or cpp_data is None:
        print("Failed to load spots. Check sidecar files.")
        return

    print(f"--- RMS Analysis for Bundle {args.bundle} ---")
    
    # 1. Relative Centroid RMS (Py vs CPP)
    # Match by (fiber, wave)
    py_spots_set = set(zip(py_data['fiber'], np.round(py_data['wave'], 3)))
    cpp_spots_set = set(zip(cpp_data['fiber'], np.round(cpp_data['wave'], 3)))
    
    common = py_spots_set.intersection(cpp_spots_set)
    print(f"Common spots found: {len(common)} / {min(len(py_spots_set), len(cpp_spots_set))}")

    rel_dist = []
    
    # Build lookups for fast access
    py_lookup = {(f, w): (x, y) for (f, w), (x, y) in zip(zip(py_data['fiber'], np.round(py_data['wave'], 3)), zip(py_data['x'], py_data['y']))}
    cpp_lookup = {(f, w): (x, y) for (f, w), (x, y) in zip(zip(cpp_data['fiber'], np.round(cpp_data['wave'], 3)), zip(cpp_data['x'], cpp_data['y']))}
    
    for f, w in common:
        px, py = py_lookup[(f, w)]
        cx, cy = cpp_lookup[(f, w)]
        rel_dist.append(np.sqrt((px - cx)**2 + (py - cy)**2))

    print(f"Relative Centroid RMS (Py vs CPP): {calculate_rms(rel_dist):.6f} px")

    # 2. Wavelength Accuracy (Absolute RMS)
    print("Wavelength fit evaluation is currently disabled (missing trace domain in FITS).")
    print("Please provide wmin/wmax for the fiber traces to enable this.")

if __name__ == "__main__":
    main()

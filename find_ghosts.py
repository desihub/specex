import numpy as np

def load_spots(path):
    spots = set()
    with open(path, 'r') as f:
        for line in f:
            parts = line.strip().split(',')
            if len(parts) == 4:
                # Use rounded values to avoid precision issues in set comparison
                # Fiber is int, wave is float, xc/yc are float.
                # C++ uses setprecision(15).
                fiber = parts[0]
                wave = float(parts[1])
                xc = float(parts[2])
                yc = float(parts[3])
                spots.add((int(fiber), round(wave, 6), round(xc, 6), round(yc, 6)))
    return spots

py_spots = load_spots('/pscratch/sd/c/cdwarner/specex/pyfit-psf-z8-00344649_05.pyspots.txt')
cpp_spots = load_spots('/pscratch/sd/c/cdwarner/specex/cpp-fit-psf-z8-00344649_05.cppspots_pass4.txt')

ghosts = py_spots - cpp_spots

print(f"Python spots: {len(py_spots)}")
print(f"C++ spots: {len(cpp_spots)}")
print(f"Ghost spots: {len(ghosts)}")

# Sort ghosts by fiber for analysis
sorted_ghosts = sorted(list(ghosts))

print("\nFiber,Wave,XC,YC")
for g in sorted_ghosts:
    print(f"{g[0]},{g[1]:.6f},{g[2]:.6f},{g[3]:.6f}")

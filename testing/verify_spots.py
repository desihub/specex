import numpy as np
import argparse
from specex.io import load_python_psf, read_preproc, read_lamp_lines
from specex.fitter import get_bundle_spots

def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("-a", "--arc", required=True)
    parser.add_argument("--in-psf", required=True)
    parser.add_argument("--lamp-lines", required=True)
    parser.add_argument("--bundle", type=int, required=True)
    parser.add_argument("--broken-fibers", type=str, default="")
    args = parser.parse_args()

    print(f"--- Spot Selection Diagnostics for Bundle {args.bundle} ---")
    
    # We need a dummy opts object for load_python_psf and read_preproc
    class Opts:
        def __init__(self, arc, psf):
            self.arc_image_filename = arc
            self.input_psf_filename = psf
            
    opts = Opts(args.arc, args.in_psf)
    
    # Load data
    ddata = read_preproc(opts)
    image = ddata['image'].T
    weight = ddata['ivar'].T
    
    psf = load_python_psf(args.in_psf, opts)
    lamp_lines = read_lamp_lines(args.lamp_lines)
    
    # Setup bundle fiber range
    bundle = psf.params_of_bundles[args.bundle]
    fmin, fmax = bundle.fiber_min, bundle.fiber_max
    print(f"Analyzing fibers {fmin} to {fmax}...")

    # Determine broken fibers
    broken_fibers = [int(f) for f in args.broken_fibers.split(",") if f.strip()] if args.broken_fibers else []

    # Call the actual selection logic
    spots = get_bundle_spots(
        psf, fmin, fmax, lamp_lines, 
        image=image, weight=weight, 
        broken_fibers=broken_fibers
    )

    print(f"\nFINAL RESULT: {len(spots)} spots selected.")
    
    # For debugging, let's print some stats about the spots
    if spots:
        waves = np.array([s['wave'] for s in spots])
        print(f"Wavelength range: {waves.min():.2f} - {waves.max():.2f} A")
        print(f"Avg S/N: {np.mean([s['snr'] for s in spots]):.2f}")

if __name__ == "__main__":
    main()

#!/bin/bash

# Exit immediately if a command exits with a non-zero status
set -e

# Source the environment setup
if [ -f "env_setup.sh" ]; then
    echo "Sourcing env_setup.sh..."
    source env_setup.sh
else
    echo "Warning: env_setup.sh not found in the current directory."
fi

# Default values for optional arguments
LAMP_LINES_C="${SPECEX}/py/specex/data/specex_linelist_desi.txt"
LAMP_LINES_PY="/global/cfs/cdirs/desi/users/cdwarner/code/specex/py/specex/da/ta/specex_linelist_desi.txt"
FIRST_BUNDLE=5
LAST_BUNDLE=5
FIRST_FIBER=125
LAST_FIBER=149
BROKEN_FIBERS="473,474"
DEG_WAVE=3

# Usage helper
usage() {
    echo "Usage: $0 -a <preproc_fits> -i <in_psf_fits> -o <outdir> [options]"
    echo ""
    echo "Required:"
    echo "  -a, --preproc        Path to preproc FITS file"
    echo "  -i, --in-psf         Path to input PSF FITS file"
    echo "  -o, --outdir         Output directory (e.g., \$SCRATCH)"
    echo ""
    echo "Options:"
    echo "  --lamp-lines-c       Lamp lines file for C++ (Default: \$SPECEX/.../specex_linelist_desi.txt)"
    echo "  --lamp-lines-py      Lamp lines file for Python (Default: /global/cfs/.../specex_linelist_desi.txt)"
    echo "  --first-bundle       Default: 5"
    echo "  --last-bundle        Default: 5"
    echo "  --first-fiber        Default: 125"
    echo "  --last-fiber         Default: 149"
    echo "  --broken-fibers      Default: 473,474"
    echo "  --deg-wave           Legendre degree wave (Default: 3)"
    exit 1
}

# Parse command line arguments
PARSED_ARGUMENTS=$(getopt -a -n run_psf_fits -o a:i:o: --long preproc:,in-psf:,outdir:,lamp-lines-c:,lamp-lines-py:,first-bundle:,last-bundle:,first-fiber:,last-fiber:,broken-fibers:,deg-wave: -- "$@")
VALID_ARGUMENTS=$?
if [ "$VALID_ARGUMENTS" != "0" ]; then
    usage
fi

eval set -- "$PARSED_ARGUMENTS"

while :
do
  case "$1" in
    -a | --preproc)        PREPROC="$2"      ; shift 2 ;;
    -i | --in-psf)         IN_PSF="$2"       ; shift 2 ;;
    -o | --outdir)         OUTDIR="$2"       ; shift 2 ;;
    --lamp-lines-c)       LAMP_LINES_C="$2" ; shift 2 ;;
    --lamp-lines-py)      LAMP_LINES_PY="$2"; shift 2 ;;
    --first-bundle)       FIRST_BUNDLE="$2" ; shift 2 ;;
    --last-bundle)        LAST_BUNDLE="$2"  ; shift 2 ;;
    --first-fiber)        FIRST_FIBER="$2"  ; shift 2 ;;
    --last-fiber)         LAST_FIBER="$2"   ; shift 2 ;;
    --broken-fibers)      BROKEN_FIBERS="$2"; shift 2 ;;
    --deg-wave)           DEG_WAVE="$2"     ; shift 2 ;;
    --) shift; break ;;
    *) echo "Unexpected option: $1"; usage ;;
  esac
done

# Validation
if [ -z "$PREPROC" ] || [ -z "$IN_PSF" ] || [ -z "$OUTDIR" ]; then
    echo "Error: Missing required arguments (-a, -i, or -o)."
    usage
fi

# Dynamically construct the file tag (e.g., extracts 'z8-00344649' from the input string)
# If the filename format shifts, this falls back safely to a default tag string.
if [[ $(basename "$PREPROC") =~ preproc-(.*)\.fits\.gz ]]; then
    FILE_TAG="${BASH_REMATCH[1]}"
else
    FILE_TAG="extracted-psf"
fi

# Construct output filenames
OUT_PSF_C="${OUTDIR}/fit-psf-${FILE_TAG}_05.fits"
OUT_PSF_PY="${OUTDIR}/pyfit-psf-${FILE_TAG}_05.fits"

echo "========================================="
echo "Running C++ desi_psf_fit..."
echo "========================================="
desi_psf_fit \
  -a "$PREPROC" \
  --in-psf "$IN_PSF" \
  --lamp-lines "$LAMP_LINES_C" \
  --out-psf "$OUT_PSF_C" \
  --first-bundle "$FIRST_BUNDLE" \
  --last-bundle "$LAST_BUNDLE" \
  --first-fiber "$FIRST_FIBER" \
  --last-fiber "$LAST_FIBER" \
  --legendre-deg-wave "$DEG_WAVE" \
  --fit-continuum \
  --broken-fibers "$BROKEN_FIBERS"

echo "========================================="
echo "Running Python specex.specex (GPU)..."
echo "========================================="
python -m specex.specex \
  -a "$PREPROC" \
  --in-psf "$IN_PSF" \
  --lamp-lines "$LAMP_LINES_PY" \
  --out-psf "$OUT_PSF_PY" \
  --first-bundle "$FIRST_BUNDLE" \
  --last-bundle "$LAST_BUNDLE" \
  --first-fiber "$FIRST_FIBER" \
  --last-fiber "$LAST_FIBER" \
  --legendre-deg-wave "$DEG_WAVE" \
  --fit-continuum \
  --broken-fibers "$BROKEN_FIBERS" \
  --gpu 4

echo "Both fits completed successfully."

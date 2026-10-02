#!/bin/bash

# Setup environment for specex on Perlmutter
# This script should be sourced: source env_setup.sh

# Load necessary modules if on a compute node
if command -v module &> /dev/null; then
    module load cudatoolkit
fi

# Set PYTHONPATH using absolute paths, derived from this script's own
# location so the same env_setup.sh works whether sourced from the
# Perlmutter CFS checkout or a local clone on a personal machine.
BASE_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
export PYTHONPATH=$BASE_DIR/py:$BASE_DIR/build:$PYTHONPATH

export NVLIBS=$(python -c "
import glob, os
site = [p for p in __import__('sys').path if 'site-packages' in p][0]
dirs = glob.glob(os.path.join(site, 'nvidia', '*', 'lib'))
print(':'.join(dirs))
")

export LD_LIBRARY_PATH=$NVLIBS:$LD_LIBRARY_PATH

echo "Environment setup complete."
echo "PYTHONPATH set to: $PYTHONPATH"
echo "LD_LIBRARY_PATH set to: $LD_LIBRARY_PATH"

# Check for GPU availability
if command -v nvidia-smi &> /dev/null; then
    echo "GPU detected:"
    nvidia-smi --list-gpus
else
    echo "No GPU detected (Login node?)"
fi

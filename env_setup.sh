#!/bin/bash

# Setup environment for specex on Perlmutter
# This script should be sourced: source env_setup.sh

# Load necessary modules if on a compute node
if command -v module &> /dev/null; then
    module load cudatoolkit
fi

# Set PYTHONPATH using absolute paths
BASE_DIR=/global/cfs/cdirs/desicollab/users/cdwarner/code/specex
export PYTHONPATH=$BASE_DIR/py:$BASE_DIR/build:$PYTHONPATH

echo "Environment setup complete."
echo "PYTHONPATH set to: $PYTHONPATH"

# Check for GPU availability
if command -v nvidia-smi &> /dev/null; then
    echo "GPU detected:"
    nvidia-smi --list-gpus
else
    echo "No GPU detected (Login node?)"
fi

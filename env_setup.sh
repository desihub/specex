#!/bin/bash

# Setup environment for specex on Perlmutter
# This script should be sourced: source env_setup.sh

# Load necessary modules if on a compute node or if modules are available
if command -v module &> /dev/null; then
    module load cudatoolkit
    module load pytorch # Often contains useful GPU libs, or load specific ones
    # Add any other specific DESI/Perlmutter modules here
fi

# Set PYTHONPATH to include current build and py directories
export PYTHONPATH=$(pwd)/py:$(pwd)/build:$PYTHONPATH

echo "Environment setup complete."
echo "PYTHONPATH set to: $PYTHONPATH"

# Check for GPU availability
if command -v nvidia-smi &> /dev/null; then
    echo "GPU detected:"
    nvidia-smi --list-gpus
else
    echo "No GPU detected (Login node?)"
fi

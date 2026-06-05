import os
import sys

# 1. MPI setup and GPU environment variables MUST come before ANY other imports
from mpi4py import MPI
comm = MPI.COMM_WORLD
rank = comm.Get_rank()
size = comm.Get_size()

# Assign GPU based on rank (assume 4 GPUs per node)
gpu_id = rank % 4
os.environ["CUDA_VISIBLE_DEVICES"] = str(gpu_id)

# Prevent JAX from pre-allocating all VRAM
os.environ["XLA_PYTHON_CLIENT_PREALLOCATE"] = "false"

import jax
import jax.numpy as jnp

print(f"Rank {rank}/{size}: Assigned GPU {gpu_id}. JAX devices: {jax.devices()}")

# Verify work
x = jnp.ones((1000, 1000))
y = jnp.dot(x, x)
print(f"Rank {rank}: Computation successful on {y.device}")

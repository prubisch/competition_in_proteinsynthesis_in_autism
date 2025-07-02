#!/bin/bash

# Parameters for test run
L='2'
N_AMINO=1000
N_RIBO=100
N_MRNA=50
L1=10
L2=25

#for simulations run on cluster
#source /appl/MCell/Blender-2.93-CellBlender/blender_2.93_environment


python3.9 run_sim.py "$L" "$N_MRNA" "$N_RIBO" "$N_AMINO" "$L1" "$L2"
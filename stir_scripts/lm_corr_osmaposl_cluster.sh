#!/bin/bash

cd ${PBS_O_WORKDIR}
module load stir/180202-d0c9261
./run_root_GATE.sh
correct_projdata correct_projdata.par
time mpirun OSMAPOSL OSMAPOSL.par


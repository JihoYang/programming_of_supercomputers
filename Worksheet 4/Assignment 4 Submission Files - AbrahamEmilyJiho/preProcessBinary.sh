#!/bin/bash 
#@ wall_clock_limit = 00:20:00
#@ job_name = ws4-mpiIO-intel-createBinary-sb
#@ job_type = MPICH
#@ output = cannon_64_$(jobid).out
#@ error = cannon_64_$(jobid).out
#@ class = test
#@ node = 4
#@ total_tasks = 64
#@ node_usage = not_shared
#@ energy_policy_tag = cannon
#@ minimize_time_to_solution = yes
#@ notify_user = jiho.yang@tum.de
#@ notification = never
#@ island_count = 1
#@ queue

. /etc/profile
. /etc/profile.d/modules.sh

export LANG=en_US.utf8
export LC_ALL=en_US.utf8

module unload mpi.intel
module unload mpi.ibm
module load mpi.intel/2018

mpiexec -n 64 ./generateBinaryInput cannon_matrices/64x64-1.in cannon_matrices/64x64-2.in >> preProcessResults.txt
mpiexec -n 64 ./generateBinaryInput cannon_matrices/128x128-1.in cannon_matrices/128x128-2.in >> preProcessResults.txt
mpiexec -n 64 ./generateBinaryInput cannon_matrices/256x256-1.in cannon_matrices/256x256-2.in >> preProcessResults.txt
mpiexec -n 64 ./generateBinaryInput cannon_matrices/512x512-1.in cannon_matrices/512x512-2.in >> preProcessResults.txt
mpiexec -n 64 ./generateBinaryInput cannon_matrices/1024x1024-1.in cannon_matrices/1024x1024-2.in >> preProcessResults.txt
mpiexec -n 64 ./generateBinaryInput cannon_matrices/2048x2048-1.in cannon_matrices/2048x2048-2.in >> preProcessResults.txt
mpiexec -n 64 ./generateBinaryInput cannon_matrices/4096x4096-1.in cannon_matrices/4096x4096-2.in >> preProcessResults.txt

mv a_* cannon_matrices_binary
mv b_* cannon_matrices_binary


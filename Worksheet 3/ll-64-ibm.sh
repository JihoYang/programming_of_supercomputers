#!/bin/bash 
#@ wall_clock_limit = 00:20:00
#@ job_name = pos-cannon-mpi-ibm
#@ job_type = Parallel
#@ output = cannon_64_$(jobid).out
#@ error = cannon_64_$(jobid).out
#@ class = test
#@ node = 4
#@ total_tasks = 64
#@ node_usage = not_shared
#@ energy_policy_tag = cannon
#@ minimize_time_to_solution = yes
#@ notification = never
#@ island_count = 1
#@ queue

. /etc/profile
. /etc/profile.d/modules.sh

module load mpi.intel
module load mpi.ibm

date > results.txt

for ((i=0; i<15; i++))
do 
    mpiexec -n 64 ./cannon cannon_matrices/64x64-1.in cannon_matrices/64x64-2.in >> results.txt
    date >> results.txt
done

for ((i=0; i<15; i++))
do 
    mpiexec -n 64 ./cannon cannon_matrices/128x128-1.in cannon_matrices/128x128-2.in >> results.txt
    date >> results.txt
done

for ((i=0; i<15; i++))
do 
    mpiexec -n 64 ./cannon cannon_matrices/256x256-1.in cannon_matrices/256x256-2.in >> results.txt
    date >> results.txt
done

for ((i=0; i<15; i++))
do 
    mpiexec -n 64 ./cannon cannon_matrices/512x512-1.in cannon_matrices/512x512-2.in >> results.txt
    date >> results.txt
done

for ((i=0; i<15; i++))
do 
    mpiexec -n 64 ./cannon cannon_matrices/1024x1024-1.in cannon_matrices/1024x1024-2.in >> results.txt
    date >> results.txt
done

for ((i=0; i<15; i++))
do 
    mpiexec -n 64 ./cannon cannon_matrices/2048x2048-1.in cannon_matrices/2048x2048-2.in >> results.txt
    date >> results.txt
done

for ((i=0; i<15; i++))
do 
    mpiexec -n 64 ./cannon cannon_matrices/4096x4096-1.in cannon_matrices/4096x4096-2.in >> results.txt
    date >> results.txt
done

./collectResults results.txt > statistics.txt

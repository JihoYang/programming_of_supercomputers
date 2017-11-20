#!/bin/bash
#@ wall_clock_limit = 02:00:00
#@ job_name = pos-lulesh-flag-haswell
#@ job_type = MPICH
#@ class = test
#@ output = pos_lulesh_flag_haswell_$(jobid).out
#@ error = pos_lulesh_flag_haswell_$(jobid).out
#@ node = 1
#@ total_tasks = 28
#@ node_usage = not_shared
#@ energy_policy_tag = lulesh
#@ minimize_time_to_solution = yes
#@ island_count = 1
#@ notify_user = jiho.yang@tum.de
#@ notification = always
#@ queue
. /etc/profile
. /etc/profile.d/modules.sh

module list

# Array of flags
array='-march=native -xHost -unroll -ipo'
X=4
filename=flags$X

# With flags
for N in $array
do
	for j in {1..5}
	do

		echo ========================================================================================================================================= >> $filename
		echo Current flag : $N >> $filename
		echo ========================================================================================================================================= >> $filename

		./lulesh_$N >> $filename

	done
done

rm lulesh_*

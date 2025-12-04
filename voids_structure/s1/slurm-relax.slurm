#!/bin/bash

#SBATCH --job-name=V9
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=8
#SBATCH --gres=gpu:1
#SBATCH --time=120:00:00
#SBATCH --partition=xhhgnormal02
##SBATCH --mem=5G

# for local source install
source /work/home/acl7i3im69/soft/deepmd-kit/bin/activate
source /work/home/acl7i3im69/scripts/dp_env.sh

conda activate deepmd

export OMP_NUM_THREADS=2
export TF_INTER_OP_PARALLELISM_THREADS=4

T=200
alatt_ori=1.882
scale_ori=0.54
delta=0.0025
i=22
alatt=$(awk "BEGIN {print $alatt_ori*($scale_ori+$i*$delta)}")

step=1000
nfile=200

iseed=8713
lmp -i in.run-nvt -v alatt ${alatt} -v i ${i} -v T ${T} -v iseed ${iseed}

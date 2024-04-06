#!/bin/bash
#SBATCH -c 8
#SBATCH --mem=64GB
#SBATCH -J AGG_lite
#SBATCH -o AGG_lite.%j.out
#SBATCH -e AGG_lite.%j.err
#SBATCH -p intel-sc3,amd-ep2
#SBATCH -q normal
#SBATCH --mail-type=ALL
#SBATCH --mail-user=liuyifang@westlake.edu.cn

export PATH=/storage/liuxiaodongLab/liuyifang/app/cellranger-7.2.0:$PATH

cellranger aggr --id=AGG_lite --csv=AGG_aggr.csv

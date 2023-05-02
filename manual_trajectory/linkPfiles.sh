#!/bin/bash

#Pfile_dir=/scratch/${USER}/LAGRANTO_runs/out_ARTofMELT_20200417
Pfile_dir=/scratch/sr9/testdata/ec.oper/20230428_00/cdf/

echo ${Pfile_dir}

ln -s ${Pfile_dir}/P* .
ln -s ${Pfile_dir}/ml_cst .

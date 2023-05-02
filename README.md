# TrajectoryScripts
this repository includes trajectory scripts for computing and plotting forecast trajectories 

The atos branch is prepared for running on Atos hpc/ecs and requires access to the lagranto installation of Moritz Pickl


Add the follwonf to profile or bashrc

###
# Lagranto settings
export model="ecmwf"
export mode="ive"
export DYN_TOOLS=/home/nemp/programs/eth_tools
export LAGRANTO=${DYN_TOOLS}/lagranto.${model}-${mode}
export PATH=$PATH:${DYN_TOOLS}/bin
export PATH=$PATH:${LAGRANTO}/bin
export PATH=$PATH:/home/nemp/programs/bin

LD_LIBRARY_PATH=/usr/local/lib
LD_LIBRARY_PATH=$LD_LIBRARY_PATH:"/home/nemp/programs/eth_tools/local/lib"
LD_LIBRARY_PATH=$LD_LIBRARY_PATH:"/usr/local/apps/netcdf4/4.7.4/GNU/8.4/lib"

export LD_LIBRARY_PATH

###














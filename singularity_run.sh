#!/usr/bin/env bash        

# Singularity run raser image     
# Author SHI Xin <shixin@ihep.ac.cn>  
# Created [2021-06-03 Thu 10:05] 


echo "Singularity run RASER ..."
#singularity shell --bind /raser -e raser.simg
#singularity shell raser.simg
singularity shell raser.simg -c "/bin/bash -rcfile singularity_bashrc.sh"


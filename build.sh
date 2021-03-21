#!/usr/bin/env bash

# Main driver to build RASER 
# Author SHI Xin <shixin@ihep.ac.cn>
# Created [2021-03-21 Sun 11:00]

usage() {
    printf "NAME\n\tbuild.sh - Main driver to build RASER\n"
    printf "\nSYNOPSIS\n"
    printf "\n\t%-5s\n" "./build.sh [OPTION]" 
    printf "\nOPTIONS\n" 
    printf "\n\t%-5s  %-40s\n"  "1"  "build raser docker image" 
    printf "\n\t%-5s  %-40s\n"  "2"  "build KDetSim"
    printf "\n\n" 
}

if [[ $# -eq 0 ]]; then
    usage
    echo "Please enter your option: "
    read option
else
    option=$1    
fi


case $option in 
    1) echo "Building raser docker image..."
       docker build -t raser .  
       ;;
    2) echo "Building KDetSim ..."
        cd KDetSim
        mkdir obj 
        make -f makeLinuxMacRoot6 KDetSimLinux   
       ;;
esac






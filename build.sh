#!/usr/bin/env bash

# Main driver to build RASER 
# Author SHI Xin <shixin@ihep.ac.cn>
# Created [2021-03-21 Sun 11:00]

usage() {
    printf "NAME\n\tbuild.sh - Main driver to build RASER\n"
    printf "\nSYNOPSIS\n"
    printf "\n\t%-5s\n" "./build.sh [OPTION]" 
    printf "\nOPTIONS\n" 
    printf "\n\t%-5s  %-40s\n"  "docker"  "build docker" 
    printf "\n\t%-5s  %-40s\n"  "raser"  "build raser"
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
    docker) echo "Building docker image ..."
       docker build -t raser .  
       ;;
    raser) echo "Building raser binary ... "
       make    
       ;;
    3) echo "Building binary file ...on Windows"
       docker run --rm -it -h raser -e DISPLAY=192.168.237.32:0  --mount type=bind,source=C:\\Users\\dell\\raser,target=/home/physicist raser bash  
       #make    
       ;;
esac






#!/usr/bin/env bash

# Main driver to build RASER 
# Author SHI Xin <shixin@ihep.ac.cn>
# Created [2021-03-21 Sun 11:00]

usage() {
    printf "NAME\n\tbuild.sh - Main driver to build RASER\n"
    printf "\nSYNOPSIS\n"
    printf "\n\t%-5s\n" "./build.sh [OPTION]" 
    printf "\nOPTIONS\n" 
    printf "\n\t%-5s  %-40s\n"  "1"  "build docker image" 
    printf "\n\t%-5s  %-40s\n"  "2"  "build binary file"
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
    1) echo "Building docker image ..."
       docker build -t raser .  
       ;;
    2) echo "Building binary file ... on Mac"
       docker run --rm -it -h raser  -v /tmp/.X11-unix:/tmp/.X11-unix -e DISPLAY=$ip:0 --mount type=bind,source=$HOME/raser,target=/home/physicist raser bash  
       make    
       ;;
    3) echo "Building binary file ...on Windows"
       docker run --rm -it -h raser -e DISPLAY=10.0.75.1:0  --mount type=bind,source=C:\\Users\\hello\\raser,target=/home/physicist raser bash  
       make    
       ;;
esac






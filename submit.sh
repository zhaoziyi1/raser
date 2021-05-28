#!/bin/bash 
#Raser run step
#FEnics solver the electric field
#The electrics get the waveform 
#Get the time resolution

option=$1
Time=2021_5_28
case $option in 
0.1) echo "raser 2D SiC basic information"
	python3 python/raser.py 2D output
	;;

0.2) echo "raser 2D scan for time resolution"
	mkdir out/fenics_${Time}/ -p
	rm out/fenics_${Time}/*
	python3 python/raser.py 2D_scan out/fenics_${Time}/
	;;

0.3) echo "data processing"
	mkdir result/sic_${Time}_result/ -p
	rm result/sic_${Time}_result/*
	python3 python/raser_read_sample.py out/fenics_${Time}/ result/sic_${Time}_result/ V
	;;

0.4) echo "get time resolution"
	mkdir result/root -p
	python3 python/add_noise_raser.py result/sic_${Time}_result/ result/root/sic_${Time}_result.root
	esac
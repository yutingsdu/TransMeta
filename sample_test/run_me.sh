#!/bin/bash


## Please first dowload the demo data from https://sourceforge.net/projects/transmeta/files/DemoData/

# !!!!!!!!!!!!!! IMPORTANT !!!!!!!!!!!!!!!!!!

#export LD_LIBRARY_PATH=/home/yuting/yuting/Software/local/boost/lib:$LD_LIBRARY_PATH
export LD_LIBRARY_PATH=/home/yuting/yuting/Software/local/boost_1_60/lib:$LD_LIBRARY_PATH

../TransMeta -B bamlist -s first -o transmeta_oudir -p 2

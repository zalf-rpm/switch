#!/bin/bash

export MONICA_PARAMETERS=/home/berg/GitHub/monica-parameters 
/home/berg/GitHub/monica/_cmake_linux_release/monica-zmq-proxy -pps -f 6666 -b 6677 &
/home/berg/GitHub/monica/_cmake_linux_release/monica-zmq-proxy -pps -f 7788 -b 7777 &
for i in `seq $2 $3`
do 
	#echo $i
	/home/berg/GitHub/monica/_cmake_linux_release/monica-zmq-server -ci -i tcp://localhost:6677 -co -o tcp://localhost:7788 &
done


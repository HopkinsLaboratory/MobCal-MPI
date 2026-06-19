#!/bin/bash

#Finds number of processors in the system
cat /proc/cpuinfo | grep "processor" | wc -l

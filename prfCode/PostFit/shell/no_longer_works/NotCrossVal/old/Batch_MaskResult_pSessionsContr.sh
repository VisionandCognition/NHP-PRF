#!/bin/bash

for pv in 0.1 0.2 0.5; do
	for rth in 1 2 4 5 10; do
		./MaskResult_pSessionsContr.sh $pv $rth danny 
		./MaskResult_pSessionsContr.sh $pv $rth eddy 
	done
done

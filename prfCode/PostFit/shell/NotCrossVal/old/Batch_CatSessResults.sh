#!/bin/bash
for rth in 0 1 2 4 5 10; do
	./CatSessResults.sh $rth danny 
	./CatSessResults.sh $rth eddy
done 
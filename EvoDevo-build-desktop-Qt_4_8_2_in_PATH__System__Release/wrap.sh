#!/bin/bash

#declare an array of the migration rate values:
array=(1e-5 1e-4 5e-4 7e-4 1e-3 2e-3  5e-3 7e-3 1e-2 2e-2  3e-2  4e-2  5e-2 6e-2 7e-2  8e-2 1e-1 2e-1 5e-1)


	echo "Migration rate	Location1:1 	Location 2:	Between location 1&2" >> foo.txt
for i in "${array[@]}";
do
	
	./EvoDevo $i >> foo.txt
done


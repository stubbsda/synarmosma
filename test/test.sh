#!/bin/bash

modules="Rational Polynomial Lattice Poset Graph Pseudograph Directed_Graph Group Integer_Matrix"
np=0

for case in $modules
do
	echo "Testing the class $case..."
	./unit_test $case
	if [ $? -eq 0 ]; 
	then 
		echo "Passed!" 
		let np++
	else 
		echo "Failed!" 
	fi
done 
echo "Passed $np of 9 tests"


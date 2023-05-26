#!/usr/bin/env csh

seqkit fx2tab --length --name --header-line $1 > ../res/lengths_double.txt

while IFS="	" read -r plas lengths
do
	echo "The length is ${lengths} bp"
	half=$((lengths/2))
	echo "The half length is ${half} bp"
	quat=$((half/2))
	pos1=$((half-quat))
	pos2=$((half+quat))
	echo "The positions are 1: ${pos1} and 2: ${pos2} bp"
	sam_command=$(samtools coverage --ff 0 -H -r ${plas}:${pos1}-${pos2} $2 >> $3)
	echo "$sam_command"
	samtools index $2
	${sam_command}
done < ../res/lengths_double.txt
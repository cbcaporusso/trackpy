#!/usr/bin/env bash

filepath=$1
time=$2
DT=$3

CURR_DIR=$(pwd)

echo "Computing displacement field from Trj..."

cd "$filepath" || exit
mkdir -p Dspl_dt_${DT}

file="Trj/xyz.dump.${time}"
if [ ! -f Dspl_dt_${DT}/xyz.dump.${time}.dspl ]; then 
	if [ $((time%DT)) -eq 0 ]; then
		time_next=$((time+DT))
		#echo $file Trj/xyz.dump.${time_next}
		awk -f $CURR_DIR/displacement.awk -v dt=$DT $file Trj/xyz.dump.${time_next} > Dspl_dt_${DT}/xyz.dump.${time}.dspl
	fi
fi


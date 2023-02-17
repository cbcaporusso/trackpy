#!/bin/bash

FILE_PATH=$1
TIME=$2

curr_dir=`pwd` 
vor_dir="/home/ub35/ub35873/voro++-0.4.6/src" ### questa è la dir con l'eseguibile voro++ ###

echo $FILE_PATH $TIME
cd "$FILE_PATH" || exit 
pwd	

###################
# rm -r local_hexatic
###################

mkdir -p local_hexatic
first_xyz=$(ls -v Trj/xyz.dump.* | head -1)

if [ ! -f $first_xyz ];then
	echo "Configuration not found"
fi

n_beads=$(awk < $first_xyz '{if(NR==4) print $1}')
min_box=$(awk < $first_xyz '{if(NR==6) print $1}')
max_box=$(awk < $first_xyz '{if(NR==6) print $2}')

sizeX=$(awk < $first_xyz '{if(NR==6) print $2-$1}')
sizeY=$(awk < $first_xyz '{if(NR==7) print $2-$1}')

echo "Lx=${sizeX} Ly=${sizeY}" 

for conf in $(ls -rt Trj/xyz.dump.${TIME})
do
	conf=$(basename $conf)
	### analysis ###
	time=$(echo ${conf} | sed 's/xy.frame.//g' | sed 's/xy.log.frame.//' | sed 's/xyz.log.dump.//g' | sed 's/xyz.dump.//')
	
	if [ ! -f local_hexatic/${conf}.hexatic ] ; then
		echo " >  ${conf}.hexatic"

		awk < Trj/${conf} '{if(NR>=10) print $1,$2,$3,'0.0'}' > conf

		${vor_dir}/voro++ -c "%i %x %y %v %n" -o -px -py  0 $sizeX 0 $sizeY -0.5 0.5 conf

		awk < conf.vol -v L=$size -f $curr_dir/elaborate.awk > local_hexatic/${conf}.hexatic
		#rm elaborate.awk conf conf.vol

	fi

done



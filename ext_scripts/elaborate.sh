#!/bin/bash

SCRIPT_PATH=$1
FILE_PATH=$2
TIME=$3

curr_dir=`pwd` 

# check if voro++ is installed
if ! command -v voro++ &> /dev/null
then
	echo "voro++ could not be found"
	echo "Please install voro++ and try again"
	exit
fi

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

		voro++ -c "%i %x %y %v %n" -o -px -py  0 $sizeX 0 $sizeY -0.5 0.5 conf

		awk < conf.vol -v L=$size -f $SCRIPT_PATH/elaborate.awk > local_hexatic/${conf}.hexatic
		#rm elaborate.awk conf conf.vol

	fi

done



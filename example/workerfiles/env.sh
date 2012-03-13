

export LD_LIBRARY_PATH=/afs/crc.nd.edu/x86_64_linux/Modules/tcltk/current/lib/:/afs/crc.nd.edu/x86_64_linux/Modules/tcltk/current/lib/tclx8.4/:$LD_LIBRARY_PATH


echo '++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++'
ldd /afs/crc.nd.edu/x86_64_linux/Modules/3.2.6/bin/modulecmd
echo '^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^'


echo '++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++'
env
echo '^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^'

MODULESHOME=/afs/crc.nd.edu/x86_64_linux/Modules/3.2.6
function module()
{
    eval $($MODULESHOME/bin/modulecmd sh $*)
}

MODULEFILES=(/opt/crc/Modules/modulefiles /afs/nd.edu/user37/ccl/software/modulefiles /afs/crc.nd.edu/user/c/cabdulwa/Public/modulefiles /afs/crc.nd.edu/user/i/izaguirr/Public/modulefiles)
for m in ${MODULEFILES[@]}; do
	echo "Using module system in $m"
	module use $m
done

modules=(epd gromacs/4.5.3 msmbuilder/env2)

for m in ${modules[@]}; do
	echo "Loaded module $m"
	module load $m
done

export OMP_NUM_THREADS=1



set -o errexit
set -o verbose

CONF_IN=structure.pdb
CONF_OUT=structure2.pdb
ASSIGNMENT=cell2.dat
DESIRED_FILES="$CONF_OUT $ASSIGNMENT"
RESULTFILE=results.tar
CLEANUP="traj* *.tpr"

### disable gmx automatic backups.
export GMX_MAXBACKUP=-1
### access the topologies sent to the worker
export GMXLIB=$PWD/gmxtopologies
### assume each worker is allocated one processor
NPROCS=1

puts() {
	echo "================================================================================"
	echo "[worker] $@"
}

prelude() {
        puts "Prelude: current environment"
	env
	puts "Prelude: uname -a"
	uname -a
	puts "Prelude: bash options"
	set -o
	echo
}

prepare-filenames() {
	structure=`ls $CONF_IN.*`
	id="${structure##*.}"
	mv -f "$CONF_IN.$id" $CONF_IN
	#mv -f "$WALKER.$id" $WALKER
}

check-initial() {
	puts "Initial file listing"
	ls
	echo
}

run-md() {
	puts "Running simulation"
	python simulate.py -s 100000
	echo
}

assign() {
    echo "Running assignment"
}

check-result() {
        puts "Generated files:"
	ls
	puts "Checking if result files ($DESIRED_FILES) exist"
	ls $CONF_OUT $ASSIGNMENT
	echo
}

package() {
	puts "Packaging results"
	tar cvf "$RESULTFILE.$id" $CONF_OUT $ASSIGNMENT $WALKER
	ls "$RESULTFILE.$id"
	echo
}

cleanup() {
	puts "Cleaning up"
	#rm -rv $CLEANUP
	echo
}

MODULESHOME=/afs/crc.nd.edu/x86_64_linux/Modules/3.2.6
module() {
	eval $(MODULESHOME/bin/modulecmd sh $*)
}

MODULEFILES=(/afs/crc.nd.edu/user/j/jkinniso/Public/modulefiles)
for m in ${MODULEFILES[@]}; do
	echo "Using modulefiles at $m"
	module use $m
done

modules=(anaconda3/2.3)
for m in ${modules[@]}; do
	echo "Loading module $m"
	module load $m
done

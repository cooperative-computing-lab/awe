
set -o errexit
set -o verbose

CONF_IN=structure.pdb
CONF_OUT=structure2.pdb
ASSIGNMENT=cell2.dat
DESIRED_FILES="$CONF_OUT $ASSIGNMENT"
RESULTFILE=results.tar
WALKER=walker.pkl
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

check-initial() {
	puts "Initial file listing"
	ls
	echo
}

run-md() {
	puts "Running simulation"
	echo $GMXLIB
	./pdb2gmx -f $CONF_IN -ff amber96 -water none
	./grompp -f sim.mdp
	./mdrun -s topol.tpr -c $CONF_OUT -nt $NPROCS
	echo
}

assign() {
	puts "Assigning trajectory"
	./awe-assign cells.dat CellIndices.dat traj.xtc StructureIndices.dat $ASSIGNMENT
	echo
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
	tar cvf $RESULTFILE $CONF_OUT $ASSIGNMENT $WALKER
	ls $RESULTFILE
	echo
}

cleanup() {
	puts "Cleaning up"
	rm -rv $CLEANUP
	echo
}


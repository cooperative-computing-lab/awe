

CONF_IN=structure.pdb
CONF_OUT=structure2.pdb
ASSIGNMENT=cell2.dat
RESULTFILE=results.tar
WALKER=walker.pkl
CLEANUP="Data output Trajectories ProjectInfo.h5 frame* *.tpr"

source ~/.bash_modules
module load gromacs/4.5.3
module load msmbuilder/lcls/2.1.1
export OMP_NUM_THREADS=1

puts() {
	echo "================================================================================"
	echo "[worker] $@"
}

check-initial() {
	puts "Initial file listing"
	ls
	echo
}

run-md() {
	puts "Running simulation"
	pdb2gmx -f $CONF_IN -ff amber96 -water none
	grompp -f sim.mdp
	mdrun -s topol.tpr -c $CONF_OUT -deffnm frame0
	echo
}

setup-msmbuilder() {
	puts "Prepping for MSMBuilder"
	mkdir -p output/traj
	mv *.xtc output/traj
	mkdir -v Data
	cp -v Gens.lh5 Data
	echo
}

assign() {
	puts "Assigning trajectory"
	ConvertDataToHDF.py -s state0.pdb -I output
	Assign.py
	ConvertAssignToText.py
	tail -1 Data/discrete.traj > $ASSIGNMENT
	echo
}

check-result() {
	puts "Checking if result files exist"
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


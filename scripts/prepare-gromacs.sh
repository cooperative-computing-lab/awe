#!/usr/bin/env bash

ME=$(basename $0)
PREFIX=awedata

GMXBINARIES=(pdb2gmx grompp mdrun)

usage() {
	cat <<EOF
USAGE: $ME [prefix]

  prefix: where the output is stored (default: awedata)

SUMMARY:

  Scan the environment for gromacs and prepare binaries and
  forcefields for AWE use. Run this on the master node.

DESCRIPTION:

  AWE runs gromacs binaries (pdb2gmx, grompp) on the workers. In order
  for these commands to work, they require the forcefield parameters
  to be transferred as well.

  $ME scans the environemt for gromacs, then populates the <prefix>
  directory with the necessary files.

EOF
}

die() {
	echo "$1"
	exit $2
}

check-args() {
	case $1 in

		-h|--help|-?)
			usage
			exit 1
			;;

		*)
			[ ! -z $1 ] && PREFIX="$1"
			;;

	esac
}


check-gmx() {
	for b in ${GMXBINARIES[@]}; do
		which $b >/dev/null 2>&1
		[ $? -eq 0 ] || die "Cannot find $b" 2
	done
}

abspath() {
	case `uname` in
		Linux)
			readlink -f $@
			;;
		*)
			echo "Cannot determin `abspath` on $(uname)"
			exit 2
			;;
	esac
}

find-gmxroot() {
	abspath $(which pdb2gmx)/../..
}

find-gmxtop() {
	if [ ! -z $GMXLIB ]; then
		echo $GMXLIB
	else
		local top=$(find-gmxroot)/share/gromacs/top
		[ file "$top" ] || die "Cannot find gromacs topology directory" 2
		echo $top
	fi
}

copy-gmx() {
	local prefix=$1
	local path=$prefix/binaries/`uname`-`arch`
	local gmxtop=$(find-gmxtop)
	local mygmxtop=$prefix/gmxtopologies

	[ ! -d $path ] && mkdir -vp $path

	for b in ${GMXBINARIES[@]}; do
		cp -v `which $b` $path
	done

	echo "\`$gmxtop\` -> \`$mygmxtop\`"
	rsync -r $gmxtop/ $mygmxtop

}


check-args $@
check-gmx
copy-gmx $PREFIX

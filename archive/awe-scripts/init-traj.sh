#!/usr/bin/env bash


set -o errexit

PARAMS=("8 50 150" "8 50 50" "8 50 10" "8 50 5" "8 50 0.5" "8 20 150" "8 20 50" "8 20 10" "8 20 5" "8 20 0.5" "17 10 150" "17 10 50" "17 10 10" "17 10 5" "17 10 0.5")


compute-rediag-freq() {
    dt=$1
    diagfreq_ps=$2
    python -c "print $diagfreq_ps * 10.**3 / $dt"
}


compute-output-freq() {
    # 2 ps
    dt=$1
    python -c "print 2 * 10.**3 / $dt"
}

compute-traj-output-freq() {
	# 20 ps
	dt=$1
	python -c "print 20 * 10.**3 / $dt"
}

compute-numsteps() {
    # 100 ns
    dt=$1
    python -c "print 100 * 10.**6 / $dt"
}


home="$1"
gendir="$2"
run=$3
clone=$4

me=$(basename $0)


echo "[$me] $home $gendir $run $clone"

[ ! -d $gendir ] &&
mkdir -vp $gendir

cp -rv $home/datafiles-init-all/* $gendir

simconf=$gendir/protomol.conf
cp $home/template-protomol.conf $simconf

read -ra params <<< ${PARAMS[$run]}
modes=${params[0]}
dt=${params[1]}
rediag=$(compute-rediag-freq $dt ${params[2]})
numsteps=$(compute-numsteps $dt)
outputfreq=$(compute-output-freq $dt)
trajoutputfreq=$(compute-traj-output-freq $dt)

echo "[$me] Params: ${params[@]}"

echo "[$me] Setting REGEX_OUTPUTFREQ: $outputfreq"
sed -i "s:REGEX_OUTPUTFREQ:$outputfreq:" $simconf

echo "[$me] Setting REGEX_TRAJOUTPUTFREQ: $trajoutputfreq"
sed -i "s:REGEX_TRAJOUTPUTFREQ:$trajoutputfreq:" $simconf

echo "[$me] Setting REGEX_REDIAGFREQ: $rediag"
sed -i "s:REGEX_REDIAGFREQ:$rediag:" $simconf

echo "[$me] Setting REGEX_NUMMODES: $modes"
sed -i "s:REGEX_NUMMODES:$modes:" $simconf

echo "[$me] Setting REGEX_TIMESTEP: $dt"
sed -i "s:REGEX_TIMESTEP:$dt:" $simconf

echo "[$me] Setting REGEX_NUMSTEPS: $numsteps"
sed -i "s:REGEX_NUMSTEPS:$numsteps:" $simconf



echo "[$me] Setting REGEX_SEEDME_*"

for i in 1 2; do
    rand=$RANDOM
    sed -i "s:REGEX_SEEDME_$i:$rand:" $simconf
done


mappings=$home/CONFIGS.mapping
if [ $run -eq 0 ]; then
	[ -f $mappings ] && rm $mappings
fi
message="RUN $run modes $modes dt $dt rediag $rediag"
echo "[$me] Prepared: $message"
echo $message  >> $mappings

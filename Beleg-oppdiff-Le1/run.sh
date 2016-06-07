#!/bin/bash
# This script calls ULF for different velocities and extracts Tmax, the scalar diffusion rate Chi_stoich and the strain rate k

#velocitiesnotworking=( '16.0', '8.0' '4.0' '3.0', '2.7' '2.65' '1.75' '1.25' '0.30' '0.35')
velocities=( '0.0125' '0.05' '0.10' '0.15' '0.20' '0.25' '0.5' '0.75' '1.0' '1.5' '2.5')

if grep -q '#include smooke-changes.ulf' oppdifJet_ct_CH4Air_smooke.ulf; then
    echo 'Already include-path added in ulf settings.'
fi

for curVel in ${velocities[*]}; do
	echo "====== Run for U_LOWER=-U_UPPER=$curVel ======"
	mv './oppdifJet_ct_CH4Air_smooke_final.ulf' './oppdifJet_ct_CH4Air_smooke_final.ulf.bak'

	# write out current velocity and chemical setup (chemical setup actually doesn't change)
	cat > smooke-changes.ulf << EOF2
#define U_LOWER  $curVel; // velocity of the Fluid expelled from the upper jet
#define U_UPPER -$curVel; // velocity for other jet
EOF2

    ulf -f oppdifJet_ct_CH4Air_smooke.ulf | tee "v_${curVel}.log"

    while test ! -f ./oppdifJet_ct_CH4Air_smooke_final.ulf; do
		sleep 0.1
	done

#    k=$(printf "print ($j+$j)/0.1.\n" | python)

	mkdir -p './results'
    mv 'oppdifJet_ct_CH4Air_smooke_final.ulf' "./results/v_${curVel}.ulf"

done

echo '#velocity of both nozzles in m/s' > results.dat
for curVel in ${velocities[*]}; do
	echo $curVel >> results.dat

# python plotScurve.py

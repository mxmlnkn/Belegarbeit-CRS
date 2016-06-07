#!/bin/bash
chist=( '0.01' '0.1' '0.25' '0.4918' '0.8800' '1.2686' '1.6383' '3.7859' '1.9816' )
#chist=( '3.7859' )
# vel=( '0.05'   '0.10'   '0.15'   '0.20'   '0.25'   '0.5'    ) # corresponding velocites for oppdifjet

if grep -q '#include smooke-changes.ulf' flameletLe1_CH4_air.ulf; then
    echo 'Already include-path added in ulf settings.'
fi

for curChi in ${chist[*]}; do
	echo "====== Run for chist=$curChi ======"
	printf "#define CHI_STOIC ${curChi};\n" > smooke-changes.ulf
	mv './flameletLe1Template_final.ulf' './flameletLe1Template_final.ulf.bak'
    ulf -f flameletLe1_CH4_air.ulf | tee "chist_${curChi}.log"
    while test ! -f ./flameletLe1Template_final.ulf; do
		sleep 0.1
	done
	mkdir -p './results'
    mv 'flameletLe1Template_final.ulf' "./results/chist_${curChi}.ulf"
done

#!/bin/bash
bash -c "cd 'Beleg-flamelet-Le1'  && python plotTests.py" & pid1=$!
bash -c "cd 'Beleg-flamelet-Levar'&& python ../Beleg-flamelet-Le1/plotTests.py" & pid2=$!
bash -c "cd 'Beleg-oppdiff-Le1'   && python plotTests.py" & pid3=$!
bash -c "cd 'Beleg-oppdiff-Levar' && python ../Beleg-oppdiff-Le1/plotTests.py" & pid4=$!
wait $pid1
wait $pid2
wait $pid3
wait $pid4
pdflatex beleg1.tex
bibtex beleg1.aux
pdflatex beleg1.tex
pdflatex beleg1.tex

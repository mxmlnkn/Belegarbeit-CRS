###################################  
# show flamelet solution for Le=1
###################################  

outputFile = 'flameletLe1_1.000.ulf'
RIFFile = './RIF/Le1/FLUT.dat'
# eps output
set terminal postscript eps size 12.0cm,8.0cm enhanced color \
    font 'Helvetica,20' linewidth 2

# plot temperature
set output 'flameletLe1_T.eps'
set xlabel 'Z'
set ylabel 'T in K'
plot outputFile u (column('Z')):(column('T')) w l title 'temperature' ,\
RIFFile u 6:24 w p t 'temperature RIF'

# plot enthalpy
set output 'flameletLe1_h.eps'
set xlabel 'Z'
set ylabel 'h in J/kg'
plot outputFile u (column('Z')):(column('hMean')) w l title 'enthalpy' ,\
RIFFile u 6:25 w p t 'enthalpy RIF'
# plot several species

# plot O2 
set output 'flameletLe1_O2.eps'
set xlabel 'Z'
set ylabel 'Y_{O_2}'
plot outputFile u (column('Z')):(column('O2')) w l title 'Y_{O_2}' ,\
RIFFile u 6:8 w p t 'Y_{O_2,RIF}'

# plot OH
set output 'flameletLe1_OH.eps'
set xlabel 'Z'
set ylabel 'Y_{OH}'
plot outputFile u (column('Z')):(column('OH')) w l title 'Y_{OH}' ,\
RIFFile u 6:18 w p t 'Y_{OH,RIF}'


# plot H2O
set output 'flameletLe1_H2O.eps'
set xlabel 'Z'
set ylabel 'Y_{H_2O}'
plot outputFile u (column('Z')):(column('H2O')) w l title 'Y_{H_2O}' ,\
RIFFile u 6:9 w p t 'Y_{H_2O,RIF}'


# plot CO2
set output 'flameletLe1_CO2.eps'
set xlabel 'Z'
set ylabel 'Y_{CO_2}'
plot outputFile u (column('Z')):(column('CO2')) w l title 'Y_{CO_2}' ,\
RIFFile u 6:11 w p t 'Y_{CO_2,RIF}'


# plot H2
set output 'flameletLe1_H2.eps'
set xlabel 'Z'
set ylabel 'Y_{H_2}'
plot outputFile u (column('Z')):(column('H2')) w l title 'Y_{H_2}' ,\
RIFFile u 6:21 w p t 'Y_{H_2,RIF}'


# plot CO
set output 'flameletLe1_CO.eps'
set xlabel 'Z'
set ylabel 'Y_{CO}'
plot outputFile u (column('Z')):(column('CO')) w l title 'Y_{CO}' ,\
RIFFile u 6:17 w p t 'Y_{CO,RIF}'


# plot CH4
set output 'flameletLe1_CH4.eps'
set xlabel 'Z'
set ylabel 'Y_{CH_4}'
plot outputFile u (column('Z')):(column('CH4')) w l title 'Y_{CH_4}' ,\
RIFFile u 6:7 w p t 'Y_{CH_4,RIF}'


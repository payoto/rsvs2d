#!/bin/bash 
# transforms grid files to the required .xyz format for mesh deformation
DIR="$( cd "$( dirname "${BASH_SOURCE[0]}" )" && pwd )"

cd $DIR
# -----------------------------#
# Section 1 writing volume.xyz #
# -----------------------------#
in="griduns"
line=$(head -n 1 $in | sed 's/\r$//' | sed 's/$\r//') 
numedge=$(echo -n $line | awk '{print $2}' | sed 's/\r$//' | sed 's/$\r//') 
numpoints=$(echo -n $line | awk '{print $3}' | sed 's/\r$//') 

echo Number of Edges:  $numedge 
echo Number of Points: $numpoints 


startl=$(($numedge+2)) 
endl=$(($numedge+1+$numpoints)) 
echo $numpoints > volume.xyz
awk -v awk_strt=$startl -v awk_end=$endl 'NR >= awk_strt {print $2,$3} ; NR == awk_end {exit}' $in | sed 's/\r$//' | sed 's/$/ 0.0/' >> volume.xyz 

touch symplane.xyz
numSym=$(wc -l symplane.xyz | awk '{print $1}')
echo Symmetry plane is $numSym Points
numDisp=$(head -n 1 displacements.xyz | sed 's/\r$//' | sed 's/$\r//') 
numTot=$(($numSym+$numDisp))

echo $numTot > displacements2.xyz
echo $numTot > surface2.xyz
awk -v awk_strt=2 'NR >= awk_strt {print $1,$2,$3} ' displacements.xyz | sed 's/\r$//' >> displacements2.xyz 
awk -v awk_strt=2 'NR >= awk_strt {print $1,$2,$3} ' surface.xyz | sed 's/\r$//' >> surface2.xyz 
awk '{print $1,$2,$3} ' symplane.xyz | sed 's/\r$//' >> surface2.xyz 
awk '{print 0.0,0.0,0.0} ' symplane.xyz | sed 's/\r$//' >> displacements2.xyz 
cp -p displacements.xyz displacements_orig.xyz
cp -p surface.xyz surface_orig.xyz
cp -p displacements2.xyz displacements.xyz
cp -p surface2.xyz surface.xyz

# -------------------------------------#
# Section 2 Call Mesh deformation code #
# -------------------------------------#
./meshprep.exe
./meshdef.exe volume.xyz.meshdef displacements.xyz volume2.xyz

# ---------------------------#
# Section 3 writing new grid #
# ---------------------------#
topblock=$(($startl-1))
head -n $topblock $in > griduns2

awk -v awk_end=$(($numpoints+1)) 'NR >= 2 {print NR-1,$1,$2} ; NR == awk_end {exit}' volume2.xyz >>${in}2
awk -v awk_end=$endl 'NR > awk_end ' $in  >> ${in}2

mv $in ${in}_old
mv ${in}2 ${in}

rm volume2.xyz
rm volume.xyz
rm volume.xyz.meshdef

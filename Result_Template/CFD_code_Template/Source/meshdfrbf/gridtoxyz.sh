#!/bin/bash 
# transforms grid files to the required .xyz format for mesh deformation

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

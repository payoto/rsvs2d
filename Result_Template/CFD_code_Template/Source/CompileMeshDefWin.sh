#!/bin/bash

DIR="$( cd "$( dirname "${BASH_SOURCE[0]}" )" && pwd )"

cd "$DIR/meshdfrbf"

mv makefilewin makefile

mkdir obj/
mkdir mod/
mkdir bin/

make

cp -rp bin/meshdef.exe ../meshdef.exe
cp -rp bin/meshprep.exe ../meshprep.exe
cp -rp x64/*.dll ../
chmod 755 gridtoxyz.sh
cp -rp gridtoxyz.sh ../gridtoxyz.sh

cd ..
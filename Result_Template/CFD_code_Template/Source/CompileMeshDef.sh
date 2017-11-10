#!/bin/bash

DIR="$( cd "$( dirname "${BASH_SOURCE[0]}" )" && pwd )"

cd "$DIR/meshdfrbf"

mv makefilelinux makefile

mkdir obj/
mkdir mod/
mkdir bin/

make

cp -rp bin/meshdef ../meshdef.exe
cp -rp bin/meshprep ../meshprep.exe
chmod 755 gridtoxyz.sh
cp -rp gridtoxyz.sh ../gridtoxyz.sh

cd ..
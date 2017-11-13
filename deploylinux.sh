#!/bin/bash


cd $1

unzip *.zip

git reset --hard HEAD

chmod 755 ./Result_Template/CFD_code_Template/Source/CompileAll.sh
./Result_Template/CFD_code_Template/Source/CompileAll.sh

options="-nosplash -nodesktop -noFigureWindows"
matlab $options -r PrepareLinux

cd ..


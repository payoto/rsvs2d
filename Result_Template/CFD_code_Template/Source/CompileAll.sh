#!/bin/bash

DIR="$( cd "$( dirname "${BASH_SOURCE[0]}" )" && pwd )"

cd $DIR

gfortran -O3 -o cartcell.exe CartCellv29.f90
gfortran -O3 -freal-8-real-16 -o eulerflowuns.exe EulerFlowUns14.f90
gfortran -O3 -o postproc.exe pproc11_AP.f90

cp -rp cartcell.exe ../supersonic_biplane/cartcell.exe
cp -rp eulerflowuns.exe ../supersonic_biplane/eulerflowuns.exe
cp -rp postproc.exe ../supersonic_biplane/postproc.exe
chmod 755 ../supersonic_biplane/RunFlow.sh
chmod 755 ../supersonic_biplane/RunOnlyFlow.sh
chmod 755 ../supersonic_biplane/RunPost.sh

cp -rp cartcell.exe ../supersonic_ogive/cartcell.exe
cp -rp eulerflowuns.exe ../supersonic_ogive/eulerflowuns.exe
cp -rp postproc.exe ../supersonic_ogive/postproc.exe
chmod 755 ../supersonic_ogive/RunFlow.sh
chmod 755 ../supersonic_ogive/RunOnlyFlow.sh
chmod 755 ../supersonic_ogive/RunPost.sh

cp -rp cartcell.exe ../transonic/cartcell.exe
cp -rp eulerflowuns.exe ../transonic/eulerflowuns.exe
cp -rp postproc.exe ../transonic/postproc.exe
chmod 755 ../transonic/RunFlow.sh
chmod 755 ../transonic/RunOnlyFlow.sh
chmod 755 ../transonic/RunPost.sh
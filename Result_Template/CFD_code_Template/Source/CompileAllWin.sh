#!/bin/bash

DIR="$( cd "$( dirname "${BASH_SOURCE[0]}" )" && pwd )"
echo $DIR
cd "$DIR"

gfortran -O3 -o cartcell.exe CartCellv31.f90
gfortran -O3 -o eulerflowuns.exe EulerFlowUns42_AP.f90
gfortran -O3 -o postproc.exe pproc12_AP.f90
gfortran -O3 -o showmesh.exe showmesh2_AP.f90
gcc -O3  meshsym_1.c -o meshsym.exe
gcc -O3  meshvol_1.c -o meshvol.exe
gcc -O3 "../../../MEX_Function_Directory/C_Development/Triangle2PLT/triangle2plt.c" -o triangle2plt.exe


chmod 755 CompileMeshDef.sh
./CompileMeshDefWin.sh


for foo in `find .. -maxdepth 1 -mindepth 1 -type d -not -name Source`
do
	cp -p *.exe $foo
	cp -p *.dll $foo
	cp -p gridtoxyz.sh $foo
	cp -p cmdfiles/*.sh $foo
	cp -p cmdfiles/*.bat $foo
	cp -p cmdfiles/bashtimeout $foo
done

chmod 755 ../*/bashtimeout
chmod 755 ../*/*.sh

# cp -rp cartcell.exe ../supersonic_biplane/cartcell.exe
# cp -rp eulerflowuns.exe ../supersonic_biplane/eulerflowuns.exe
# cp -rp postproc.exe ../supersonic_biplane/postproc.exe
# cp -rp meshsym.exe ../supersonic_biplane/meshsym.exe
# chmod 755 ../supersonic_biplane/RunFlow.sh
# chmod 755 ../supersonic_biplane/RunFlowSym.sh
# chmod 755 ../supersonic_biplane/RunOnlyFlow.sh
# chmod 755 ../supersonic_biplane/RunPost.sh

# cp -rp cartcell.exe ../supersonic_ogive/cartcell.exe
# cp -rp eulerflowuns.exe ../supersonic_ogive/eulerflowuns.exe
# cp -rp postproc.exe ../supersonic_ogive/postproc.exe
# cp -rp meshsym.exe ../supersonic_ogive/meshsym.exe
# chmod 755 ../supersonic_ogive/RunFlow.sh
# chmod 755 ../supersonic_ogive/RunFlowSym.sh
# chmod 755 ../supersonic_ogive/RunOnlyFlow.sh
# chmod 755 ../supersonic_ogive/RunPost.sh

# cp -rp cartcell.exe ../supersonic_ogivecoarse/cartcell.exe
# cp -rp eulerflowuns.exe ../supersonic_ogivecoarse/eulerflowuns.exe
# cp -rp postproc.exe ../supersonic_ogivecoarse/postproc.exe
# cp -rp meshsym.exe ../supersonic_ogivecoarse/meshsym.exe
# chmod 755 ../supersonic_ogivecoarse/RunFlow.sh
# chmod 755 ../supersonic_ogivecoarse/RunFlowSym.sh
# chmod 755 ../supersonic_ogivecoarse/RunOnlyFlow.sh
# chmod 755 ../supersonic_ogivecoarse/RunPost.sh

# cp -rp cartcell.exe ../supersonic_ogivefine/cartcell.exe
# cp -rp eulerflowuns.exe ../supersonic_ogivefine/eulerflowuns.exe
# cp -rp postproc.exe ../supersonic_ogivefine/postproc.exe
# cp -rp meshsym.exe ../supersonic_ogivefine/meshsym.exe
# chmod 755 ../supersonic_ogivefine/RunFlow.sh
# chmod 755 ../supersonic_ogivefine/RunFlowSym.sh
# chmod 755 ../supersonic_ogivefine/RunOnlyFlow.sh
# chmod 755 ../supersonic_ogivefine/RunPost.sh

# cp -rp cartcell.exe ../transonic/cartcell.exe
# cp -rp eulerflowuns.exe ../transonic/eulerflowuns.exe
# cp -rp postproc.exe ../transonic/postproc.exe
# cp -rp meshsym.exe ../transonic/meshsym.exe
# chmod 755 ../transonic/RunFlow.sh
# chmod 755 ../transonic/RunFlowSym.sh
# chmod 755 ../transonic/RunOnlyFlow.sh
# chmod 755 ../transonic/RunPost.sh

# cp -rp cartcell.exe ../transonicfine/cartcell.exe
# cp -rp eulerflowuns.exe ../transonicfine/eulerflowuns.exe
# cp -rp postproc.exe ../transonicfine/postproc.exe
# cp -rp meshsym.exe ../transonicfine/meshsym.exe
# chmod 755 ../transonicfine/RunFlow.sh
# chmod 755 ../transonicfine/RunFlowSym.sh
# chmod 755 ../transonicfine/RunOnlyFlow.sh
# chmod 755 ../transonicfine/RunPost.sh


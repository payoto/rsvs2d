cd %~dp0
gfortran -O3 -o cartcell.exe CartCellv31.f90
gfortran -O3 -o eulerflowuns.exe EulerFlowUns41_AP.f90
gfortran -O3 -o postproc.exe pproc12_AP.f90
gfortran -O3 -o showmesh.exe showmesh2_AP.f90
gcc -O3 -o meshsym.exe meshsym_1.c

copy "cartcell.exe" "..\supersonic_biplane\cartcell.exe"
copy "eulerflowuns.exe" "..\supersonic_biplane\eulerflowuns.exe"
copy "postproc.exe" "..\supersonic_biplane\postproc.exe"
copy "meshsym.exe" "..\supersonic_biplane\meshsym.exe"

copy "cartcell.exe" "..\supersonic_ogive\cartcell.exe"
copy eulerflowuns.exe ..\supersonic_ogive\eulerflowuns.exe
copy postproc.exe ..\supersonic_ogive\postproc.exe
copy meshsym.exe ..\supersonic_ogive\meshsym.exe

copy "cartcell.exe" "..\supersonic_ogivecoarse\cartcell.exe"
copy eulerflowuns.exe ..\supersonic_ogivecoarse\eulerflowuns.exe
copy postproc.exe ..\supersonic_ogivecoarse\postproc.exe
copy meshsym.exe ..\supersonic_ogivecoarse\meshsym.exe

copy "cartcell.exe" "..\supersonic_ogivefine\cartcell.exe"
copy eulerflowuns.exe ..\supersonic_ogivefine\eulerflowuns.exe
copy postproc.exe ..\supersonic_ogivefine\postproc.exe
copy meshsym.exe ..\supersonic_ogivefine\meshsym.exe

copy "cartcell.exe" "..\transonic\cartcell.exe"
copy eulerflowuns.exe ..\transonic\eulerflowuns.exe
copy postproc.exe ..\transonic\postproc.exe
copy meshsym.exe ..\transonic\meshsym.exe

copy "cartcell.exe" "..\transonicfine\cartcell.exe"
copy eulerflowuns.exe ..\transonicfine\eulerflowuns.exe
copy postproc.exe ..\transonicfine\postproc.exe
copy meshsym.exe ..\transonicfine\meshsym.exe
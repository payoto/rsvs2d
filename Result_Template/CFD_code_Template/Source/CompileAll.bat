cd %~dp0
gfortran -O3 -freal-8-real-16 -o cartcell.exe CartCellv29.f90
gfortran -O3 -freal-8-real-16 -o eulerflowuns.exe EulerFlowUns39.f90
gfortran -O3 -o postproc.exe pproc11_AP.f90
gfortran -O3 -o showmesh.exe showmesh1_AP.f90
gcc -O3 -o meshsym.exe meshsym_1.c

copy "cartcell.exe" "..\supersonic_biplane\cartcell.exe"
copy "eulerflowuns.exe" "..\supersonic_biplane\eulerflowuns.exe"
copy "postproc.exe" "..\supersonic_biplane\postproc.exe"
copy "meshsym.exe" "..\supersonic_biplane\meshsym.exe"

copy "cartcell.exe" "..\supersonic_ogive\cartcell.exe"
copy eulerflowuns.exe ..\supersonic_ogive\eulerflowuns.exe
copy postproc.exe ..\supersonic_ogive\postproc.exe
copy meshsym.exe ..\supersonic_ogive\meshsym.exe

copy "cartcell.exe" "..\transonic\cartcell.exe"
copy eulerflowuns.exe ..\transonic\eulerflowuns.exe
copy postproc.exe ..\transonic\postproc.exe
copy meshsym.exe ..\transonic\meshsym.exe

copy "cartcell.exe" "..\transonicfine\cartcell.exe"
copy eulerflowuns.exe ..\transonicfine\eulerflowuns.exe
copy postproc.exe ..\transonicfine\postproc.exe
copy meshsym.exe ..\transonicfine\meshsym.exe
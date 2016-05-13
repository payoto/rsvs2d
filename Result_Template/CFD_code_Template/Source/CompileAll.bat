cd %~dp0
gfortran -O3 -o cartcell.exe CartCellv29.f90
gfortran -O3 -freal-8-real-16 -o eulerflowuns.exe EulerFlowUns14.f90
gfortran -O3 -o postproc.exe pproc10_AP.f90

copy "cartcell.exe" "..\supersonic_biplane\cartcell.exe"
copy "eulerflowuns.exe" "..\supersonic_biplane\eulerflowuns.exe"
copy "postproc.exe" "..\supersonic_biplane\postproc.exe"

copy "cartcell.exe" "..\supersonic_ogive\cartcell.exe"
copy eulerflowuns.exe ..\supersonic_ogive\eulerflowuns.exe
copy postproc.exe ..\supersonic_ogive\postproc.exe
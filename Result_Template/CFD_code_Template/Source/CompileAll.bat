cd %~dp0
gfortran -O3 -o cartcell.exe CartCellv27.f90
gfortran -O3 -freal-8-real-16 -o eulerflowuns.exe EulerFlowUns14.f90
gfortran -O3 -o postproc.exe pproc8.f90
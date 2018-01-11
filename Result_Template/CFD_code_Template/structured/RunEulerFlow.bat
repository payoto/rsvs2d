cd %~dp0
touch res_hist.dat
cat res_hist.dat >> full_hist.dat
eulerflowuns.exe

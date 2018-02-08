#!/bin/bash

DIR="$( cd "$( dirname "${BASH_SOURCE[0]}" )" && pwd )"

cd $DIR
echo $$ > actpid

./eulerflowuns.exe

touch res_hist.dat
touch full_hist.dat
NFHIST=$(wc -l full_hist.dat | awk '{print $1}')
awk -v NFH=$NFHIST '/^[[:space:]]*[0-9]/  {print NR+NFH,$2,$0}' res_hist.dat >> full_hist.dat


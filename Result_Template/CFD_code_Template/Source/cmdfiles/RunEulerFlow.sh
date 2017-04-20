#!/bin/bash

DIR="$( cd "$( dirname "${BASH_SOURCE[0]}" )" && pwd )"

cd $DIR
touch res_hist.dat
cat res_hist.dat >> full_hist.dat
./eulerflowuns.exe

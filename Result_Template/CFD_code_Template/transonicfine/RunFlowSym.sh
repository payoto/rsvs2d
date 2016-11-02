#!/bin/bash

DIR="$( cd "$( dirname "${BASH_SOURCE[0]}" )" && pwd )"

cd $DIR

cartcell.exe
meshsym.exe
mv griduns griduns_full
mv griduns2 griduns
eulerflowuns.exe

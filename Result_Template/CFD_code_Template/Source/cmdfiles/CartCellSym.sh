#!/bin/bash

DIR="$( cd "$( dirname "${BASH_SOURCE[0]}" )" && pwd )"

cd $DIR

bashtimeout ./meshsym.exe
mv griduns griduns_full
mv griduns2 griduns

#!/bin/bash

DIR="$( cd "$( dirname "${BASH_SOURCE[0]}" )" && pwd )"

cd $DIR


head -n 11 settings > newsettings
echo "0 1"  >> newsettings
tail -n 3 settings >> newsettings

cp newsettings settings

./eulerflowuns.exe

#!/bin/bash

DIR="$( cd "$( dirname "${BASH_SOURCE[0]}" )" && pwd )"

cd $DIR
echo $$ > actpid

bashtimeout triangle -pq33aDenV boundtriangle
triangle2plt boundtriangle.1 griduns

mv boundtriangle.1.griduns griduns
rm boundtriangle.1.*
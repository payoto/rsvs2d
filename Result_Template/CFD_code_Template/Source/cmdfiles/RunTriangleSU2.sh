#!/bin/bash

DIR="$( cd "$( dirname "${BASH_SOURCE[0]}" )" && pwd )"

cd "$DIR"
echo $$ > actpid
bashtimeout triangle ${1:--pq33aDenV} boundtriangle
triangle2plt boundtriangle.1 su2

mv boundtriangle.1.su2 triangularmesh.su2
rm boundtriangle.1.*
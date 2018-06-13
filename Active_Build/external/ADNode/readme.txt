Potential Flow Panel Method with Reverse Mode Automatic Differentiation
The potential flow program has been modified to allow use with an
object-oriented implementation of reverse mode AD using operator overloading.

Example usage is given in inv_pressureAD.m and example.m which evaluate
an inverse pressure design objective function and compares gradients from
AD with those from finite differences.

The AD implementation is provided through ADNode.m which is available 
from: https://github.com/gaika/madiff and distributed under the 
GNU General Public License; this license permits the distribution of both
the original and modified versions of the code.
See ADNode_LICENSE.txt for the full license terms.
See ADNode.m for the details of modifications made.

Summary of root files:
 - ADNode.m              Automatic differentiation object class
 - ADNode_LICENSE.txt    Original license for ADNode.m (GNU GPL)
 - aerofoils.mat         Matlab variables with aerofoil coordinates
 - example.m             Example script to compare with finite differences
 - inv_pressureAD.m      Inverse pressure design objective function with AD
 - inv_pressureFD.m                                                 with finite differences
 - potflowsolver.m       Modified potential flow panel method



Laurence Kedward, March 2017
# README #
Last updated 13/11/2017


For this code to work necessary programs:  
+	MATLAB installed (2015a or later) - including parallel toolbox  
+	c compiler compatible with MATLAB for the compilation of mex files  
+	Standalone c compiler for the compilation of console programs  
+	fortran (90+) compiler for compilation of flow solvers   
+	Cygwin needed for compilation of some flow solver features(Mesh motion) on windows with the "make" utility  

	
## What is this repository for? ##

This repository is the MATLAB implementation of the 2D R-Snake Volume of Solid (RSVS) parameterisation.
It includes an optimisation framework designed for to explore generic optimisation problems using this parameterisation.
A number of objective functions are included

## How do I get set up? ##

To get started:
(Linux - First Time) 
	Run: `./Result_Template/Source/CompileAll.sh`
		This compiles all the separate utilities (mesher, flow solver, mesh deformation)
		requires fortran and c compilers
	Launch Matlab from the source directory
	in matlab:
		>> InitialiseSnakeFlow
		>> Compile_Mex
		>> ExecInclude

		
(Linux) If the distro contains all the necessary program and utilities this is designed to deply the code from a zipped git repo

	mv deploylinux.sh ../
	cd ..
	deploylinux.sh <source directory name>
		
(Windows) 
	Run: `./Result_Template/Source/CompileAllWin.sh`
	
This compiles all the separate utilities (mesher, flow solver, mesh deformation)
	requires fortran and c compilers. Requires cygwin with the make tool.
	Launch Matlab from the source directory
	in matlab:
	
		>> InitialiseSnakeFlow
		>> Compile_Mex
		>> ExecInclude

(Any - At the start of a new MATLAB session) 
	in matlab:
	
		>> InitialiseSnakeFlow
		>> ExecInclude
		
Tests to run:
	Main('SnakesFoilVVSmall4')
	% Tests the Snaking process on a multi-topology case
	ExecuteOptimisation('Test_Rosenbrock')
	% Tests the optimisation framework on an analytical function
	
These tests will save results in:

	../results/Standard_Execution/<Archive_YYYY_MM>/Day_YY_MM_DD/
	
and
	
	../results/Optimisation/<Archive_YYYY_MM>/Day_YY_MM_DD/
	
respectively. All result files and folders are time stamped such that a named sort will return the files in chronological order of creation.
These files and their location are entered into a file  at:
	../results/<archive_name>/Index_<archive_name>.txt

## Outline of the code ##

The code is centred around an optimisation framework executed through the function `ExecuteOptimisation`
This framework was developped to exploit the RSVS which is implemented in function `Snakes` and executed through function `Main`


ExecuteOptimisation can accept a minimum of 1 input and up to 3.
	In normal execution only the first input is required. It must be the name of a function stored as a string.
	This function must return a valid parameter structure in the form of those generated by StructOptimParam.
	In practice new parameter configurations are added as local functions to StructOptimParam.m
	A example is below:
	
		function [paroptim]=<function_name_callable_string>()
			% With parameter structure with default settings
			[paroptim]=DefaultOptim(); 
			% Apply a Standard modification to default settings: for example change all the 
			% settings required to change the optimiser for aero, differential evolution cases:
			paroptim=CutCellObjective(paroptim);
			paroptim=OptimDE(paroptim);
			
			% Single changes to specific parameters:
			paroptim.general.maxIter=5; % maximum number of optimisation iteration
			paroptim.general.worker=2; % number of parallel worker
		end


Input 2 of ExecuteOptimisation is a cell array used to restart an optimisation. The form of this is:
{'path_to_restart_data.mat',{'optim_type',isRestart}}  
'path_to_restart_data.mat' is the path to the matlab data where the restart information is stored. The strict
		minimum  is for this to contain the result structure from the optimisation process called `optimstruct`
		(Saved in OptimRes_<...>.mat in the result directory)  
'optim_type' is the type of optimisation that was used as a string usually either 'conjgrad' or 'DE'  
isRestart stores whether a gradient based optimisation should be restarted from a Line search (0) or a gradient calculation (1)
	in general mod(numel(optimstruct)+1,2) will give you the right answer for gradient based optimisation. Always 1 for DE.  


## Contribution guidelines ##

To add new features a few things must be considered (roughly in order of importance):
+ Maintaining the git archive's integrity and minimising merge conflicts for yourself and other users.
+ Execution by other users.
+ Simplicity of deployment.
+ Maintaining consistent coding style.

### Managing git ###

To manage the git repo's integrity new features should be developed in new folders contained in `./Active_Build/`
Make sure that names of functions are logical and express what they do. If large numbers of files are needed, make sure each contain information
about the creator, the feature they are part of and the intended purpose.

Use a separate git branch for development of new experimental features that may break the main code (and not just the new feature).

Large new features should add their own parameter default to `StructOptimParam.m` 

### Managing git ###

## What does this ACTUALLY do and who do I talk to?##

For more information about what the code does (i.e. the science of it)  
[Restricted Snakes: a Flexible Topology Parameterisation Method for Aerodynamic Optimisation](https://arc.aiaa.org/doi/pdf/10.2514/6.2017-1410)  
[Mixing and Refinement of Design Variables for Geometry and Topology Optimization in Aerodynamics](https://arc.aiaa.org/doi/pdfplus/10.2514/6.2017-3577)  

(Also available on research gate)

Alexandre Payot - a.payot@bristol.ac.uk  
[ResearchGate profile](https://www.researchgate.net/profile/Alexandre_Payot)  
[Google Scholar profile](https://scholar.google.co.uk/citations?user=JX_AmkwAAAAJ&hl=en)  
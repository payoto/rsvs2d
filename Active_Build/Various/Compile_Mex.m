%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%              University of Bristol
%      Department of Aerospace Engineering
%                     2015
%
%          Subdivision of Surfaces
%      for Aerodynamic shape parametrisation
%          - Compile Mex code -
%             Alexandre Payot
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Compile Grid Initialisation
mex  .\MEX_Function_Directory\MEX_Sources\gridgen\GridInit_MEX.c...
    -outdir .\MEX_Function_Directory\MEX_Executables\gridgen 
%%
%mex .\MEX_Function_Directory\MEX_Sources\gridgen\GridInit_MEX_OLD.c...
%    -outdir .\MEX_Function_Directory\MEX_Executables\gridgen
%% Compile grid refinement

mex  .\MEX_Function_Directory\MEX_Sources\gridgen\GridRefine_MEX.c...
    -outdir .\MEX_Function_Directory\MEX_Executables\gridgen 

%% Compile FindObjNum_MEx

mex .\MEX_Function_Directory\MEX_Sources\findobjnum\FindObjNum_MEX.c...
    -outdir .\MEX_Function_Directory\MEX_Executables\findobjnum
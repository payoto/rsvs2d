% Parameter File 
% 13-May-2016 23:01:40 
% 2016-05-13T230140_Test_CG_Aero 

param.general.optimCase = ['Test_CG_Aero']; 
param.general.paramCase = ['optimSupersonic']; 
param.general.optimMethod = ['conjgradls']; 
param.general.desVarRange = [0,1]; 
param.general.nPop = [6]; 
param.general.startPop = ['randuniformsharp']; 
param.general.maxIter = [50]; 
param.general.worker = [8]; 
param.general.objectiveName = ['CutCellFlow']; 
param.general.direction = ['min']; 
param.general.defaultVal = [1000]; 
param.general.knownOptim = [0.14609]; 
param.general.restartSource = ''; 
param.general.symType = ['horz']; 
param.general.nDesVar = [0]; 
param.general.symDesVarList = [ ]; 
param.general.notDesInd = [ ]; 
param.general.iterGap = [1]; 
param.optim.DE.diffAmplification = [0.5]; 
param.optim.DE.xOverRatio = [0.5]; 
param.optim.DE.geneType = ['single']; 
param.optim.CG.diffStepSize = [0.001]; 
param.optim.CG.varOverflow = ['truncate']; 
param.optim.CG.varActive = ['all']; 
param.optim.CG.lineSearch = [false]; 
param.optim.CG.validVol = [0.3]; 
param.spline.splineCase = ['aerosnake']; 
param.spline.domain = ['normalizeX']; 
param.obj.flow.CFDfolder = ['C:\Users\ap1949\Local Documents\PhD\Development Work\Snakes Volume Parametrisation\source\Result_Template\CFD_code_Template\supersonic_ogive']; 
param.obj.flow.stoponerror = [false]; 
param.obj.flow.targConv = [-6]; 
param.obj.flow.lengthConvTest = [100]; 
param.obj.flow.restartIter = [1000]; 
param.obj.flow.maxRestart = [5]; 
param.obj.flow.nMach = [2]; 
param.constraint.desVarConstr = {['MeanVolFrac']}; 
param.constraint.desVarVal = {[0.4]}; 
param.constraint.resConstr = {['AeroResidualBarrier']}; 
param.constraint.resVal = {[-3.5,-0.5]}; 

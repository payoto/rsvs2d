% Parameter File 
% 09-Jun-2016 16:54:39 
% 2016-06-09T165439_Test_MultiTopo_M2_CG 

param.general.optimCase = ['Test_MultiTopo_M2_CG']; 
param.general.paramCase = ['optimSupersonicMultiTopo']; 
param.general.optimMethod = ['conjgrad']; 
param.general.desVarRange = [0,1]; 
param.general.nPop = [60]; 
param.general.startPop = ['halfuniformsharp']; 
param.general.maxIter = [10]; 
param.general.worker = [4]; 
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
param.general.desvarconnec = [ ]; 
param.optim.DE.diffAmplification = [0.5]; 
param.optim.DE.xOverRatio = [0.5]; 
param.optim.DE.geneType = ['single']; 
param.optim.CG.diffStepSize = [0.001,-0.001]; 
param.optim.CG.varOverflow = ['truncate']; 
param.optim.CG.varActive = ['border']; 
param.optim.CG.lineSearch = [false]; 
param.optim.CG.validVol = [0.1]; 
param.optim.CG.nLineSearch = [12]; 
param.spline.splineCase = ['aerosnake']; 
param.spline.domain = ['normalizeX']; 
param.obj.flow.CFDfolder = ['C:\Users\ap1949\Local Documents\PhD\Development Work\Snakes Volume Parametrisation\source\Result_Template\CFD_code_Template\supersonic_ogive']; 
param.obj.flow.stoponerror = [false]; 
param.obj.flow.targConv = [-6]; 
param.obj.flow.lengthConvTest = [100]; 
param.obj.flow.restartIter = [1000]; 
param.obj.flow.maxRestart = [5]; 
param.obj.flow.nMach = [2]; 
param.constraint.desVarConstr = {['MinSumVolFrac']}; 
param.constraint.desVarVal = {[3.8]}; 
param.constraint.resConstr = {['AeroResidualBarrier']}; 
param.constraint.resVal = {[-0.5,0.5]}; 

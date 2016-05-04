% Parameter File 
% 28-Apr-2016 19:34:13 
% 2016-04-28T193413_Test_MultiTopo_M2_D 

param.general.optimCase = ['Test_MultiTopo_M2_D']; 
param.general.paramCase = ['optimSupersonicMultiTopo']; 
param.general.optimMethod = ['DEtan']; 
param.general.desVarRange = [0,1]; 
param.general.nDesVar = [0]; 
param.general.nPop = [48]; 
param.general.startPop = ['initbusemann']; 
param.general.maxIter = [20]; 
param.general.worker = [8]; 
param.general.objectiveName = ['CutCellFlow']; 
param.general.direction = ['min']; 
param.general.defaultVal = [1000]; 
param.general.knownOptim = [0.14609]; 
param.general.restartSource = ''; 
param.general.symType = ['horz']; 
param.general.symDesVarList = [ ]; 
param.optim.DE.diffAmplification = [0.5]; 
param.optim.DE.xOverRatio = [0.5]; 
param.optim.DE.geneType = ['horz']; 
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
param.constraint.resConstr = {['AeroResidual']}; 
param.constraint.resVal = {[-3.5]}; 

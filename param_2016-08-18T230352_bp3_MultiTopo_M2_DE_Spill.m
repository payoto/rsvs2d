% Parameter File 
% 18-Aug-2016 23:03:52 
% 2016-08-18T230352_bp3_MultiTopo_M2_DE_Spill 
param.general.optimCase = ['bp3_MultiTopo_M2_DE_Spill']; 
param.general.paramCase = ['optimSupersonicMultiTopo']; 
param.general.optimMethod = ['DE']; 
param.general.desVarRange = [0,1]; 
param.general.nPop = [100]; 
param.general.startPop = ['initbusemann']; 
param.general.specificFillName = ['24DVaverage']; 
param.general.maxIter = [100]; 
param.general.worker = [12]; 
param.general.objectiveName = ['CutCellFlow']; 
param.general.direction = ['min']; 
param.general.defaultVal = [1000]; 
param.general.knownOptim = [0.146088675]; 
param.general.symType = ['horz']; 
param.general.nDesVar = [0]; 
param.general.symDesVarList = [ ]; 
param.general.notDesInd = [ ]; 
param.general.iterGap = [1]; 
param.general.desvarconnec = [ ]; 
param.general.restartSource = {'',''}; 
param.general.isRestart = [false]; 
param.general.varOverflow = ['spill']; 
param.optim.DE.diffAmplification = [0.5]; 
param.optim.DE.xOverRatio = [0.5]; 
param.optim.DE.geneType = ['horz']; 
param.optim.CG.diffStepSize = [0.01,-0.01]; 
param.optim.CG.varActive = ['all']; 
param.optim.CG.borderActivation = [0.15]; 
param.optim.CG.lineSearch = [false]; 
param.optim.CG.validVol = [0.5]; 
param.optim.CG.openVol = [0.1]; 
param.optim.CG.nLineSearch = [12]; 
param.spline.splineCase = ['aerosnake']; 
param.spline.domain = ['normalizeX']; 
param.obj.flow.CFDfolder = ['/panfs/panasas01/aero/ap1949/SnakVolParam/source\Result_Template\CFD_code_Template\supersonic_ogive']; 
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
param.constraint.initConstr = {[' ']}; 
param.constraint.initVal = {[' ']}; 

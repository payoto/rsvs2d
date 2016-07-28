% Parameter File 
% 27-Jul-2016 23:20:56 
% 2016-07-27T232056_desk_Aero_CG_20L_po 
param.general.optimCase = ['desk_Aero_CG_20L_po']; 
param.general.paramCase = ['optimSupersonic_Long']; 
param.general.optimMethod = ['conjgrad']; 
param.general.desVarRange = [0,1]; 
param.general.nPop = [12]; 
param.general.startPop = ['halfuniformsharp']; 
param.general.maxIter = [10]; 
param.general.worker = [4]; 
param.general.objectiveName = ['CutCellFlow']; 
param.general.direction = ['min']; 
param.general.defaultVal = [1000]; 
param.general.knownOptim = [0.146088675]; 
param.general.symType = ['horz']; 
param.general.nDesVar = [0]; 
param.general.symDesVarList = [ ]; 
param.general.notDesInd = [ ]; 
param.general.iterGap = [2]; 
param.general.desvarconnec = [ ]; 
param.general.restartSource = {'',''}; 
param.general.isRestart = [false]; 
param.general.varOverflow = ['spill']; 
param.optim.DE.diffAmplification = [0.5]; 
param.optim.DE.xOverRatio = [0.5]; 
param.optim.DE.geneType = ['single']; 
param.optim.CG.diffStepSize = [0.001,-0.001]; 
param.optim.CG.varActive = ['snaksensiv']; 
param.optim.CG.borderActivation = [0.15]; 
param.optim.CG.lineSearch = [false]; 
param.optim.CG.validVol = [0.3]; 
param.optim.CG.openVol = [0.1]; 
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
param.constraint.desVarConstr = {['MeanVolFrac']}; 
param.constraint.desVarVal = {[0.3]}; 
param.constraint.resConstr = {['AeroResidualBarrier']}; 
param.constraint.resVal = {[-0.5,0.5]}; 
param.constraint.initConstr = {[' ']}; 
param.constraint.initVal = {[' ']}; 

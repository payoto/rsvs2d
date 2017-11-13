% Parameter File 
% 23-Jun-2016 17:19:27 
% 2016-06-23T171927_Full_Aero_CG_component 

param.general.optimCase = ['Full_Aero_CG_component']; 
param.general.paramCase = ['SupersonicComponent']; 
param.general.optimMethod = ['conjgrad']; 
param.general.desVarRange = [0,1]; 
param.general.nPop = [12]; 
param.general.startPop = ['innerbound']; 
param.general.maxIter = [20]; 
param.general.worker = [4]; 
param.general.objectiveName = ['CutCellFlow']; 
param.general.direction = ['min']; 
param.general.defaultVal = [1000]; 
param.general.knownOptim = [0.14609]; 
param.general.symType = ['none']; 
param.general.nDesVar = [0]; 
param.general.symDesVarList = [ ]; 
param.general.notDesInd = [ ]; 
param.general.iterGap = [2]; 
param.general.desvarconnec = [ ]; 
param.general.restartSource = {'',''}; 
param.general.isRestart = [false]; 
param.optim.DE.diffAmplification = [0.5]; 
param.optim.DE.xOverRatio = [0.5]; 
param.optim.DE.geneType = ['single']; 
param.optim.CG.diffStepSize = [0.01,-0.01]; 
param.optim.CG.varOverflow = ['truncate']; 
param.optim.CG.varActive = ['wideborder']; 
param.optim.CG.borderActivation = [0.15]; 
param.optim.CG.lineSearch = [false]; 
param.optim.CG.validVol = [0.5]; 
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
param.constraint.desVarConstr = {[' ']}; 
param.constraint.desVarVal = {[ ]}; 
param.constraint.resConstr = {['AeroResidualBarrier']}; 
param.constraint.resVal = {[-0.5,0.5]}; 
param.constraint.initConstr = {['LocalVolFrac_image']}; 
param.constraint.initVal = {{['.\Active_Build\ConstraintFiles\smile_5b12.png'],['min']}}; 
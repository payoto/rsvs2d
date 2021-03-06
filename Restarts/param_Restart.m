% Parameter File 
% 14-Apr-2016 11:53:38 
% Restart 

param.general.passDomBounds = [-1,1;-0.5,0.5]; 
param.general.passGridSteps = [3]; 
param.general.refineSteps = [4]; 
param.general.passPadding = [1]; 
param.general.typDat = ['optimRand']; 
param.general.typeBound = ['snaxel']; 
param.general.subdivType = ['area']; 
param.general.loadLogical = [false]; 
param.general.useSnakes = [true]; 
param.general.execTest = [false]; 
param.general.boundstr = {['boundaryis0'],['solidnotIn0'],['0bound']}; 
param.general.restart = [true]; 
param.general.case = ['optimTest']; 
param.results.archiveName = ['Standard_Execution']; 
param.results.resultRoot = ['C:\Users\ap1949\Local Documents\PhD\Development Work\Snakes Volume Parametrisation\source\..\results\']; 
param.results.noteFiles = {['CurrentBuild']}; 
param.results.tags = {['snakes'],['optimisation']}; 
param.plotting.isCheckRes = [true]; 
param.plotting.plotInterval = [0]; 
param.plotting.makeMov = [false]; 
param.plotting.debugPlot = [0]; 
param.snakes.step.snakesSteps = [24]; 
param.snakes.step.mergeTopo = [false]; 
param.snakes.step.maxStep = [0.4]; 
param.snakes.step.maxDt = [0.25]; 
param.snakes.step.convLevel = [1e-08]; 
param.snakes.step.arrivalTolerance = [0.0001]; 
param.snakes.step.subStep = [1]; 
param.snakes.step.snakesMinSteps = [5]; 
param.snakes.step.snakData = ['all']; 
param.snakes.step.snakesConsole = [true]; 
param.snakes.step.stepType = ['indiv']; 
param.snakes.step.vSwitch = [1e-15]; 
param.snakes.step.dtRatio = [5]; 
param.snakes.step.convCheckRate = [5]; 
param.snakes.step.convCheckRange = [15]; 
param.snakes.step.convDistance = [200]; 
param.snakes.refine.refineGrid = [4]; 
param.snakes.refine.typeRefine = ['all']; 
param.snakes.force.maxForceVel = [1]; 
param.snakes.force.bendingVelInfluence = [0]; 
param.snakes.force.tensVelInfluence = [1]; 
param.snakes.force.maxVelRatio = [1]; 
param.snakes.force.dampBase = [1]; 
param.snakes.force.dampSides = [0]; 
param.snakes.force.vectorMagAveraging = [true]; 
param.snakes.force.lengthEpsilon = [1e-05]; 
param.snakes.force.typeSmear = ['length']; 
param.snakes.force.isLast = [false]; 
param.snakes.force.velType = ['default']; 
param.snakes.force.vel.Type = {['default']}; 
param.snakes.force.vel.ChangeStep = [0]; 
param.snakes.force.vel.ChangeConv = [10]; 
param.snakes.force.vel.ChangeTrigger = ['none']; 
param.optiminit.cellLevels = [8,2]; 
param.optiminit.refineCellLvl = [0]; 
param.optiminit.defaultfill = [0.5]; 
param.optiminit.defaultCorner = [0.001]; 
param.optiminit.corneractive = [false]; 

% Parameter File 
% 23-Sep-2016 02:10:59 
% 2016-09-23T021059_desk_0012cos_parametrisation 
param.general.passDomBounds = [-1.2,1.2;-0.2,0.2]; 
param.general.passGridSteps = [3]; 
param.general.refineSteps = [2]; 
param.general.passPadding = [1]; 
param.general.typDat = ['optimInit']; 
param.general.typeBound = ['snaxel']; 
param.general.subdivType = ['chaikin']; 
param.general.loadLogical = [false]; 
param.general.useSnakes = [true]; 
param.general.execTest = [false]; 
param.general.boundstr = {['boundaryis0'],['solidnotIn0'],['0bound']}; 
param.general.restart = [true]; 
param.general.case = ['optimInverseDesign_Lcos']; 
param.results.archiveName = ['Optimisation']; 
param.results.resultRoot = ['C:\Users\ap1949\Local Documents\PhD\Development Work\Snakes Volume Parametrisation\source\..\results\']; 
param.results.noteFiles = {['CurrentBuild']}; 
param.results.tags = {['snakes'],['optimisation']}; 
param.plotting.isCheckRes = [true]; 
param.plotting.plotInterval = [0]; 
param.plotting.makeMov = [false]; 
param.plotting.debugPlot = [0]; 
param.plotting.checkSensitivities = [false]; 
param.snakes.step.snakesSteps = [200]; 
param.snakes.step.mergeTopo = [true]; 
param.snakes.step.maxStep = [0.4]; 
param.snakes.step.maxDt = [0.5]; 
param.snakes.step.convLevel = [1e-08]; 
param.snakes.step.arrivalTolerance = [0.1]; 
param.snakes.step.subStep = [1]; 
param.snakes.step.snakesMinSteps = [5]; 
param.snakes.step.snakData = ['light']; 
param.snakes.step.snakesConsole = [false]; 
param.snakes.step.stepType = ['indiv']; 
param.snakes.step.vSwitch = [1e-15]; 
param.snakes.step.dtRatio = [5]; 
param.snakes.step.snaxInitPos = [1e-05]; 
param.snakes.step.convCheckRate = [100]; 
param.snakes.step.convCheckRange = [15]; 
param.snakes.step.convDistance = [500]; 
param.snakes.step.fillLooseStep = [5]; 
param.snakes.step.fillLooseCut = [0.001]; 
param.snakes.refine.refineGrid = [5,1]; 
param.snakes.refine.gridDistrib = ['cosX1']; 
param.snakes.refine.typeRefine = ['all']; 
param.snakes.refine.LEShrink = [false]; 
param.snakes.refine.TEShrink = [true]; 
param.snakes.refine.edgeFinish = ['sharpen']; 
param.snakes.refine.resampleSnak = [true]; 
param.snakes.refine.axisRatio = [1.5]; 
param.snakes.refine.typeCorner = ['global']; 
param.snakes.force.maxForceVel = [1]; 
param.snakes.force.bendingVelInfluence = [0]; 
param.snakes.force.tensVelInfluence = [1]; 
param.snakes.force.maxVelRatio = [1]; 
param.snakes.force.dampBase = [1]; 
param.snakes.force.dampSides = [0]; 
param.snakes.force.vectorMagAveraging = [true]; 
param.snakes.force.lengthEpsilon = [1e-06]; 
param.snakes.force.typeSmear = ['length']; 
param.snakes.force.isLast = [false]; 
param.snakes.force.velType = ['default']; 
param.snakes.force.vel.Type = {['default']}; 
param.snakes.force.vel.ChangeStep = [0]; 
param.snakes.force.vel.ChangeConv = [10]; 
param.snakes.force.vel.ChangeTrigger = ['none']; 
param.optiminit.cellLevels = [22,2]; 
param.optiminit.refineCellLvl = [0]; 
param.optiminit.defaultfill = [0.5]; 
param.optiminit.defaultCorner = [1e-05]; 
param.optiminit.corneractive = [false]; 
param.optiminit.modeSmoothType = ['peaksmooth']; 
param.optiminit.modeSmoothNum = [6]; 

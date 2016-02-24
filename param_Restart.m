% Parameter File 
% 24-Feb-2016 12:38:03 
% Restart 

param.general.passDomBounds = [-1,1;-1,1]; 
param.general.passGridSteps = [3]; 
param.general.refineSteps = [1]; 
param.general.passPadding = [1]; 
param.general.typDat = ['testsmooth3_1']; 
param.general.typeBound = ['snaxel']; 
param.general.loadLogical = [false]; 
param.general.useSnakes = [true]; 
param.general.execTest = [false]; 
param.general.boundstr = {['boundaryis0'],['solidnotIn0'],['0bound']}; 
param.general.restart = [true]; 
param.general.case = ['Snakestestsmooth3_1']; 
param.results.archiveName = ['Standard_Execution']; 
param.results.resultRoot = ['C:\Users\ap1949\Local Documents\PhD\Development Work\Snakes Volume Parametrisation\source\..\results\']; 
param.results.noteFiles = {['CurrentBuild'],['OptimSQP']}; 
param.results.tags = {['snakes'],['Opimisation'],['SQP'],['Profile Length']}; 
param.plotting.isCheckRes = [true]; 
param.plotting.plotInterval = [0]; 
param.plotting.makeMov = [false]; 
param.plotting.debugPlot = [0]; 
param.snakes.step.snakesSteps = [200]; 
param.snakes.step.mergeTopo = [true]; 
param.snakes.step.maxStep = [0.9]; 
param.snakes.step.maxDt = [1]; 
param.snakes.step.convLevel = [1e-08]; 
param.snakes.step.arrivalTolerance = [0.01]; 
param.snakes.step.subStep = [1]; 
param.snakes.step.restartFlag = [false]; 
param.snakes.step.snakesMinSteps = [5]; 
param.snakes.refine.refineGrid = [8]; 
param.snakes.refine.typeRefine = ['all']; 
param.snakes.force.maxForceVel = [1]; 
param.snakes.force.bendingVelInfluence = [0]; 
param.snakes.force.tensVelInfluence = [1]; 
param.snakes.force.maxVelRatio = [1]; 
param.snakes.force.dampBase = [1]; 
param.snakes.force.dampSides = [0]; 
param.snakes.force.vectorMagAveraging = [true]; 
param.snakes.force.lengthEpsilon = [1e-05]; 
param.snakes.force.velType = ['default']; 

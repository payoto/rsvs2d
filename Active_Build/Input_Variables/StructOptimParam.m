%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%              University of Bristol
%      Department of Aerospace Engineering
%                     2015
%
%       Optimisation Using
%          Parametric Snakes for
%      for Aerodynamic shape
%         parametrisation
%             Alexandre Payot
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


function [paroptim]=StructOptimParam(caseStr)
    % Main function that allows changes
    
    [paroptim]=eval(caseStr);
    paroptim.general.optimCase=regexprep(regexprep(...
        regexprep(caseStr,'(\(|\)|,)','_'),'\.','_'),'''','');
    
end



%% Deafult Optimisation Inputs

function [paroptim]=DefaultOptim()
    
    paroptim.general=DefaultOptimGeneral();
    [paroptim.optim.DE]=DefaultOptimDE();
    [paroptim.optim.CG]=DefaultOptimCG();
    [paroptim.spline]=DefaultOptimSpline();
    paroptim.obj.flow=DefaultCutCell_Flow();
    [paroptim.obj.invdes]=DefaultInversedesign();
    [paroptim.constraint]=DefaultConstraint();
    
    paroptim.optim.supportOptim=[];
    paroptim.structdat=GetStructureData(paroptim);
    varExtract={'paramCase'};
    [paramCase]=ExtractVariables(varExtract,paroptim);
    paroptim.parametrisation=structInputVar(paramCase);
    
    paroptim.initparam=DefaultSnakeInit(paroptim.parametrisation);
end

function [paroptimgeneral]=DefaultOptimGeneral()
    
    paroptimgeneral.optimCase='optimDefault';
    paroptimgeneral.paramCase='optimDefault';
    paroptimgeneral.optimMethod='DE';
    paroptimgeneral.desVarRange=[0,1];
    paroptimgeneral.nPop=6;
    paroptimgeneral.startPop='rand';
    paroptimgeneral.specificFillName='24DVaverage';
    paroptimgeneral.maxIter=5;
    paroptimgeneral.worker=6; % Max 4 on this computer
    paroptimgeneral.objectiveName='LengthArea'; % 'InverseDesign' 'CutCellFlow'
    paroptimgeneral.direction='max';
    paroptimgeneral.defaultVal=-1e3;
    paroptimgeneral.knownOptim=[0.146088675]; %#ok<*NBRAK>
    paroptimgeneral.useSnake=true;
    
    paroptimgeneral.symType='none'; % 'horz'
    paroptimgeneral.nDesVar=[0];
    paroptimgeneral.symDesVarList=[];
    paroptimgeneral.notDesInd=[];
    paroptimgeneral.iterGap=1;
    paroptimgeneral.desvarconnec=[]; % Structure assigned later
    
    paroptimgeneral.refineOptim=[];
    
    paroptimgeneral.restartSource={'',''};
    paroptimgeneral.isRestart=false;
    paroptimgeneral.varOverflow='spill'; % 'truncate' 'spill'
    paroptimgeneral.initInterp={};
end

function [paroptimDE]=DefaultOptimDE()
    
    paroptimDE.diffAmplification=0.5; %[0,2]
    paroptimDE.xOverRatio=0.5;
    paroptimDE.geneType='single'; % 'horz' 'vert'
end

function [paroptimoptimCG]=DefaultOptimCG()
    
    paroptimoptimCG.diffStepSize=[1e-3,-1e-3]; %[0,2]
    paroptimoptimCG.minDiffStep=1e-6;
    paroptimoptimCG.maxDiffStep=0;
    paroptimoptimCG.varActive='all'; % 'all' 'border' 'wideborder' 'snaksensiv'
    paroptimoptimCG.sensCalc='snake'; % 'analytical'
    paroptimoptimCG.stepAlgo='conjgrad'; %'BFGS'
    paroptimoptimCG.sensAnalyticalType='raw'; % 'raw' 'smooth'
    paroptimoptimCG.borderActivation=0.15;
    paroptimoptimCG.lineSearch=false;
    paroptimoptimCG.lineSearchType='backbisection';
    paroptimoptimCG.validVol=0; % Interval of validity of the derivatives
    paroptimoptimCG.minVol=1e-4;
    paroptimoptimCG.startVol=0.5;
    paroptimoptimCG.openVol=0.1;
    paroptimoptimCG.nLineSearch=12;
    paroptimoptimCG.wolfeC1=1e-4;
    paroptimoptimCG.wolfeC2=0.1;
    paroptimoptimCG.gradScaleType='none';
    paroptimoptimCG.gradScale=[];
end

function paroptimobjflow=DefaultCutCell_Flow()
    
    paroptimobjflow.CFDfolder=[cd,'\Result_Template\CFD_code_Template\supersonic_ogive'];
    paroptimobjflow.stoponerror=true;
    paroptimobjflow.targConv=-6;
    paroptimobjflow.lengthConvTest=100;
    paroptimobjflow.restartIter=1000;
    paroptimobjflow.startIterFlow=4000;
    paroptimobjflow.maxRestart=5;
    paroptimobjflow.nMach=2;
    paroptimobjflow.isSymFlow=false;
    
end

function [paroptimspline]=DefaultOptimSpline()
    
    
    paroptimspline.splineCase='aerosnake';
    paroptimspline.resampleSnak=false;
    
end

function [paroptimobjinvdes]=DefaultInversedesign()
    
    
    paroptimobjinvdes.aeroClass='NACA'; % 'NACA' 'UIUC'
    paroptimobjinvdes.aeroName='4412';
    paroptimobjinvdes.profileComp='distance'; % 'distance' or 'area'
    
    
end

function [paroptimconstraint]=DefaultConstraint()
    
    paroptimconstraint.desVarConstr={' '};
    paroptimconstraint.desVarVal={[0]};
    paroptimconstraint.resConstr={' '};
    paroptimconstraint.resVal={[0]};
    paroptimconstraint.initConstr={' '};
    paroptimconstraint.initVal={' '};
    
end

function [paraminit]=DefaultSnakeInit(paraminit)
    
    paraminit.snakes.step.snakesSteps=20;
    paraminit.general.restart=false;
end

function paroptim=ModifySnakesParam(paroptim,paramCase)
    
    paroptim.general.paramCase=paramCase;
    paroptim.parametrisation=structInputVar(paramCase);
    paroptim.initparam=DefaultSnakeInit(paroptim.parametrisation);
    
end

%% Standard Modifications

% Optimisation Method

function [paroptim]=OptimDE(paroptim)
    paroptim.general.optimMethod='DEtan';
    paroptim.optim.DE.diffAmplification=0.5; %[0,2]
    paroptim.optim.DE.xOverRatio=0.5;
    paroptim.optim.DE.geneType='single'; % 'single' 'horz' 'vert'
    
end

function [paroptim]=OptimDE_weak(paroptim)
    paroptim.general.optimMethod='DEtan';
    paroptim.optim.DE.diffAmplification=1; %[0,2]
    paroptim.optim.DE.xOverRatio=0.2;
    paroptim.optim.DE.geneType='single'; % 'single' 'horz' 'vert'
    
end

function [paroptim]=OptimDE_horiz(paroptim)
    paroptim.general.optimMethod='DEtan';
    paroptim.optim.DE.diffAmplification=0.5; %[0,2]
    paroptim.optim.DE.xOverRatio=0.5;
    paroptim.optim.DE.geneType='horz'; % 'single' 'horz' 'vert'
end

function [paroptim]=OptimCG(paroptim)
    
    paroptim.general.optimMethod='conjgrad';
    paroptim.general.stepAlgo='conjgrad';
    paroptim.general.varOverflow='spill';
    paroptim.general.iterGap=2;
end

function [paroptim]=OptimBFGS(paroptim)
    
    paroptim.general.optimMethod='conjgrad';
    paroptim.optim.CG.stepAlgo='BFGS';
    paroptim.general.varOverflow='spill';
    paroptim.general.iterGap=2;
end

% Constraints
function [paroptim]=MeanVolumeConstraint(paroptim)
    
    paroptim.constraint.desVarConstr={'MeanVolFrac'};
    paroptim.constraint.desVarVal={0.4};
    paroptim.constraint.resConstr={'AeroResidualBarrier'};
    paroptim.constraint.resVal={[-0.5,0.5]};
    
end

function [paroptim]=ValVolumeConstraint(paroptim)
    
    paroptim.constraint.desVarConstr={'ValVolFrac'};
    paroptim.constraint.desVarVal={0.1};
    paroptim.constraint.resConstr={'AeroResidualBarrier'};
    paroptim.constraint.resVal={[-0.5,0.5]};
    
end

function [paroptim]=MeanVolumeConstraint_30(paroptim)
    
    paroptim.constraint.desVarConstr={'MeanVolFrac'};
    paroptim.constraint.desVarVal={0.3};
    paroptim.constraint.resConstr={'AeroResidualBarrier'};
    paroptim.constraint.resVal={[-0.5,0.5]};
    
end

function [paroptim]=SumVolumeConstraint(paroptim)
    
    paroptim.constraint.desVarConstr={'MinSumVolFrac'};
    paroptim.constraint.desVarVal={3.8};
    paroptim.constraint.resConstr={'AeroResidualBarrier'};
    paroptim.constraint.resVal={[-0.5,0.5]};
    
end

function [paroptim]=SumVolumeConstraint2(paroptim)
    
    paroptim.constraint.desVarConstr={'MinSumVolFrac'};
    paroptim.constraint.desVarVal={4};
    paroptim.constraint.resConstr={' '};
    paroptim.constraint.resVal={[]};
    
end

function [paroptim]=ConstraintArea(paroptim)
    
    paroptim.constraint.resConstr={' '};
    paroptim.constraint.resVal={[]};
    
end

function [paroptim]=LocalVolumeConstraint(paroptim)
    
    paroptim.constraint.desVarConstr={' '};
    paroptim.constraint.desVarVal={[]};
    paroptim.constraint.initConstr={'LocalVolFrac_image'};
    paroptim.constraint.initVal={{'.\Active_Build\ConstraintFiles\smile_5b12.png','min'}};
    paroptim.constraint.resConstr={'AeroResidualBarrier'};
    paroptim.constraint.resVal={[-0.5,0.5]};
    
end

function [paroptim]=NoConstraint(paroptim)
    paroptim.constraint.initConstr={};
    paroptim.constraint.initVal={};
    paroptim.constraint.desVarConstr={};
    paroptim.constraint.desVarVal={};
    paroptim.constraint.resConstr={};
    paroptim.constraint.resVal={};
    
end

function [paroptim]=NACA0012Constraint(paroptim)
    
    paroptim.constraint.desVarConstr={'Naca0012'};
    paroptim.constraint.desVarVal={[]};
    paroptim.constraint.resConstr={'AeroResidualBarrier'};
    paroptim.constraint.resVal={[-0.5,0.5]};
    
end

function [paroptim]=SnaxVolResConstraint(paroptim)
    paroptim.constraint.initConstr={};
    paroptim.constraint.initVal={};
    paroptim.constraint.desVarConstr={};
    paroptim.constraint.desVarVal={};
    paroptim.constraint.resConstr={'SnakResBarrier'};
    paroptim.constraint.resVal={1};
    
end

% Objectives
function paroptim=CutCellObjective(paroptim)
    
    paroptim.general.objectiveName='CutCellFlow';
    paroptim.general.direction='min';
    paroptim.general.defaultVal=1e3;
end

function paroptim=LengthAreaObjective(paroptim)
    
    paroptim.general.objectiveName='LengthArea';
    paroptim.general.direction='max';
    paroptim.general.defaultVal=-1e3;
end

function paroptim=InvDesObjective(paroptim)
    
    paroptim.general.objectiveName='InverseDesign';
    paroptim.general.direction='min';
    paroptim.general.defaultVal=1000;
    
end

% Run Sizes
function paroptim=FullOpt_bp3(paroptim)
    
    paroptim.general.nPop=48;
    paroptim.general.maxIter=150;
    paroptim.general.worker=12;
    
end

function paroptim=FullOpt_bp2(paroptim) %#ok<*DEFNU>
    
    paroptim.general.nPop=48;
    paroptim.general.maxIter=100;
    paroptim.general.worker=8;
    
end

function paroptim=FullOpt_Desktop(paroptim)
    
    paroptim.general.nPop=32;
    paroptim.general.maxIter=50;
    paroptim.general.worker=8;
    
end

function paroptim=Test_bp3(paroptim)
    paroptim.general.nPop=12;
    paroptim.general.maxIter=3;
    paroptim.general.worker=12;
    
end

function paroptim=Test_bp2(paroptim)
    paroptim.general.nPop=8;
    paroptim.general.maxIter=3;
    paroptim.general.worker=8;
    
end

function paroptim=Test_Desktop(paroptim)
    
    paroptim.general.nPop=6;
    paroptim.general.maxIter=3;
    paroptim.general.worker=6;
    
end

%% Standard Blocks (No Iter, pop, worker included)

function [paroptim]=CG_Aero()
    
    [paroptim]=DefaultOptim();
    
    paroptim=ModifySnakesParam(paroptim,'optimSupersonic');
    
    [paroptim]=MeanVolumeConstraint(paroptim);
    paroptim=CutCellObjective(paroptim);
    
    [paroptim]=OptimCG(paroptim);
    paroptim.optim.CG.varActive='snaksensiv';
    paroptim.parametrisation.optiminit.modeSmoothType='peaksmooth'; % 'peaksmooth' 'polysmooth';
    paroptim.parametrisation.optiminit.modeSmoothNum=4;
    paroptim.parametrisation.general.subdivType='chaikin';
    paroptim.parametrisation.snakes.refine.axisRatio=1;
    paroptim.general.symType='horz'; % 'horz'
    paroptim.optim.CG.diffStepSize=[1e-4,-1e-4];
    paroptim.optim.CG.minDiffStep=1e-6;
    paroptim.optim.CG.validVol=0.2;
end

function [paroptim]=CG_Area()
    
    [paroptim]=DefaultOptim();
    
    paroptim=ModifySnakesParam(paroptim,'optimSupersonic');
    
    [paroptim]=SumVolumeConstraint2(paroptim);
    
    [paroptim]=OptimCG(paroptim);
    paroptim.parametrisation.general.subdivType='chaikin';
    paroptim.parametrisation.snakes.refine.axisRatio=1;
    paroptim.general.symType='horz'; % 'horz'
    
end

function [paroptim]=DE_Aero()
    
    % Load
    [paroptim]=DefaultOptim();
    
    % Standard modifications
    paroptim=ModifySnakesParam(paroptim,'optimSupersonic');
    [paroptim]=MeanVolumeConstraint(paroptim);
    paroptim=CutCellObjective(paroptim);
    [paroptim]=OptimDE(paroptim);
    
    % Additional modifications
    paroptim.general.startPop='rand';
    paroptim.parametrisation.general.subdivType='chaikin';
    paroptim.parametrisation.snakes.refine.axisRatio=1;
    paroptim.general.symType='horz'; % 'horz'
    
end

function [paroptim]=CG_NACA0012()
    
    [paroptim]=DefaultOptim();
    
    paroptim=ModifySnakesParam(paroptim,'optimNACA0012');
    
    [paroptim]=NACA0012Constraint(paroptim);
    paroptim=CutCellObjective(paroptim);
    
    [paroptim]=OptimCG(paroptim);
    
    paroptim.general.startPop='NACA0012';
    paroptim.optim.CG.diffStepSize=[1e-6,-1e-6];
    paroptim.optim.CG.minDiffStep=1e-7;
    paroptim.optim.CG.maxDiffStep=1e-5;
    paroptim.optim.CG.varActive='snaksensiv';
    paroptim.optim.CG.validVol=0.2;
    
    paroptim.parametrisation.optiminit.modeSmoothType='peaksmooth'; % 'peaksmooth' 'polysmooth';
    paroptim.parametrisation.optiminit.modeSmoothNum=5;
    
    paroptim.obj.flow.nMach=0.85;
    paroptim.obj.flow.CFDfolder=[cd,'\Result_Template\CFD_code_Template\transonic'];
    paroptim.obj.flow.stoponerror=false;
    paroptim.obj.flow.targConv=-6;
    paroptim.obj.flow.lengthConvTest=100;
    paroptim.obj.flow.startIterFlow=10000;
    paroptim.obj.flow.restartIter=5000;
    paroptim.obj.flow.maxRestart=10;
    paroptim.obj.flow.isSymFlow=true;
    
    paroptim.parametrisation.general.subdivType='chaikin';
    paroptim.parametrisation.snakes.refine.axisRatio=1;
    paroptim.general.symType='horz'; % 'horz'
    
end

function [paroptim]=MultiTopo_DEhoriz()
    
    % Root param
    [paroptim]=DefaultOptim();
    % Standard Modifications
    paroptim=Test_Desktop(paroptim);
    paroptim=ModifySnakesParam(paroptim,'optimSupersonicMultiTopo');
    [paroptim]=SumVolumeConstraint(paroptim);
    paroptim=CutCellObjective(paroptim);
    [paroptim]=OptimDE_horiz(paroptim);
    
    paroptim.general.symType='horz'; % 'horz'
    paroptim.general.startPop='initbusemann';
    
    
    
end

function [paroptim]=MultiTopo_CGhoriz()
    
    [paroptim]=DefaultOptim();
    paroptim=ModifySnakesParam(paroptim,'optimSupersonicMultiTopo');
    [paroptim]=SumVolumeConstraint(paroptim);
    paroptim=CutCellObjective(paroptim);
    [paroptim]=OptimCG(paroptim);
    
    paroptim.optim.CG.varActive='border'; % 'wideborder'
    paroptim.parametrisation.general.subdivType='chaikin';
    paroptim.general.symType='horz'; % 'horz'
    
end

function [paroptim]=Component_CG()
    
    [paroptim]=DefaultOptim();
    paroptim=ModifySnakesParam(paroptim,'SupersonicComponent');
    [paroptim]=LocalVolumeConstraint(paroptim);
    paroptim=CutCellObjective(paroptim);
    [paroptim]=OptimCG(paroptim);
    paroptim.general.startPop='innerbound';
    paroptim.optim.CG.varActive='wideborder'; % 'wideborder'
    paroptim.parametrisation.general.subdivType='chaikin';
    paroptim.general.symType='none'; % 'horz'
    
end

function [paroptim]=Component_DE()
    
    [paroptim]=DefaultOptim();
    paroptim=ModifySnakesParam(paroptim,'SupersonicComponent');
    [paroptim]=LocalVolumeConstraint(paroptim);
    paroptim=CutCellObjective(paroptim);
    [paroptim]=OptimDE_horiz(paroptim);
    
    paroptim.general.startPop='initaeroshell';
    paroptim.general.varOverflow='spill'; % 'truncate' 'spill'
    paroptim.general.optimMethod='DE';
    paroptim.parametrisation.general.subdivType='chaikin';
    paroptim.general.symType='none'; % 'horz'
    paroptim.general.varOverflow='spill'; % 'truncate' 'spill'
    
end

function [paroptim]=Inverse_CG()
    
    [paroptim]=DefaultOptim();
    paroptim=ModifySnakesParam(paroptim,'optimInverseDesign');
    [paroptim]=SnaxVolResConstraint(paroptim);
    paroptim=InvDesObjective(paroptim);
    [paroptim]=OptimCG(paroptim);
    paroptim.general.startPop='halfuniformsharp';
    paroptim.optim.CG.varActive='snaksensiv';
    paroptim.optim.CG.sensCalc='snake';
    paroptim.parametrisation.optiminit.modeSmoothType='peaksmooth'; % 'peaksmooth' 'polysmooth';
    paroptim.parametrisation.optiminit.modeSmoothNum=6;
    paroptim.optim.CG.diffStepSize=[1e-3,-1e-3];
    paroptim.optim.CG.validVol=0.2;
    
    paroptim.spline.splineCase='inversedesign2';
    paroptim.spline.resampleSnak=true;
    
    paroptim.parametrisation.general.subdivType='chaikin';
    paroptim.general.symType='none'; % 'horz'
    
end

function [paroptim]=Inverse_Bulk()
    
    [paroptim]=DefaultOptim();
    paroptim=ModifySnakesParam(paroptim,'optimInverseDesign');
    %[paroptim]=SnaxVolResConstraint(paroptim);
    paroptim=InvDesObjective(paroptim);
    [paroptim]=OptimDE(paroptim);
    paroptim.general.objectiveName='InverseDesignBulk';
    
    paroptim.spline.splineCase='inversedesign2';
    paroptim.spline.resampleSnak=true;
    
    paroptim.parametrisation.general.subdivType='chaikin';
    paroptim.general.symType='none'; % 'horz'
    
end
%% Callable functions

function [paroptim]=TestOptim()
    
    [paroptim]=DefaultOptim();
    
    paroptim=ModifySnakesParam(paroptim,'optimSupersonic');
    
    [paroptim]=SumVolumeConstraint2(paroptim);
    
    [paroptim]=OptimCG(paroptim);
    paroptim.general.optimMethod='conjgrad';
    paroptim.parametrisation.general.subdivType='chaikin';
    paroptim.parametrisation.snakes.refine.axisRatio=1;
    paroptim.general.symType='horz'; % 'horz'
    
end

function [paroptim]=StandardOptim()
    
    [paroptim]=DefaultOptim();
    
end

function [paroptim]=Test_CG_desktop()
    
    [paroptim]=DefaultOptim();
    [paroptim]=OptimCG(paroptim);
    paroptim=Test_Desktop(paroptim);
    paroptim.parametrisation.general.subdivType='chaikin';
    paroptim.general.nPop=9;
    paroptim.general.symType='horz'; % 'horz'
    
    paroptim.parametrisation.snakes.refine.axisRatio=0.5;
    paroptim.general.maxIter=20;
end

function [paroptim]=Test_CG_Aero()
    
    [paroptim]=CG_Aero();
    
    paroptim.parametrisation.snakes.refine.axisRatio=2.2;
    
    paroptim.general.nPop=12;
    paroptim.general.maxIter=50;
    paroptim.general.worker=4;
end

function [paroptim]=Test_CG_Area()
    
    [paroptim]=CG_Area();
    
    paroptim.parametrisation.snakes.refine.axisRatio=1;
    
    paroptim.general.nPop=12;
    paroptim.general.maxIter=4;
    paroptim.general.worker=4;
end

function [paroptim]=TestParOptim_desktop()
    
    [paroptim]=DefaultOptim();
    
    paroptim=Test_Desktop(paroptim);
    paroptim.parametrisation.general.subdivType='chaikin';
end

function [paroptim]=TestParOptimSym_desktop()
    
    [paroptim]=DefaultOptim();
    paroptim=Test_Desktop(paroptim);
    
    paroptim.parametrisation.general.subdivType='chaikin';
    paroptim.general.symType='horz'; % 'horz'
end

function [paroptim]=Test_MultiTopo_M2_CG()
    
    [paroptim]=MultiTopo_CGhoriz();
    paroptim=Test_Desktop(paroptim);
    
    paroptim.general.nPop=60;
    paroptim.general.maxIter=10;
    
    paroptim.general.worker=12;
    
end

function [paroptim]=Test_MultiTopo_Init_D()
    
    [paroptim]=MultiTopo_DEhoriz();
    
    paroptim=Test_Desktop(paroptim);
    [paroptim]=SumVolumeConstraint2(paroptim);
    paroptim=LengthAreaObjective(paroptim);
    
    paroptim.general.nPop=12;
    paroptim.general.maxIter=3;
    paroptim.general.worker=6;
end

function [paroptim]=Test_MultiTopo_Init_D2()
    
    [paroptim]=MultiTopo_DEhoriz();
    paroptim=Test_Desktop(paroptim);
    
    paroptim.general.nPop=32;
    paroptim.general.maxIter=3;
    
    paroptim.general.worker=8;
end

function [paroptim]=TestParOptimAero_desktop()
    
    % Root param
    [paroptim]=DefaultOptim();
    % Standard Modifications
    paroptim=Test_Desktop(paroptim);
    paroptim=ModifySnakesParam(paroptim,'optimSupersonic');
    [paroptim]=MeanVolumeConstraint(paroptim);
    paroptim=CutCellObjective(paroptim);
    
    paroptim.general.optimMethod='DEtan';
    paroptim.general.symType='horz'; % 'horz'
    
end

function [paroptim]=TestParOptim_HPC()
    
    [paroptim]=DefaultOptim();
    paroptim=Test_bp3(paroptim);
    
end

function [paroptim]=HPC_LengthArea()
    
    [paroptim]=DefaultOptim();
    paroptim=LengthAreaObjective(paroptim);
    paroptim=FullOpt_bp3(paroptim);
    
    
    paroptim.general.knownOptim=[0.2*(8+pi)/(8+2*pi)];
    
end

function [paroptim]=Test_Init()
    
    % Root param
    [paroptim]=DefaultOptim();
    % Standard Modifications
    paroptim=Test_Desktop(paroptim);
    paroptim=ModifySnakesParam(paroptim,'TestInitOptim');
    paroptim.constraint.resConstr={' '};
    paroptim.constraint.resVal={[]};
    
    paroptim.constraint.initConstr={'LocalVolFrac_image'};
    paroptim.constraint.initVal={{'.\Active_Build\ConstraintFiles\inline_5b12.png','min'}};
    paroptim=LengthAreaObjective(paroptim);
    [paroptim]=OptimDE_horiz(paroptim);
    
    paroptim.general.symType='none'; % 'horz'
    paroptim.general.startPop='initaeroshell';
    
    paroptim.general.varOverflow='spill'; % 'truncate' 'spill'
    paroptim.general.optimMethod='DE';
    paroptim.general.nPop=100;
    paroptim.general.maxIter=2;
    paroptim.general.worker=12;
end

%% Test cases for length Area

function [paroptim]=Test_smoothCG_Area()
    
    [paroptim]=DefaultOptim();
    
    paroptim=ModifySnakesParam(paroptim,'optimSupersonic');
    
    [paroptim]=MeanVolumeConstraint(paroptim);
    [paroptim]=ConstraintArea(paroptim);
    
    [paroptim]=OptimCG(paroptim);
    paroptim.parametrisation.general.subdivType='chaikin';
    paroptim.parametrisation.snakes.refine.axisRatio=1;
    paroptim.general.symType='horz'; % 'horz'
    
    paroptim.optim.CG.varActive='all';
    paroptim.parametrisation.snakes.refine.axisRatio=1;
    
    paroptim.general.nPop=12;
    paroptim.general.maxIter=4;
    paroptim.general.worker=4;
end

function [paroptim]=Test_polypeakCG_outmissile_Area()
    
    [paroptim]=DefaultOptim();
    
    paroptim=ModifySnakesParam(paroptim,'SupersonicComponent');
    
    [paroptim]=LocalVolumeConstraint(paroptim);
    [paroptim]=ConstraintArea(paroptim);
    paroptim.constraint.initVal={{'.\Active_Build\ConstraintFiles\missile_5b12.png','min'}};
    
    
    [paroptim]=OptimCG(paroptim);
    paroptim.general.startPop='outerbound';
    paroptim.parametrisation.general.subdivType='chaikin';
    paroptim.parametrisation.snakes.refine.axisRatio=1;
    paroptim.general.symType='none'; % 'horz'
    
    paroptim.optim.CG.varActive='snaksensiv';
    paroptim.parametrisation.optiminit.modeSmoothType='polypeaksmooth'; % 'peaksmooth' 'polysmooth';
    paroptim.parametrisation.optiminit.modeSmoothNum=4;
    
    paroptim.parametrisation.snakes.refine.axisRatio=1;
    
    paroptim.general.nPop=12;
    paroptim.general.maxIter=4;
    paroptim.general.worker=4;
end

function [paroptim]=Test_peakCG_outmissile_Area()
    
    [paroptim]=DefaultOptim();
    
    paroptim=ModifySnakesParam(paroptim,'SupersonicComponent');
    
    [paroptim]=LocalVolumeConstraint(paroptim);
    [paroptim]=ConstraintArea(paroptim);
    paroptim.constraint.initVal={{'.\Active_Build\ConstraintFiles\missile_5b12.png','min'}};
    
    
    [paroptim]=OptimCG(paroptim);
    paroptim.general.startPop='outerbound';
    paroptim.parametrisation.general.subdivType='chaikin';
    paroptim.parametrisation.snakes.refine.axisRatio=1;
    paroptim.general.symType='none'; % 'horz'
    
    paroptim.optim.CG.varActive='snaksensiv';
    paroptim.parametrisation.optiminit.modeSmoothType='peaksmooth'; % 'peaksmooth' 'polysmooth';
    paroptim.parametrisation.optiminit.modeSmoothNum=4;
    
    paroptim.parametrisation.snakes.refine.axisRatio=1;
    
    paroptim.general.nPop=12;
    paroptim.general.maxIter=4;
    paroptim.general.worker=4;
end

function [paroptim]=Test_polyCG_outmissile_Area()
    
    [paroptim]=DefaultOptim();
    
    paroptim=ModifySnakesParam(paroptim,'SupersonicComponent');
    
    [paroptim]=LocalVolumeConstraint(paroptim);
    [paroptim]=ConstraintArea(paroptim);
    paroptim.constraint.initVal={{'.\Active_Build\ConstraintFiles\missile_5b12.png','min'}};
    
    
    [paroptim]=OptimCG(paroptim);
    paroptim.general.startPop='outerbound';
    paroptim.parametrisation.general.subdivType='chaikin';
    paroptim.parametrisation.snakes.refine.axisRatio=1;
    paroptim.general.symType='none'; % 'horz'
    
    paroptim.optim.CG.varActive='snaksensiv';
    paroptim.parametrisation.optiminit.modeSmoothType='polysmooth'; % 'peaksmooth' 'polysmooth';
    paroptim.parametrisation.optiminit.modeSmoothNum=4;
    
    paroptim.parametrisation.snakes.refine.axisRatio=1;
    
    paroptim.general.nPop=12;
    paroptim.general.maxIter=4;
    paroptim.general.worker=4;
end

function [paroptim]=SmoothCG_outmis_Aero()
    
    [paroptim]=DefaultOptim();
    
    paroptim=ModifySnakesParam(paroptim,'SupersonicComponent');
    
    [paroptim]=LocalVolumeConstraint(paroptim);
    paroptim.constraint.initVal={{'.\Active_Build\ConstraintFiles\missile2_5b12.png','min'}};
    
    
    paroptim=CutCellObjective(paroptim);
    
    [paroptim]=OptimCG(paroptim);
    paroptim.general.startPop='outerbound';
    paroptim.parametrisation.general.subdivType='chaikin';
    paroptim.parametrisation.snakes.refine.axisRatio=1;
    paroptim.general.symType='none'; % 'horz'
    
    paroptim.optim.CG.varActive='snaksensiv';
    
    paroptim.parametrisation.optiminit.modeSmoothType='polypeaksmooth'; % 'peaksmooth' 'polysmooth';
    paroptim.parametrisation.optiminit.modeSmoothNum=6;
    
    paroptim.parametrisation.snakes.refine.axisRatio=1;
    
    paroptim.general.nPop=12;
    paroptim.general.maxIter=20;
    paroptim.general.worker=4;
end

%% Analytical test cases

function [paroptim]=Test_Rosenbrock()
    
    [paroptim]=DefaultOptim();
    [paroptim]=NoConstraint(paroptim);
    
    paroptim.optim.CG.diffStepSize=[1e-3,-1e-3];
    paroptim.optim.CG.minDiffStep=1e-8;
    paroptim.optim.CG.varActive='all';
    paroptim.optim.CG.validVol=1;
    paroptim.optim.CG.stepAlgo='conjgrad';
    
    paroptim.general.startPop='Rosen';
    paroptim.general.direction='min';
    paroptim.general.optimMethod='conjgrad';
    paroptim.general.symType='none';
    paroptim.general.objectiveName='Rosenbrock';
    paroptim.general.useSnake=false;
    paroptim.general.nDesVar=4*2;
    paroptim.general.varOverflow='truncate';
    paroptim.general.desVarRange=[-2.048, 2.048];
    paroptim.general.nPop=25;
    paroptim.general.maxIter=100;
    paroptim.general.worker=4;
    paroptim.general.knownOptim=[0 1 1];
end

function [paroptim]=Test_Rosenbrock_BFGS()
    [paroptim]=Test_Rosenbrock();
    
    paroptim.optim.CG.stepAlgo='BFGS';
    
end
%% refinement

function paroptim=test_refine1()
    [paroptim]=Inverse_CG();
    paroptim=ModifySnakesParam(paroptim,'optimInverseDesign_cosref');
    paroptim.obj.invdes.aeroName='0012';
    
    paroptim.optim.CG.lineSearchType='backbisection';
    paroptim.general.refineOptim=[2 1 20; 2 1 20];
    %paroptim.optim.CG.varActive='all';
    paroptim.general.startPop='halfuniformthin';
    paroptim.general.nPop=12;
    paroptim.general.maxIter=4;
    paroptim.general.worker=4;
end
function paroptim=test_refine2()
    [paroptim]=Inverse_CG();
    paroptim=ModifySnakesParam(paroptim,'optimInverseDesign_cosref');
    paroptim.obj.invdes.aeroName='2212';
    paroptim.optim.CG.lineSearchType='backbisection';
    paroptim.general.refineOptim=[2 1 20; 2 1 20];
    %paroptim.optim.CG.varActive='all';
    paroptim.general.startPop='halfuniformthin';
    paroptim.general.nPop=12;
    paroptim.general.maxIter=20;
    paroptim.general.worker=4;
end
function paroptim=test_refine3()
    [paroptim]=Inverse_CG();
    paroptim=ModifySnakesParam(paroptim,'optimInverseDesign_cosref');
    paroptim.obj.invdes.aeroName='0012';
    
    paroptim.optim.CG.lineSearchType='polydistrib';
    paroptim.general.refineOptim=[2 1 20; 2 1 20];
    %paroptim.optim.CG.varActive='all';
    paroptim.general.startPop='halfuniformthin';
    paroptim.general.nPop=12;
    paroptim.general.maxIter=20;
    paroptim.general.worker=4;
end
function paroptim=test_refine4()
    [paroptim]=Inverse_CG();
    paroptim=ModifySnakesParam(paroptim,'optimInverseDesign_cosref');
    paroptim.obj.invdes.aeroName='2212';
    paroptim.optim.CG.lineSearchType='polydistrib';
    paroptim.general.refineOptim=[2 1 20; 2 1 20];
    %paroptim.optim.CG.varActive='all';
    paroptim.general.startPop='halfuniformthin';
    paroptim.general.nPop=12;
    paroptim.general.maxIter=20;
    paroptim.general.worker=4;
end
function paroptim=test_refine5()
    [paroptim]=Inverse_CG();
    paroptim=ModifySnakesParam(paroptim,'optimInverseDesign_cosref');
    paroptim.obj.invdes.aeroName='4412';
    paroptim.optim.CG.lineSearchType='backbisection';
    paroptim.general.refineOptim=[2 1 20; 2 1 20];
    %paroptim.optim.CG.varActive='all';
    paroptim.general.startPop='halfuniformthin';
    paroptim.general.nPop=12;
    paroptim.general.maxIter=20;
    paroptim.general.worker=4;
end
function paroptim=test_refine6()
    [paroptim]=Inverse_CG();
    paroptim=ModifySnakesParam(paroptim,'optimInverseDesign_cosref');
    paroptim.obj.invdes.aeroName='4412';
    
    paroptim.optim.CG.lineSearchType='polydistrib';
    paroptim.general.refineOptim=[2 1 20; 2 1 20];
    %paroptim.optim.CG.varActive='all';
    paroptim.general.startPop='halfuniformthin';
    paroptim.general.nPop=12;
    paroptim.general.maxIter=20;
    paroptim.general.worker=4;
end

function paroptim=test_NACA0012()
    
    [paroptim]=CG_NACA0012();
    
    paroptim.general.refineOptim=[2 1; 2 1];
    %paroptim.optim.CG.varActive='all';
    paroptim.general.nPop=12;
    paroptim.general.maxIter=6;
    paroptim.general.worker=4;
end

function paroptim=bp3_NACA0012()
    
    [paroptim]=CG_NACA0012();
    
    paroptim.general.refineOptim=[2 1 100; 2 1 100];
    %paroptim.optim.CG.varActive='all';
    paroptim.general.nPop=12;
    paroptim.general.maxIter=60;
    paroptim.general.worker=12;
end

function paroptim=bp3_NACA0012S()
    
    [paroptim]=CG_NACA0012();
    paroptim=ModifySnakesParam(paroptim,'optimNACA0012Sc');
    paroptim.general.refineOptim=[2 1 100; 2 1 100; 2 1 100];
    %paroptim.optim.CG.varActive='all';
    paroptim.general.nPop=12;
    paroptim.general.maxIter=100;
    paroptim.general.worker=12;
end

function paroptim=bp3_NACA0012Su()
    
    [paroptim]=CG_NACA0012();
    paroptim=ModifySnakesParam(paroptim,'optimNACA0012Su');
    paroptim.general.refineOptim=[2 1 100; 2 1 100; 2 1 100];
    %paroptim.optim.CG.varActive='all';
    paroptim.general.nPop=100;
    paroptim.general.maxIter=60;
    paroptim.general.worker=12;
end

function paroptim=bp3_NACA0012L()
    
    [paroptim]=CG_NACA0012();
    
    paroptim=ModifySnakesParam(paroptim,'optimNACA0012L');
    paroptim.general.refineOptim=[2 1 60];
    %paroptim.optim.CG.varActive='all';
    paroptim.general.nPop=12;
    paroptim.general.maxIter=60;
    paroptim.general.worker=12;
end

function paroptim=bp3_NACA0012NoRef()
    
    [paroptim]=CG_NACA0012();
    
    paroptim=ModifySnakesParam(paroptim,'optimNACA0012L');
    paroptim.general.refineOptim=zeros([0 2]);
    %paroptim.optim.CG.varActive='all';
    paroptim.general.nPop=12;
    paroptim.general.maxIter=120;
    paroptim.general.worker=12;
end

function [paroptim]=Refine10to40()
    [paroptim]=CG_Aero();
    paroptim=ModifySnakesParam(paroptim,'optimSupersonicCos');
    paroptim.parametrisation.snakes.refine.axisRatio=1.25;
    paroptim.general.startPop='halfuniformsharp';
    
    paroptim.optim.CG.diffStepSize=[1e-3,-1e-3];
    paroptim.optim.CG.varActive='snaksensiv';
    paroptim.optim.CG.validVol=0.3;
    paroptim.parametrisation.optiminit.modeSmoothType='peaksmooth'; % 'peaksmooth' 'polysmooth';
    paroptim.parametrisation.optiminit.modeSmoothNum=4;
    paroptim.general.refineOptim=[2 1; 2 1];
    paroptim.general.nPop=12;
    paroptim.general.maxIter=32;
    paroptim.general.worker=12;
end

function [paroptim]=Refine10to40_val()
    [paroptim]=CG_Aero();
    [paroptim]=ValVolumeConstraint(paroptim);
    paroptim=ModifySnakesParam(paroptim,'optimSupersonicCos');
    paroptim.parametrisation.snakes.refine.axisRatio=1.25;
    paroptim.general.startPop='halfuniformsharp';
    
    paroptim.optim.CG.diffStepSize=[1e-3,-1e-3];
    paroptim.optim.CG.varActive='snaksensiv';
    paroptim.optim.CG.validVol=0.3;
    paroptim.parametrisation.optiminit.modeSmoothType='peaksmooth'; % 'peaksmooth' 'polysmooth';
    paroptim.parametrisation.optiminit.modeSmoothNum=4;
    paroptim.general.refineOptim=[2 1 36; 2 1 36];
    paroptim.general.nPop=12;
    paroptim.general.maxIter=32;
    paroptim.general.worker=12;
end

function [paroptim]=Refine40()
    [paroptim]=CG_Aero();
    paroptim=ModifySnakesParam(paroptim,'optimSupersonicCos');
    paroptim.parametrisation.optiminit.cellLevels=[42,2];
    paroptim.parametrisation.general.passDomBounds=MakeCartesianGridBoundsInactE(paroptim.parametrisation.optiminit.cellLevels);
    paroptim.parametrisation.snakes.refine.axisRatio=1.25*4;
    paroptim.parametrisation.optiminit.modeSmoothType='peaksmooth'; % 'peaksmooth' 'polysmooth';
    paroptim.parametrisation.optiminit.modeSmoothNum=4;
    paroptim.initparam=DefaultSnakeInit(paroptim.parametrisation);
    
    paroptim.general.startPop='halfuniformsharp';
    paroptim.optim.CG.diffStepSize=[1e-3,-1e-3];
    paroptim.optim.CG.varActive='snaksensiv';
    paroptim.optim.CG.validVol=0.3;
    paroptim.general.refineOptim=[0];
    paroptim.general.nPop=12;
    paroptim.general.maxIter=32;
    paroptim.general.worker=12;
end

function paroptim=bp3_refsweep_cv02212_interp()
    paroptim=bp3_refsweep_cv00012();
    paroptim.obj.invdes.aeroName='2212';
    paroptim=ModifySnakesParam(paroptim,'optimInverseDesign_cv');
    
    paroptim.general.startPop='NACA0012';
end

% NACA 0012 sweep

function paroptim=bp3_NACA0012_sweep()
    
    [paroptim]=CG_NACA0012();
    paroptim=ModifySnakesParam(paroptim,'optimNACA0012Sc');
    paroptim.general.refineOptim=[2 1 100; 2 1 100; 2 1 100];
    %paroptim.optim.CG.varActive='all';
    paroptim.general.nPop=12;
    paroptim.general.maxIter=100;
    paroptim.general.worker=8;
end

function paroptim=NACA0012sweep_Sc()
    
    paroptim=bp3_NACA0012_sweep();
    paroptim=ModifySnakesParam(paroptim,'optimNACA0012Scv');
    paroptim.general.refineOptim=[2 1 100; 2 1 100; 2 1 100];
    paroptim.parametrisation.general.subdivType='chaikinNaca0012';
end
function paroptim=NACA0012sweep_Su()
    
    paroptim=bp3_NACA0012_sweep();
    paroptim=ModifySnakesParam(paroptim,'optimNACA0012Suv');
    paroptim.general.refineOptim=[2 1 100; 2 1 100; 2 1 100];
    paroptim.parametrisation.general.subdivType='chaikinNaca0012';
end
function paroptim=NACA0012sweep_Nc()
    
    paroptim=bp3_NACA0012_sweep();
    paroptim=ModifySnakesParam(paroptim,'optimNACA0012Ncv');
    paroptim.general.refineOptim=[2 1 100; 2 1 100];
    paroptim.parametrisation.general.subdivType='chaikinNaca0012';
end
function paroptim=NACA0012sweep_Nu()
    
    paroptim=bp3_NACA0012_sweep();
    paroptim=ModifySnakesParam(paroptim,'optimNACA0012Nuv');
    paroptim.general.refineOptim=[2 1 100; 2 1 100];
    paroptim.parametrisation.general.subdivType='chaikinNaca0012';
end
function paroptim=NACA0012sweep_Lc()
    
    paroptim=bp3_NACA0012_sweep();
    paroptim=ModifySnakesParam(paroptim,'optimNACA0012Lcv');
    paroptim.general.refineOptim=[2 1 100];
    paroptim.parametrisation.general.subdivType='chaikinNaca0012';
end
function paroptim=NACA0012sweep_Lu()
    
    paroptim=bp3_NACA0012_sweep();
    paroptim=ModifySnakesParam(paroptim,'optimNACA0012Luv');
    paroptim.general.refineOptim=[2 1 100];
    paroptim.parametrisation.general.subdivType='chaikinNaca0012';
end
function paroptim=NACA0012sweep_Vc()
    
    paroptim=bp3_NACA0012_sweep();
    paroptim=ModifySnakesParam(paroptim,'optimNACA0012Vcv');
    paroptim.general.refineOptim=[0];
    paroptim.parametrisation.general.subdivType='chaikinNaca0012';
end
function paroptim=NACA0012sweep_Vu()
    
    paroptim=bp3_NACA0012_sweep();
    paroptim=ModifySnakesParam(paroptim,'optimNACA0012Vuv');
    paroptim.general.refineOptim=[0];
    paroptim.parametrisation.general.subdivType='chaikinNaca0012';
end

function paroptim=bp3_NACA0012_sweep_BFGS()
    
    [paroptim]=CG_NACA0012();
    [paroptim]=OptimBFGS(paroptim);
    paroptim=ModifySnakesParam(paroptim,'optimNACA0012Scv');
    paroptim.general.refineOptim=[2 1 100; 2 1 100; 2 1 100];
    %paroptim.optim.CG.varActive='all';
    paroptim.general.nPop=12;
    paroptim.general.maxIter=100;
    paroptim.general.worker=8;
end

function paroptim=NACA0012sweepBFGS_Sc()
    
    paroptim=bp3_NACA0012_sweep_BFGS();
    paroptim=ModifySnakesParam(paroptim,'optimNACA0012Scv');
    paroptim.general.refineOptim=[2 1 100; 2 1 100; 2 1 100];
    paroptim.parametrisation.general.subdivType='chaikinNaca0012';
end
function paroptim=NACA0012sweepBFGS_Su()
    
    paroptim=bp3_NACA0012_sweep_BFGS();
    paroptim=ModifySnakesParam(paroptim,'optimNACA0012Suv');
    paroptim.general.refineOptim=[2 1 100; 2 1 100; 2 1 100];
    paroptim.parametrisation.general.subdivType='chaikinNaca0012';
end
function paroptim=NACA0012sweepBFGS_Nc()
    
    paroptim=bp3_NACA0012_sweep_BFGS();
    paroptim=ModifySnakesParam(paroptim,'optimNACA0012Ncv');
    paroptim.general.refineOptim=[2 1 100; 2 1 100];
    paroptim.parametrisation.general.subdivType='chaikinNaca0012';
end
function paroptim=NACA0012sweepBFGS_Nu()
    
    paroptim=bp3_NACA0012_sweep_BFGS();
    paroptim=ModifySnakesParam(paroptim,'optimNACA0012Nuv');
    paroptim.general.refineOptim=[2 1 100; 2 1 100];
    paroptim.parametrisation.general.subdivType='chaikinNaca0012';
end
function paroptim=NACA0012sweepBFGS_Lc()
    
    paroptim=bp3_NACA0012_sweep_BFGS();
    paroptim=ModifySnakesParam(paroptim,'optimNACA0012Lcv');
    paroptim.general.refineOptim=[2 1 100];
    paroptim.parametrisation.general.subdivType='chaikinNaca0012';
end
function paroptim=NACA0012sweepBFGS_Lu()
    
    paroptim=bp3_NACA0012_sweep_BFGS();
    paroptim=ModifySnakesParam(paroptim,'optimNACA0012Luv');
    paroptim.general.refineOptim=[2 1 100];
    paroptim.parametrisation.general.subdivType='chaikinNaca0012';
end
function paroptim=NACA0012sweepBFGS_Vc()
    
    paroptim=bp3_NACA0012_sweep_BFGS();
    paroptim=ModifySnakesParam(paroptim,'optimNACA0012Vcv');
    paroptim.general.refineOptim=[0];
    paroptim.parametrisation.general.subdivType='chaikinNaca0012';
end
function paroptim=NACA0012sweepBFGS_Vu()
    
    paroptim=bp3_NACA0012_sweep_BFGS();
    paroptim=ModifySnakesParam(paroptim,'optimNACA0012Vuv');
    paroptim.general.refineOptim=[0];
    paroptim.parametrisation.general.subdivType='chaikinNaca0012';
end

function paroptim=RestartNACA0012sweep_Lc()
    
    paroptim=bp3_NACA0012_sweep();
    paroptim=ModifySnakesParam(paroptim,'optimNACA0012Lcv');
    paroptim.general.refineOptim=[0];
    paroptim.parametrisation.general.subdivType='chaikinNaca0012';
    paroptim.general.worker=12;
end

function paroptim=RestartNACA0012sweep_Nc()
    
    paroptim=bp3_NACA0012_sweep();
    paroptim=ModifySnakesParam(paroptim,'optimNACA0012Ncv');
    paroptim.general.refineOptim=[0];
    paroptim.optim.CG.diffStepSize=[1e-6,-1e-6]; %[0,2]
    paroptim.optim.CG.minDiffStep=1e-7;
    paroptim.optim.CG.validVol=0.0025;
    paroptim.parametrisation.general.subdivType='chaikinNaca0012';
    paroptim.general.worker=12;
end

function paroptim=RestartNACA0012sweep_Nc2()
    
    paroptim=bp3_NACA0012_sweep();
    paroptim=ModifySnakesParam(paroptim,'optimNACA0012Ncv');
    paroptim.obj.flow.CFDfolder=[cd,'\Result_Template\CFD_code_Template\transonicfine'];
    paroptim.general.refineOptim=[0];
    paroptim.optim.CG.diffStepSize=[1e-6,-1e-6]; %[0,2]
    paroptim.optim.CG.minDiffStep=1e-7;
    paroptim.optim.CG.validVol=0.0025;
    paroptim.parametrisation.general.subdivType='chaikinNaca0012';
    paroptim.general.worker=12;
end
% Area M2 sweep

function paroptim=bp3_AreaM2sweep()
    
    [paroptim]=CG_Aero();
    [paroptim]=ValVolumeConstraint(paroptim);
    paroptim.general.startPop='halfuniformsharp';
    paroptim=ModifySnakesParam(paroptim,'optimNACA0012Sc');
    paroptim.parametrisation.snakes.refine.LEShrink=true;
    paroptim.general.refineOptim=[2 1 100; 2 1 100; 2 1 100];
    %paroptim.optim.CG.varActive='all';
    paroptim.general.nPop=12;
    paroptim.general.maxIter=100;
    paroptim.general.worker=8;
end

function paroptim=AreaM2sweep_Sc()
    
    paroptim=bp3_AreaM2sweep();
    paroptim=ModifySnakesParam(paroptim,'optimNACA0012Sc');
    paroptim.parametrisation.snakes.refine.LEShrink=true;
    paroptim.general.refineOptim=[2 1 100; 2 1 100; 2 1 100];
end
function paroptim=AreaM2sweep_Su()
    
    paroptim=bp3_AreaM2sweep();
    paroptim=ModifySnakesParam(paroptim,'optimNACA0012Su');
    paroptim.parametrisation.snakes.refine.LEShrink=true;
    paroptim.general.refineOptim=[2 1 100; 2 1 100; 2 1 100];
end
function paroptim=AreaM2sweep_Nc()
    
    paroptim=bp3_AreaM2sweep();
    paroptim=ModifySnakesParam(paroptim,'optimNACA0012Nc');
    paroptim.parametrisation.snakes.refine.LEShrink=true;
    paroptim.general.refineOptim=[2 1 100; 2 1 100];
end
function paroptim=AreaM2sweep_Nu()
    
    paroptim=bp3_AreaM2sweep();
    paroptim=ModifySnakesParam(paroptim,'optimNACA0012Nu');
    paroptim.parametrisation.snakes.refine.LEShrink=true;
    paroptim.general.refineOptim=[2 1 100; 2 1 100];
end
function paroptim=AreaM2sweep_Lc()
    
    paroptim=bp3_AreaM2sweep();
    paroptim=ModifySnakesParam(paroptim,'optimNACA0012Lc');
    paroptim.parametrisation.snakes.refine.LEShrink=true;
    paroptim.general.refineOptim=[2 1 100];
end
function paroptim=AreaM2sweep_Lu()
    
    paroptim=bp3_AreaM2sweep();
    paroptim=ModifySnakesParam(paroptim,'optimNACA0012Lu');
    paroptim.parametrisation.snakes.refine.LEShrink=true;
    paroptim.general.refineOptim=[2 1 100];
end
function paroptim=AreaM2sweep_Vc()
    
    paroptim=bp3_AreaM2sweep();
    paroptim=ModifySnakesParam(paroptim,'optimNACA0012Vc');
    paroptim.parametrisation.snakes.refine.LEShrink=true;
    paroptim.general.refineOptim=[0];
end
function paroptim=AreaM2sweep_Vu()
    
    paroptim=bp3_AreaM2sweep();
    paroptim=ModifySnakesParam(paroptim,'optimNACA0012Vu');
    paroptim.parametrisation.snakes.refine.LEShrink=true;
    paroptim.general.refineOptim=[0];
end


function paroptim=bp3_AreaM2_sweep_BFGS()
    
    [paroptim]=CG_Aero();
    [paroptim]=OptimBFGS(paroptim);
    [paroptim]=ValVolumeConstraint(paroptim);
    paroptim.general.startPop='halfuniformsharp';
    paroptim=ModifySnakesParam(paroptim,'optimNACA0012Sc');
    paroptim.parametrisation.snakes.refine.LEShrink=true;
    paroptim.general.refineOptim=[2 1 100; 2 1 100; 2 1 100];
    %paroptim.optim.CG.varActive='all';
    paroptim.general.nPop=12;
    paroptim.general.maxIter=100;
    paroptim.general.worker=8;
end

function paroptim=AreaM2sweepBFGS_Sc()
    
    paroptim=bp3_AreaM2_sweep_BFGS();
    paroptim=ModifySnakesParam(paroptim,'optimNACA0012Sc');
    paroptim.parametrisation.snakes.refine.LEShrink=true;
    paroptim.general.refineOptim=[2 1 100; 2 1 100; 2 1 100];
end
function paroptim=AreaM2sweepBFGS_Su()
    
    paroptim=bp3_AreaM2_sweep_BFGS();
    paroptim=ModifySnakesParam(paroptim,'optimNACA0012Su');
    paroptim.parametrisation.snakes.refine.LEShrink=true;
    paroptim.general.refineOptim=[2 1 100; 2 1 100; 2 1 100];
end
function paroptim=AreaM2sweepBFGS_Nc()
    
    paroptim=bp3_AreaM2_sweep_BFGS();
    paroptim=ModifySnakesParam(paroptim,'optimNACA0012Nc');
    paroptim.parametrisation.snakes.refine.LEShrink=true;
    paroptim.general.refineOptim=[2 1 100; 2 1 100];
end
function paroptim=AreaM2sweepBFGS_Nu()
    
    paroptim=bp3_AreaM2_sweep_BFGS();
    paroptim=ModifySnakesParam(paroptim,'optimNACA0012Nu');
    paroptim.parametrisation.snakes.refine.LEShrink=true;
    paroptim.general.refineOptim=[2 1 100; 2 1 100];
end
function paroptim=AreaM2sweepBFGS_Lc()
    
    paroptim=bp3_AreaM2_sweep_BFGS();
    paroptim=ModifySnakesParam(paroptim,'optimNACA0012Lc');
    paroptim.parametrisation.snakes.refine.LEShrink=true;
    paroptim.general.refineOptim=[2 1 100];
end
function paroptim=AreaM2sweepBFGS_Lu()
    
    paroptim=bp3_AreaM2_sweep_BFGS();
    paroptim=ModifySnakesParam(paroptim,'optimNACA0012Lu');
    paroptim.parametrisation.snakes.refine.LEShrink=true;
    paroptim.general.refineOptim=[2 1 100];
end
function paroptim=AreaM2sweepBFGS_Vc()
    
    paroptim=bp3_AreaM2_sweep_BFGS();
    paroptim=ModifySnakesParam(paroptim,'optimNACA0012Vc');
    paroptim.parametrisation.snakes.refine.LEShrink=true;
    paroptim.general.refineOptim=[0];
end
function paroptim=AreaM2sweepBFGS_Vu()
    
    paroptim=bp3_AreaM2_sweep_BFGS();
    paroptim=ModifySnakesParam(paroptim,'optimNACA0012Vu');
    paroptim.parametrisation.snakes.refine.LEShrink=true;
    paroptim.general.refineOptim=[0];
end

% Test derivatives

function paroptim=TestDerivDesktop(e)
    
    paroptim=AreaM2sweep_Sc();
    paroptim.optim.CG.diffStepSize=[e,-e]; %[0,2]
    paroptim.optim.CG.minDiffStep=e;
    paroptim.general.maxIter=6;
    paroptim.general.refineOptim=0;
    paroptim.general.startPop='loadshape';
    paroptim.general.specificFillName='.\Active_Build\Input_Variables\Parabola.mat';
    paroptim.optim.CG.gradScaleType='none';
    paroptim.obj.flow.CFDfolder=[cd,'\Result_Template\CFD_code_Template\supersonic_ogivecoarse'];
    paroptim.general.worker=4;
end
function paroptim=TestDerivRosen(e)
    
    paroptim=Test_Rosenbrock();
    paroptim.optim.CG.diffStepSize=[e,-e]; %[0,2]
    paroptim.optim.CG.minDiffStep=e;
    
    paroptim.general.refineOptim=0;
    paroptim.general.maxIter=100;
    paroptim.general.worker=4;
end

function paroptim=TestDerivBFGS(e)
    
    paroptim=AreaM2sweepBFGS_Lc();
    paroptim.optim.CG.diffStepSize=[e,-e]; %[0,2]
    paroptim.optim.CG.minDiffStep=e;
    paroptim.general.maxIter=6;
    paroptim.general.refineOptim=0;
end
function paroptim=TestDeriv(e)
    
    paroptim=AreaM2sweep_Lc();
    paroptim.optim.CG.diffStepSize=[e,-e]; %[0,2]
    paroptim.optim.CG.minDiffStep=e;
    paroptim.general.maxIter=6;
    paroptim.general.refineOptim=0;
    
end
function paroptim=BuseDerivBFGS(e)
    
    paroptim=AreaM2sweepBFGS_Lc();
    paroptim.general.startPop='specificfill';
    paroptim.general.specificFillName='testgrad_busemann';
    paroptim.optim.CG.diffStepSize=[e,-e]; %[0,2]
    paroptim.optim.CG.minDiffStep=e;
    paroptim.general.maxIter=6;
    paroptim.general.refineOptim=0;
end
function paroptim=BuseDeriv(e)
    
    paroptim=AreaM2sweep_Lc();
    paroptim.general.startPop='specificfill';
    paroptim.general.specificFillName='testgrad_busemann';
    paroptim.optim.CG.diffStepSize=[e,-e]; %[0,2]
    paroptim.optim.CG.minDiffStep=e;
    paroptim.general.maxIter=6;
    paroptim.general.refineOptim=0;
    
end
function paroptim=BuseDerivCoarseBFGS(e)
    
    paroptim=AreaM2sweepBFGS_Lc();
    paroptim.general.startPop='specificfill';
    paroptim.general.specificFillName='testgrad_busemann';
    paroptim.obj.flow.CFDfolder=[cd,'\Result_Template\CFD_code_Template\supersonic_ogivecoarse'];
    paroptim.optim.CG.diffStepSize=[e,-e]; %[0,2]
    paroptim.optim.CG.minDiffStep=e;
    paroptim.general.maxIter=6;
    paroptim.general.refineOptim=0;
end
function paroptim=BuseDerivCoarse(e)
    
    paroptim=AreaM2sweep_Lc();
    paroptim.general.startPop='specificfill';
    paroptim.general.specificFillName='testgrad_busemann';
%     paroptim.general.startPop='loadshape';
%     paroptim.general.specificFillName='.\Active_Build\Input_Variables\Parabola.mat';
%     paroptim.optim.CG.gradScaleType='none';
    paroptim.obj.flow.CFDfolder=[cd,'\Result_Template\CFD_code_Template\supersonic_ogivecoarse'];
    paroptim.optim.CG.diffStepSize=[e,-e]; %[0,2]
    paroptim.optim.CG.minDiffStep=e;
    paroptim.general.maxIter=6;
    paroptim.general.refineOptim=0;
    
end
function paroptim=ParaDerivBFGS(e)
    
    paroptim=AreaM2sweepBFGS_Lc();
    paroptim.general.startPop='loadshape';
    paroptim.general.specificFillName='.\Active_Build\Input_Variables\Parabola.mat';
    paroptim.optim.CG.diffStepSize=[e,-e]; %[0,2]
    paroptim.optim.CG.minDiffStep=e;
    paroptim.general.maxIter=6;
    paroptim.general.refineOptim=0;
end
function paroptim=ParaDeriv(e)
    
    paroptim=AreaM2sweep_Lc();
    paroptim.general.startPop='loadshape';
    paroptim.general.specificFillName='.\Active_Build\Input_Variables\Parabola.mat';
    paroptim.optim.CG.diffStepSize=[e,-e]; %[0,2]
    paroptim.optim.CG.minDiffStep=e;
    paroptim.general.maxIter=6;
    paroptim.general.refineOptim=0;
    
end
function paroptim=ParaDerivUniBFGS(e)
    
    paroptim=AreaM2sweepBFGS_Lu();
    paroptim.general.startPop='loadshape';
    paroptim.general.specificFillName='.\Active_Build\Input_Variables\Parabola.mat';
    paroptim.optim.CG.diffStepSize=[e,-e]; %[0,2]
    paroptim.optim.CG.minDiffStep=e;
    paroptim.general.maxIter=6;
    paroptim.general.refineOptim=0;
end
function paroptim=ParaDerivUni(e)
    
    paroptim=AreaM2sweep_Lu();
    paroptim.general.startPop='loadshape';
    paroptim.general.specificFillName='.\Active_Build\Input_Variables\Parabola.mat';
    paroptim.optim.CG.diffStepSize=[e,-e]; %[0,2]
    paroptim.optim.CG.minDiffStep=e;
    paroptim.general.maxIter=6;
    paroptim.general.refineOptim=0;
    
end
function paroptim=ParaDerivAllBFGS(e)
    
    paroptim=AreaM2sweepBFGS_Lc();
    paroptim.general.startPop='loadshape';
    paroptim.optim.CG.varActive='all';
    paroptim.general.specificFillName='.\Active_Build\Input_Variables\Parabola.mat';
    paroptim.optim.CG.diffStepSize=[e,-e]; %[0,2]
    paroptim.optim.CG.minDiffStep=e;
    paroptim.general.maxIter=6;
    paroptim.general.refineOptim=0;
end
function paroptim=ParaDerivAll(e)
    
    paroptim=AreaM2sweep_Lc();
    paroptim.general.startPop='loadshape';
    paroptim.optim.CG.varActive='all';
    paroptim.general.specificFillName='.\Active_Build\Input_Variables\Parabola.mat';
    paroptim.optim.CG.diffStepSize=[e,-e]; %[0,2]
    paroptim.optim.CG.minDiffStep=e;
    paroptim.general.maxIter=6;
    paroptim.general.refineOptim=0;
    
end

function paroptim=LongDerivBFGS(e)
    
    paroptim=AreaM2sweepBFGS_Nc();
    %paroptim.general.startPop='loadshape';
    %paroptim.general.specificFillName='.\Active_Build\Input_Variables\Parabola.mat';
    paroptim.optim.CG.diffStepSize=[e,-e]; %[0,2]
    paroptim.optim.CG.minDiffStep=e;
    paroptim.general.maxIter=100;
    paroptim.general.refineOptim=0;
end
function paroptim=LongDeriv(e)
    
    paroptim=AreaM2sweep_Nc();
    %paroptim.general.startPop='loadshape';
    %paroptim.general.specificFillName='.\Active_Build\Input_Variables\Parabola.mat';
    paroptim.optim.CG.diffStepSize=[e,-e]; %[0,2]
    paroptim.optim.CG.minDiffStep=e;
    paroptim.general.maxIter=100;
    paroptim.general.refineOptim=0;
    
end

function paroptim=VolSweep(e)
    
    paroptim=AreaM2sweep_Nc();
    %paroptim.general.startPop='loadshape';
    %paroptim.general.specificFillName='.\Active_Build\Input_Variables\Parabola.mat';
    paroptim.optim.CG.diffStepSize=[1e-6,-1e-6]; %[0,2
    paroptim.constraint.desVarVal={e};
    paroptim.parametrisation.snakes.refine.axisRatio = min(10*e*1.5,1); 
    paroptim.optim.CG.minDiffStep=1e-6;
    paroptim.general.maxIter=100;
    paroptim.general.refineOptim=0;
    
end
% Bp3 Inverse design sweep


function paroptim=desk_refsweep_cv00012()
    [paroptim]=Inverse_CG();
    paroptim=ModifySnakesParam(paroptim,'optimInverseDesign_cv');
    paroptim.obj.invdes.aeroName='0012';
    
    paroptim.optim.CG.lineSearchType='backbisection';
    paroptim.general.refineOptim=[2 1 4; 2 1 4];
    %paroptim.optim.CG.varActive='all';
    paroptim.general.startPop='halfuniformthin';
    paroptim.general.nPop=12;
    paroptim.general.maxIter=4;
    paroptim.general.worker=4;
end


function paroptim=bp3_refsweep_cv00012()
    [paroptim]=Inverse_CG();
    paroptim=ModifySnakesParam(paroptim,'optimInverseDesign_cv');
    paroptim.obj.invdes.aeroName='0012';
    
    paroptim.optim.CG.lineSearchType='backbisection';
    paroptim.general.refineOptim=[2 1 100; 2 1 100];
    %paroptim.optim.CG.varActive='all';
    paroptim.general.startPop='halfuniformthin';
    paroptim.general.nPop=12;
    paroptim.general.maxIter=100;
    paroptim.general.worker=8;
end
function paroptim=bp3_refsweep_cv02212()
    paroptim=bp3_refsweep_cv00012();
    paroptim.obj.invdes.aeroName='2212';
    paroptim=ModifySnakesParam(paroptim,'optimInverseDesign_cv');
    
end
function paroptim=bp3_refsweep_cu00012()
    paroptim=bp3_refsweep_cv00012();
    paroptim.obj.invdes.aeroName='0012';
    paroptim=ModifySnakesParam(paroptim,'optimInverseDesign_cu');
    paroptim.general.refineOptim=[2 2 100; 2 2 100];
    
    
end
function paroptim=bp3_refsweep_cu02212()
    paroptim=bp3_refsweep_cu00012();
    paroptim.obj.invdes.aeroName='2212';
    paroptim=ModifySnakesParam(paroptim,'optimInverseDesign_cu');
    
end
function paroptim=bp3_refsweep_uv02212()
    paroptim=bp3_refsweep_cv00012();
    paroptim.obj.invdes.aeroName='2212';
    paroptim=ModifySnakesParam(paroptim,'optimInverseDesign_uv');
    
end
function paroptim=bp3_refsweep_uv00012()
    paroptim=bp3_refsweep_cv00012();
    paroptim.obj.invdes.aeroName='0012';
    paroptim=ModifySnakesParam(paroptim,'optimInverseDesign_uv');
    
    
end
function paroptim=bp3_refsweep_uu02212()
    paroptim=bp3_refsweep_cu00012();
    paroptim.obj.invdes.aeroName='2212';
    paroptim=ModifySnakesParam(paroptim,'optimInverseDesign_uu');
    
end
function paroptim=bp3_refsweep_uu00012()
    paroptim=bp3_refsweep_cu00012();
    paroptim.obj.invdes.aeroName='0012';
    paroptim=ModifySnakesParam(paroptim,'optimInverseDesign_uu');
end


function paroptim=bp3_refsweep_cv10012()
    [paroptim]=Inverse_CG();
    paroptim=ModifySnakesParam(paroptim,'optimInverseDesign_cv1');
    paroptim.obj.invdes.aeroName='0012';
    
    paroptim.optim.CG.lineSearchType='backbisection';
    paroptim.general.refineOptim=[2 1 100];
    %paroptim.optim.CG.varActive='all';
    paroptim.general.startPop='halfuniformthin';
    paroptim.general.nPop=12;
    paroptim.general.maxIter=100;
    paroptim.general.worker=8;
end
function paroptim=bp3_refsweep_cv12212()
    paroptim=bp3_refsweep_cv10012();
    paroptim.obj.invdes.aeroName='2212';
    paroptim=ModifySnakesParam(paroptim,'optimInverseDesign_cv1');
    
end
function paroptim=bp3_refsweep_cu10012()
    paroptim=bp3_refsweep_cv10012();
    paroptim.obj.invdes.aeroName='0012';
    paroptim=ModifySnakesParam(paroptim,'optimInverseDesign_cu1');
    paroptim.general.refineOptim=[2 2 100];
    paroptim.general.startPop='outerbound';
    
end
function paroptim=bp3_refsweep_cu12212()
    paroptim=bp3_refsweep_cu10012();
    paroptim.obj.invdes.aeroName='2212';
    paroptim=ModifySnakesParam(paroptim,'optimInverseDesign_cu1');
    
end
function paroptim=bp3_refsweep_uv12212()
    paroptim=bp3_refsweep_cv10012();
    paroptim.obj.invdes.aeroName='2212';
    paroptim=ModifySnakesParam(paroptim,'optimInverseDesign_uv1');
    
end
function paroptim=bp3_refsweep_uv10012()
    paroptim=bp3_refsweep_cv10012();
    paroptim.obj.invdes.aeroName='0012';
    paroptim=ModifySnakesParam(paroptim,'optimInverseDesign_uv1');
    
    
end
function paroptim=bp3_refsweep_uu12212()
    paroptim=bp3_refsweep_cu10012();
    paroptim.obj.invdes.aeroName='2212';
    paroptim=ModifySnakesParam(paroptim,'optimInverseDesign_uu1');
    
end
function paroptim=bp3_refsweep_uu10012()
    paroptim=bp3_refsweep_cu10012();
    paroptim.obj.invdes.aeroName='0012';
    paroptim=ModifySnakesParam(paroptim,'optimInverseDesign_uu1');
    
end


function paroptim=bp3_refsweep_cv20012()
    [paroptim]=Inverse_CG();
    paroptim=ModifySnakesParam(paroptim,'optimInverseDesign_cv2');
    paroptim.obj.invdes.aeroName='0012';
    
    paroptim.optim.CG.lineSearchType='backbisection';
    paroptim.general.refineOptim=[0];
    %paroptim.optim.CG.varActive='all';
    paroptim.general.startPop='halfuniformthin';
    paroptim.general.nPop=12;
    paroptim.general.maxIter=100;
    paroptim.general.worker=8;
end
function paroptim=bp3_refsweep_cv22212()
    paroptim=bp3_refsweep_cv20012();
    paroptim.obj.invdes.aeroName='2212';
    paroptim=ModifySnakesParam(paroptim,'optimInverseDesign_cv2');
    
end
function paroptim=bp3_refsweep_cu20012()
    paroptim=bp3_refsweep_cv20012();
    paroptim.obj.invdes.aeroName='0012';
    paroptim=ModifySnakesParam(paroptim,'optimInverseDesign_cu2');
    paroptim.general.refineOptim=[0];
    paroptim.general.startPop='outerbound';
    
end
function paroptim=bp3_refsweep_cu22212()
    paroptim=bp3_refsweep_cu20012();
    paroptim.obj.invdes.aeroName='2212';
    paroptim=ModifySnakesParam(paroptim,'optimInverseDesign_cu2');
    
end
function paroptim=bp3_refsweep_uv22212()
    paroptim=bp3_refsweep_cv20012();
    paroptim.obj.invdes.aeroName='2212';
    paroptim=ModifySnakesParam(paroptim,'optimInverseDesign_uv2');
    
end
function paroptim=bp3_refsweep_uv20012()
    paroptim=bp3_refsweep_cv20012();
    paroptim.obj.invdes.aeroName='0012';
    paroptim=ModifySnakesParam(paroptim,'optimInverseDesign_uv2');
    
    
end
function paroptim=bp3_refsweep_uu22212()
    paroptim=bp3_refsweep_cu20012();
    paroptim.obj.invdes.aeroName='2212';
    paroptim=ModifySnakesParam(paroptim,'optimInverseDesign_uu2');
    
end
function paroptim=bp3_refsweep_uu20012()
    paroptim=bp3_refsweep_cu20012();
    paroptim.obj.invdes.aeroName='0012';
    paroptim=ModifySnakesParam(paroptim,'optimInverseDesign_uu2');
    
end

% Test

function paroptim=AreaM2sweep_ScTest()
    
    paroptim=bp3_AreaM2sweep();
    paroptim=ModifySnakesParam(paroptim,'optimNACA0012Sc');
    paroptim.parametrisation.snakes.refine.LEShrink=true;
    paroptim.general.refineOptim=[2 1 4; 2 1 4];
    paroptim.general.maxIter=0;
    paroptim.general.worker=4;
end
function paroptim=AreaM2sweep_ScTest2()
    
    paroptim=bp3_AreaM2sweep();
    paroptim=ModifySnakesParam(paroptim,'optimNACA0012Sc');
    paroptim.parametrisation.snakes.refine.LEShrink=true;
    paroptim.general.refineOptim=[2 1 4; 2 1 4];
    paroptim.general.maxIter=4;
    paroptim.general.worker=12;
end

%% christmas cases

% Busemann sweep
function [paroptim]=areabusesweep(e)
    
    [paroptim]=MultiTopo_DEhoriz();
    [paroptim]=SumVolumeConstraint(paroptim);
    paroptim=Test_Desktop(paroptim);
    paroptim.general.varOverflow='spill'; % 'truncate' 'spill'
    paroptim.general.optimMethod='DE';
    
    paroptim.general.varOverflow='spill'; % 'truncate' 'spill'
    paroptim.general.nPop=100;
    paroptim.general.maxIter=200;
    paroptim.general.worker=8;
    
    paroptim=ModifySnakesParam(paroptim,'optimSupersonicMultiTopo');
    paroptim.parametrisation.snakes.refine.axisRatio =e*10; % min(10*e*1.5,1); 
    paroptim.constraint.desVarVal={e};
    paroptim.constraint.desVarConstr={'MinValVolFrac'};
end

% Inverse Design refinement
function paroptim=refsweep(gridCase,airfoil,lvl)
    [paroptim]=Inverse_CG();
    
    paroptim.optim.CG.lineSearchType='backbisection';
    
    %paroptim.optim.CG.varActive='all';
    paroptim.general.startPop='halfuniformthin';
    paroptim.general.nPop=12;
    paroptim.general.maxIter=100;
    paroptim.general.worker=8; 
    paroptim.optim.CG.diffStepSize=[1e-6,-1e-6]; %[0,2]
    paroptim.optim.CG.minDiffStep=1e-7;
    paroptim.optim.CG.validVol=0.2;
    
    paroptim.obj.invdes.aeroName=airfoil;
    
    paroptim=ModifySnakesParam(paroptim,['optimInverseDesign']);
    
    switch gridCase(1)
        case 'c'
            paroptim.parametrisation.snakes.refine.gridDistrib='cosX01';
        case 'u'
            paroptim.parametrisation.snakes.refine.gridDistrib='none';
    end
    
    switch gridCase(2)
        case 'v'
            paroptim.parametrisation.snakes.refine.axisRatio=2^lvl;
            paroptim.parametrisation.optiminit.cellLevels=[(6*2^lvl+2),2];
            paroptim.parametrisation.snakes.refine.refineGrid=[4 1];
            paroptim.general.refineOptim=[2 1 100; 2 1 100];
        case 'u'
            paroptim.parametrisation.snakes.refine.axisRatio=1;
            paroptim.parametrisation.optiminit.cellLevels=[(6*2^lvl+2),2^(lvl+1)];
            paroptim.parametrisation.snakes.refine.refineGrid=[4 4];
            paroptim.general.refineOptim=[2 2 100; 2 2 100];
    end
    
    paroptim.general.refineOptim=paroptim.general.refineOptim(lvl+1:end,:);
    
    paroptim.parametrisation.general.passDomBounds=...
        MakeCartesianGridBoundsInactE(paroptim.parametrisation.optiminit.cellLevels);
    
    
    paroptim.initparam=DefaultSnakeInit(paroptim.parametrisation);
end

% NACA0012 refine

function [paroptim]=testRefineNACA()
    paroptim=NACA0012Sweep('c',0,'BFGS');
    paroptim.general.maxIter=0;
    
end

function paroptim=NACA0012Sweep(gridCase,lvl,optimiser)
    
    paroptim=bp3_NACA0012_sweep();
    paroptim=ModifySnakesParam(paroptim,'optimNACA0012');
    paroptim.obj.flow.CFDfolder=[cd,'\Result_Template\CFD_code_Template\transonicfine'];
    paroptim.general.refineOptim=[0];
    paroptim.optim.CG.diffStepSize=[1e-5,-1e-5]; %[0,2]
    paroptim.optim.CG.minDiffStep=1e-7;
    paroptim.optim.CG.validVol=0.2;
    paroptim.parametrisation.general.subdivType='chaikinNaca0012';
    paroptim.general.worker=8;
    
    paroptim=ModifySnakesParam(paroptim,['optimNACA0012']);
    
    switch gridCase(1)
        case 'c'
            paroptim.parametrisation.snakes.refine.gridDistrib='cosX01';
        case 'u'
            paroptim.parametrisation.snakes.refine.gridDistrib='none';
    end
    
    paroptim.optim.CG.stepAlgo=optimiser;
    % only horizontal refinement
    paroptim.parametrisation.snakes.refine.axisRatio=2*2^lvl;
    paroptim.parametrisation.optiminit.cellLevels=[(6*2^lvl+2),2];
    paroptim.parametrisation.snakes.refine.refineGrid=[8 1];
    paroptim.general.refineOptim=[2 1 100; 2 1 100];
        
    
    paroptim.general.refineOptim=paroptim.general.refineOptim(lvl+1:end,:);
    
    paroptim.parametrisation.general.passDomBounds=...
        MakeCartesianGridBoundsInactE(paroptim.parametrisation.optiminit.cellLevels);
    
    
    paroptim.initparam=DefaultSnakeInit(paroptim.parametrisation);
    
end

% Supersonic refine
function paroptim=volsweeprefine(e,gridCase,lvl)
    % e= [0.01:0.12]
    % gridCase={c u}
    % lvl=[0:2]
    
    paroptim=AreaM2sweep_Nc();
    paroptim.general.startPop='loadshape';
    paroptim.general.specificFillName='.\Active_Build\Input_Variables\Parabola.mat';
    paroptim.optim.CG.diffStepSize=[1e-6,-1e-6]; %[0,2
    paroptim.constraint.desVarVal={e};
    paroptim.optim.CG.minDiffStep=1e-7;
    paroptim.general.maxIter=100;
    paroptim.general.worker=8;
    paroptim.general.refineOptim=0;
    
    paroptim=ModifySnakesParam(paroptim,['optimNACA0012Nc']);
    
    switch gridCase(1)
        case 'c'
            paroptim.parametrisation.snakes.refine.gridDistrib='cosX01';
        case 'u'
            paroptim.parametrisation.snakes.refine.gridDistrib='none';
    end
    
    % only horizontal refinement
    %paroptim.parametrisation.snakes.refine.axisRatio=2^lvl;
    paroptim.parametrisation.snakes.refine.axisRatio = min(2*e*1.5,1); 
    paroptim.parametrisation.optiminit.cellLevels=[(8*2^lvl+2),2];
    paroptim.parametrisation.snakes.refine.refineGrid=[8 1];
    paroptim.general.refineOptim=[2 1 100; 2 1 100];
    paroptim.parametrisation.snakes.refine.LEShrink=true;
    paroptim.parametrisation.optiminit.defaultCorner=1e-6;
        
    
    paroptim.general.refineOptim=paroptim.general.refineOptim(lvl+1:end,:);
    
    paroptim.parametrisation.general.passDomBounds=...
        MakeCartesianGridBoundsInactE(paroptim.parametrisation.optiminit.cellLevels);
    
    
    paroptim.initparam=DefaultSnakeInit(paroptim.parametrisation);
end

%% Inverse Design Cases


function paroptim=test_invdes()
    [paroptim]=Inverse_CG();
    
    paroptim=ModifySnakesParam(paroptim,'optimInverseDesign');
    paroptim.obj.invdes.aeroName='0012';
    %paroptim.optim.CG.varActive='all';
    paroptim.general.startPop='halfuniform';
    paroptim.general.nPop=12;
    paroptim.general.maxIter=2;
    paroptim.general.worker=4;
end

function paroptim=desk_0012cos()
    [paroptim]=Inverse_CG();
    
    paroptim.optim.CG.validVol=0.3;
    paroptim.spline.splineCase='inversedesign3';
    paroptim.optim.CG.sensCalc='analytical';
    paroptim.optim.CG.sensAnalyticalType='raw';
    paroptim.optim.CG.diffStepSize=[1e-4,-1e-4];
    paroptim=ModifySnakesParam(paroptim,'optimInverseDesign_Lcos');
    paroptim.obj.invdes.aeroName='0012';
    %paroptim.optim.CG.varActive='all';
    paroptim.general.startPop='halfuniform';
    paroptim.general.nPop=12;
    paroptim.general.maxIter=40;
    paroptim.general.worker=4;
end

function paroptim=desk_4412cos()
    
    paroptim=desk_0012cos();
    paroptim.obj.invdes.aeroName='4412';
    %paroptim.optim.CG.varActive='all';
    paroptim.parametrisation.snakes.refine.axisRatio=2;
    
end 

function paroptim=desk_2212cos()
    
    paroptim=desk_0012cos();
    paroptim.obj.invdes.aeroName='2212';
    %paroptim.optim.CG.varActive='all';
    paroptim.parametrisation.snakes.refine.axisRatio=2;
    
end 

function paroptim=bp3_0012cos()
    [paroptim]=Inverse_CG();
    
    paroptim.optim.CG.diffStepSize=[1e-3,-1e-3];
    paroptim.optim.CG.validVol=0.3;
    paroptim.spline.splineCase='inversedesign';
    
    paroptim=ModifySnakesParam(paroptim,'optimInverseDesign_Lcos');
    paroptim.obj.invdes.aeroName='0012';
    %paroptim.optim.CG.varActive='all';
    paroptim.general.startPop='halfuniform';
    paroptim.general.nPop=12;
    paroptim.general.maxIter=40;
    paroptim.general.worker=8;
    
    paroptim.parametrisation.snakes.refine.axisRatio=2;
end

function paroptim=bp3_4412cos()
    
    paroptim=bp3_0012cos();
    paroptim.obj.invdes.aeroName='4412';
    %paroptim.optim.CG.varActive='all';
    
end 

function paroptim=bp3_2212cos()
    
    paroptim=bp3_0012cos();
    paroptim.obj.invdes.aeroName='2212';
    %paroptim.optim.CG.varActive='all';
    
end 

function paroptim=bp32_0012cos()
    [paroptim]=Inverse_CG();
    
    
    paroptim.general.refineOptim=[2 1; 2 1];
    paroptim.optim.CG.diffStepSize=[1e-3,-1e-3];
    paroptim.optim.CG.validVol=0.3;
    paroptim.spline.splineCase='inversedesign3';
    paroptim.optim.CG.varActive='snaksensiv';
    paroptim=ModifySnakesParam(paroptim,'optimInverseDesign_cosref');
    paroptim.obj.invdes.aeroName='0012';
    %paroptim.optim.CG.varActive='all';
    paroptim.general.startPop='halfuniform';
    paroptim.general.nPop=12;
    paroptim.general.maxIter=10;
    paroptim.general.worker=8;
    paroptim.parametrisation.snakes.refine.axisRatio=1;
    
end

function paroptim=bp32_4412cos()
    
    paroptim=bp32_0012cos();
    paroptim.obj.invdes.aeroName='4412';
    %paroptim.optim.CG.varActive='all';
    
end 

function paroptim=bp32_2212cos()
    
    paroptim=bp32_0012cos();
    paroptim.obj.invdes.aeroName='2212';
    %paroptim.optim.CG.varActive='all';
    
end 

function paroptim=bp3_invdes()
    [paroptim]=Inverse_CG();
    
    paroptim.obj.invdes.aeroName='0012';
    
    paroptim.general.nPop=12;
    paroptim.general.maxIter=16;
    paroptim.general.worker=12;
end

function paroptim=bp3_invdes_L()
     paroptim=bp3_invdes();
     
    paroptim.general.startPop='halfuniform';
     paroptim=ModifySnakesParam(paroptim,'optimInverseDesign_L');
end

% Bulk inverse design


function paroptim=bulkNacaInvDes()
    [paroptim]=Inverse_Bulk();
    
    paroptim=ModifySnakesParam(paroptim,'optimInverseDesign_bulk');
    paroptim.obj.invdes.aeroName='';
    %paroptim.optim.CG.varActive='all';
    paroptim.general.startPop='NACAmulti';
    t1={'1','2','3','4'};
    t2={'1','2','3','4','5','6','7'};
    t3={'06','08','10','12','16'};
    initInterp{1}='0012';
    for ii=1:numel(t1)
        for jj=1:numel(t1)
            for kk=1:numel(t1)
                initInterp{end+1}=[t1{ii},t2{jj},t3{kk}];
            end
        end
    end
    initInterp=initInterp([1:2,50,59]);
    %initInterp={'0012','4412','1108'};
    %initInterp={'4412'};
    paroptim.general.initInterp=initInterp;
    paroptim.general.nPop=12;
    paroptim.general.maxIter=1;
    paroptim.general.worker=4;
    paroptim.spline.resampleSnak=false;
    
end


function paroptim=bulkOtherInvDes()
    [paroptim]=Inverse_Bulk();
    
    paroptim=ModifySnakesParam(paroptim,'optimInverseDesign_bulk');
    paroptim.obj.invdes.aeroName='';
    %paroptim.optim.CG.varActive='all';
    paroptim.general.startPop='AEROmulti';
    initInterp={'RAE 2822 AIRFOIL','ONERA NACA CAMBRE AIRFOIL'};
    paroptim.obj.invdes.aeroClass='lib';
    paroptim.general.initInterp=initInterp;
    paroptim.general.nPop=12;
    paroptim.general.maxIter=1;
    paroptim.general.worker=4;
    paroptim.spline.resampleSnak=false;
    
end
function paroptim=bulkNacaInvDes2()
    [paroptim]=Inverse_Bulk();
    
    paroptim=ModifySnakesParam(paroptim,'optimInverseDesign_bulk2');
    paroptim.obj.invdes.aeroName='';
    %paroptim.optim.CG.varActive='all';
    paroptim.general.startPop='NACAmulti';
    t1={'1','2','3','4'};
    t2={'1','2','3','4','5','6','7'};
    t3={'08','10','12','16','18'};
    initInterp{1}='0012';
    for ii=1:numel(t1)
        for jj=1:numel(t1)
            for kk=1:numel(t1)
                initInterp{end+1}=[t1{ii},t2{jj},t3{kk}];
            end
        end
    end
    paroptim.general.initInterp=initInterp;
    paroptim.general.nPop=12;
    paroptim.general.maxIter=1;
    paroptim.general.worker=4;
end


%% Full Aero Optimisations

% Desktop
function [paroptim]=FullSupersonicOptim_Desktop()
    
    [paroptim]=TestParOptimAero_desktop();
    
    paroptim.general.nPop=32;
    paroptim.general.maxIter=45;
    paroptim.general.worker=8;
    
end

function [paroptim]=FullSupersonicOptimSym_Desktop()
    
    [paroptim]=TestParOptimAero_desktop();
    
    paroptim.general.nPop=32;
    paroptim.general.maxIter=45;
    paroptim.general.worker=8;
    paroptim.general.optimMethod='DEtan';
    paroptim.general.knownOptim=8.82356428E-03;
    
    paroptim.general.symType='horz'; % 'horz'
    paroptim.parametrisation.snakes.refine.axisRatio=0.4;
    
end

function [paroptim]=Full_MultiTopo_M2_D()
    
    [paroptim]=MultiTopo_DEhoriz();
    
    paroptim.general.nPop=56;
    paroptim.general.maxIter=150;
    paroptim.general.worker=8;
    
    
end

function [paroptim]=Full_Aero_CG_20()
    [paroptim]=CG_Aero();
    [paroptim]=LocalVolumeConstraint(paroptim);
    paroptim.parametrisation.snakes.refine.axisRatio=2.5;
    
    paroptim.general.nPop=12;
    paroptim.general.maxIter=40;
    paroptim.general.worker=8;
end

function [paroptim]=BP3_Aero_CG_20_Long()
    
    [paroptim]=CG_Aero();
    [paroptim]=MeanVolumeConstraint_30(paroptim);
    
    paroptim=ModifySnakesParam(paroptim,'optimSupersonic_Long');
    paroptim.parametrisation.snakes.refine.axisRatio=2.8333;
    
    paroptim.general.nPop=12;
    paroptim.general.maxIter=100;
    paroptim.general.worker=12;
end

function [paroptim]=Full_Aero_CG_10()
    [paroptim]=CG_Aero();
    
    paroptim.parametrisation.snakes.refine.axisRatio=1.25;
    
    paroptim.general.nPop=12;
    paroptim.general.maxIter=40;
    paroptim.general.worker=8;
end

function [paroptim]=Full_Aero_CG_05()
    [paroptim]=CG_Aero();
    
    paroptim.parametrisation.snakes.refine.axisRatio=0.6;
    
    paroptim.general.nPop=12;
    paroptim.general.maxIter=40;
    paroptim.general.worker=8;
end

function [paroptim]=Full_Aero_CG_smile()
    
    [paroptim]=Component_CG();
    paroptim.parametrisation.snakes.refine.axisRatio=0.5;
    
    paroptim.general.nPop=12;
    paroptim.general.maxIter=20;
    paroptim.general.worker=4;
end

% bp3
function [paroptim]=FullSupersonicOptim_HPC()
    
    [paroptim]=TestParOptimAero_desktop();
    
    paroptim.general.nPop=48;
    paroptim.general.maxIter=100;
    paroptim.general.worker=12;
    
end

function [paroptim]=FullSupersonicOptimSym_bp3_025()
    
    [paroptim]=TestParOptimAero_desktop();
    
    paroptim=FullOpt_bp3(paroptim);
    paroptim.general.optimMethod='DEtan';
    paroptim.general.knownOptim=8.82356428E-03;
    
    paroptim.general.symType='horz'; % 'horz'
    paroptim.parametrisation.snakes.refine.axisRatio=0.5;
    
end

function [paroptim]=FullSupersonicOptimSym_bp3_05()
    
    [paroptim]=TestParOptimAero_desktop();
    paroptim=FullOpt_bp3(paroptim);
    
    paroptim.general.optimMethod='DEtan';
    paroptim.general.knownOptim=8.82356428E-03;
    
    paroptim.general.symType='horz'; % 'horz'
    paroptim.parametrisation.snakes.refine.axisRatio=1;
    
end

function [paroptim]=FullSupersonicOptimSym_bp3_1()
    
    [paroptim]=TestParOptimAero_desktop();
    
    paroptim=FullOpt_bp3(paroptim);
    paroptim.general.optimMethod='DEtan';
    paroptim.general.knownOptim=8.82356428E-03;
    
    paroptim.general.symType='horz'; % 'horz'
    paroptim.parametrisation.snakes.refine.axisRatio=2;
    
end

function [paroptim]=bp3_Aero_CG_20()
    [paroptim]=CG_Aero();
    
    paroptim.parametrisation.snakes.refine.axisRatio=2.5;
    
    paroptim.general.nPop=12;
    paroptim.general.maxIter=100;
    paroptim.general.worker=12;
end

function [paroptim]=bp3_Aero_CG_10()
    [paroptim]=CG_Aero();
    
    paroptim.parametrisation.snakes.refine.axisRatio=1.25;
    
    paroptim.general.nPop=12;
    paroptim.general.maxIter=100;
    paroptim.general.worker=12;
end

function [paroptim]=bp3_Aero_CG_05()
    [paroptim]=CG_Aero();
    
    paroptim.parametrisation.snakes.refine.axisRatio=0.6;
    
    paroptim.general.nPop=12;
    paroptim.general.maxIter=100;
    paroptim.general.worker=12;
end

function [paroptim]=bp3_MultiTopo_M2_CG()
    
    [paroptim]=MultiTopo_CGhoriz();
    paroptim=Test_Desktop(paroptim);
    %paroptim.general.restartSource={'bp3BorderOptim','DE'};
    paroptim.general.nPop=60;
    paroptim.general.maxIter=40;
    
    paroptim.general.worker=12;
    
end

function [paroptim]=bp3_MultiTopo_M2_DE_Spill()
    
    [paroptim]=MultiTopo_DEhoriz();
    paroptim=Test_Desktop(paroptim);
    paroptim.general.varOverflow='spill'; % 'truncate' 'spill'
    paroptim.general.optimMethod='DE';
    
    paroptim.general.varOverflow='spill'; % 'truncate' 'spill'
    paroptim.general.nPop=100;
    paroptim.general.maxIter=100;
    
    paroptim.general.worker=12;
    
end

function [paroptim]=bp3_MultiTopo_M2_DE_SpillT()
    
    [paroptim]=bp3_MultiTopo_M2_DE_Spill();
    paroptim.general.nPop=100;
    paroptim.general.maxIter=4;
    
    paroptim.general.worker=12;
    
end

function [paroptim]=bp3_MultiTopo_M2_DE_SpillTT()
    
    [paroptim]=bp3_MultiTopo_M2_DE_Spill();
    paroptim.general.nPop=100;
    paroptim.general.maxIter=1;
    
    paroptim.general.worker=12;
    
end

function [paroptim]=bp3_MultiTopo_M2_CG_wide()
    
    [paroptim]=MultiTopo_CGhoriz();
    paroptim.optim.CG.varActive='wideborder';
    paroptim=Test_Desktop(paroptim);
    paroptim.general.restartSource={'bp3BorderOptim','DE'};
    paroptim.general.nPop=60;
    paroptim.general.maxIter=40;
    
    paroptim.general.worker=12;
    
end

function [paroptim]=Desk_CG_missile_in()
    
    [paroptim]=Component_CG();
    paroptim.parametrisation.snakes.refine.axisRatio=0.5;
    paroptim.constraint.initVal={{'.\Active_Build\ConstraintFiles\missile_5b12.png','min'}};
    paroptim.general.nPop=12;
    paroptim.general.maxIter=50;
    paroptim.general.worker=4;
end

function [paroptim]=Desk_DE_missile_horz()
    
    [paroptim]=Component_DE();
    paroptim.parametrisation.snakes.refine.axisRatio=0.5;
    
    paroptim.constraint.initVal={{'.\Active_Build\ConstraintFiles\missile_5b12.png','min'}};
    paroptim.general.nPop=50;
    paroptim.general.maxIter=50;
    paroptim.general.worker=4;
end

function [paroptim]=Desk_DE_smile_horz()
    
    [paroptim]=Component_DE();
    paroptim.parametrisation.snakes.refine.axisRatio=0.5;
    
    paroptim.general.nPop=50;
    paroptim.general.maxIter=50;
    paroptim.general.worker=4;
end

function [paroptim]=Desk_CG_missile_out()
    
    [paroptim]=Component_CG();
    paroptim.parametrisation.snakes.refine.axisRatio=0.5;
    paroptim.constraint.initVal={{'.\Active_Build\ConstraintFiles\missile_5b12.png','min'}};
    paroptim.general.startPop='outerbound';
    
    paroptim.optim.CG.varActive='border'; % 'wideborder'
    
    paroptim.general.varOverflow='spill';
    paroptim.general.nPop=12;
    paroptim.general.maxIter=50;
    paroptim.general.worker=4;
end

function [paroptim]=Desk_CG_inline_out()
    
    [paroptim]=Component_CG();
    paroptim.parametrisation.snakes.refine.axisRatio=0.5;
    paroptim.constraint.initVal={{'.\Active_Build\ConstraintFiles\inline_5b12.png','min'}};
    paroptim.general.startPop='outerbound';
    
    paroptim.optim.CG.varActive='border'; % 'wideborder'
    
    paroptim.general.varOverflow='spill';
    paroptim.general.nPop=12;
    paroptim.general.maxIter=50;
    paroptim.general.worker=4;
end

% Test pop size effect on multi topo optim


function [paroptim]=bp3_MT_pop_100()
    
    [paroptim]=MultiTopo_DEhoriz();
    paroptim=Test_Desktop(paroptim);
    paroptim.general.varOverflow='spill'; % 'truncate' 'spill'
    paroptim.general.optimMethod='DE';
    
    paroptim.general.varOverflow='spill'; % 'truncate' 'spill'
    paroptim.general.nPop=100;
    paroptim.general.maxIter=200;
    
    paroptim.general.worker=12;
    
end

function [paroptim]=bp3_MT_pop_75()
    
    [paroptim]=bp3_MT_pop_100();
    
    paroptim.general.nPop=75;
    
end

function [paroptim]=bp3_MT_pop_25()
    
    [paroptim]=bp3_MT_pop_100();
    
    paroptim.general.nPop=25;
    
end

function [paroptim]=bp3_MT_pop_50()
    
    [paroptim]=bp3_MT_pop_100();
    
    paroptim.general.nPop=50;
    
end

function [paroptim]=bp3_MT_pop_125()
    
    [paroptim]=bp3_MT_pop_100();
    
    paroptim.general.nPop=125;
    
end

function [paroptim]=bp3_MT_pop_150()
    
    [paroptim]=bp3_MT_pop_100();
    
    paroptim.general.nPop=150;
    
end

%% Smoothed mode

% desktop

function [paroptim]=desk_Aero_CG_10_pk()
    
    [paroptim]=CG_Aero();
    
    paroptim.parametrisation.snakes.refine.axisRatio=1;
    paroptim.optim.CG.diffStepSize=[1e-3,-1e-3];
    paroptim.optim.CG.varActive='snaksensiv';
    paroptim.optim.CG.validVol=0.3;
    paroptim.parametrisation.optiminit.modeSmoothType='peaksmooth'; % 'peaksmooth' 'polysmooth';
    paroptim.parametrisation.optiminit.modeSmoothNum=4;
    paroptim.general.nPop=12;
    paroptim.general.maxIter=18;
    paroptim.general.worker=4;
end

function [paroptim]=desk_Aero_CG_10_po()
    [paroptim]=CG_Aero();
    
    paroptim.parametrisation.snakes.refine.axisRatio=1;
    paroptim.optim.CG.diffStepSize=[1e-3,-1e-3];
    paroptim.optim.CG.varActive='snaksensiv';
    paroptim.optim.CG.validVol=0.3;
    paroptim.parametrisation.optiminit.modeSmoothType='polysmooth'; % 'peaksmooth' 'polysmooth';
    paroptim.parametrisation.optiminit.modeSmoothNum=4;
    paroptim.general.nPop=12;
    paroptim.general.maxIter=18;
    paroptim.general.worker=4;
end

function [paroptim]=desk_Aero_CG_10_n()
    [paroptim]=CG_Aero();
    
    paroptim.parametrisation.snakes.refine.axisRatio=1;
    paroptim.optim.CG.diffStepSize=[1e-3,-1e-3];
    paroptim.optim.CG.varActive='all';
    paroptim.optim.CG.validVol=0.3;
    paroptim.parametrisation.optiminit.modeSmoothType='polysmooth'; % 'peaksmooth' 'polysmooth';
    paroptim.parametrisation.optiminit.modeSmoothNum=4;
    paroptim.general.nPop=12;
    paroptim.general.maxIter=18;
    paroptim.general.worker=4;
end

function [paroptim]=desk_Aero_CG_10_popk()
    [paroptim]=CG_Aero();
    
    paroptim.parametrisation.snakes.refine.axisRatio=1;
    paroptim.optim.CG.diffStepSize=[1e-3,-1e-3];
    paroptim.optim.CG.varActive='snaksensiv';
    paroptim.optim.CG.validVol=0.3;
    paroptim.parametrisation.optiminit.modeSmoothType='polypeaksmooth'; % 'peaksmooth' 'polysmooth';
    paroptim.parametrisation.optiminit.modeSmoothNum=4;
    paroptim.general.nPop=12;
    paroptim.general.maxIter=18;
    paroptim.general.worker=4;
end

function [paroptim]=desk_outmis_pk()
    
    [paroptim]=DefaultOptim();
    
    paroptim=ModifySnakesParam(paroptim,'SupersonicComponent');
    
    [paroptim]=LocalVolumeConstraint(paroptim);
    paroptim.constraint.initVal={{'.\Active_Build\ConstraintFiles\missile2_5b12.png','min'}};
    
    
    paroptim=CutCellObjective(paroptim);
    
    [paroptim]=OptimCG(paroptim);
    paroptim.general.startPop='outerbound';
    paroptim.parametrisation.general.subdivType='chaikin';
    paroptim.parametrisation.snakes.refine.axisRatio=1;
    paroptim.general.symType='none'; % 'horz'
    
    paroptim.optim.CG.varActive='snaksensiv';
    
    paroptim.parametrisation.optiminit.modeSmoothType='peaksmooth'; % 'peaksmooth' 'polysmooth';
    paroptim.parametrisation.optiminit.modeSmoothNum=4;
    
    paroptim.parametrisation.snakes.refine.axisRatio=1;
    
    paroptim.general.nPop=12;
    paroptim.general.maxIter=10;
    paroptim.general.worker=4;
end

function [paroptim]=desk_outmis_po()
    
    [paroptim]=DefaultOptim();
    
    paroptim=ModifySnakesParam(paroptim,'SupersonicComponent');
    
    [paroptim]=LocalVolumeConstraint(paroptim);
    paroptim.constraint.initVal={{'.\Active_Build\ConstraintFiles\missile2_5b12.png','min'}};
    
    
    paroptim=CutCellObjective(paroptim);
    
    [paroptim]=OptimCG(paroptim);
    paroptim.general.startPop='outerbound';
    paroptim.parametrisation.general.subdivType='chaikin';
    paroptim.parametrisation.snakes.refine.axisRatio=1;
    paroptim.general.symType='none'; % 'horz'
    
    paroptim.optim.CG.varActive='snaksensiv';
    
    paroptim.parametrisation.optiminit.modeSmoothType='polysmooth'; % 'peaksmooth' 'polysmooth';
    paroptim.parametrisation.optiminit.modeSmoothNum=4;
    
    paroptim.parametrisation.snakes.refine.axisRatio=1;
    
    paroptim.general.nPop=12;
    paroptim.general.maxIter=10;
    paroptim.general.worker=4;
end

function [paroptim]=desk_outmis_popk()
    
    [paroptim]=DefaultOptim();
    
    paroptim=ModifySnakesParam(paroptim,'SupersonicComponent');
    
    [paroptim]=LocalVolumeConstraint(paroptim);
    paroptim.constraint.initVal={{'.\Active_Build\ConstraintFiles\missile2_5b12.png','min'}};
    
    
    paroptim=CutCellObjective(paroptim);
    
    [paroptim]=OptimCG(paroptim);
    paroptim.general.startPop='outerbound';
    paroptim.parametrisation.general.subdivType='chaikin';
    paroptim.parametrisation.snakes.refine.axisRatio=1;
    paroptim.general.symType='none'; % 'horz'
    
    paroptim.optim.CG.varActive='snaksensiv';
    
    paroptim.parametrisation.optiminit.modeSmoothType='polypeaksmooth'; % 'peaksmooth' 'polysmooth';
    paroptim.parametrisation.optiminit.modeSmoothNum=4;
    
    paroptim.parametrisation.snakes.refine.axisRatio=1;
    
    paroptim.general.nPop=12;
    paroptim.general.maxIter=10;
    paroptim.general.worker=4;
end

function [paroptim]=desk_outmis_()
    
    [paroptim]=DefaultOptim();
    
    paroptim=ModifySnakesParam(paroptim,'SupersonicComponent');
    
    [paroptim]=LocalVolumeConstraint(paroptim);
    paroptim.constraint.initVal={{'.\Active_Build\ConstraintFiles\missile2_5b12.png','min'}};
    
    
    paroptim=CutCellObjective(paroptim);
    
    [paroptim]=OptimCG(paroptim);
    paroptim.general.startPop='outerbound';
    paroptim.parametrisation.general.subdivType='chaikin';
    paroptim.parametrisation.snakes.refine.axisRatio=1;
    paroptim.general.symType='none'; % 'horz'
    
    paroptim.optim.CG.varActive='border';
    
    paroptim.parametrisation.optiminit.modeSmoothType='polysmooth'; % 'peaksmooth' 'polysmooth';
    paroptim.parametrisation.optiminit.modeSmoothNum=4;
    
    paroptim.parametrisation.snakes.refine.axisRatio=1;
    
    paroptim.general.nPop=12;
    paroptim.general.maxIter=10;
    paroptim.general.worker=4;
end

function [paroptim]=desk_Aero_CG_20L_pk()
    
    [paroptim]=CG_Aero();
    
    paroptim=ModifySnakesParam(paroptim,'optimSupersonic_Long');
    paroptim.constraint.desVarVal={0.3};
    paroptim.parametrisation.snakes.refine.axisRatio=4/3*2*2;
    paroptim.obj.flow.nMach=4;
    
    paroptim.optim.CG.diffStepSize=[1e-3,-1e-3];
    paroptim.optim.CG.varActive='snaksensiv';
    paroptim.optim.CG.validVol=0.3;
    paroptim.parametrisation.optiminit.modeSmoothType='peaksmooth'; % 'peaksmooth' 'polysmooth';
    paroptim.parametrisation.optiminit.modeSmoothNum=4;
    paroptim.general.nPop=12;
    paroptim.general.maxIter=20;
    paroptim.general.worker=4;
end

function [paroptim]=desk_Aero_CG_20L_po()
    [paroptim]=CG_Aero();
    
    paroptim=ModifySnakesParam(paroptim,'optimSupersonic_Long');
    paroptim.constraint.desVarVal={0.3};
    paroptim.parametrisation.snakes.refine.axisRatio=4/3*2*2;
    paroptim.general.startPop='randuniform';
    paroptim.obj.flow.nMach=4;
    
    paroptim.optim.CG.diffStepSize=[1e-3,-1e-3];
    paroptim.optim.CG.varActive='snaksensiv';
    paroptim.optim.CG.validVol=0.3;
    paroptim.parametrisation.optiminit.modeSmoothType='polysmooth'; % 'peaksmooth' 'polysmooth';
    paroptim.parametrisation.optiminit.modeSmoothNum=6;
    paroptim.general.nPop=12;
    paroptim.general.maxIter=10;
    paroptim.general.worker=4;
end

function [paroptim]=BP3_24to12DV_pk()
    [paroptim]=CG_Aero();
    
    paroptim.parametrisation.snakes.refine.axisRatio=1.25;
    paroptim.general.startPop='specificfill';
    paroptim.general.specificFillName='24DVaverage';
    paroptim.optim.CG.diffStepSize=[1e-4,-1e-4];
    paroptim.optim.CG.varActive='snaksensiv';
    paroptim.optim.CG.validVol=0.3;
    paroptim.parametrisation.optiminit.modeSmoothType='peaksmooth'; % 'peaksmooth' 'polysmooth';
    paroptim.parametrisation.optiminit.modeSmoothNum=4;
    
    paroptim.general.nPop=12;
    paroptim.general.maxIter=6;
    paroptim.general.worker=12;
end

% BP3
function [paroptim]=BP3_SmoothCG_outmis_Aero()
    
    [paroptim]=DefaultOptim();
    
    paroptim=ModifySnakesParam(paroptim,'SupersonicComponent');
    
    [paroptim]=LocalVolumeConstraint(paroptim);
    paroptim.constraint.initVal={{'.\Active_Build\ConstraintFiles\missile2_5b12.png','min'}};
    
    
    paroptim=CutCellObjective(paroptim);
    
    [paroptim]=OptimCG(paroptim);
    paroptim.general.startPop='outerbound';
    paroptim.parametrisation.general.subdivType='chaikin';
    paroptim.parametrisation.snakes.refine.axisRatio=1;
    paroptim.general.symType='none'; % 'horz'
    
    paroptim.optim.CG.varActive='snaksensiv';
    
    paroptim.parametrisation.optiminit.modeSmoothType='peaksmooth'; % 'peaksmooth' 'polysmooth';
    paroptim.parametrisation.optiminit.modeSmoothNum=4;
    
    paroptim.parametrisation.snakes.refine.axisRatio=1;
    
    paroptim.general.nPop=12;
    paroptim.general.maxIter=50;
    paroptim.general.worker=12;
end

    % axRatio 1.25
function [paroptim]=bp3_Aero_CG_10_smooth_peak()
    [paroptim]=CG_Aero();
    
    paroptim.parametrisation.snakes.refine.axisRatio=1.25;
    
    paroptim.optim.CG.diffStepSize=[1e-3,-1e-3];
    paroptim.optim.CG.varActive='snaksensiv';
    paroptim.optim.CG.validVol=0.3;
    paroptim.parametrisation.optiminit.modeSmoothType='peaksmooth'; % 'peaksmooth' 'polysmooth';
    
    paroptim.parametrisation.optiminit.modeSmoothNum=4;
    
    paroptim.general.nPop=12;
    paroptim.general.maxIter=100;
    paroptim.general.worker=12;
end

function [paroptim]=bp3_Aero_CG_10_smooth_poly()
    [paroptim]=CG_Aero();
    
    paroptim.parametrisation.snakes.refine.axisRatio=1.25;
    
    paroptim.optim.CG.diffStepSize=[1e-3,-1e-3];
    paroptim.optim.CG.varActive='snaksensiv';
    paroptim.optim.CG.validVol=0.3;
    paroptim.parametrisation.optiminit.modeSmoothType='polysmooth'; % 'peaksmooth' 'polysmooth';
    paroptim.parametrisation.optiminit.modeSmoothNum=4;
    paroptim.general.nPop=12;
    paroptim.general.maxIter=100;
    paroptim.general.worker=12;
end

function [paroptim]=bp3_Aero_CG_10_smooth_polypeak()
    [paroptim]=CG_Aero();
    
    paroptim.parametrisation.snakes.refine.axisRatio=1.25;
    
    
    paroptim.optim.CG.diffStepSize=[1e-3,-1e-3];
    paroptim.optim.CG.varActive='snaksensiv';
    paroptim.optim.CG.validVol=0.3;
    paroptim.parametrisation.optiminit.modeSmoothType='polypeaksmooth'; % 'peaksmooth' 'polysmooth';
    paroptim.parametrisation.optiminit.modeSmoothNum=4;
    paroptim.general.nPop=12;
    paroptim.general.maxIter=100;
    paroptim.general.worker=12;
end

function [paroptim]=bp3_Aero_CG_10_smooth_none()
    [paroptim]=CG_Aero();
    
    paroptim.parametrisation.snakes.refine.axisRatio=1.25;
    
    
    paroptim.optim.CG.varActive='all';
    paroptim.parametrisation.optiminit.modeSmoothType='peaksmooth'; % 'peaksmooth' 'polysmooth';
    paroptim.parametrisation.optiminit.modeSmoothNum=4;
    paroptim.general.nPop=12;
    paroptim.general.maxIter=100;
    paroptim.general.worker=12;
end


    % axRatio1
function [paroptim]=bp3_Aero_CG_8_smooth_peak()
    [paroptim]=CG_Aero();
    
    paroptim.parametrisation.snakes.refine.axisRatio=1;
    
    paroptim.optim.CG.diffStepSize=[1e-3,-1e-3];
    paroptim.optim.CG.varActive='snaksensiv';
    paroptim.parametrisation.optiminit.modeSmoothType='peaksmooth'; % 'peaksmooth' 'polysmooth';
    paroptim.parametrisation.optiminit.modeSmoothNum=4;
   
    paroptim.optim.CG.validVol=0.3;
    paroptim.general.nPop=12;
    paroptim.general.maxIter=80;
    paroptim.general.worker=12;
end

function [paroptim]=bp3_Aero_CG_8_smooth_poly()
    [paroptim]=bp3_Aero_CG_8_smooth_peak();
    paroptim.optim.CG.varActive='snaksensiv';
    paroptim.parametrisation.optiminit.modeSmoothType='polysmooth'; % 'peaksmooth' 'polysmooth';
    paroptim.parametrisation.optiminit.modeSmoothNum=4;
    
end

function [paroptim]=bp3_Aero_CG_8_smooth_polypeak()
    
    
    [paroptim]=bp3_Aero_CG_8_smooth_peak();
    paroptim.optim.CG.varActive='snaksensiv';
    paroptim.parametrisation.optiminit.modeSmoothType='polypeaksmooth'; % 'peaksmooth' 'polysmooth';
    paroptim.parametrisation.optiminit.modeSmoothNum=4;
  
end

function [paroptim]=bp3_Aero_CG_8_smooth_none()
    
    [paroptim]=bp3_Aero_CG_8_smooth_peak();
    paroptim.optim.CG.varActive='all';
    paroptim.parametrisation.optiminit.modeSmoothType='peaksmooth'; % 'peaksmooth' 'polysmooth';
    paroptim.parametrisation.optiminit.modeSmoothNum=4;
    
end

    % axRatio 1.25L
function [paroptim]=bp3_Aero_CG_10L_pk()
    [paroptim]=CG_Aero();
    paroptim=ModifySnakesParam(paroptim,'optimSupersonic_Long');
    paroptim.parametrisation.snakes.refine.axisRatio=1.25*2;
    paroptim.general.startPop='randuniform';
    paroptim.optim.CG.diffStepSize=[1e-3,-1e-3];
    paroptim.optim.CG.varActive='snaksensiv';
    paroptim.optim.CG.validVol=0.3;
    paroptim.parametrisation.optiminit.modeSmoothType='peaksmooth'; % 'peaksmooth' 'polysmooth';
    paroptim.parametrisation.optiminit.modeSmoothNum=4;
    
    paroptim.general.nPop=12;
    paroptim.general.maxIter=44;
    paroptim.general.worker=12;
end

function [paroptim]=bp3_Aero_CG_10L_po()
    [paroptim]=bp3_Aero_CG_10L_pk();
    
    paroptim.optim.CG.varActive='snaksensiv';
    paroptim.optim.CG.validVol=0.3;
    paroptim.parametrisation.optiminit.modeSmoothType='polysmooth'; % 'peaksmooth' 'polysmooth';
    paroptim.parametrisation.optiminit.modeSmoothNum=4;
    
end

function [paroptim]=bp3_Aero_CG_10L_popk()
    [paroptim]=bp3_Aero_CG_10L_pk();
    
    paroptim.optim.CG.varActive='snaksensiv';
    paroptim.optim.CG.validVol=0.3;
    paroptim.parametrisation.optiminit.modeSmoothType='polypeaksmooth'; % 'peaksmooth' 'polysmooth';
    paroptim.parametrisation.optiminit.modeSmoothNum=4;
    
end

function [paroptim]=bp3_Aero_CG_10L_none()
    [paroptim]=bp3_Aero_CG_10L_pk();
    
    
    paroptim.optim.CG.varActive='all';
    paroptim.parametrisation.optiminit.modeSmoothType='peaksmooth'; % 'peaksmooth' 'polysmooth';
    paroptim.parametrisation.optiminit.modeSmoothNum=4;
    
end

function [paroptim]=bp3_Aero_CG_05_smooth()
    [paroptim]=CG_Aero();
    
    paroptim.parametrisation.snakes.refine.axisRatio=0.6;
    
    paroptim.optim.CG.varActive='snaksensiv';
    paroptim.parametrisation.optiminit.modeSmoothType='peaksmooth'; % 'peaksmooth' 'polysmooth';
    paroptim.parametrisation.optiminit.modeSmoothNum=4;
    paroptim.general.nPop=12;
    paroptim.general.maxIter=100;
    paroptim.general.worker=12;
end

function [paroptim]=bp3_Aero_CG_20_smooth()
    [paroptim]=CG_Aero();
    
    paroptim=ModifySnakesParam(paroptim,'optimSupersonic_Long');
    paroptim.parametrisation.snakes.refine.axisRatio=2*19/14;
    
    paroptim.optim.CG.varActive='snaksensiv';
    paroptim.parametrisation.optiminit.modeSmoothType='peaksmooth'; % 'peaksmooth' 'polysmooth';
    paroptim.parametrisation.optiminit.modeSmoothNum=4;
    
    paroptim.general.nPop=12;
    paroptim.general.maxIter=100;
    paroptim.general.worker=12;
end

function [paroptim]=bp3_outmis_po()
    
    [paroptim]=DefaultOptim();
    
    paroptim=ModifySnakesParam(paroptim,'SupersonicComponent');
    
    [paroptim]=LocalVolumeConstraint(paroptim);
    paroptim.constraint.initVal={{'.\Active_Build\ConstraintFiles\missile2_5b12.png','min'}};
    
    
    paroptim=CutCellObjective(paroptim);
    
    [paroptim]=OptimCG(paroptim);
    paroptim.optim.CG.diffStepSize=[1e-2,-1e-2];
    paroptim.optim.CG.validVol=0.5;
    paroptim.general.startPop='outerbound';
    paroptim.parametrisation.general.subdivType='chaikin';
    paroptim.parametrisation.snakes.refine.axisRatio=1;
    paroptim.general.symType='none'; % 'horz'
    
    paroptim.optim.CG.varActive='snaksensiv';
    
    paroptim.parametrisation.optiminit.modeSmoothType='polysmooth'; % 'peaksmooth' 'polysmooth';
    paroptim.parametrisation.optiminit.modeSmoothNum=4;
    
    paroptim.parametrisation.snakes.refine.axisRatio=1;
    
    paroptim.general.nPop=12;
    paroptim.general.maxIter=30;
    paroptim.general.worker=12;
end

function [paroptim]=bp3_outmis_popk()
    
    [paroptim]=DefaultOptim();
    
    paroptim=ModifySnakesParam(paroptim,'SupersonicComponent');
    
    [paroptim]=LocalVolumeConstraint(paroptim);
    paroptim.constraint.initVal={{'.\Active_Build\ConstraintFiles\missile2_5b12.png','min'}};
    
    
    paroptim=CutCellObjective(paroptim);
    
    [paroptim]=OptimCG(paroptim);
    paroptim.optim.CG.diffStepSize=[1e-2,-1e-2];
    paroptim.optim.CG.validVol=0.5;
    paroptim.general.startPop='outerbound';
    paroptim.parametrisation.general.subdivType='chaikin';
    paroptim.parametrisation.snakes.refine.axisRatio=1;
    paroptim.general.symType='none'; % 'horz'
    
    paroptim.optim.CG.varActive='snaksensiv';
    
    paroptim.parametrisation.optiminit.modeSmoothType='polypeaksmooth'; % 'peaksmooth' 'polysmooth';
    paroptim.parametrisation.optiminit.modeSmoothNum=4;
    
    paroptim.parametrisation.snakes.refine.axisRatio=1;
    
    paroptim.general.nPop=12;
    paroptim.general.maxIter=30;
    paroptim.general.worker=12;
end

function [paroptim]=bp3_Aero_CG_12L_pk()
    [paroptim]=CG_Aero();
    paroptim=ModifySnakesParam(paroptim,'optimSupersonic_Long');
    paroptim.parametrisation.snakes.refine.axisRatio=1.25*1.25*2;
    paroptim.general.startPop='halfuniformsharp';
    paroptim.optim.CG.diffStepSize=[1e-3,-1e-3];
    paroptim.optim.CG.varActive='snaksensiv';
    paroptim.optim.CG.validVol=0.4;
    paroptim.parametrisation.optiminit.modeSmoothType='peaksmooth'; % 'peaksmooth' 'polysmooth';
    paroptim.parametrisation.optiminit.modeSmoothNum=5;
    paroptim.obj.flow.nMach=0.85;
    paroptim.obj.flow.CFDfolder=[cd,'\Result_Template\CFD_code_Template\transonic'];
    paroptim.obj.flow.stoponerror=false;
    paroptim.obj.flow.targConv=-6;
    paroptim.obj.flow.lengthConvTest=100;
    paroptim.obj.flow.restartIter=2000;
    paroptim.obj.flow.maxRestart=10;
     
    
    
    
    
    paroptim.spline.splineCase='aerosnake';
    paroptim.spline.domain='normalizeX';

    
    paroptim.general.nPop=12;
    paroptim.general.maxIter=60;
    paroptim.general.worker=12;
end

function [paroptim]=bp3_Aero_CG_12L_pk_re()
    [paroptim]=bp3_Aero_CG_12L_pk();
    paroptim.general.maxIter=16;
    
end

%% Test Local Optimum

function [paroptim]=LocOptim_12to24_12()
    [paroptim]=CG_Aero();
    paroptim=ModifySnakesParam(paroptim,'optimSupersonic');
    paroptim.parametrisation.snakes.refine.axisRatio=1.25;
    paroptim.general.startPop='initinterp';
    paroptim.general.initInterp={'InterpFunc_24DV.mat'};
    paroptim.optim.CG.diffStepSize=[1e-3,-1e-3];
    paroptim.optim.CG.varActive='snaksensiv';
    paroptim.optim.CG.validVol=0.3;
    paroptim.parametrisation.optiminit.modeSmoothType='peaksmooth'; % 'peaksmooth' 'polysmooth';
    paroptim.parametrisation.optiminit.modeSmoothNum=4;
    
    paroptim.general.nPop=12;
    paroptim.general.maxIter=30;
    paroptim.general.worker=12;
end

function [paroptim]=LocOptim_12to24_13()
    [paroptim]=LocOptim_12to24_12();
    
    nDes=13;
    
    paroptim.parametrisation.optiminit.cellLevels(1)=nDes+2;
    paroptim.parametrisation.general.passDomBounds=...
        MakeCartesianGridBoundsInactE(...
        paroptim.parametrisation.optiminit.cellLevels);
    
    paroptim.initparam=DefaultSnakeInit(paroptim.parametrisation);
    paroptim.parametrisation.snakes.refine.axisRatio=1.25*nDes/12;
    
end

function [paroptim]=LocOptim_12to24_14()
    [paroptim]=LocOptim_12to24_12();
    
    nDes=14;
    
    paroptim.parametrisation.optiminit.cellLevels(1)=nDes+2;
    paroptim.parametrisation.general.passDomBounds=...
        MakeCartesianGridBoundsInactE(...
        paroptim.parametrisation.optiminit.cellLevels);
    
    paroptim.initparam=DefaultSnakeInit(paroptim.parametrisation);
    paroptim.parametrisation.snakes.refine.axisRatio=1.25*nDes/12;
    
end

function [paroptim]=LocOptim_12to24_15()
    [paroptim]=LocOptim_12to24_12();
    
    nDes=15;
    
    paroptim.parametrisation.optiminit.cellLevels(1)=nDes+2;
    paroptim.parametrisation.general.passDomBounds=...
        MakeCartesianGridBoundsInactE(...
        paroptim.parametrisation.optiminit.cellLevels);
    
    paroptim.initparam=DefaultSnakeInit(paroptim.parametrisation);
    paroptim.parametrisation.snakes.refine.axisRatio=1.25*nDes/12;
    
end

function [paroptim]=LocOptim_12to24_16()
    [paroptim]=LocOptim_12to24_12();
    
    nDes=16;
    
    paroptim.parametrisation.optiminit.cellLevels(1)=nDes+2;
    paroptim.parametrisation.general.passDomBounds=...
        MakeCartesianGridBoundsInactE(...
        paroptim.parametrisation.optiminit.cellLevels);
    
    paroptim.initparam=DefaultSnakeInit(paroptim.parametrisation);
    paroptim.parametrisation.snakes.refine.axisRatio=1.25*nDes/12;
    
end

function [paroptim]=LocOptim_12to24_17()
    [paroptim]=LocOptim_12to24_12();
    
    nDes=17;
    
    paroptim.parametrisation.optiminit.cellLevels(1)=nDes+2;
    paroptim.parametrisation.general.passDomBounds=...
        MakeCartesianGridBoundsInactE(...
        paroptim.parametrisation.optiminit.cellLevels);
    
    paroptim.initparam=DefaultSnakeInit(paroptim.parametrisation);
    paroptim.parametrisation.snakes.refine.axisRatio=1.25*nDes/12;
    
end

function [paroptim]=LocOptim_12to24_18()
    [paroptim]=LocOptim_12to24_12();
    
    nDes=18;
    
    paroptim.parametrisation.optiminit.cellLevels(1)=nDes+2;
    paroptim.parametrisation.general.passDomBounds=...
        MakeCartesianGridBoundsInactE(...
        paroptim.parametrisation.optiminit.cellLevels);
    
    paroptim.initparam=DefaultSnakeInit(paroptim.parametrisation);
    paroptim.parametrisation.snakes.refine.axisRatio=1.25*nDes/12;
    
end

function [paroptim]=LocOptim_12to24_19()
    [paroptim]=LocOptim_12to24_12();
    
    nDes=19;
    
    paroptim.parametrisation.optiminit.cellLevels(1)=nDes+2;
    paroptim.parametrisation.general.passDomBounds=...
        MakeCartesianGridBoundsInactE(...
        paroptim.parametrisation.optiminit.cellLevels);
    
    paroptim.initparam=DefaultSnakeInit(paroptim.parametrisation);
    paroptim.parametrisation.snakes.refine.axisRatio=1.25*nDes/12;
    
end

function [paroptim]=LocOptim_12to24_20()
    [paroptim]=LocOptim_12to24_12();
    
    nDes=20;
    
    paroptim.parametrisation.optiminit.cellLevels(1)=nDes+2;
    paroptim.parametrisation.general.passDomBounds=...
        MakeCartesianGridBoundsInactE(...
        paroptim.parametrisation.optiminit.cellLevels);
    
    paroptim.initparam=DefaultSnakeInit(paroptim.parametrisation);
    paroptim.parametrisation.snakes.refine.axisRatio=1.25*nDes/12;
    
end

function [paroptim]=LocOptim_12to24_21()
    [paroptim]=LocOptim_12to24_12();
    
    nDes=21;
    
    paroptim.parametrisation.optiminit.cellLevels(1)=nDes+2;
    paroptim.parametrisation.general.passDomBounds=...
        MakeCartesianGridBoundsInactE(...
        paroptim.parametrisation.optiminit.cellLevels);
    
    paroptim.initparam=DefaultSnakeInit(paroptim.parametrisation);
    paroptim.parametrisation.snakes.refine.axisRatio=1.25*nDes/12;
    
end

function [paroptim]=LocOptim_12to24_22()
    [paroptim]=LocOptim_12to24_12();
    
    nDes=22;
    
    paroptim.parametrisation.optiminit.cellLevels(1)=nDes+2;
    paroptim.parametrisation.general.passDomBounds=...
        MakeCartesianGridBoundsInactE(...
        paroptim.parametrisation.optiminit.cellLevels);
    
    paroptim.initparam=DefaultSnakeInit(paroptim.parametrisation);
    paroptim.parametrisation.snakes.refine.axisRatio=1.25*nDes/12;
    
end

function [paroptim]=LocOptim_12to24_23()
    [paroptim]=LocOptim_12to24_12();
    
    nDes=23;
    
    paroptim.parametrisation.optiminit.cellLevels(1)=nDes+2;
    paroptim.parametrisation.general.passDomBounds=...
        MakeCartesianGridBoundsInactE(...
        paroptim.parametrisation.optiminit.cellLevels);
    
    paroptim.initparam=DefaultSnakeInit(paroptim.parametrisation);
    paroptim.parametrisation.snakes.refine.axisRatio=1.25*nDes/12;
    
end

function [paroptim]=LocOptim_12to24_24()
    [paroptim]=LocOptim_12to24_12();
    
    nDes=24;
    paroptim.general.refineOptim=[2 1; 2 1];
    paroptim.parametrisation.optiminit.cellLevels(1)=nDes+2;
    paroptim.parametrisation.general.passDomBounds=...
        MakeCartesianGridBoundsInactE(...
        paroptim.parametrisation.optiminit.cellLevels);
    
    paroptim.initparam=DefaultSnakeInit(paroptim.parametrisation);
    paroptim.parametrisation.snakes.refine.axisRatio=1.25*nDes/12;
    
end

function [paroptim]=LocOptim_12to24_2_12()
    [paroptim]=CG_Aero();
    paroptim=ModifySnakesParam(paroptim,'optimSupersonic');
    paroptim.parametrisation.snakes.refine.axisRatio=1.25;
    paroptim.general.startPop='initinterp';
    paroptim.general.initInterp={'InterpFunc_13DV.mat'};
    paroptim.optim.CG.diffStepSize=[1e-3,-1e-3];
    paroptim.optim.CG.varActive='snaksensiv';
    paroptim.optim.CG.validVol=0.3;
    paroptim.parametrisation.optiminit.modeSmoothType='peaksmooth'; % 'peaksmooth' 'polysmooth';
    paroptim.parametrisation.optiminit.modeSmoothNum=4;
    
    paroptim.general.nPop=12;
    paroptim.general.maxIter=30;
    paroptim.general.worker=12;
end

function [paroptim]=LocOptim_12to24_2_13()
    [paroptim]=LocOptim_12to24_2_12();
    
    nDes=13;
    
    paroptim.parametrisation.optiminit.cellLevels(1)=nDes+2;
    paroptim.parametrisation.general.passDomBounds=...
        MakeCartesianGridBoundsInactE(...
        paroptim.parametrisation.optiminit.cellLevels);
    
    paroptim.initparam=DefaultSnakeInit(paroptim.parametrisation);
    paroptim.parametrisation.snakes.refine.axisRatio=1.25*nDes/12;
    
end

function [paroptim]=LocOptim_12to24_2_14()
    [paroptim]=LocOptim_12to24_2_12();
    
    nDes=14;
    
    paroptim.parametrisation.optiminit.cellLevels(1)=nDes+2;
    paroptim.parametrisation.general.passDomBounds=...
        MakeCartesianGridBoundsInactE(...
        paroptim.parametrisation.optiminit.cellLevels);
    
    paroptim.initparam=DefaultSnakeInit(paroptim.parametrisation);
    paroptim.parametrisation.snakes.refine.axisRatio=1.25*nDes/12;
    
end

function [paroptim]=LocOptim_12to24_2_15()
    [paroptim]=LocOptim_12to24_2_12();
    
    nDes=15;
    
    paroptim.parametrisation.optiminit.cellLevels(1)=nDes+2;
    paroptim.parametrisation.general.passDomBounds=...
        MakeCartesianGridBoundsInactE(...
        paroptim.parametrisation.optiminit.cellLevels);
    
    paroptim.initparam=DefaultSnakeInit(paroptim.parametrisation);
    paroptim.parametrisation.snakes.refine.axisRatio=1.25*nDes/12;
    
end

function [paroptim]=LocOptim_12to24_2_16()
    [paroptim]=LocOptim_12to24_2_12();
    
    nDes=16;
    
    paroptim.parametrisation.optiminit.cellLevels(1)=nDes+2;
    paroptim.parametrisation.general.passDomBounds=...
        MakeCartesianGridBoundsInactE(...
        paroptim.parametrisation.optiminit.cellLevels);
    
    paroptim.initparam=DefaultSnakeInit(paroptim.parametrisation);
    paroptim.parametrisation.snakes.refine.axisRatio=1.25*nDes/12;
    
end

function [paroptim]=LocOptim_12to24_2_17()
    [paroptim]=LocOptim_12to24_2_12();
    
    nDes=17;
    
    paroptim.parametrisation.optiminit.cellLevels(1)=nDes+2;
    paroptim.parametrisation.general.passDomBounds=...
        MakeCartesianGridBoundsInactE(...
        paroptim.parametrisation.optiminit.cellLevels);
    
    paroptim.initparam=DefaultSnakeInit(paroptim.parametrisation);
    paroptim.parametrisation.snakes.refine.axisRatio=1.25*nDes/12;
    
end

function [paroptim]=LocOptim_12to24_2_18()
    [paroptim]=LocOptim_12to24_2_12();
    
    nDes=18;
    
    paroptim.parametrisation.optiminit.cellLevels(1)=nDes+2;
    paroptim.parametrisation.general.passDomBounds=...
        MakeCartesianGridBoundsInactE(...
        paroptim.parametrisation.optiminit.cellLevels);
    
    paroptim.initparam=DefaultSnakeInit(paroptim.parametrisation);
    paroptim.parametrisation.snakes.refine.axisRatio=1.25*nDes/12;
    
end

function [paroptim]=LocOptim_12to24_2_19()
    [paroptim]=LocOptim_12to24_2_12();
    
    nDes=19;
    
    paroptim.parametrisation.optiminit.cellLevels(1)=nDes+2;
    paroptim.parametrisation.general.passDomBounds=...
        MakeCartesianGridBoundsInactE(...
        paroptim.parametrisation.optiminit.cellLevels);
    
    paroptim.initparam=DefaultSnakeInit(paroptim.parametrisation);
    paroptim.parametrisation.snakes.refine.axisRatio=1.25*nDes/12;
    
end

function [paroptim]=LocOptim_12to24_2_20()
    [paroptim]=LocOptim_12to24_2_12();
    
    nDes=20;
    
    paroptim.parametrisation.optiminit.cellLevels(1)=nDes+2;
    paroptim.parametrisation.general.passDomBounds=...
        MakeCartesianGridBoundsInactE(...
        paroptim.parametrisation.optiminit.cellLevels);
    
    paroptim.initparam=DefaultSnakeInit(paroptim.parametrisation);
    paroptim.parametrisation.snakes.refine.axisRatio=1.25*nDes/12;
    
end

function [paroptim]=LocOptim_12to24_2_21()
    [paroptim]=LocOptim_12to24_2_12();
    
    nDes=21;
    
    paroptim.parametrisation.optiminit.cellLevels(1)=nDes+2;
    paroptim.parametrisation.general.passDomBounds=...
        MakeCartesianGridBoundsInactE(...
        paroptim.parametrisation.optiminit.cellLevels);
    
    paroptim.initparam=DefaultSnakeInit(paroptim.parametrisation);
    paroptim.parametrisation.snakes.refine.axisRatio=1.25*nDes/12;
    
end

function [paroptim]=LocOptim_12to24_2_22()
    [paroptim]=LocOptim_12to24_2_12();
    
    nDes=22;
    
    paroptim.parametrisation.optiminit.cellLevels(1)=nDes+2;
    paroptim.parametrisation.general.passDomBounds=...
        MakeCartesianGridBoundsInactE(...
        paroptim.parametrisation.optiminit.cellLevels);
    
    paroptim.initparam=DefaultSnakeInit(paroptim.parametrisation);
    paroptim.parametrisation.snakes.refine.axisRatio=1.25*nDes/12;
    
end

function [paroptim]=LocOptim_12to24_2_23()
    [paroptim]=LocOptim_12to24_2_12();
    
    nDes=23;
    
    paroptim.parametrisation.optiminit.cellLevels(1)=nDes+2;
    paroptim.parametrisation.general.passDomBounds=...
        MakeCartesianGridBoundsInactE(...
        paroptim.parametrisation.optiminit.cellLevels);
    
    paroptim.initparam=DefaultSnakeInit(paroptim.parametrisation);
    paroptim.parametrisation.snakes.refine.axisRatio=1.25*nDes/12;
    
end

function [paroptim]=LocOptim_12to24_2_24()
    [paroptim]=LocOptim_12to24_2_12();
    
    nDes=24;
    
    paroptim.parametrisation.optiminit.cellLevels(1)=nDes+2;
    paroptim.parametrisation.general.passDomBounds=...
        MakeCartesianGridBoundsInactE(...
        paroptim.parametrisation.optiminit.cellLevels);
    
    paroptim.initparam=DefaultSnakeInit(paroptim.parametrisation);
    paroptim.parametrisation.snakes.refine.axisRatio=1.25*nDes/12;
    
end

%% Refinement

%% Component
% 24/06/2016
function [paroptim]=bp3_Aero_CG_smile_in()
    
    [paroptim]=Component_CG();
    paroptim.parametrisation.snakes.refine.axisRatio=0.5;
    
    paroptim.general.nPop=12;
    paroptim.general.maxIter=70;
    paroptim.general.worker=12;
end

function [paroptim]=bp3_Aero_CG_missile_in()
    
    [paroptim]=Component_CG();
    paroptim.parametrisation.snakes.refine.axisRatio=0.5;
    paroptim.constraint.initVal={{'.\Active_Build\ConstraintFiles\missile_5b12.png','min'}};
    paroptim.general.nPop=12;
    paroptim.general.maxIter=70;
    paroptim.general.worker=12;
end

function [paroptim]=bp3_Aero_CG_smile_out()
    
    [paroptim]=Component_CG();
    paroptim.parametrisation.snakes.refine.axisRatio=0.5;
    paroptim.general.startPop='innerbound';
    paroptim.general.nPop=12;
    paroptim.general.maxIter=70;
    paroptim.general.worker=12;
end

function [paroptim]=bp3_Aero_CG_missile_out()
    
    [paroptim]=Component_CG();
    paroptim.parametrisation.snakes.refine.axisRatio=0.5;
    paroptim.constraint.initVal={{'.\Active_Build\ConstraintFiles\missile_5b12.png','min'}};
    paroptim.general.startPop='outerbound';
    paroptim.general.nPop=12;
    paroptim.general.maxIter=70;
    paroptim.general.worker=12;
end

function [paroptim]=bp3_Aero_DE_smile_horz()
    
    [paroptim]=Component_DE();
    paroptim.parametrisation.snakes.refine.axisRatio=0.5;
    
    paroptim.constraint.initVal={{'.\Active_Build\ConstraintFiles\smile_5b12.png','min'}};
    paroptim.general.nPop=100;
    paroptim.general.maxIter=100;
    paroptim.general.worker=12;
end

function [paroptim]=bp3_Aero_DE_smile()
    
    [paroptim]=Component_DE();
    [paroptim]=OptimDE(paroptim);
    paroptim.parametrisation.snakes.refine.axisRatio=0.5;
    
    paroptim.general.nPop=100;
    paroptim.general.maxIter=100;
    paroptim.general.worker=12;
end

function [paroptim]=bp3_Aero_DE_smile_weak()
    
    [paroptim]=Component_DE();
    [paroptim]=OptimDE_weak(paroptim);
    paroptim.parametrisation.snakes.refine.axisRatio=0.5;
    
    paroptim.general.nPop=100;
    paroptim.general.maxIter=100;
    paroptim.general.worker=12;
end

function [paroptim]=bp3_Aero_DE_missile_horz()
    
    [paroptim]=Component_DE();
    paroptim.parametrisation.snakes.refine.axisRatio=0.5;
    
    paroptim.constraint.initVal={{'.\Active_Build\ConstraintFiles\missile_5b12.png','min'}};
    paroptim.general.nPop=100;
    paroptim.general.maxIter=100;
    paroptim.general.worker=12;
end

function [paroptim]=bp3_Aero_DE_missile()
    
    [paroptim]=Component_DE();
    [paroptim]=OptimDE(paroptim);
    paroptim.parametrisation.snakes.refine.axisRatio=0.5;
    
    paroptim.constraint.initVal={{'.\Active_Build\ConstraintFiles\missile_5b12.png','min'}};
    paroptim.general.nPop=100;
    paroptim.general.maxIter=100;
    paroptim.general.worker=12;
end

function [paroptim]=bp3_Aero_DE_Spill_missile2()
    
    [paroptim]=Component_DE();
    paroptim.parametrisation.snakes.refine.axisRatio=0.5;
    
    paroptim.constraint.initVal={{'.\Active_Build\ConstraintFiles\missile2_5b12.png','min'}};
    paroptim.general.nPop=100;
    paroptim.general.maxIter=100;
    paroptim.general.worker=12;
end

function [paroptim]=bp3_Aero_CG_Restart_missile2()
    
    [paroptim]=Component_CG();
    paroptim.parametrisation.snakes.refine.axisRatio=0.5;
    
    paroptim.constraint.initVal={{'.\Active_Build\ConstraintFiles\missile2_5b12.png','min'}};
    
    paroptim.optim.CG.diffStepSize=[1e-3,-1e-3];
    paroptim.optim.CG.varActive='snaksensiv';
    paroptim.parametrisation.optiminit.modeSmoothType='peaksmooth'; % 'peaksmooth' 'polysmooth';
    paroptim.parametrisation.optiminit.modeSmoothNum=4;
    
    paroptim.general.nPop=12;
    paroptim.general.maxIter=40;
    paroptim.general.worker=12;
end

function [paroptim]=Desk_CG_Re_missile2()
    
    [paroptim]=Component_CG();
    paroptim.parametrisation.snakes.refine.axisRatio=0.5;
    
    paroptim.constraint.initVal={{'.\Active_Build\ConstraintFiles\missile2_5b12.png','min'}};
    
    paroptim.optim.CG.diffStepSize=[1e-3,-1e-3];
    paroptim.optim.CG.varActive='snaksensiv';
    paroptim.parametrisation.optiminit.modeSmoothType='peaksmooth'; % 'peaksmooth' 'polysmooth';
    paroptim.parametrisation.optiminit.modeSmoothNum=4;
    
    paroptim.general.nPop=12;
    paroptim.general.maxIter=40;
    paroptim.general.worker=4;
end

function [paroptim]=bp3_Aero_DE_Spill_smile()
    
    [paroptim]=Component_DE();
    paroptim.parametrisation.snakes.refine.axisRatio=0.5;
    
    paroptim.constraint.initVal={{'.\Active_Build\ConstraintFiles\smile_5b12.png','min'}};
    paroptim.general.nPop=100;
    paroptim.general.maxIter=100;
    paroptim.general.worker=12;
end

function [paroptim]=bp3_Aero_DE_missile_weak()
    
    [paroptim]=Component_DE();
    [paroptim]=OptimDE_weak(paroptim);
    paroptim.parametrisation.snakes.refine.axisRatio=0.5;
    
    paroptim.constraint.initVal={{'.\Active_Build\ConstraintFiles\missile_5b12.png','min'}};
    paroptim.general.nPop=100;
    paroptim.general.maxIter=100;
    paroptim.general.worker=12;
end
% Componetns long
function [paroptim]=bp3_missile2()
    
    [paroptim]=Component_DE();
    paroptim.parametrisation.snakes.refine.axisRatio=0.5;
    
    paroptim.constraint.initVal={{'.\Active_Build\ConstraintFiles\missile2_5b12.png','min'}};
    paroptim.general.nPop=100;
    paroptim.general.maxIter=100;
    paroptim.general.worker=12;
end

function [paroptim]=bp3_smile()
    
    [paroptim]=Component_DE();
    paroptim.parametrisation.snakes.refine.axisRatio=0.5;
    
    paroptim.constraint.initVal={{'.\Active_Build\ConstraintFiles\smile_5b12.png','min'}};
    paroptim.general.nPop=100;
    paroptim.general.maxIter=100;
    paroptim.general.worker=12;
end


% test versions

function [paroptim]=Tbp3_Aero_CG_smile_in()
    
    [paroptim]=bp3_Aero_CG_smile_in();
    
    paroptim.general.nPop=4;
    paroptim.general.maxIter=4;
    paroptim.general.worker=12;
end

function [paroptim]=Tbp3_Aero_CG_missile_in()
    [paroptim]=bp3_Aero_CG_missile_in();
    
    paroptim.general.nPop=4;
    paroptim.general.maxIter=4;
    paroptim.general.worker=12;
end

function [paroptim]=Tbp3_Aero_CG_smile_out()
    [paroptim]=bp3_Aero_CG_smile_out();
    
    paroptim.general.nPop=4;
    paroptim.general.maxIter=4;
    paroptim.general.worker=12;
end

function [paroptim]=Tbp3_Aero_CG_missile_out()
    [paroptim]=bp3_Aero_CG_missile_out();
    
    paroptim.general.nPop=4;
    paroptim.general.maxIter=4;
    paroptim.general.worker=12;
end

function [paroptim]=Tbp3_Aero_DE_smile_horz()
    [paroptim]=bp3_Aero_DE_smile_horz();
    
    paroptim.general.nPop=4;
    paroptim.general.maxIter=4;
    paroptim.general.worker=12;
end

function [paroptim]=Tbp3_Aero_DE_smile()
    [paroptim]=bp3_Aero_DE_smile();
    
    paroptim.general.nPop=4;
    paroptim.general.maxIter=4;
    paroptim.general.worker=12;
end

function [paroptim]=Tbp3_Aero_DE_smile_weak()
    [paroptim]=bp3_Aero_DE_smile_weak();
    paroptim.general.nPop=4;
    paroptim.general.maxIter=4;
    paroptim.general.worker=12;
end

function [paroptim]=Tbp3_Aero_DE_missile_horz()
    [paroptim]=bp3_Aero_DE_missile_horz();
    paroptim.general.nPop=4;
    paroptim.general.maxIter=4;
    paroptim.general.worker=12;
end

function [paroptim]=Tbp3_Aero_DE_missile()
    [paroptim]=bp3_Aero_DE_missile();
    paroptim.general.nPop=4;
    paroptim.general.maxIter=4;
    paroptim.general.worker=4;
end

function [paroptim]=Tbp3_Aero_DE_missile_weak()
    [paroptim]=bp3_Aero_DE_missile_weak();
    paroptim.general.nPop=4;
    paroptim.general.maxIter=4;
    paroptim.general.worker=12;
end

% bp2

function [paroptim]=FullSupersonicOptimSym_bp2_025()
    
    [paroptim]=TestParOptimAero_desktop();
    
    paroptim.general.nPop=32;
    paroptim.general.maxIter=50;
    paroptim.general.worker=8;
    paroptim.general.optimMethod='DEtan';
    paroptim.general.knownOptim=8.82356428E-03;
    
    paroptim.general.symType='horz'; % 'horz'
    paroptim.parametrisation.snakes.refine.axisRatio=0.25;
    
end

function [paroptim]=FullSupersonicOptimSym_bp2_05()
    
    [paroptim]=TestParOptimAero_desktop();
    
    paroptim.general.nPop=32;
    paroptim.general.maxIter=50;
    paroptim.general.worker=8;
    paroptim.general.optimMethod='DEtan';
    paroptim.general.knownOptim=8.82356428E-03;
    
    paroptim.general.symType='horz'; % 'horz'
    paroptim.parametrisation.snakes.refine.axisRatio=0.5;
    
end

function [paroptim]=FullSupersonicOptimSym_bp2_1()
    
    [paroptim]=TestParOptimAero_desktop();
    
    paroptim.general.nPop=32;
    paroptim.general.maxIter=50;
    paroptim.general.worker=8;
    paroptim.general.optimMethod='DEtan';
    paroptim.general.knownOptim=8.82356428E-03;
    
    paroptim.general.symType='horz'; % 'horz'
    paroptim.parametrisation.snakes.refine.axisRatio=1;
    
end

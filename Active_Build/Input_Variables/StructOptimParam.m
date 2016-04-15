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
    paroptim.general.optimCase=caseStr;

end

function structdat=GetStructureData(paroptim)
        
    [structdat]=ExploreStructureTree(paroptim);
    structdat.vardat.names=[structdat.vars(:).name];
    structdat.vardat.varmatch=zeros(size(structdat.vardat.names));
    for ii=1:length(structdat.vars)
        jj=regexp(structdat.vardat.names,structdat.vars(ii).name);
        structdat.vardat.varmatch(jj)=ii;
    end
    
end
%% Deafult Optimisation Inputs

function [paroptim]=DefaultOptim()
    
    paroptim.general=DefaultOptimGeneral();
    [paroptim.optim.DE]=DefaultOptimDE();
    [paroptim.spline]=DefaultOptimSpline();
    paroptim.obj.flow=DefaultCutCell_Flow();
    [paroptim.constraint]=DefaultConstraint();
    paroptim.structdat=GetStructureData(paroptim);
    
end

function [paroptimgeneral]=DefaultOptimGeneral()
    
    paroptimgeneral.optimCase='optimDefault';
    paroptimgeneral.paramCase='optimDefault';
    paroptimgeneral.optimMethod='DE';
    paroptimgeneral.desVarRange=[0,1];
    paroptimgeneral.nDesVar=[0];
    paroptimgeneral.nPop=6;
    paroptimgeneral.startPop='rand';
    paroptimgeneral.maxIter=5;
    paroptimgeneral.worker=6; % Max 4 on this computer
    paroptimgeneral.objectiveName='LengthArea';
    paroptimgeneral.direction='max';
    paroptimgeneral.knownOptim=[0.146088675];
end

function [paroptimDE]=DefaultOptimDE()
    
    paroptimDE.diffAmplification=0.5; %[0,2]
    paroptimDE.xOverRatio=0.5;
end

function paroptimobjflow=DefaultCutCell_Flow()
    
    paroptimobjflow.CFDfolder=[cd,'\Result_Template\CFD_code_Template\supersonic_ogive'];
    paroptimobjflow.stoponerror=false;
    paroptimobjflow.targConv=-6;
    paroptimobjflow.lengthConvTest=100;
    paroptimobjflow.restartIter=1000;
    paroptimobjflow.maxRestart=5;
end

function [paroptimspline]=DefaultOptimSpline()
    
    
    paroptimspline.splineCase='aerosnake';
    paroptimspline.domain='normalizeX';
end

function [paroptimconstraint]=DefaultConstraint()
    
    paroptimconstraint.desVarConstr={' '};
    paroptimconstraint.desVarVal={[0]};
    paroptimconstraint.resConstr={' '};
    paroptimconstraint.resVal={[0]};
    
end

%% Standard Modifications

function [paraminit]=ChangeSnakeInit(paraminit)
    
    paraminit.snakes.step.snakesSteps=20;
    paraminit.general.restart=false;
end

function [paroptim]=VolumeConstraint(paroptim)
    
    paroptim.constraint.desVarConstr={'MeanVolFrac'};
    paroptim.constraint.desVarVal={0.4};
    paroptim.constraint.resConstr={' '};
    paroptim.constraint.resVal={[0]};
    
end
%% Callable functions

function [paroptim]=StandardOptim()
    
    [paroptim]=DefaultOptim();
    
    varExtract={'paramCase'};
    [paramCase]=ExtractVariables(varExtract,paroptim);
    
    paroptim.parametrisation=structInputVar(paramCase);
    paroptim.initparam=ChangeSnakeInit(paroptim.parametrisation);
    
end

function [paroptim]=TestParOptim_desktop()
    
    [paroptim]=DefaultOptim();
    
    varExtract={'paramCase'};
    [paramCase]=ExtractVariables(varExtract,paroptim);
    
    paroptim.parametrisation=structInputVar(paramCase);
    
    paroptim.parametrisation.general.subdivType='chaikin';
    paroptim.general.nPop=6;
    paroptim.general.maxIter=4;
    paroptim.general.worker=6; 
    paroptim.initparam=ChangeSnakeInit(paroptim.parametrisation);
end

function [paroptim]=TestParOptimAero_desktop()
    
    
    [paroptim]=TestParOptim_desktop();
    paroptim.general.paramCase='optimSupersonic';
    varExtract={'paramCase'};
    [paramCase]=ExtractVariables(varExtract,paroptim);
    paroptim.parametrisation=structInputVar(paramCase);
    
    [paroptim]=VolumeConstraint(paroptim);
    
    paroptim.general.objectiveName='CutCellFlow';
    paroptim.general.direction='min';
    paroptim.general.optimMethod='DEtan';
     
    paroptim.parametrisation.general.subdivType='chaikin';
    paroptim.parametrisation.general.refineSteps=3;
end

function [paroptim]=FullSupersonicOptim_Desktop()
    
    [paroptim]=TestParOptimAero_desktop();
    
    paroptim.general.nPop=48;
    paroptim.general.maxIter=100;
    paroptim.general.worker=8; 
    
end



function [paroptim]=TestParOptim_HPC()
    
    [paroptim]=DefaultOptim();
    
    varExtract={'paramCase'};
    [paramCase]=ExtractVariables(varExtract,paroptim);
    
    paroptim.parametrisation=structInputVar(paramCase);
    paroptim.initparam=ChangeSnakeInit(paroptim.parametrisation);
    paroptim.general.nPop=12;
    paroptim.general.maxIter=4;
    paroptim.general.worker=12; 
end

function [paroptim]=HPC_LengthArea()
    
    [paroptim]=DefaultOptim();
    
    varExtract={'paramCase'};
    [paramCase]=ExtractVariables(varExtract,paroptim);
    
    paroptim.parametrisation=structInputVar(paramCase);
    paroptim.initparam=ChangeSnakeInit(paroptim.parametrisation);
    paroptim.general.nPop=48;
    paroptim.general.maxIter=100;
    paroptim.general.worker=12;
    paroptim.general.objectiveName='LengthArea';
    paroptim.general.direction='max';
    paroptim.general.knownOptim=[0.2*(8+pi)/(8+2*pi)];
    
end
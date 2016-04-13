%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%              University of Bristol
%      Department of Aerospace Engineering
%                     2015
%
%          Subdivision of Surfaces
%      for Aerodynamic shape parametrisation
%           Main Function Input Variables
%             Alexandre Payot
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%



function [param]=structInputVar(caseStr)
    % Main function that allows changes
    
    [param]=eval(caseStr);
    param.general.case=caseStr;
    
    [param.structdat]=GetStructureData(param);
    param.structdat.vardat.names=[param.structdat.vars(:).name];
    param.structdat.vardat.varmatch=zeros(size(param.structdat.vardat.names));
    for ii=1:length(param.structdat.vars)
        jj=regexp(param.structdat.vardat.names,param.structdat.vars(ii).name);
        param.structdat.vardat.varmatch(jj)=ii;
    end
end

function structdat=GetStructureData(param)
        
    [structdat]=ExploreStructureTree(param);
    structdat.vardat.names=[structdat.vars(:).name];
    structdat.vardat.varmatch=zeros(size(structdat.vardat.names));
    for ii=1:length(structdat.vars)
        jj=regexp(structdat.vardat.names,structdat.vars(ii).name);
        structdat.vardat.varmatch(jj)=ii;
    end
    
end

%% Default Inputs
function paramgeneral=default_general()
    
    paramgeneral.passDomBounds=[-1,1;-1,1];
    paramgeneral.passGridSteps=3; 
    paramgeneral.refineSteps=1;
    paramgeneral.passPadding=1;
    paramgeneral.typDat='vvlofoil';
    paramgeneral.typeBound='snaxel'; % 'vertex' or 'snaxel'
    paramgeneral.subdivType='chaikin';
    paramgeneral.loadLogical=false;
    paramgeneral.useSnakes=true;
    paramgeneral.execTest=false;
    paramgeneral.boundstr{1}='boundaryis0'; %'boundaryis0'
    paramgeneral.boundstr{2}='solidnotIn0';
    paramgeneral.boundstr{3}='0bound';
    paramgeneral.restart=true;
    
end

function paramplotting=default_plotting()
    
    paramplotting.isCheckRes=true;
    paramplotting.plotInterval=0;
    paramplotting.makeMov=false;
    paramplotting.debugPlot=[0];
    
end

function paramresults=default_results()
    
    paramresults.archiveName='Standard_Execution';
    paramresults.resultRoot=[cd,'\..\results\'];
    paramresults.noteFiles={'CurrentBuild'};
    paramresults.tags={'snakes'};
    
end

function paramsnakesstep=default_snakes_step()
    
    paramsnakesstep.snakesSteps=100;
    paramsnakesstep.mergeTopo=true;
    paramsnakesstep.maxStep=0.4;
    paramsnakesstep.maxDt=0.5;
    paramsnakesstep.convLevel=10^-8;
    paramsnakesstep.arrivalTolerance=2e-1;
    paramsnakesstep.subStep=1;
    paramsnakesstep.snakesMinSteps=5;
    paramsnakesstep.snakData='all';
    paramsnakesstep.snakesConsole=true;
    paramsnakesstep.stepType='indiv'; % 'strict' 'bounded' 'indiv' 'mixed'
    paramsnakesstep.vSwitch=1e-5;
    paramsnakesstep.dtRatio=5;
    
    paramsnakesstep.convCheckRate=5;
    paramsnakesstep.convCheckRange=10;
    paramsnakesstep.convDistance=50;
end

function paramsnakesrefine=default_snakes_refine()
    
    paramsnakesrefine.refineGrid=4;
    paramsnakesrefine.typeRefine='grey';
    
end

function paramsnakesforce=default_snakes_force()
    
    paramsnakesforce.maxForceVel=1;
    paramsnakesforce.bendingVelInfluence=0;
    paramsnakesforce.tensVelInfluence=1;
    paramsnakesforce.maxVelRatio=1;
    paramsnakesforce.dampBase=1;
    paramsnakesforce.dampSides=0;
    paramsnakesforce.vectorMagAveraging=true;
    paramsnakesforce.lengthEpsilon=1e-5;
    paramsnakesforce.typeSmear='length';
    paramsnakesforce.isLast=false;
    
    paramsnakesforce.velType='default';
    paramsnakesforce.vel.Type={'default'};
    paramsnakesforce.vel.ChangeStep=[0];
    paramsnakesforce.vel.ChangeConv=[10];
    paramsnakesforce.vel.ChangeTrigger='none';
end

function paramsnakes=default_snakes()
    
    paramsnakes.step=default_snakes_step();
    paramsnakes.refine= default_snakes_refine();
    paramsnakes.force=default_snakes_force();
end

function paramoptiminit=default_optimInit()
    
    paramoptiminit.cellLevels=[8,2];
    paramoptiminit.refineCellLvl=[0];
    paramoptiminit.defaultfill=0.5;
    paramoptiminit.defaultCorner=1e-3;
    paramoptiminit.corneractive=false;
end

%% Standard Parameter sub sections

function param=OptimConvergence(param)
    

    param.results.noteFiles={'CurrentBuild','OptimSQP'};
    param.results.tags={'snakes','Opimisation','SQP','Profile Length'};
    
end

function param=AvoidLocalOptim(param)
    

    param.snakes.force.lengthEpsilon=1e-5;
    param.snakes.force.typeSmear='length';
    param.snakes.step.arrivalTolerance=10e-2;
    
end

function param=SmoothFEDynamic(param)
    

    param.results.noteFiles={'CurrentBuild','FESmoothing'};
    param.results.tags={'snakes','Dynamic','Curvature'};
    
end

function param=DualOptimSmoothing(param)
    
    
    param.snakes.force.velType='default';
    param.snakes.force.vel.Type={'default','velAreaOnly'};
    param.snakes.force.vel.ChangeStep=[0, 100];
    param.snakes.force.vel.ChangeConv=[10,1e-6];
    param.snakes.force.vel.ChangeTrigger='both'; 
    
    
end

function param=LinOptimSmoothing(param)
    param.snakes.force.vel.ChangeStep=[0];
    param.snakes.force.vel.ChangeConv=[10];
    param.snakes.force.vel.ChangeTrigger='none'; 
    param.snakes.force.velType='velMinLin';
    param.snakes.force.vel.Type={'velMinLin'};
    
end

%% Callable functions

% Default
function [param]=DefaultCase()
    
    % Load defaults
    param.general=default_general();
    param.results=default_results();
    param.plotting=default_plotting();
    param.snakes=default_snakes();
    param.optiminit=default_optimInit();
   

end

% Smoothness tests

function [param]=Snakestestsmooth1()
    
    [param]=DefaultCase();
    param=OptimConvergence(param);
    %param=LinOptimSmoothing(param);
    param=AvoidLocalOptim(param);
    param.general.typDat='testsmooth1';
    
    param.snakes.step.snakesSteps=200;
    param.snakes.refine.refineGrid=8;
    param.snakes.refine.typeRefine='all';
    
end

function [param]=Snakestestsmooth1_2()
    
    [param]=DefaultCase();
    param=OptimConvergence(param);
    param=AvoidLocalOptim(param);
    %param=LinOptimSmoothing(param);
    param.general.typDat='testsmooth1_2';
    
    param.snakes.step.snakesSteps=200;
    param.snakes.refine.refineGrid=6;
    param.snakes.refine.typeRefine='all';
    param.plotting.debugPlot=[0];
    
end

function [param]=Snakestestsmooth2()
    
    
    [param]=DefaultCase();
    param=OptimConvergence(param);
    
    param.general.typDat='testsmooth2';

    param.snakes.step.snakesSteps=100;
    param.snakes.refine.refineGrid=4;
    param.snakes.refine.typeRefine='all';
    
    param.general.passDomBounds=[-1,1;-0.75,0.75];
end

function [param]=Snakestestsmooth3()
    
    
    [param]=DefaultCase();
    param=OptimConvergence(param);
    
    param.general.typDat='testsmooth3';

    param.snakes.step.snakesSteps=200;
    param.snakes.refine.refineGrid=8;
    param.snakes.refine.typeRefine='all';
    
end

function [param]=Snakestestsmooth3_1()
    
    
    [param]=DefaultCase();
    param=OptimConvergence(param);
    
    param.general.typDat='testsmooth3_1';

    param.snakes.step.snakesSteps=200;
    param.snakes.refine.refineGrid=4;
    param.snakes.refine.typeRefine='grey';
    
    param.general.passDomBounds=[-1,1;-0.6,0.6];
    param=AvoidLocalOptim(param);
end

function [param]=Snakestestsmooth4()
    
    
    [param]=DefaultCase();
    param=OptimConvergence(param);
    
    param.general.typDat='testsmooth4';

    param.snakes.step.snakesSteps=200;
    param.snakes.refine.refineGrid=8;
    param.snakes.refine.typeRefine='all';
    
end

% Optimisation Cases

function [param]=optimDefault()
    
    [param]=DefaultCase();
    
    param=OptimConvergence(param);
    param=AvoidLocalOptim(param);
    
    param.general.typDat='optimInit';
    param.general.restart=true;
    param.general.refineSteps=4;
    param.general.subdivType='area';
    
    param.snakes.refine.refineGrid=4;
    param.snakes.refine.typeRefine='all';
    
    param.snakes.step.mergeTopo=false;
    param.snakes.step.snakesSteps=120;
    param.snakes.step.snakData='light';
    param.snakes.step.snakesConsole=false;
    
    param.results.archiveName='Optimisation';
    param.results.resultRoot=[cd,'\..\results\'];
    param.results.noteFiles={'CurrentBuild'};
    param.results.tags={'snakes','optimisation'};
    
    
    sizeRatio=param.optiminit.cellLevels(1,:)+2;
    sizeRatio=sizeRatio(2)/sizeRatio(1);
    param.general.passDomBounds(2,:)=param.general.passDomBounds(2,:)*sizeRatio;
    
    
end

function [param]=optimTest()
    
    [param]=optimDefault();
    
    param.general.restart=false;
    param.snakes.step.snakData='all';
    param.snakes.step.snakesConsole=true;
    param.snakes.step.snakesSteps=50;
    param.general.typDat='optimRand';
    param.results.archiveName='Standard_Execution';
    
    param.general.passDomBounds=[-1,1;-0.5,0.5];
    %param=DualOptimSmoothing(param);
    
end

% Small Shapes

function [param]=SnakesFoilVVSmall()
    
    [param]=DefaultCase();
    param=OptimConvergence(param);
    param=AvoidLocalOptim(param);
    
    param.general.typDat='vvlofoil';
    param.snakes.step.snakesSteps=150;
    param.snakes.refine.refineGrid=4;
    param.snakes.refine.typeRefine='grey';
    param.general.passDomBounds=[-1,1;-0.5,0.5];
    param.general.refineSteps=5;
    param.snakes.step.mergeTopo=false;
    param.snakes.step.convLevel=10^-8;
end

function [param]=Supersonic()
    
    [param]=DefaultCase();
    param=OptimConvergence(param);
    param=AvoidLocalOptim(param);
    
    param.general.typDat='supersonic';
    param.snakes.step.snakesSteps=150;
    param.snakes.refine.refineGrid=4;
    param.snakes.refine.typeRefine='grey';
    param.general.passDomBounds=[-1,1;-0.5,0.5];
    param.general.refineSteps=5;
    param.snakes.step.mergeTopo=false;
    
end

function [param]=Line()
    
    [param]=DefaultCase();
    param=OptimConvergence(param);
    param=AvoidLocalOptim(param);
    
    param.general.typDat='line';
    param.snakes.step.snakesSteps=150;
    param.snakes.refine.refineGrid=4;
    param.snakes.refine.typeRefine='grey';
    param.general.passDomBounds=[-1,1;-0.5,0.5];
    param.general.refineSteps=5;
    param.snakes.step.mergeTopo=false;
    
    param.general.subdivType='area';
end

function [param]=SnakesFoilVVSmall2()
    
    [param]=DefaultCase();
    param=OptimConvergence(param);
    param=AvoidLocalOptim(param);
    
    param.general.typDat='vvlofoil2';
    param.snakes.step.snakesSteps=150;
    param.snakes.refine.refineGrid=4;
    param.snakes.refine.typeRefine='all';
    param.general.passDomBounds=[-1,1;-0.6,0.6];
    param.general.refineSteps=5;
    param.snakes.step.mergeTopo=true;
    param.plotting.debugPlot=[1];
    
end

function [param]=SnakesFoilVVSmall3()
    
    [param]=DefaultCase();
    param=OptimConvergence(param);
    param=AvoidLocalOptim(param);
    
    param.general.typDat='vvlofoil3';
    param.snakes.step.snakesSteps=150;
    param.snakes.refine.refineGrid=4;
    param.snakes.refine.typeRefine='all';
    param.general.passDomBounds=[-1,1;-0.4,0.4];
    param.general.refineSteps=5;
    param.snakes.step.mergeTopo=true;
    param.plotting.debugPlot=[1];
    
end

function [param]=SnakesFoilVVSmall4()
    
    [param]=DefaultCase();
    param=OptimConvergence(param);
    param=AvoidLocalOptim(param);
    
    param.general.typDat='vvlofoil4';
    param.snakes.step.snakesSteps=150;
    param.snakes.refine.refineGrid=4;
    param.snakes.refine.typeRefine='all';
    param.general.passDomBounds=[-1,1;-0.25,0.25];
    param.general.refineSteps=5;
    param.snakes.step.mergeTopo=true;
    param.plotting.debugPlot=[0];
    
end

function [param]=WeirdShape()
    
    % Load defaults
    param.general=default_general();
    param.results=default_results();
    param.plotting=default_plotting();
    param.snakes=default_snakes();
    
    param=OptimConvergence(param);
    param=DualOptimSmoothing(param);
    
    param.snakes.refine.refineGrid=4;
    param.general.typDat='low5shape';
    param.snakes.refine.typeRefine='grey';
    param.snakes.step.snakesSteps=100;
   

end

function [param]=WeirdShape2()
    
    % Load defaults
    param.general=default_general();
    param.results=default_results();
    param.plotting=default_plotting();
    param.snakes=default_snakes();
    
    %param=OptimConvergence(param);
    param=DualOptimSmoothing(param);
    
    param=AvoidLocalOptim(param);
    param.snakes.refine.refineGrid=4;
    param.general.typDat='low5shape2';
    param.snakes.refine.typeRefine='all';
    param.snakes.step.snakesSteps=500;
   
    

end

%% Not updated yet
function [param]=SnakesFoilVSmall()
    
    passDomBounds=[-1,1;-1,1];
     % number of steps in design domain
    passGridSteps=3; 
     % number of refining steps
    refineSteps=2;
    passPadding=1;
    
    typDat='vlofoil';
    typeBound='snaxel'; % 'vertex' or 'snaxel'
    loadLogical=false;
    useSnakes=true;
    isCheckRes=true;
    snakesPlotInterval=0;
    snakesSteps=50;
    refineGrid=2;
    typeRefine='grey';
    execTest=false;
    makeMov=false;
    
end

function [param]=SnakesdoubleBody()
    
    passDomBounds=[-1,1;-1,1];
     % number of steps in design domain
    passGridSteps=3; 
     % number of refining steps
    refineSteps=0;
    passPadding=1;
    
    typDat='doubleBody';
    typeBound='snaxel'; % 'vertex' or 'snaxel'
    loadLogical=false;
    useSnakes=true;
    isCheckRes=true;
    snakesPlotInterval=1;
    snakesSteps=100;
    refineGrid=2;
    typeRefine='grey';
    execTest=false;
    makeMov=false;
    
end

function [param]=RandSmall()
    
    passDomBounds=[-1,1;-1,1];
     % number of steps in design domain
    passGridSteps=10; 
     % number of refining steps
    refineSteps=4;
    passPadding=1;
    
    typDat='rand';
    typeBound='vertex'; % 'vertex' or 'snaxel'
    loadLogical=false;
    useSnakes=false;
    isCheckRes=true;
    snakesPlotInterval=1;
    snakesSteps=2;
    refineGrid=2;
    typeRefine='all';
    execTest=false;
    makeMov=false;
    
end

function [param]=SquareSnakes()
    
    passDomBounds=[-1,1;-1,1];
     % number of steps in design domain
    passGridSteps=10; 
     % number of refining steps
    refineSteps=4;
    passPadding=1;
    
    typDat='square';
    typeBound='snaxel'; % 'vertex' or 'snaxel'
    loadLogical=false;
    useSnakes=true;
    isCheckRes=true;
    snakesPlotInterval=5;
    snakesSteps=50;
    
    refineGrid=2;
    typeRefine='all';
    execTest=false;
    makeMov=false;
    
    
end

function [param]=UniNoPlot()
    
    passDomBounds=[-1.2,1.2;-1.2,1.2];
     % number of steps in design domain
    passGridSteps=10; 
     % number of refining steps
    refineSteps=2;
    passPadding=1;
    
    typDat='uni';
    typeBound='snaxel'; % 'vertex' or 'snaxel'
    loadLogical=false;
    useSnakes=true;
    isCheckRes=false;
    snakesPlotInterval=0;
    snakesSteps=20;
    
    refineGrid=2;
    typeRefine='grey';
    execTest=false;
    makeMov=false;
    
    
end

function [param]=UniLo2InitPlot() %#ok<*DEFNU>
    
    passDomBounds=[-2,2;-1,1];
     % number of steps in design domain
    passGridSteps=10; 
     % number of refining steps
    refineSteps=2;
    passPadding=1;
    
    typDat='uniLo2';
    typeBound='snaxel'; % 'vertex' or 'snaxel'
    loadLogical=false;
    useSnakes=true;
    isCheckRes=true;
    snakesPlotInterval=0;
    snakesSteps=20;
    
    refineGrid=2;
    typeRefine='grey';
    execTest=false;
    makeMov=false;
    
    
end

function [param]=WedgInitPlot() %#ok<*DEFNU>
    
    passDomBounds=[-2,2;-1,1];
     % number of steps in design domain
    passGridSteps=10; 
     % number of refining steps
    refineSteps=2;
    passPadding=1;
    
    typDat='wedg';
    typeBound='snaxel'; % 'vertex' or 'snaxel'
    loadLogical=false;
    useSnakes=true;
    isCheckRes=true;
    snakesPlotInterval=0;
    snakesSteps=40;
    
    refineGrid=2;
    typeRefine='grey';
    execTest=false;
    makeMov=false;
    
    
end

function [param]=WedgPlot() %#ok<*DEFNU>
    
    passDomBounds=[-2,2;-1,1];
     % number of steps in design domain
    passGridSteps=10; 
     % number of refining steps
    refineSteps=2;
    passPadding=1;
    
    typDat='wedg';
    typeBound='snaxel'; % 'vertex' or 'snaxel'
    loadLogical=false;
    useSnakes=true;
    isCheckRes=true;
    snakesPlotInterval=4;
    snakesSteps=20;
    
    refineGrid=2;
    typeRefine='grey';
    execTest=false;
    makeMov=false;
    
    
end

function [param]=WeirdShapeNoPlot()
    
    passDomBounds=[-1.2,1.2;-1.2,1.2];
     % number of steps in design domain
    passGridSteps=10; 
     % number of refining steps
    refineSteps=3;
    passPadding=1;
    
    typDat='low5shape';
    typeBound='snaxel'; % 'vertex' or 'snaxel'
    loadLogical=false;
    useSnakes=true;
    isCheckRes=true;
    snakesPlotInterval=0;
    snakesSteps=100;
    
    refineGrid=4;
    typeRefine='all';
    execTest=false;
    makeMov=false;
    
end







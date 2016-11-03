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
    paramplotting.checkSensitivities=false;
    
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
    paramsnakesstep.maxStep=0.3;
    paramsnakesstep.maxDt=0.5;
    paramsnakesstep.convLevel=10^-8;
    paramsnakesstep.arrivalTolerance=1e-2;
    paramsnakesstep.subStep=1;
    paramsnakesstep.snakesMinSteps=5;
    paramsnakesstep.snakData='all';
    paramsnakesstep.snakesConsole=true;
    paramsnakesstep.stepType='indiv'; % 'strict' 'bounded' 'indiv' 'mixed'
    paramsnakesstep.vSwitch=1e-15;
    paramsnakesstep.dtRatio=5;
    paramsnakesstep.snaxInitPos=1e-4;
    
    paramsnakesstep.convCheckRate=100;
    paramsnakesstep.convCheckRange=15;
    paramsnakesstep.convDistance=500;
    
    paramsnakesstep.fillLooseStep=5;
    paramsnakesstep.fillLooseCut=1e-3;
    
end

function paramsnakesrefine=default_snakes_refine()
    
    paramsnakesrefine.refineGrid=4;
    paramsnakesrefine.gridDistrib='none'; % 'cosX1'
    paramsnakesrefine.typeRefine='grey';
    paramsnakesrefine.LEShrink=0;
    paramsnakesrefine.TEShrink=0.008175297200000/2;
    paramsnakesrefine.edgeFinish='none'; % 'none' 'sharpen' 'shrink'
    paramsnakesrefine.resampleSnak=false;
    paramsnakesrefine.axisRatio=1; % y/x ratio
    paramsnakesrefine.typeCorner='global'; % y/x ratio
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
    paramsnakesforce.distEpsilon=1e-5;
    paramsnakesforce.dirEpsilon=1e-5;
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
    paramoptiminit.defaultCorner=1e-4;
    paramoptiminit.corneractive=false;
    paramoptiminit.modeSmoothType='peaksmooth';
    paramoptiminit.modeSmoothNum=6;
end

%% Standard Parameter sub sections

function param=OptimConvergence(param)
    

    param.results.noteFiles={'CurrentBuild','OptimSQP'};
    param.results.tags={'snakes','Opimisation','SQP','Profile Length'};
    
end

function param=AvoidLocalOptim(param)
    

    param.snakes.force.lengthEpsilon=1e-5;
    param.snakes.force.typeSmear='length';
    param.snakes.step.arrivalTolerance=1e-1;
    param.snakes.step.snaxInitPos=1*param.snakes.force.lengthEpsilon;
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

%%  tests cases

function [param]=Snakestestsmooth1()
    
    [param]=DefaultCase();
    param=OptimConvergence(param);
    %param=LinOptimSmoothing(param);
    param.general.typDat='testsmooth1';
    
    param.snakes.force.lengthEpsilon=1e-5;
    param.snakes.force.typeSmear='length';
    param.snakes.step.arrivalTolerance=10e-3;
    param.snakes.step.snaxInitPos=10*param.snakes.force.lengthEpsilon;
    
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
    param=AvoidLocalOptim(param);
    param.general.typDat='testsmooth3_1';

    param.snakes.step.snakesSteps=200;
    param.snakes.refine.refineGrid=8;
    param.snakes.refine.typeRefine='all';
    
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

function [param]=Snakestestsmooth4_ref()
    
    
    [param]=DefaultCase();
    param=OptimConvergence(param);
    
    param.general.typDat='testsmooth4';

    param.snakes.step.snakesSteps=50;
    param.snakes.refine.refineGrid=8;
    param.snakes.refine.typeRefine='all';
    
end

function [param]=ModeAnalysis()
    
    [param]=DefaultCase();
    param=OptimConvergence(param);
    param=AvoidLocalOptim(param);
    
    param.general.refineSteps=4;
    param.general.subdivType='chaikin';
    param.general.typDat='basisanalysis';
    param.snakes.step.snakesSteps=150;
    param.snakes.refine.refineGrid=8;
    param.snakes.refine.typeRefine='all';
    param.general.passDomBounds=[-1,1;-0.4444,0.4444];
    param.general.refineSteps=4;
    param.snakes.step.mergeTopo=true;
     param.snakes.refine.axisRatio=1;
    param.snakes.refine.TEShrink=false;
    param.snakes.refine.LEShrink=false;
    param.snakes.refine.edgeFinish='none';
    param.plotting.checkSensitivities=true;
    
end

function [param]=ModeAnalysis2body()
    
    [param]=DefaultCase();
    param=OptimConvergence(param);
    param=AvoidLocalOptim(param);
    
    param.general.refineSteps=4;
    param.general.subdivType='chaikin';
    param.general.typDat='2bodymode';
    param.snakes.step.snakesSteps=150;
    param.snakes.refine.refineGrid=4;
    param.snakes.refine.typeRefine='all';
    param.general.passDomBounds=[-1,1;-0.4444,0.4444];
    param.general.refineSteps=4;
    param.snakes.step.mergeTopo=true;
     param.snakes.refine.axisRatio=1;
    param.snakes.refine.TEShrink=false;
    param.snakes.refine.LEShrink=false;
    param.snakes.refine.edgeFinish='none';
    param.plotting.checkSensitivities=true;
    
end

function [param]=FlatPlate()
    
    [param]=DefaultCase();
    param=OptimConvergence(param);
    param=AvoidLocalOptim(param);
    
    param.general.refineSteps=4;
    param.general.subdivType='chaikin';
    param.general.typDat='flatplate';
    param.snakes.step.snakesSteps=150;
    param.snakes.refine.refineGrid=4;
    param.snakes.refine.typeRefine='all';
    param.general.passDomBounds=[-1,1;-0.1818,0.1818];
    param.general.refineSteps=2;
    param.snakes.step.mergeTopo=true;
     param.snakes.refine.axisRatio=1;
    param.snakes.refine.TEShrink=false;
    param.snakes.refine.LEShrink=false;
    param.snakes.refine.edgeFinish='none';
    param.plotting.checkSensitivities=false;
    
end

function [param]=Missile2In()
    
    [param]=DefaultCase();
    param=OptimConvergence(param);
    param=AvoidLocalOptim(param);
    
    param.general.refineSteps=4;
    param.general.subdivType='chaikin';
    param.general.typDat='missile2_input';
    param.snakes.step.snakesSteps=150;
    param.snakes.refine.refineGrid=4;
    param.snakes.refine.typeRefine='all';
    param.general.passDomBounds=[-1,1;-2/3,2/3];
    param.general.refineSteps=2;
    param.snakes.step.mergeTopo=true;
     param.snakes.refine.axisRatio=1;
    param.snakes.refine.TEShrink=false;
    param.snakes.refine.LEShrink=false;
    param.snakes.refine.edgeFinish='none';
    param.plotting.checkSensitivities=false;
    
end

%% Optimisation Cases

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
    param.snakes.refine.LEShrink=0.008175297200000/2;
    param.snakes.refine.TEShrink=0.008175297200000/2;
    param.snakes.refine.resampleSnak=true;
    param.snakes.refine.axisRatio=0.25;
    
    param.snakes.step.mergeTopo=false;
    param.snakes.step.snakesSteps=50;
    param.snakes.step.snakData='light';
    param.snakes.step.snakesConsole=false;
    
    param.results.archiveName='Optimisation';
    param.results.resultRoot=[cd,'\..\results\'];
    param.results.noteFiles={'CurrentBuild'};
    param.results.tags={'snakes','optimisation'};
    
    param.snakes.refine.axisRatio=0.25;
    
    sizeRatio=param.optiminit.cellLevels(1,:)+2;
    sizeRatio=sizeRatio(2)/sizeRatio(1);
    param.general.passDomBounds(2,:)=param.general.passDomBounds(2,:)*sizeRatio;
    
    
end

function [param]=OptimNoRestart()
    
    [param]=optimDefault();
    param.general.restart=false;
    param.snakes.step.snakesSteps=3;
end

function [param]=optimTest()
    
    [param]=optimDefault();
    
    param.general.restart=false;
    param.snakes.step.snakData='all';
    param.snakes.step.snakesConsole=true;
    param.snakes.step.snakesSteps=150;
    param.general.typDat='optimRand';
    param.results.archiveName='Standard_Execution';
    
    param.general.passDomBounds=[-1,1;-0.5,0.5];
    %param=DualOptimSmoothing(param);
    
end

function [param]=optimSupersonic()
   
    [param]=DefaultCase();
    
    param=OptimConvergence(param);
    param=AvoidLocalOptim(param);
    
    param.general.typDat='optimInit';
    param.general.restart=true;
    param.general.refineSteps=4;
    param.general.subdivType='area';
    
    param.snakes.refine.refineGrid=4;
    param.snakes.refine.typeRefine='all';
    param.snakes.refine.LEShrink=true;
    param.snakes.refine.TEShrink=true;
    param.snakes.refine.resampleSnak=true;
    param.snakes.refine.axisRatio=0.25;
    
    param.snakes.step.mergeTopo=true;
    param.snakes.step.snakesSteps=100;
    param.snakes.step.snakData='light';
    param.snakes.step.snakesConsole=false;
    
    param.results.archiveName='Optimisation';
    param.results.resultRoot=[cd,'\..\results\'];
    param.results.noteFiles={'CurrentBuild'};
    param.results.tags={'snakes','optimisation'};
    
    param.snakes.refine.axisRatio=0.25;
    
    param.optiminit.cellLevels=[12,2];
    param.general.passDomBounds=MakeCartesianGridBoundsInactE(param.optiminit.cellLevels);
    
    param.general.subdivType='chaikin';
    param.snakes.refine.TEShrink=true;
    param.snakes.refine.LEShrink=true;
    param.snakes.refine.edgeFinish='sharpen';
    param.snakes.refine.resampleSnak=false;
    param.general.refineSteps=3;
    param.optiminit.corneractive=false;
end

function [param]=optimSupersonicCos()
    [param]=optimSupersonic();
    
    param.snakes.refine.axisRatio=1;
    param.snakes.refine.gridDistrib='cosX1';
    param.optiminit.cellLevels=[12,2];
    param.general.passDomBounds=MakeCartesianGridBoundsInactE(param.optiminit.cellLevels);
end

function [param]=optimSupersonic_Long()
   
    [param]=optimSupersonic();
    
    param.snakes.refine.axisRatio=2.8333;
    param.optiminit.cellLevels=[26,2];
    param.general.passDomBounds=MakeCartesianGridBoundsInactE(param.optiminit.cellLevels);
    
end

function [param]=optimSupersonicMultiTopo()
    
    [param]=DefaultCase();
    
    param=OptimConvergence(param);
    param=AvoidLocalOptim(param);
    
    param.general.typDat='optimInit';
    param.general.restart=false;
    param.general.refineSteps=3;
    param.general.subdivType='chaikin';
    
    param.snakes.refine.refineGrid=4;
    param.snakes.refine.typeRefine='all';
    param.snakes.refine.LEShrink=true;
    param.snakes.refine.TEShrink=true;
    param.snakes.refine.edgeFinish='sharpen';
    param.snakes.refine.resampleSnak=false;
    param.snakes.refine.axisRatio=0.17;
    
    param.snakes.step.mergeTopo=true;
    param.snakes.step.snakesSteps=150;
    param.snakes.step.snakData='light';
    param.snakes.step.snakesConsole=false;
    param.snakes.step.maxStep=0.2;
    param.snakes.step.maxDt=0.5;
    param.snakes.step.fillLooseStep=5;
    param.snakes.step.fillLooseCut=1e-3;
    param.snakes.step.stepType='indiv';
     
    param.results.archiveName='Optimisation';
    param.results.resultRoot=[cd,'\..\results\'];
    param.results.noteFiles={'CurrentBuild'};
    param.results.tags={'snakes','optimisation'};
    
    param.optiminit.corneractive=true;
    param.optiminit.cellLevels=[6,10];
    param.general.passDomBounds=MakeCartesianGridBounds(param.optiminit.cellLevels);
    
end

function [param]=SupersonicComponent()
    
    [param]=DefaultCase();
    
    param=OptimConvergence(param);
    param=AvoidLocalOptim(param);
    
    param.general.typDat='optimInit';
    param.general.restart=false;
    param.general.refineSteps=3;
    param.general.subdivType='chaikin';
    
    
    param.optiminit.modeSmoothType='polypeaksmooth'; % 'peaksmooth';
    
    
    param.snakes.refine.refineGrid=4;
    param.snakes.refine.typeRefine='all';
    param.snakes.refine.LEShrink=true;
    param.snakes.refine.TEShrink=true;
    param.snakes.refine.typeCorner='global';
    param.snakes.refine.edgeFinish='sharpen';
    param.snakes.refine.resampleSnak=false;
    param.snakes.refine.axisRatio=1;
    
    param.snakes.step.mergeTopo=true;
    param.snakes.step.snakesSteps=150;
    param.snakes.step.snakData='light';
    param.snakes.step.snakesConsole=false;
    param.snakes.step.maxStep=0.2;
    param.snakes.step.maxDt=0.5;
    param.snakes.step.fillLooseStep=5;
    param.snakes.step.fillLooseCut=1e-3;
    param.snakes.step.stepType='indiv';
     
    param.results.archiveName='Optimisation';
    param.results.resultRoot=[cd,'\..\results\'];
    param.results.noteFiles={'CurrentBuild'};
    param.results.tags={'snakes','optimisation'};
    
    param.optiminit.corneractive=true;
    param.optiminit.cellLevels=[1,1];
    param.general.passDomBounds=MakeCartesianGridBounds(param.optiminit.cellLevels);
    
end

function [param]=TestInitOptim()
    
    [param]=DefaultCase();
    
    param=OptimConvergence(param);
    param=AvoidLocalOptim(param);
    
    param.general.typDat='optimInit';
    param.general.restart=false;
    param.general.refineSteps=3;
    param.general.subdivType='chaikin';
    
    param.snakes.refine.refineGrid=4;
    param.snakes.refine.typeRefine='all';
    param.snakes.refine.LEShrink=true;
    param.snakes.refine.TEShrink=true;
    param.snakes.refine.typeCorner='global';
    param.snakes.refine.edgeFinish='sharpen';
    param.snakes.refine.resampleSnak=false;
    param.snakes.refine.axisRatio=1;
    
    param.snakes.step.mergeTopo=true;
    param.snakes.step.snakesSteps=20;
    param.snakes.step.snakData='light';
    param.snakes.step.snakesConsole=false;
    param.snakes.step.maxStep=0.2;
    param.snakes.step.maxDt=0.5;
    param.snakes.step.fillLooseStep=5;
    param.snakes.step.fillLooseCut=1e-3;
    param.snakes.step.stepType='indiv';
     
    param.results.archiveName='Optimisation';
    param.results.resultRoot=[cd,'\..\results\'];
    param.results.noteFiles={'CurrentBuild'};
    param.results.tags={'snakes','optimisation'};
    
    param.optiminit.corneractive=true;
    param.optiminit.cellLevels=[1,1];
    param.general.passDomBounds=MakeCartesianGridBounds(param.optiminit.cellLevels);
    
end

function [param]=optimInverseDesign()
   
    [param]=DefaultCase();
    
    param=OptimConvergence(param);
    param=AvoidLocalOptim(param);
    
    param.general.typDat='optimInit';
    param.general.restart=true;
    param.general.refineSteps=2;
    
    param.snakes.refine.refineGrid=[4 4];
    param.snakes.refine.typeRefine='all';
    
    param.snakes.step.mergeTopo=true;
    param.snakes.step.snakesSteps=200;
    param.snakes.step.snakData='light';
    param.snakes.step.snakesConsole=false;
    
    param.results.archiveName='Optimisation';
    param.results.resultRoot=[cd,'\..\results\'];
    param.results.noteFiles={'CurrentBuild'};
    param.results.tags={'snakes','optimisation'};
    
    param.snakes.refine.axisRatio=1;
    
    param.optiminit.cellLevels=[10,2];
    param.optiminit.defaultCorner=1e-4;
    param.general.passDomBounds=MakeCartesianGridBoundsInactE(param.optiminit.cellLevels);
    
    param.general.subdivType='chaikin';
    param.snakes.refine.TEShrink=true;
    param.snakes.refine.LEShrink=false;
    param.snakes.refine.edgeFinish='sharpen';
    param.snakes.refine.resampleSnak=true;
    param.general.refineSteps=2;
    param.optiminit.corneractive=false;
end

function [param]=optimInverseDesign_L()
   
    [param]=optimInverseDesign();
    
    param.snakes.refine.axisRatio=1.5;
    
    param.optiminit.cellLevels=[22,2];
    param.general.passDomBounds=MakeCartesianGridBoundsInactE(param.optiminit.cellLevels);
end

function [param]=optimInverseDesign_Lcos()
   
    [param]=optimInverseDesign();
    
    param.snakes.refine.axisRatio=1.5;
    
    param.snakes.refine.refineGrid=[5 1];
    
    param.snakes.refine.gridDistrib='cosX1';
    param.optiminit.cellLevels=[22,2];
    param.general.passDomBounds=MakeCartesianGridBoundsInactE(param.optiminit.cellLevels);
end

function [param]=optimInverseDesign_cosref()
   
    [param]=optimInverseDesign();
    
    param.snakes.refine.axisRatio=1;
    
    param.snakes.refine.refineGrid=[4 4];
    
    param.snakes.refine.gridDistrib='cosX1';
    param.optiminit.cellLevels=[10,2];
    param.general.passDomBounds=MakeCartesianGridBoundsInactE(param.optiminit.cellLevels);
end

function [param]=optimNACA0012()
   
    [param]=DefaultCase();
    
    param=OptimConvergence(param);
    param=AvoidLocalOptim(param);
    
    param.general.typDat='optimInit';
    param.general.restart=true;
    param.general.refineSteps=4;
    param.general.subdivType='area';
    
    param.snakes.refine.refineGrid=4;
    param.snakes.refine.typeRefine='all';
    param.snakes.refine.LEShrink=true;
    param.snakes.refine.TEShrink=true;
    param.snakes.refine.resampleSnak=true;
    
    
    param.snakes.step.mergeTopo=true;
    param.snakes.step.snakesSteps=100;
    param.snakes.step.snakData='light';
    param.snakes.step.snakesConsole=false;
    
    param.results.archiveName='Optimisation';
    param.results.resultRoot=[cd,'\..\results\'];
    param.results.noteFiles={'CurrentBuild'};
    param.results.tags={'snakes','optimisation'};
    
    param.snakes.refine.axisRatio=1;
    param.snakes.refine.gridDistrib='cosX01';
    param.optiminit.cellLevels=[12,2];
    param.general.passDomBounds=MakeCartesianGridBoundsInactE(param.optiminit.cellLevels);
    
    param.general.subdivType='chaikin';
    param.snakes.refine.TEShrink=true;
    param.snakes.refine.LEShrink=false;
    param.snakes.refine.edgeFinish='sharpen';
    param.snakes.refine.resampleSnak=false;
    param.general.refineSteps=3;
    param.optiminit.corneractive=false;
end

function [param]=optimNACA0012L()
    [param]=optimNACA0012();
    
    param.snakes.refine.axisRatio=2;
    param.snakes.refine.gridDistrib='cosX1';
    param.optiminit.cellLevels=[22,2];
    param.general.passDomBounds=MakeCartesianGridBoundsInactE(param.optiminit.cellLevels);
end


% NACA0012

function [param]=optimNACA0012Sc()
    [param]=optimNACA0012();
    
    param.snakes.refine.axisRatio=0.5;
    param.snakes.refine.gridDistrib='cosX01';
    param.optiminit.cellLevels=[7,2];
    param.general.passDomBounds=MakeCartesianGridBoundsInactE(param.optiminit.cellLevels);
end
function [param]=optimNACA0012Su()
    [param]=optimNACA0012();
    
    param.snakes.refine.axisRatio=0.5;
    param.snakes.refine.gridDistrib='none';
    param.optiminit.cellLevels=[7,2];
    param.general.passDomBounds=MakeCartesianGridBoundsInactE(param.optiminit.cellLevels);
end
function [param]=optimNACA0012Nc()
    [param]=optimNACA0012();
    
    param.snakes.refine.axisRatio=1;
    param.snakes.refine.gridDistrib='cosX01';
    param.optiminit.cellLevels=[12,2];
    param.general.passDomBounds=MakeCartesianGridBoundsInactE(param.optiminit.cellLevels);
end
function [param]=optimNACA0012Nu()
    [param]=optimNACA0012();
    
    param.snakes.refine.axisRatio=1;
    param.snakes.refine.gridDistrib='none';
    param.optiminit.cellLevels=[12,2];
    param.general.passDomBounds=MakeCartesianGridBoundsInactE(param.optiminit.cellLevels);
end
function [param]=optimNACA0012Lc()
    [param]=optimNACA0012();
    
    param.snakes.refine.axisRatio=2;
    param.snakes.refine.gridDistrib='cosX01';
    param.optiminit.cellLevels=[22,2];
    param.general.passDomBounds=MakeCartesianGridBoundsInactE(param.optiminit.cellLevels);
end
function [param]=optimNACA0012Lu()
    [param]=optimNACA0012();
    
    param.snakes.refine.axisRatio=2;
    param.snakes.refine.gridDistrib='none';
    param.optiminit.cellLevels=[22, 2];
    param.general.passDomBounds=MakeCartesianGridBoundsInactE(param.optiminit.cellLevels);
end
function [param]=optimNACA0012Vc()
    [param]=optimNACA0012();
    
    param.snakes.refine.axisRatio=4;
    param.snakes.refine.gridDistrib='cosX01';
    param.optiminit.cellLevels=[42,2];
    param.general.passDomBounds=MakeCartesianGridBoundsInactE(param.optiminit.cellLevels);
end
function [param]=optimNACA0012Vu()
    [param]=optimNACA0012();
    
    param.snakes.refine.axisRatio=4;
    param.snakes.refine.gridDistrib='none';
    param.optiminit.cellLevels=[42, 2];
    param.general.passDomBounds=MakeCartesianGridBoundsInactE(param.optiminit.cellLevels);
end

% BP3 inverse design testing procedure

function [param]=optimInverseDesign_cv()
   
    [param]=optimInverseDesign();
    
    param.snakes.refine.axisRatio=1;
    
    param.snakes.refine.refineGrid=[4 1];
    
    param.snakes.refine.gridDistrib='cosX1';
    param.optiminit.cellLevels=[10,2];
    param.general.passDomBounds=MakeCartesianGridBoundsInactE(param.optiminit.cellLevels);
end
function [param]=optimInverseDesign_cu()
   
    [param]=optimInverseDesign();
    
    param.snakes.refine.axisRatio=1;
    
    param.snakes.refine.refineGrid=[4 4];
    
    param.snakes.refine.gridDistrib='cosX1';
    param.optiminit.cellLevels=[10,2];
    param.general.passDomBounds=MakeCartesianGridBoundsInactE(param.optiminit.cellLevels);
end
function [param]=optimInverseDesign_uv()
   
    [param]=optimInverseDesign();
    
    param.snakes.refine.axisRatio=1;
    
    param.snakes.refine.refineGrid=[4 1];
    
    param.snakes.refine.gridDistrib='none';
    param.optiminit.cellLevels=[10,2];
    param.general.passDomBounds=MakeCartesianGridBoundsInactE(param.optiminit.cellLevels);
end
function [param]=optimInverseDesign_uu()
   
    [param]=optimInverseDesign();
    
    param.snakes.refine.axisRatio=1;
    
    param.snakes.refine.refineGrid=[4 4];
    
    param.snakes.refine.gridDistrib='none';
    param.optiminit.cellLevels=[10,2];
    param.general.passDomBounds=MakeCartesianGridBoundsInactE(param.optiminit.cellLevels);
end

% Bp3 inverse design reference grids

function [param]=optimInverseDesign_cv1()
   
    [param]=optimInverseDesign();
    
    param.snakes.refine.axisRatio=2;
    
    param.snakes.refine.refineGrid=[4 1];
    
    param.snakes.refine.gridDistrib='cosX1';
    param.optiminit.cellLevels=[18,2];
    param.general.passDomBounds=MakeCartesianGridBoundsInactE(param.optiminit.cellLevels);
end
function [param]=optimInverseDesign_cu1()
   
    [param]=optimInverseDesign();
    
    param.snakes.refine.axisRatio=1;
    
    param.snakes.refine.refineGrid=[4 4];
    
    param.snakes.refine.gridDistrib='cosX1';
    param.optiminit.cellLevels=[18,4];
    param.general.passDomBounds=MakeCartesianGridBoundsInactE(param.optiminit.cellLevels);
end
function [param]=optimInverseDesign_uv1()
   
    [param]=optimInverseDesign();
    
    param.snakes.refine.axisRatio=2;
    
    param.snakes.refine.refineGrid=[4 1];
    
    param.snakes.refine.gridDistrib='none';
    param.optiminit.cellLevels=[18,2];
    param.general.passDomBounds=MakeCartesianGridBoundsInactE(param.optiminit.cellLevels);
end
function [param]=optimInverseDesign_uu1()
   
    [param]=optimInverseDesign();
    
    param.snakes.refine.axisRatio=1;
    
    param.snakes.refine.refineGrid=[4 4];
    
    param.snakes.refine.gridDistrib='none';
    param.optiminit.cellLevels=[18,4];
    param.general.passDomBounds=MakeCartesianGridBoundsInactE(param.optiminit.cellLevels);
end
function [param]=optimInverseDesign_cv2()
   
    [param]=optimInverseDesign();
    
    param.snakes.refine.axisRatio=4;
    
    param.snakes.refine.refineGrid=[4 1];
    
    param.snakes.refine.gridDistrib='cosX1';
    param.optiminit.cellLevels=[34,2];
    param.general.passDomBounds=MakeCartesianGridBoundsInactE(param.optiminit.cellLevels);
end
function [param]=optimInverseDesign_cu2()
   
    [param]=optimInverseDesign();
    
    param.snakes.refine.axisRatio=1;
    
    param.snakes.refine.refineGrid=[4 4];
    
    param.snakes.refine.gridDistrib='cosX1';
    param.optiminit.cellLevels=[34,8];
    param.general.passDomBounds=MakeCartesianGridBoundsInactE(param.optiminit.cellLevels);
end
function [param]=optimInverseDesign_uv2()
   
    [param]=optimInverseDesign();
    
    param.snakes.refine.axisRatio=4;
    
    param.snakes.refine.refineGrid=[4 1];
    
    param.snakes.refine.gridDistrib='none';
    param.optiminit.cellLevels=[34,2];
    param.general.passDomBounds=MakeCartesianGridBoundsInactE(param.optiminit.cellLevels);
end
function [param]=optimInverseDesign_uu2()
   
    [param]=optimInverseDesign();
    
    param.snakes.refine.axisRatio=1;
    
    param.snakes.refine.refineGrid=[4 4];
    
    param.snakes.refine.gridDistrib='none';
    param.optiminit.cellLevels=[34,8];
    param.general.passDomBounds=MakeCartesianGridBoundsInactE(param.optiminit.cellLevels);
end

%% Surrogate modelling Cases

function [param]=surrogateDefault()
    
    [param]=DefaultCase();
    
    param=OptimConvergence(param);
    param=AvoidLocalOptim(param);
    
    param.general.typDat='optimInit';
    param.general.restart=true;
    param.general.refineSteps=4;
    param.general.subdivType='chaikin';
    
    param.snakes.refine.refineGrid=6;
    param.snakes.refine.typeRefine='all';
    param.snakes.refine.LEShrink=true;
    param.snakes.refine.TEShrink=true;
    param.snakes.refine.resampleSnak=false;
    
    param.snakes.step.mergeTopo=true;
    param.snakes.step.snakesSteps=100;
    param.snakes.step.snakData='light';
    param.snakes.step.snakesConsole=false;
    
    param.results.archiveName='Optimisation';
    param.results.resultRoot=[cd,'\..\results\'];
    param.results.noteFiles={'CurrentBuild'};
    param.results.tags={'snakes','optimisation'};
    
    param.snakes.refine.axisRatio=0.5;
    
    param.optiminit.refineCellLvl=[0];
    param.optiminit.defaultfill=0.5;
    param.optiminit.defaultCorner=0.25;
    param.optiminit.corneractive=true;
    param.optiminit.cellLevels=[7,1];
    sizeRatio=param.optiminit.cellLevels(1,:)+2;
    sizeRatio=sizeRatio(2)/sizeRatio(1);
    param.general.passDomBounds(2,:)=param.general.passDomBounds(2,:)*sizeRatio;
    
    
end

%% Small Shapes

function [param]=SnakesFoilVVSmall()
    
    [param]=DefaultCase();
    param=OptimConvergence(param);
    param=AvoidLocalOptim(param);
    
    param.general.typDat='vvlofoil_extra';
    param.snakes.step.snakesSteps=150;
    
    
    param.snakes.refine.refineGrid=4;
    param.snakes.refine.typeRefine='all';
    param.general.passDomBounds=[-1,1;-0.4,0.4];
    param.general.refineSteps=5;
    param.snakes.step.mergeTopo=true;
    param.snakes.step.convLevel=10^-8;
    param.snakes.refine.TEShrink=true;
    param.snakes.refine.LEShrink=false;
    param.snakes.refine.edgeFinish='sharpen';
end

function [param]=ManualRefine()
    
    [param]=DefaultCase();
    param=OptimConvergence(param);
    param=AvoidLocalOptim(param);
    
    param.general.typDat='ManualRefine2';
    param.snakes.step.snakesSteps=150;
    
    param.snakes.force.typeSmear='dir';
    param.snakes.step.arrivalTolerance=1e-5;
    param.snakes.force.lengthEpsilon=1e-5;
    param.snakes.force.distEpsilon=0;
    param.snakes.force.dirEpsilon=1e-5;
    param.snakes.step.snaxInitPos=1e-5;
    
    param.general.loadLogical=true;
    param.snakes.refine.refineGrid=[4 4];
    param.snakes.refine.typeRefine='all';
    param.general.passDomBounds=[-1,1;-0.4,0.4];
    param.general.refineSteps=5;
    param.snakes.step.mergeTopo=true;
    param.snakes.step.convLevel=10^-8;
    param.snakes.refine.TEShrink=true;
    param.snakes.refine.LEShrink=false;
    param.snakes.refine.edgeFinish='sharpen';
end

function [param]=testRefinement()
    
    [param]=DefaultCase();
    param=OptimConvergence(param);
    param=AvoidLocalOptim(param);
    
    param.general.typDat='testsmooth1';
    param.snakes.step.snakesSteps=100;
    param.snakes.refine.refineGrid=[3 2];
    param.snakes.refine.typeRefine='all';
    param.general.passDomBounds=[-1,1;-1,1];
    param.general.refineSteps=5;
    param.snakes.step.mergeTopo=true;
    param.snakes.step.convLevel=10^-8;
    param.snakes.refine.TEShrink=true;
    param.snakes.refine.LEShrink=false;
    param.snakes.refine.edgeFinish='none';
end

function [param]=SnakNaca0012()
    
    [param]=DefaultCase();
    param=OptimConvergence(param);
    param=AvoidLocalOptim(param);
    param.general.typDat='naca0012cos';
    param.snakes.step.snakesSteps=150;
    param.snakes.refine.refineGrid=[4 4];
    param.snakes.refine.gridDistrib='cosX1';
    param.snakes.refine.typeRefine='all';
    param.general.passDomBounds=[-1.4,1.4;-0.4,0.4];
    param.general.refineSteps=5;
    param.snakes.step.mergeTopo=false;
    param.snakes.step.convLevel=10^-8;
    param.snakes.refine.TEShrink=true;
    param.snakes.refine.LEShrink=false;
    param.snakes.refine.edgeFinish='sharpen';
end

function [param]=SnakesFoilVVSmall_ref()
    
    [param]=DefaultCase();
    param=OptimConvergence(param);
    param=AvoidLocalOptim(param);
    
    param.general.typDat='vvlofoil';
    param.snakes.step.snakesSteps=150;
    param.snakes.refine.refineGrid=8;
    param.snakes.refine.typeRefine='all';
    param.general.passDomBounds=[-1,1;-0.4,0.4];
    param.general.refineSteps=5;
    param.snakes.step.mergeTopo=false;
    param.snakes.step.convLevel=10^-8;
    param.snakes.refine.TEShrink=true;
    param.snakes.refine.LEShrink=false;
    param.snakes.refine.edgeFinish='sharpen';
end

function [param]=Supersonic()
    
    [param]=DefaultCase();
    param=OptimConvergence(param);
    param=AvoidLocalOptim(param);
    
    param.general.typDat='supersonic';
    param.snakes.step.snakesSteps=150;
    param.snakes.refine.refineGrid=4;
    param.snakes.refine.typeRefine='grey';
    param.general.passDomBounds=[-1,1;-0.4,0.4];
    param.general.refineSteps=4;
    param.snakes.step.mergeTopo=false;
     param.snakes.refine.axisRatio=0.25;
    param.snakes.refine.TEShrink=true;
    param.snakes.refine.LEShrink=true;
    param.snakes.refine.edgeFinish='sharpen';
end

function [param]=Klunker()
    
    [param]=DefaultCase();
    param=OptimConvergence(param);
    param=AvoidLocalOptim(param);
    
    param.general.typDat='klunker';
    param.snakes.step.snakesSteps=150;
    param.snakes.refine.refineGrid=4;
    param.snakes.refine.typeRefine='all';
    param.general.refineSteps=4;
    param.snakes.step.mergeTopo=false;
     param.snakes.refine.axisRatio=2.2;
    param.snakes.refine.TEShrink=true;
    param.snakes.refine.LEShrink=true;
    param.snakes.refine.edgeFinish='sharpen';
    
    param.optiminit.cellLevels=[13,2];
    sizeRatio=param.optiminit.cellLevels(1,:)+2;
    sizeRatio=sizeRatio(2)/sizeRatio(1);
    param.general.passDomBounds(2,:)=param.general.passDomBounds(2,:)*sizeRatio;
    
end

function [param]=Klunker2()
    
    [param]=DefaultCase();
    param=OptimConvergence(param);
    param=AvoidLocalOptim(param);
    
    param.general.typDat='klunker2';
    param.snakes.step.snakesSteps=150;
    param.snakes.refine.refineGrid=2;
    param.snakes.refine.typeRefine='all';
    param.general.refineSteps=4;
    param.snakes.step.mergeTopo=false;
     param.snakes.refine.axisRatio=1.0;
    param.snakes.refine.TEShrink=true;
    param.snakes.refine.LEShrink=true;
    param.snakes.refine.edgeFinish='sharpen';
    
    param.optiminit.cellLevels=[8,4];
    sizeRatio=param.optiminit.cellLevels(1,:)+2;
    sizeRatio=sizeRatio(2)/sizeRatio(1);
    param.general.passDomBounds(2,:)=param.general.passDomBounds(2,:)*sizeRatio;
    
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
    param.snakes.refine.refineGrid=[3 2];
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
    param.snakes.refine.typeRefine='all';
    param.snakes.step.snakesSteps=100;
   
    param.general.boundstr{1}='boundaryis1'; %'boundaryis0'
    param.general.boundstr{2}='solidisIn1';
    param.general.boundstr{3}='1bound';

end

function [param]=InitTest_in()
    
    % Load defaults
    [param]=InitTest_out();
   
    param.general.boundstr{1}='boundaryis1'; %'boundaryis0'
    param.general.boundstr{2}='solidisIn1';
    param.general.boundstr{3}='1bound';

end

function [param]=InitTest_out()
    
    % Load defaults
    param.general=default_general();
    param.results=default_results();
    param.plotting=default_plotting();
    param.snakes=default_snakes();
    
    param=OptimConvergence(param);
    param=AvoidLocalOptim(param);
    
    param.snakes.refine.refineGrid=4;
    param.general.typDat='inittest';
    param.snakes.refine.typeRefine='all';
    param.snakes.step.snakesSteps=100;
   
%     param.general.boundstr{1}='boundaryis1'; %'boundaryis0'
%     param.general.boundstr{2}='solidisIn1';
%     param.general.boundstr{3}='1bound';

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

function [param]=Square()
    
    [param]=DefaultCase();
    param=OptimConvergence(param);
    %param=LinOptimSmoothing(param);
    param=AvoidLocalOptim(param);
    param.general.typDat='sqare';
    
    param.snakes.step.snakesSteps=150;
    param.snakes.refine.refineGrid=4;
    param.snakes.refine.typeRefine='all';
    
    param.snakes.step.maxStep=0.4;
    param.snakes.step.maxDt=0.5;
    
    param.snakes.step.mergeTopo=true;
    
    param.snakes.step.fillLooseStep=40;
    param.snakes.step.fillLooseCut=0.90;
end

function [param]=HalfWedge()
    
    [param]=DefaultCase();
    param=OptimConvergence(param);
    %param=LinOptimSmoothing(param);
    param=AvoidLocalOptim(param);
    param.general.typDat='halfwedge';
    
    param.snakes.step.snakesSteps=250;
    param.snakes.refine.refineGrid=4;
    param.snakes.refine.typeRefine='all';
    
    param.snakes.step.maxStep=0.4;
    param.snakes.step.maxDt=0.5;
    
    param.snakes.step.mergeTopo=true;
    
    param.snakes.step.fillLooseStep=40;
    param.snakes.step.fillLooseCut=0.90;
    
    param.snakes.refine.TEShrink=true;
    param.snakes.refine.LEShrink=true;
    
    param.snakes.refine.edgeFinish='sharpen';
end

function [param]=BuzmanBiplane()
    
    [param]=HalfWedge();
    
    param.general.typDat='buzmanbiplane';
    param.snakes.refine.edgeFinish='sharpen';
    param.snakes.refine.axisRatio=0.135;
    param.optiminit.cellLevels=[4,9];
    sizeRatio=param.optiminit.cellLevels(1,:)+2;
    sizeRatio=sizeRatio(2)/sizeRatio(1);
    param.general.passDomBounds(2,:)=param.general.passDomBounds(2,:)*sizeRatio;
    
end

function [param]=BuzmanBiplane2()
    
    [param]=HalfWedge();
    
    param.general.typDat='buzmanbiplane2';
    param.snakes.refine.edgeFinish='sharpen';
    param.snakes.refine.axisRatio=0.05;
    param.optiminit.cellLevels=[4,14];
    sizeRatio=param.optiminit.cellLevels(1,:)+2;
    sizeRatio=sizeRatio(2)/sizeRatio(1);
    param.general.passDomBounds(2,:)=param.general.passDomBounds(2,:)*sizeRatio;
    
end

function [param]=BuzmanBiplane3()
    
    [param]=HalfWedge();
    
    param.general.typDat='buzmanbiplane3';
    param.snakes.refine.edgeFinish='sharpen';
    param.snakes.refine.axisRatio=0.2;
    param.optiminit.cellLevels=[6,9];
    sizeRatio=param.optiminit.cellLevels(1,:)+2;
    sizeRatio=sizeRatio(2)/sizeRatio(1);
    param.general.passDomBounds(2,:)=param.general.passDomBounds(2,:)*sizeRatio;
    param.snakes.step.fillLooseStep=30;
    param.snakes.step.fillLooseCut=0.5;
    param.snakes.step.convCheckRate=25;
    param.snakes.step.convCheckRange=15;
    param.snakes.step.convDistance=50;
end

function [param]=BuzmanBiplane4()
    
    [param]=HalfWedge();
    
    param.general.typDat='buzmanbiplane4';
    param.snakes.refine.edgeFinish='sharpen';
    param.snakes.refine.axisRatio=0.23;
    param.optiminit.cellLevels=[6,9];
    sizeRatio=param.optiminit.cellLevels(1,:)+2;
    sizeRatio=sizeRatio(2)/sizeRatio(1);
    param.general.passDomBounds(2,:)=param.general.passDomBounds(2,:)*sizeRatio;
    param.snakes.step.fillLooseStep=30;
    param.snakes.step.fillLooseCut=1e-4;
    
end

function [param]=Donught()
    
    [param]=HalfWedge();
    
    param.general.typDat='donught';
    param.snakes.refine.edgeFinish='none';
    param.snakes.refine.axisRatio=1;
    param.optiminit.cellLevels=[3,3];
    sizeRatio=param.optiminit.cellLevels(1,:)+2;
    sizeRatio=sizeRatio(2)/sizeRatio(1);
    param.general.passDomBounds(2,:)=param.general.passDomBounds(2,:)*sizeRatio;
    param.snakes.step.fillLooseStep=20;
    param.snakes.step.fillLooseCut=0.5;
    
end

function [param]=Donught2()
    
    [param]=HalfWedge();
    
    param.general.typDat='donught2';
    param.snakes.refine.edgeFinish='none';
    param.snakes.refine.axisRatio=1;
    param.optiminit.cellLevels=[3,3];
    sizeRatio=param.optiminit.cellLevels(1,:)+2;
    sizeRatio=sizeRatio(2)/sizeRatio(1);
    param.general.passDomBounds(2,:)=param.general.passDomBounds(2,:)*sizeRatio;
    param.snakes.step.fillLooseStep=20;
    param.snakes.step.fillLooseCut=0.5;
    
end

function [param]=TestInit()
    
    [param]=HalfWedge();
    
    param.snakes.step.snakesSteps=150;
    param.general.typDat='testinit';
    param.snakes.refine.edgeFinish='sharpen';
    param.snakes.refine.axisRatio=1;
    param.snakes.step.fillLooseStep=20;
    param.snakes.step.fillLooseCut=1e-3;
    param.general.passDomBounds=param.general.passDomBounds*1.4;
    param.optiminit.cellLevels=[9,5];
    sizeRatio=param.optiminit.cellLevels(1,:)+2;
    sizeRatio=sizeRatio(1)/sizeRatio(2);
    param.general.passDomBounds(1,:)=param.general.passDomBounds(1,:)*sizeRatio;
end

%% validation cases
% Validation cases need:
%   CurrentValidation()
%   param.snakes.refine.refineGrid
%   param.snakes.step.snakesSteps
function [param]=CurrentValidation()
    
    [param]=DefaultCase();
    
    % Note Files
    param.results.noteFiles={'CurrentBuild'};
    param.results.tags={'snakes','Opimisation','VALIDATION','SQP','Profile Length'};
    param.results.archiveName='ParamValidation';
    
    % Local optimum avoidance params
    param.snakes.step.mergeTopo=true;

    param.snakes.force.typeSmear='dir';
    param.snakes.step.arrivalTolerance=1e-1;
    param.snakes.force.lengthEpsilon=1e-5;
    param.snakes.force.distEpsilon=0;
    param.snakes.force.dirEpsilon=1e-4;
    param.snakes.step.snaxInitPos=2e-4;
    param.snakes.step.convCheckRate=100;
    param.snakes.step.convCheckRange=15;
    param.snakes.step.convDistance=500;
    param.snakes.step.fillLooseStep=5;
    param.snakes.step.fillLooseCut=1e-3;
    param.snakes.step.maxStep=0.3;
    param.snakes.step.maxDt=0.5;
    
    % Default stepping params for validation (some cases might need more)
    param.snakes.step.snakesSteps=100;
    param.snakes.refine.refineGrid=4;
    param.snakes.refine.typeRefine='all';
    
    % Need to sort out the domain sizes to be always square?
end

function [param]=val_Snakestestsmooth1()
    
    [param]=CurrentValidation();
    
    param.general.typDat='testsmooth1';
    
    
    param.snakes.step.snakesSteps=100;
    param.snakes.refine.refineGrid=8;
    
end

function [param]=val_Snakestestsmooth1_2()
    
    [param]=CurrentValidation();
    
    param.general.typDat='testsmooth1_2';
    
    param.snakes.step.snakesSteps=100;
    param.snakes.refine.refineGrid=8;
    
end

function [param]=val_Snakestestsmooth2()
    
    
    [param]=CurrentValidation();
    
    param.general.typDat='testsmooth2';

    param.snakes.step.snakesSteps=100;
    param.snakes.refine.refineGrid=4;
    param.snakes.refine.typeRefine='all';
    
    param.general.passDomBounds=[-1,1;-0.75,0.75];
end

function [param]=val_Snakestestsmooth3()
    
    
    [param]=CurrentValidation();
    
    param.general.typDat='testsmooth3';

    param.snakes.step.snakesSteps=200;
    param.snakes.refine.refineGrid=8;
    param.snakes.refine.typeRefine='all';
    
end

function [param]=val_Snakestestsmooth3_1()
    
    
    [param]=CurrentValidation();
    
    param.general.typDat='testsmooth3_1';

    param.snakes.step.snakesSteps=200;
    param.snakes.refine.refineGrid=8;
    param.snakes.refine.typeRefine='all';
    
    param.general.passDomBounds=[-1,1;-0.6,0.6];
    param=AvoidLocalOptim(param);
end

function [param]=val_SnakesFoilVVSmall()
    
    [param]=CurrentValidation();
    
    param.general.typDat='vvlofoil';
    param.snakes.step.snakesSteps=150;
    param.snakes.refine.refineGrid=4;
    param.snakes.refine.typeRefine='all';
    param.general.passDomBounds=[-1,1;-0.4,0.4];
    
    param.snakes.step.mergeTopo=true;
    param.snakes.step.convLevel=10^-8;
    param.snakes.refine.TEShrink=true;
    param.snakes.refine.LEShrink=false;
    param.snakes.refine.edgeFinish='sharpen';
end

function [param]=val_SnakesFoilVVSmall4()
    
    [param]=CurrentValidation();
    
    param.general.typDat='vvlofoil4';
    param.snakes.step.snakesSteps=150;
    param.snakes.refine.refineGrid=4;
    param.snakes.refine.typeRefine='all';
    param.general.passDomBounds=[-1,1;-0.25,0.25];
    param.general.refineSteps=3;
    param.snakes.step.mergeTopo=true;
    
end

function [param]=val_BuzmanBiplane3()
    
    [param]=CurrentValidation();
    
    param.general.typDat='buzmanbiplane3';
    
    param.optiminit.cellLevels=[6,9];
    sizeRatio=param.optiminit.cellLevels(1,:)+2;
    sizeRatio=sizeRatio(2)/sizeRatio(1);
    param.general.passDomBounds(2,:)=param.general.passDomBounds(2,:)*sizeRatio;
    
end

function [param]=val_Donught()
    
    [param]=CurrentValidation();
    
    param.general.typDat='donught';
    param.snakes.refine.edgeFinish='none';
    param.snakes.refine.axisRatio=1;
    param.optiminit.cellLevels=[3,3];
    sizeRatio=param.optiminit.cellLevels(1,:)+2;
    sizeRatio=sizeRatio(2)/sizeRatio(1);
    param.general.passDomBounds(2,:)=param.general.passDomBounds(2,:)*sizeRatio;
    param.snakes.step.fillLooseStep=20;
    param.snakes.step.fillLooseCut=0.5;
    
end

function [param]=val_Donught2()
    
    [param]=CurrentValidation();
    
    param.general.typDat='donught2';
    param.snakes.refine.edgeFinish='none';
    
    param.snakes.step.fillLooseStep=20;
    param.snakes.step.fillLooseCut=0.5;
    
end

function [param]=val_WeirdShapeIn()
    
    [param]=CurrentValidation();
    
    param.general.typDat='low5shape';
    
   
    param.general.boundstr{1}='boundaryis1'; %'boundaryis0'
    param.general.boundstr{2}='solidisIn1';
    param.general.boundstr{3}='1bound';

end

function [param]=val_WeirdShapeOut()
    
    [param]=CurrentValidation();
    
    param.general.typDat='low5shape';
   
    
    param.general.typDat='low5shape';

end

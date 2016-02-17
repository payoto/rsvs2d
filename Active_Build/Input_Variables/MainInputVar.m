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
%
%
%
%
% Valid MainInpuVar case strings are:
%    ->  SnakesFoilSmall
%    ->  RandSmall
%    ->  SquareSnakes
%    ->
%    ->
%    ->


function [passDomBounds,passGridSteps,refineSteps,passPadding...
        ,typDat,typeBound,loadLogical,useSnakes,isCheckRes,snakesPlotInterval...
        ,snakesSteps,typeRefine,execTest,refineGrid,makeMov]=MainInputVar(caseStr)
    % Main function that allows changes
    
    [passDomBounds,passGridSteps,refineSteps,passPadding...
        ,typDat,typeBound,loadLogical,useSnakes,isCheckRes,snakesPlotInterval...
        ,snakesSteps,typeRefine,execTest,refineGrid,makeMov]=eval(caseStr);
    

end


function [passDomBounds,passGridSteps,refineSteps,passPadding...
        ,typDat,typeBound,loadLogical,useSnakes,isCheckRes,snakesPlotInterval...
        ,snakesSteps,typeRefine,execTest,refineGrid,makeMov]=SnakesFoilSmall()
    
    passDomBounds=[-1,1;-1,1];
     % number of steps in design domain
    passGridSteps=3; 
     % number of refining steps
    refineSteps=1;
    passPadding=1;
    
    typDat='foillo';
    typeBound='snaxel'; % 'vertex' or 'snaxel'
    loadLogical=false;
    useSnakes=true;
    isCheckRes=true;
    snakesPlotInterval=1;
    snakesSteps=20;
    refineGrid=2;
    typeRefine='grey';
    execTest=false;
    makeMov=true;
    
end


function [passDomBounds,passGridSteps,refineSteps,passPadding...
        ,typDat,typeBound,loadLogical,useSnakes,isCheckRes,snakesPlotInterval...
        ,snakesSteps,typeRefine,execTest,refineGrid,makeMov]=SnakesFoilVSmall()
    
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

function [passDomBounds,passGridSteps,refineSteps,passPadding...
        ,typDat,typeBound,loadLogical,useSnakes,isCheckRes,snakesPlotInterval...
        ,snakesSteps,typeRefine,execTest,refineGrid,makeMov]=SnakesdoubleBody()
    
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

function [passDomBounds,passGridSteps,refineSteps,passPadding...
        ,typDat,typeBound,loadLogical,useSnakes,isCheckRes,snakesPlotInterval...
        ,snakesSteps,typeRefine,execTest,refineGrid,makeMov]=SnakesFoilVVSmall()
    
    passDomBounds=[-1,1;-1,1];
     % number of steps in design domain
    passGridSteps=3; 
     % number of refining steps
    refineSteps=2;
    passPadding=1;
    
    typDat='vvlofoil';
    typeBound='snaxel'; % 'vertex' or 'snaxel'
    loadLogical=false;
    useSnakes=true;
    isCheckRes=true;
    snakesPlotInterval=0;
    snakesSteps=100;
    refineGrid=2;
    typeRefine='grey';
    execTest=false;
    makeMov=false;
    
end

function [passDomBounds,passGridSteps,refineSteps,passPadding...
        ,typDat,typeBound,loadLogical,useSnakes,isCheckRes,snakesPlotInterval...
        ,snakesSteps,typeRefine,execTest,refineGrid,makeMov]=Snakestestsmooth1()
    
    passDomBounds=[-1,1;-1,1];
     % number of steps in design domain
    passGridSteps=3; 
     % number of refining steps
    refineSteps=2;
    passPadding=1;
    
    typDat='testsmooth1';
    typeBound='snaxel'; % 'vertex' or 'snaxel'
    loadLogical=false;
    useSnakes=true;
    isCheckRes=true;
    snakesPlotInterval=0;
    snakesSteps=15;
    refineGrid=4;
    typeRefine='all';
    execTest=false;
    makeMov=false;
    
end

function [passDomBounds,passGridSteps,refineSteps,passPadding...
        ,typDat,typeBound,loadLogical,useSnakes,isCheckRes,snakesPlotInterval...
        ,snakesSteps,typeRefine,execTest,refineGrid,makeMov]=Snakestestsmooth4()
    
    passDomBounds=[-1,1;-1,1];
     % number of steps in design domain
    passGridSteps=3; 
     % number of refining steps
    refineSteps=2;
    passPadding=1;
    
    typDat='testsmooth4';
    typeBound='snaxel'; % 'vertex' or 'snaxel'
    loadLogical=false;
    useSnakes=true;
    isCheckRes=true;
    snakesPlotInterval=0;
    snakesSteps=200;
    refineGrid=4;
    typeRefine='grey';
    execTest=false;
    makeMov=false;
    
end

function [passDomBounds,passGridSteps,refineSteps,passPadding...
        ,typDat,typeBound,loadLogical,useSnakes,isCheckRes,snakesPlotInterval...
        ,snakesSteps,typeRefine,execTest,refineGrid,makeMov]=RandSmall()
    
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

function [passDomBounds,passGridSteps,refineSteps,passPadding...
        ,typDat,typeBound,loadLogical,useSnakes,isCheckRes,snakesPlotInterval...
        ,snakesSteps,typeRefine,execTest,refineGrid,makeMov]=SquareSnakes()
    
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

function [passDomBounds,passGridSteps,refineSteps,passPadding...
        ,typDat,typeBound,loadLogical,useSnakes,isCheckRes,snakesPlotInterval...
        ,snakesSteps,typeRefine,execTest,refineGrid,makeMov]=WeirdShape()
    
    passDomBounds=[-1.2,1.2;-1.2,1.2];
     % number of steps in design domain
    passGridSteps=10; 
     % number of refining steps
    refineSteps=2;
    passPadding=1;
    
    typDat='low5shape';
    typeBound='snaxel'; % 'vertex' or 'snaxel'
    loadLogical=false;
    useSnakes=true;
    isCheckRes=false;
    snakesPlotInterval=0;
    snakesSteps=200;
    
    refineGrid=4;
    typeRefine='all';
    execTest=false;
    makeMov=false;
    
    
end

function [passDomBounds,passGridSteps,refineSteps,passPadding...
        ,typDat,typeBound,loadLogical,useSnakes,isCheckRes,snakesPlotInterval...
        ,snakesSteps,typeRefine,execTest,refineGrid,makeMov]=UniNoPlot()
    
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

function [passDomBounds,passGridSteps,refineSteps,passPadding...
        ,typDat,typeBound,loadLogical,useSnakes,isCheckRes,snakesPlotInterval...
        ,snakesSteps,typeRefine,execTest,refineGrid,makeMov]=UniLo2InitPlot() %#ok<*DEFNU>
    
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


function [passDomBounds,passGridSteps,refineSteps,passPadding...
        ,typDat,typeBound,loadLogical,useSnakes,isCheckRes,snakesPlotInterval...
        ,snakesSteps,typeRefine,execTest,refineGrid,makeMov]=WedgInitPlot() %#ok<*DEFNU>
    
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


function [passDomBounds,passGridSteps,refineSteps,passPadding...
        ,typDat,typeBound,loadLogical,useSnakes,isCheckRes,snakesPlotInterval...
        ,snakesSteps,typeRefine,execTest,refineGrid,makeMov]=WedgPlot() %#ok<*DEFNU>
    
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


function [passDomBounds,passGridSteps,refineSteps,passPadding...
        ,typDat,typeBound,loadLogical,useSnakes,isCheckRes,snakesPlotInterval...
        ,snakesSteps,typeRefine,execTest,refineGrid,makeMov]=WeirdShapeNoPlot()
    
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

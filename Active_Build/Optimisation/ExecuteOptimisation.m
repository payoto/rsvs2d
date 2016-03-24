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

function []=ExecuteOptimisation(caseStr)
    close all
    clc
    
    [paramoptim,outinfo,unstrGrid,baseGrid,gridrefined,connectstructinfo,unstrRef,restartsnake]...
        =InitialiseOptimisation(caseStr);
    
    
    
    % Specify starting population
    
    %% Star optimisation Loop
    %
    % Assign design variables to grid
    
    % Compute Shape using snakes
    
    % Evaluate Objective Function
    
    % create new population
    
    %% Finish Optimisation
    
    
    newFill=ones([paramoptim.general.nDesVar,1]);
    [newGrid,newRefGrid,newrestart]=ReFillGrids(baseGrid,gridrefined,restartsnake,connectstructinfo,newFill);
    diary off
end

%%  Optimisation Operation Blocks

function [paramoptim,outinfo,iterstruct,unstrGrid,baseGrid,gridrefined,connectstructinfo,unstrRef,restartsnake]...
        =InitialiseOptimisation(caseStr)
    
    procStr='INITIALISE OPTIMISATION PROCESS';
    [tStart]=PrintStart(procStr,1);
    % Initialise Workspace
    
    diaryFile=[cd,'\Result_Template\Latest_Diary.log'];
    fidDiary=fopen(diaryFile,'w');
    fclose(fidDiary);
    diary(diaryFile);
    
    include_EdgeInformation
    include_SnakeParam
    include_EdgeInformation
    include_Utilities
    include_PostProcessing
    include_Mex_Wrapper
    
    
    
    % Initialise Optimisation
    % Get Parametrisation parameters
    paramoptim=StructOptimParam(caseStr);
    [outinfo]=OptimisationOutput('init',paramoptim);
    % Initialise Grid
    [unstrGrid,baseGrid,gridrefined,connectstructinfo,unstrRef,loop]...
        =GridInitAndRefine(paramoptim.parametrisation);
    
    paramoptim.general.nDesVar=sum([baseGrid.cell(:).isactive]);
    [iterstruct]=InitialiseIterationStruct(paramoptim,paramoptim.general.nDesVar);
    [~,~,~,~,restartsnake]=ExecuteSnakes(gridrefined,loop,...
        baseGrid,connectstructinfo,paramoptim.initparam,outinfo,0,0);
    
    
    
    [~]=PrintEnd(procStr,1,tStart);
end

function [population]=PerformIteration(paramoptim,outinfo,nIter,population,gridrefined,restartsnake,...
        baseGrid,connectstructinfo)
    procStr=['ITERATION ',int2str(nIter)];
    [tStart]=PrintStart(procStr,1);
    
    varExtract={'nPop'};
    [nPop]=ExtractVariables(varExtract,param);
    
    for ii=1:nPop
        
        currentMember=population(ii).fill;
        [newGrid,newRefGrid,newrestartsnake]=ReFillGrids(baseGrid,gridrefined,restartsnake,connectstructinfo,currentMember);
        
        [~,~,snakSave,loop,~]=ExecuteSnakes(newRefGrid,newrestartsnake,...
            newGrid,connectstructinfo,paramoptim.parametrisation,outinfo,nIter,1);

    end
    
    [~]=PrintEnd(procStr,1,tStart);
end

%% Parametrisation Interface

function [unstructured,unstructReshape,gridrefined,connectstructinfo,unstructuredrefined,loop]...
        =GridInitAndRefine(param)
    % Executes the Grid Initialisation process
    procStr='INITIAL GRID OPERATIONS';
    [tStart]=PrintStart(procStr,2);
    
    
    [unstructured,~,unstructReshape]=...
        GridInitialisationV2(param);
    [gridrefined,connectstructinfo,unstructuredrefined,loop]=...
        GridRefinement(unstructReshape,param);
    
    
    [~]=PrintEnd(procStr,2,tStart);
    
end

function [snaxel,snakposition,snakSave,looprestart,restartsnake]...
        =ExecuteSnakes(gridrefined,looprestart,baseGrid,connectstructinfo,param,outinfo,nIter,nProf)
    % Executes the snakes edge detection process
    %
    procStr='SNAKE PROCESS';
    [tStart]=PrintStart(procStr,2);
    
    varExtract={'refineSteps'};
    [refineSteps]=ExtractVariables(varExtract,param);
    
    
    [snaxel,snakposition,snakSave,loopsnaxel,restartsnake]=Snakes(gridrefined,looprestart,...
        baseGrid,connectstructinfo,param);

    if length(loopsnaxel)==length(looprestart)
        for ii=1:length(loopsnaxel)
            looprestart(ii).snaxel=loopsnaxel(ii).snaxel;
        end
    else
        looprestart=loopsnaxel;
    end
    
    
    looprestart=SubdivisionSurface_Snakes(looprestart,refineSteps,param);
    
    OptimisationOutput('profile',param,outinfo,nIter,nProf,looprestart,restartsnake,snakSave);
    [~]=PrintEnd(procStr,2,tStart);
end

function [newGrid,newRefGrid,newRestart]=ReFillGrids(baseGrid,refinedGrid,restartsnake,connectstructinfo,newFill)
    
    activeCell=logical([baseGrid.cell(:).isactive]);
    activeInd=[baseGrid.cell((activeCell)).index];
    
    connecInd=[connectstructinfo.cell(:).old];
    activConnecSub=FindObjNum([],activeInd,connecInd);
    activeCellSub=find(activeCell);
    refCellInd=[refinedGrid.cell(:).index];
    
    cellCentreInd=[restartsnake.cellCentredGrid(:).index];
    
    if numel(newFill)~=numel(activeCellSub)
        error('Fill and Active Set do not match in size')
    end
    
    newGrid=baseGrid;
    newRefGrid=refinedGrid;
    newRestart=restartsnake;
    for ii=1:length(activeCellSub)
        newGrid.cell(activeCellSub(ii)).fill=newFill(ii);
        newRestart.volfracconnec.cell(activeCellSub(ii)).targetfill=newFill(ii);
        newSub=FindObjNum([],[connectstructinfo.cell(activConnecSub(ii)).new],refCellInd);
        [newRefGrid.cell(newSub).fill]=deal(newFill(ii));
        newSub=FindObjNum([],[connectstructinfo.cell(activConnecSub(ii)).new],cellCentreInd);
        [newRestart.cellCentredGrid(newSub).fill]=deal(newFill(ii));
    end
    
end

%% Optimisation Specific Operations

function []=InitialisePopulation()
    
    
end

function [iterstruct]=InitialiseIterationStruct(paroptim,nDesVar)
    
    varExtract={'nPop','nIter'};
    [nPop,nIter]=ExtractVariables(varExtract,paroptim);
    
    [valFill{1:nPop}]=deal(zeros([1,nDesVar]));
    
    population=struct('fill',valFill,'location','','objective',[],'contraint',true);
    
    [iterstruct(1:nIter).population]=deal(population);
    
end

%% Print to screen functions

function [tStart]=PrintStart(procStr,lvl)
    
    procStart=[procStr,' start'];
    tStart=now;
    switch lvl
        case 1
            disp('  ')
            disp('--------------------------------------------------------------------------------------------')
            disp('--------------------------------------------------------------------------------------------')
            disp(procStart)
            disp(datestr(tStart,0))
            disp('--------------------------------------------------------------------------------------------')
            disp('  ')
            
        case 2
            disp('  ')
            disp('-----------------------')
            disp(procStart)
            
        case 3
            disp('----------')
            disp(procStart)
    end
    
end

function [tElapsed]=PrintEnd(procStr,lvl,tStart)
    
    procStart=[procStr,' end'];
    tEnd=now;
    tElapsed=tEnd-tStart;
    switch lvl
        case 1
            disp('  ')
            disp('--------------------------------------------------------------------------------------------')
            disp(['    Time Elapsed:',datestr(tElapsed,'HH:MM:SS:FFF')]);
            disp(procStart)
            disp('--------------------------------------------------------------------------------------------')
            disp('--------------------------------------------------------------------------------------------')
            disp('  ')
            
        case 2
            
            
            disp(['    Time Elapsed:',datestr(tElapsed,'HH:MM:SS:FFF')]);
            disp(procStart)
            disp('-----------------------')
            disp('  ')
        case 3
            disp('----------')
            
    end
    
end

%% Subdivision process

function [loop]=SubdivisionSurface_Snakes(loop,refineSteps,param)
    % Function taking in a closed loop of vertices and applying the subdivision
    % process
    % typBOund is te type of boundary that is expected, it can either be the
    % string 'vertex' (default) or the string 'snaxel' to show that the snaxel
    % process has been used
     varExtract={'typeBound','subdivType'};
    [typeBound,subdivType]=ExtractVariables(varExtract,param);
    if ~exist('typeBound','var'), typeBound='vertex'; end
    
    for ii=1:length(loop)
        startPoints=loop(ii).(typeBound).coord;
        loop(ii).isccw=CCWLoop(startPoints);
        newPoints=SubDivision(startPoints,refineSteps,subdivType);
        %newPoints=SubSurfBSpline(startPoints(1:end-2,:),refineSteps);
        loop(ii).subdivision=newPoints;
    end
end

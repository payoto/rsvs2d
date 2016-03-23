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
    
    [paramoptim,unstrGrid,baseGrid,gridrefined,connectstructinfo,unstrRef,restartsnake]...
        =InitialiseOptimisation(caseStr);
    
    [~,~,snakSave,loop,~]=ExecuteSnakes(unstrRef,restartsnake,...
        unstrGrid,connectstructinfo,paramoptim.parametrisation);
    
    % Specify starting population
    
    %% Star optimisation Loop
    %
    % Assign design variables to grid
    
    % Compute Shape using snakes
    
    % Evaluate Objective Function
    
    % create new population
    
    %% Finish Optimisation
    
    
    newFill=ones([paramoptim.general.nDesVar,1]);
    [newGrid,newRefGrid]=ReFillGrids(baseGrid,gridrefined,connectstructinfo,newFill);
    diary off
end

%%  Optimisation Operation Blocks

function [paramoptim,unstrGrid,baseGrid,gridrefined,connectstructinfo,unstrRef,restartsnake]...
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
    % Initialise Grid
    [unstrGrid,baseGrid,gridrefined,connectstructinfo,unstrRef,loop]...
        =GridInitAndRefine(paramoptim.parametrisation);
    
    paramoptim.general.nDesVar=sum([baseGrid.cell(:).isactive]);
    
    [~,~,~,~,restartsnake]=ExecuteSnakes(unstrRef,loop,...
        unstrGrid,connectstructinfo,paramoptim.initparam);
    
    
    
    
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

function [snaxel,snakposition,snakSave,looprestart,restartsnake]=ExecuteSnakes(unstrRef,looprestart,...
        unstrGrid,connectstructinfo,param)
    % Executes the snakes edge detection process
    %
    procStr='SNAKE PROCESS';
    [tStart]=PrintStart(procStr,2);
    
    [snaxel,snakposition,snakSave,loopsnaxel,restartsnake]=Snakes(unstrRef,looprestart,...
        unstrGrid,connectstructinfo,param);

    if length(loopsnaxel)==length(looprestart)
        for ii=1:length(loopsnaxel)
            looprestart(ii).snaxel=loopsnaxel(ii).snaxel;
        end
    else
        looprestart=loopsnaxel;
    end
    
    
    
    [~]=PrintEnd(procStr,2,tStart);
end

function [newGrid,newRefGrid]=ReFillGrids(baseGrid,refinedGrid,connectstructinfo,newFill)
    
    activeCell=logical([baseGrid.cell(:).isactive]);
    activeInd=[baseGrid.cell((activeCell)).index];
    
    connecInd=[connectstructinfo.cell(:).old];
    activConnecSub=FindObjNum([],activeInd,connecInd);
    activeCellSub=find(activeCell);
    refCellInd=[refinedGrid.cell(:).index];
    
    if numel(newFill)~=numel(activeCellSub)
        error('Fill and Active Set do not match in size')
    end
    
    newGrid=baseGrid;
    newRefGrid=refinedGrid;
    
    for ii=1:length(activeCellSub)
        newGrid.cell(activeCellSub(ii)).fill=newFill(ii);
        
        newSub=FindObjNum([],[connectstructinfo.cell(activConnecSub(ii)).new],refCellInd);
        [newRefGrid.cell(newSub).fill]=deal(newFill(ii));
    end
    
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




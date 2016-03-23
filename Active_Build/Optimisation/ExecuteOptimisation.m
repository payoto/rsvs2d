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
    %% Initialise Workspace
    close all
    clc
    
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
    
    
    
    %% Initialise Optimisation
    % Get Parametrisation parameters
    paramoptim=StructOptimParam(caseStr);
    % Initialise Grid
    [unstructured,unstructReshape,gridrefined,connectstructinfo,unstructuredrefined,loop]...
        =GridInitAndRefine(paramoptim.parametrisation);
    
    % Specify starting population
    %% Star optimisation Loop
    
    % Assign design variables to grid
    
    % Compute Shape using snakes
    
    % Evaluate Objective Function
    
    % create new population
    
    %% Finish Optimisation
    
    
    
    
    diary off
end


%% Parametrisation

function [unstructured,unstructReshape,gridrefined,connectstructinfo,unstructuredrefined,loop]...
        =GridInitAndRefine(param)
    % Executes the Grid Initialisation process
    disp('  ')
    disp('----------------------------------------------------------------------')
    disp('INITIAL GRID OPERATIONS')
    t1=now;
    [unstructured,~,unstructReshape]=...
        GridInitialisationV2(param);
    [gridrefined,connectstructinfo,unstructuredrefined,loop]=...
        GridRefinement(unstructReshape,param);
    t2=now;
    disp(['    Time taken:',datestr(t2-t1,'HH:MM:SS:FFF')]);
    disp('GENERATION PROCESS end')
    
    disp('----------------------------------------------------------------------')
    disp('  ')
    
end



function []=ReFillGrids(baseGrid,refineGrid,connectstructinfo,newFill)
    
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
    newRefGrid=refineGrid;
    
    for ii=1:length(activeCellSub)
        newGrid.cell(activeCellSub(ii)).fill=newFill(ii);
        
        newSub=FindObjNum([],connectstructinfo.cell(activConnecSub).old,refCellInd);
        [newRefGrid.cell(newSub).fill]=deal(newFill(ii));
    end
    
end











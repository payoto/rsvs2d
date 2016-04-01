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

function [iterstruct]=ExecuteOptimisation(caseStr)
    close all
    clc
    procStr2=['OPTIMISATION - ',caseStr];
    [tStartOpt]=PrintStart(procStr2,0);
    %clusterObj=parcluster('OptimSnakes');
    [paramoptim,outinfo,iterstruct,unstrGrid,baseGrid,gridrefined,connectstructinfo,unstrRef,restartsnake]...
        =InitialiseOptimisation(caseStr);
    varExtract={'maxIter'};
    [maxIter]=ExtractVariables(varExtract,paramoptim);
    
    
    % Specify starting population
    
    % Start optimisation Loop
    for nIter=1:maxIter
        % Assign design variables to grid
        procStr=['ITERATION ',int2str(nIter)];
        [tStart]=PrintStart(procStr,1);
        % Compute Shape using snakes
        [iterstruct(nIter).population]=PerformIteration(paramoptim,outinfo,nIter,iterstruct(nIter).population,gridrefined,restartsnake,...
            baseGrid,connectstructinfo);
        % Evaluate Objective Function
        [iterstruct]=GenerateNewPop(paramoptim,iterstruct,nIter);
        % create new population
        [~]=PrintEnd(procStr,1,tStart);
    end
    %% Finish Optimisation
    iterstruct(end)=[];
    OptimisationOutput('final',paramoptim,outinfo,iterstruct);
    
    [~]=PrintEnd(procStr2,0,tStartOpt);
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
    
    % Start Parallel Pool
    
    if numel(gcp('nocreate'))==0
        poolName=parallel.importProfile('ExportOptimSnakes.settings');
        clusterObj=parcluster(poolName);
        clusterObj.NumWorkers=paramoptim.general.worker;
        saveProfile(clusterObj);
        parpool(poolName)
    end
   
    paramoptim.general.nDesVar=sum([baseGrid.cell(:).isactive]);
    [iterstruct]=InitialiseIterationStruct(paramoptim,paramoptim.general.nDesVar);
    [iterstruct]=InitialisePopulation(paramoptim,iterstruct);
    
    [~,~,~,~,restartsnake]=ExecuteSnakes_Optim(gridrefined,loop,...
        baseGrid,connectstructinfo,paramoptim.initparam,paramoptim.spline,outinfo,0,0);
    
    [~]=PrintEnd(procStr,1,tStart);
end

function [population]=PerformIteration(paramoptim,outinfo,nIter,population,gridrefined,restartsnake,...
        baseGrid,connectstructinfo)
    
    
    varExtract={'nPop','objectiveName'};
    [nPop,objectiveName]=ExtractVariables(varExtract,paramoptim);
    
    paramsnake=paramoptim.parametrisation;
    paramspline=paramoptim.spline;
    
    parfor ii=1:nPop
        
        currentMember=population(ii).fill;
        [newGrid,newRefGrid,newrestartsnake]=ReFillGrids(baseGrid,gridrefined,restartsnake,connectstructinfo,currentMember);
        
        [~,~,snakSave,loop,~,outTemp]=ExecuteSnakes_Optim(newRefGrid,newrestartsnake,...
            newGrid,connectstructinfo,paramsnake,paramspline,outinfo,nIter,ii);
        population(ii).location=outTemp.dirprofile;
        population(ii).objective=EvaluateObjective(objectiveName,paramoptim,population(ii),loop);
    end
    
    [outinfo]=OptimisationOutput('iteration',paramoptim,nIter,outinfo,population);
    
end

function [iterstruct]=GenerateNewPop(paramoptim,iterstruct,nIter)
    procStr=['Generate New Population'];
    [tStart]=PrintStart(procStr,2);
    
    varExtract={'nPop',};
    [nPop]=ExtractVariables(varExtract,paramoptim);
    
    [newPop,iterstruct(nIter).population]=OptimisationMethod(paramoptim,iterstruct(nIter).population,...
        iterstruct(max([nIter-1,1])).population);
    
    for ii=1:nPop
        iterstruct(nIter+1).population(ii).fill=newPop(ii,:);
    end
    
    [~]=PrintEnd(procStr,2,tStart);
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

function [iterstruct]=InitialisePopulation(paroptim,iterstruct)
    
    varExtract={'nDesVar','nPop','startPop'};
    [nDesVar,nPop,startPop]=ExtractVariables(varExtract,paroptim);
    
    switch startPop
        case 'rand'
           origPop=rand([nPop,nDesVar]);
        case 'randuniform'
            
           origPop=rand([nPop,1])*ones([1 nDesVar]);
    end
    
    
    for ii=1:nPop
        iterstruct(1).population(ii).fill=origPop(ii,:);
    end
end

function [iterstruct]=InitialiseIterationStruct(paroptim,nDesVar)
    
    varExtract={'nPop','maxIter'};
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
        case 0
            disp('  ')
            disp('--------------------------------------------------------------------------------------------')
            disp('********************************************************************************************')
            disp('--------------------------------------------------------------------------------------------')
            disp(procStart)
            disp(datestr(tStart,0))
            disp('--------------------------------------------------------------------------------------------')
            disp('  ')
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
        case 0
            disp('  ')
            disp('--------------------------------------------------------------------------------------------')
            disp(procStart)
            disp(['    Time Elapsed:',datestr(tElapsed,'HH:MM:SS:FFF')]);
            disp(datestr(tStart,0))
            disp('--------------------------------------------------------------------------------------------')
            disp('********************************************************************************************')
            disp('--------------------------------------------------------------------------------------------')
            disp('  ')
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

%% Objective Function

function [objValue]=EvaluateObjective(objectiveName,paramoptim,member,loop)
    
    objValue=[];
    objValue=eval([objectiveName,'(paramoptim,member,loop);']);
    
    
end

function [objValue]=LengthArea(paramoptim,member,loop)
    
    points=loop.snaxel.coord(1:end-1,:);
    [A]=abs(CalculatePolyArea(points));
    vec=points([end,1:end-1],:)-points;
    L=sum(sqrt(sum(vec.^2,2)));
    
    objValue=A/L;
    
end

function [A]=CalculatePolyArea(points)
    
    pointsVec=points';
    pointsVec=pointsVec(:);
    plot(points(:,1),points(:,2));
    n=length(points(:,1));
    centreMat=eye(2*n);
    centreMat=(centreMat+centreMat(:,[end-1:end,1:end-2]))*0.5;
    
    [rotDif]=[0 -1 0 1; 1 0 -1 0];
    normMat=zeros(2*n);
    for ii=1:n-1
        normMat((2*(ii-1)+1):(2*(ii-1)+2),(2*(ii-1)+1):(2*(ii-1)+4))=rotDif;
    end
    ii=n;
    normMat((2*(ii-1)+1):(2*(ii-1)+2),(2*(ii-1)+1):(2*(ii-1)+2))=rotDif(:,1:2);
    normMat((2*(ii-1)+1):(2*(ii-1)+2),1:2)=rotDif(:,3:4);
    A=0.5*(normMat*pointsVec)'*(centreMat*pointsVec);
    
end

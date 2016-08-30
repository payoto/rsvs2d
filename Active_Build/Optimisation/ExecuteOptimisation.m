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

function [iterstruct,outinfo]=ExecuteOptimisation(caseStr,restartFromPop)
    close all
    clc
    procStr2=['OPTIMISATION - ',caseStr];
    [tStartOpt]=PrintStart(procStr2,0);
    %clusterObj=parcluster('OptimSnakes');
    [paramoptim,outinfo,iterstruct,unstrGrid,baseGrid,gridrefined,...
        connectstructinfo,unstrRef,restartsnake]...
        =InitialiseOptimisation(caseStr);
    
    varExtract={'maxIter','restartSource'};
    [maxIter,restartSource]=ExtractVariables(varExtract,paramoptim);
    startIter=1;
    firstValidIter=1;
    % Restart
    inNFlag=nargin;
    if inNFlag==2 || ~isempty(restartSource{1})
        if inNFlag==2
            restartSource=restartFromPop;
        end
        [iterstruct,startIter,maxIter,paramoptim,firstValidIter]=RestartOptions(paramoptim,inNFlag,...
            restartSource,maxIter,iterstruct);
    end
    
    % Specify starting population
    
    % Start optimisation Loop
    for nIter=startIter:maxIter
        % Assign design variables to grid
        procStr=['ITERATION ',int2str(nIter)];
        [tStart]=PrintStart(procStr,1);
        % Compute Shape using snakes
        [iterstruct(nIter).population]=PerformIteration(paramoptim,outinfo,nIter,iterstruct(nIter).population,gridrefined,restartsnake,...
            baseGrid,connectstructinfo);
        % Evaluate Objective Function
        [iterstruct,paramoptim]=GenerateNewPop(paramoptim,iterstruct,nIter,firstValidIter);
        % create new population
        [~]=PrintEnd(procStr,1,tStart);
    end
    %% Finish Optimisation
    iterstruct(end)=[];
    [~]=PrintEnd(procStr2,0,tStartOpt);
    pause(0.01)
    diary off
    OptimisationOutput('final',paramoptim,outinfo,iterstruct);
end

function [iterstruct,startIter,maxIter,paramoptim,firstValidIter]=RestartOptions(paramoptim,inNFlag,...
        restartSource,maxIter,iterstruct)
    
    
    
    load(restartSource{1})
    
    startIter=length(optimstruct);
    maxIter=startIter+maxIter;
    iterstruct=[optimstruct,iterstruct];
    
    [iterstruct,paramoptim,firstValidIter]=GenerateRestartPop(paramoptim,iterstruct,startIter,restartSource{2});
    
    paramoptim.general.restartSource=restartSource;
    startIter=startIter+1;
    
    
    
    
end

function [iterstruct,paroptim,firstValidIter]=GenerateRestartPop(paroptim,iterstruct,startIter,precOptMeth)
    
    varExtract={'optimMethod','direction','iterGap','nPop'};
    [optimMethod,direction,iterGap,nPop]=ExtractVariables(varExtract,paroptim);
    
    firstValidIter=startIter+1;
    [isGradient]=CheckIfGradient(optimMethod);
    
    if isGradient
        switch direction
            case 'min'
                [~,rootInd]=min([iterstruct(startIter).population(:).objective]);
            case 'max'
                [~,rootInd]=max([iterstruct(startIter).population(:).objective]);
        end
        
        %         [origPop,nPop,deltas]=InitialiseGradientBased(...
        %             iterstruct(startIter).population(rootInd).fill,paroptim);
        %         paroptim.general.nPop=nPop;
        [isGradientPrec]=CheckIfGradient(precOptMeth{1});
        
        if isGradientPrec
            paroptim.optim.CG.lineSearch=precOptMeth{2};
            [origPop,~,paroptim,deltas]...
                =OptimisationMethod(paroptim,iterstruct(startIter).population,...
                iterstruct(max([startIter-iterGap,1])).population);
            firstValidIter=1;
        else
            
            paroptim.optim.CG.lineSearch=true;
            paroptim.general.isRestart=true;
            [origPop,~,paroptim,deltas]...
                =OptimisationMethod(paroptim,iterstruct(startIter).population,...
                iterstruct(max([startIter-iterGap,startIter])).population);
        end
        
        
        varExtract={'nPop'};
        [nPop]=ExtractVariables(varExtract,paroptim);
        [iterstruct(startIter+1).population]=GeneratePopulationStruct(paroptim);
        for ii=1:nPop
            iterstruct(startIter+1).population(ii).fill=origPop(ii,:);
            iterstruct(startIter+1).population(ii).optimdat.var=deltas{ii}(1,:);
            iterstruct(startIter+1).population(ii).optimdat.value=deltas{ii}(2,:);
            
        end
        iterstruct(startIter+1).population=ApplySymmetry(paroptim,...
            iterstruct(startIter+1).population);
    else
        nPop
        length(iterstruct(startIter).population)
        if nPop==length(iterstruct(startIter).population)
            firstValidIter=1;
        end
        [iterstruct]=GenerateNewPop(paroptim,iterstruct,startIter,startIter);
        
    end
    
end

%%  Optimisation Operation Blocks

function [paramoptim,outinfo,iterstruct,unstrGrid,baseGrid,gridrefined,...
        connectstructinfo,unstrRef,restartsnake]=InitialiseOptimisation(caseStr)
    
    procStr='INITIALISE OPTIMISATION PROCESS';
    [tStart]=PrintStart(procStr,1);
    % Initialise Workspace
    include_EdgeInformation
    include_SnakeParam
    include_EdgeInformation
    include_Utilities
    include_PostProcessing
    include_Mex_Wrapper
    include_Optimisation
    
    diaryFile=[cd,'\Result_Template\Latest_Diary.log'];
    diaryFile=MakePathCompliant(diaryFile);
    fidDiary=fopen(diaryFile,'w');
    fclose(fidDiary);
    diary(diaryFile);
    
    
    % Initialise Optimisation
    % Get Parametrisation parameters
    paramoptim=StructOptimParam(caseStr);
    [outinfo]=OptimisationOutput('init',paramoptim);
    [~,paramoptim]=ConstraintMethod('init',paramoptim,[]);
    % Initialise Grid
    [unstrGrid,baseGrid,gridrefined,connectstructinfo,unstrRef,loop]...
        =GridInitAndRefine(paramoptim.parametrisation);
    [desvarconnec,~,~]=ExtractVolumeCellConnectivity(baseGrid);
    [paramoptim.general.desvarconnec]=...
        ExtractDesignVariableConnectivity(baseGrid,desvarconnec);
    % Start Parallel Pool
    
    if numel(gcp('nocreate'))==0
        comStr=computer;
        if strcmp(comStr(1:2),'PC')
            poolName=parallel.importProfile('ExportOptimSnakes.settings');
        else
            poolName=parallel.importProfile('ExportOptimSnakesLinux.settings');
        end
        clusterObj=parcluster(poolName);
        clusterObj.NumWorkers=paramoptim.general.worker;
        saveProfile(clusterObj);
        parpool(poolName)
    end
    
    [paramoptim]=OptimisationParametersModif(paramoptim,baseGrid);
    [iterstruct,paramoptim]=InitialisePopulation(paramoptim);
    
    iterstruct(1).population=ApplySymmetry(paramoptim,iterstruct(1).population);
    
    [~,~,~,~,restartsnake]=ExecuteSnakes_Optim(gridrefined,loop,...
        baseGrid,connectstructinfo,paramoptim.initparam,...
        paramoptim.spline,outinfo,0,0,0);
    
    [outinfo]=OptimisationOutput('iteration',...
        paramoptim,0,outinfo,iterstruct(1),{});
    
    [~]=PrintEnd(procStr,1,tStart);
end

function [population]=PerformIteration(paramoptim,outinfo,nIter,population,...
        gridrefined,restartsnake,baseGrid,connectstructinfo)
    
    
    varExtract={'nPop','objectiveName','defaultVal','lineSearch'};
    [nPop,objectiveName,defaultVal,lineSearch]=ExtractVariables(varExtract,paramoptim);
    

    if (~CheckSnakeSensitivityAlgorithm(paramoptim)) || lineSearch
        [population,supportstruct,captureErrors]=IterateNoSensitivity(paramoptim,outinfo,nIter,population,...
            gridrefined,restartsnake,baseGrid,connectstructinfo);
    else
        [population,supportstruct,captureErrors]=IterateSensitivity(paramoptim,outinfo,nIter,population,...
            gridrefined,restartsnake,baseGrid,connectstructinfo);
    end
    nPop=numel(population);
    [paramoptim]=SetVariables({'nPop'},{nPop},paramoptim);
    
    [population,captureErrors]=ParallelObjectiveCalc...
        (objectiveName,paramoptim,population,supportstruct,captureErrors);
    
    [population]=ConstraintMethod('Res',paramoptim,population);
    population=EnforceConstraintViolation(population,defaultVal);
    
    [outinfo]=OptimisationOutput('iteration',paramoptim,nIter,outinfo,population,captureErrors);
    
end


function [population,captureErrors]=ParallelObjectiveCalc...
        (objectiveName,paramoptim,population,supportstruct,captureErrors)
    
    nPop=numel(population);
    
    parfor ii=1:nPop %
        try
            [population(ii).objective,additional]=...
                EvaluateObjective(objectiveName,paramoptim,population(ii),...
                supportstruct(ii).loop);
            fieldsAdd=fieldnames(additional);
            for jj=1:numel(fieldsAdd)
                population(ii).additional.(fieldsAdd{jj})=additional.(fieldsAdd{jj});
            end
        catch MEexception
            % Error Capture
            population(ii).constraint=false;
            population(ii).exception=[population(ii).exception,'error: ',MEexception.identifier];
            captureErrors{ii}=[captureErrors{ii},MEexception.getReport];
        end
    end
    
end

%% Normal Iteration

function [population,supportstruct,captureErrors]=IterateNoSensitivity(paramoptim,outinfo,nIter,population,...
        gridrefined,restartsnake,baseGrid,connectstructinfo)
    
    varExtract={'nPop'};
    [nPop]=ExtractVariables(varExtract,paramoptim);
    
    paramsnake=paramoptim.parametrisation;
    paramspline=paramoptim.spline;
    [population]=ConstraintMethod('DesVar',paramoptim,population);
    
    [captureErrors{1:nPop}]=deal('');
    supportstruct=repmat(struct('loop',[]),[1,nPop]);
    parfor ii=1:nPop
        %for ii=1:nPop
        
        currentMember=population(ii).fill;
        [newGrid,newRefGrid,newrestartsnake]=ReFillGrids(baseGrid,gridrefined,restartsnake,connectstructinfo,currentMember);
        try
            % Normal Execution
            [population(ii),supportstruct(ii)]=NormalExecutionIteration(population(ii),newRefGrid,newrestartsnake,...
                newGrid,connectstructinfo,paramsnake,paramspline,outinfo,nIter,ii,paramoptim);
            
            
        catch MEexception
            % Error Capture
            population(ii).constraint=false;
            population(ii).exception=['error: ',MEexception.identifier];
            captureErrors{ii}=MEexception.getReport;
        end
    end
    
end

%% Gradient Iteration
function [population,supportstruct,captureErrors]=IterateSensitivity(paramoptim,outinfo,nIter,population,...
        gridrefined,restartsnake,baseGrid,connectstructinfo)
    
    
    paramsnake=paramoptim.parametrisation;
    paramspline=paramoptim.spline;
    
    [population]=ConstraintMethod('DesVar',paramoptim,population);
    
    [population,supportstruct,restartsnake,paramsnake,paramoptim,captureErrors]=ComputeRootSensitivityPop...
        (paramsnake,paramspline,paramoptim,population,baseGrid,gridrefined,...
        restartsnake,connectstructinfo,outinfo,nIter);
    
    population=ApplySymmetry(paramoptim,population);
    [population]=ConstraintMethod('DesVar',paramoptim,population);
    
    varExtract={'nPop'};
    [nPop]=ExtractVariables(varExtract,paramoptim);
    [captureErrors{2:nPop}]=deal('');
    
    supportstruct=[supportstruct,repmat(struct('loop',[]),[1,nPop-1])];
    
    parfor ii=1:nPop-1
        %for ii=1:nPop
        currentMember=population(ii+1).fill;
        [newGrid,newRefGrid,newrestartsnake]=ReFillGrids(baseGrid,gridrefined,restartsnake,connectstructinfo,currentMember);
        try
            % Normal Execution
            [population(ii+1),supportstruct(ii+1)]=NormalExecutionIteration(population(ii+1),newRefGrid,newrestartsnake,...
                newGrid,connectstructinfo,paramsnake,paramspline,outinfo,nIter,ii+1,paramoptim);
            
        catch MEexception
            % Error Capture
            population(ii+1).constraint=false;
            population(ii+1).exception=['error: ',MEexception.identifier];
            captureErrors{ii+1}=MEexception.getReport;
        end
    end
    
end

function [newpopulation,supportstruct,restartsnake,paramsnake,paramoptim,captureErrors]=...
        ComputeRootSensitivityPop(paramsnake,paramspline,paramoptim,...
        population,baseGrid,gridrefined,restartsnake,connectstructinfo,outinfo,nIter)
    
    captureErrors{1}='';
    try
        % Compute root member profile
        currentMember=population(1).fill;
        [newGrid,newRefGrid,newrestartsnake]=ReFillGrids(baseGrid,gridrefined,restartsnake,connectstructinfo,currentMember);
        
        [population(1),supportstruct,restartsnake]=NormalExecutionIteration(population(1),newRefGrid,...
            newrestartsnake,newGrid,connectstructinfo,paramsnake,paramspline...
            ,outinfo,nIter,1,paramoptim);
        
        % Compute sensitivity
        newfillstruct=SnakesSensitivity(newRefGrid,restartsnake,...
            newGrid,paramsnake,paramoptim,currentMember);
        
        nPop=numel(newfillstruct)+1;
        [paramoptim]=SetVariables({'nPop'},{nPop},paramoptim);
        [newpopulation]=GeneratePopulationStruct(paramoptim);
        newpopulation(1)=population(1);
        for ii=2:nPop
            newpopulation(ii).fill=newfillstruct(ii-1).fill;
            newpopulation(ii).optimdat.var=newfillstruct(ii-1).optimdat.var;
            newpopulation(ii).optimdat.value=newfillstruct(ii-1).optimdat.val;
        end
        varExtract={'restart'};
        [paramsnake]=SetVariables(varExtract,{true},paramsnake);
    catch MEexception
        population(1).constraint=false;
        population(1).exception=['error: ',MEexception.identifier];
        captureErrors{1}=MEexception.getReport;
        newpopulation=population;
        supportstruct.loop=[];
        warning('Sensitivity Extraction has failed, basis will not be smoothed.')
    end
    
    
end

function [population,supportstruct,restartsnake]=NormalExecutionIteration(population,newRefGrid,...
        newrestartsnake,newGrid,connectstructinfo,paramsnake,paramspline...
        ,outinfo,nIter,ii,paramoptim)
    
    varExtract={'restart','boundstr'};
    [isRestart,boundstr]=ExtractVariables(varExtract,paramsnake);
    varExtract={'nPop'};
    [nPop]=ExtractVariables(varExtract,paramoptim);
    
    if ~isRestart
        [newRefGrid]=EdgePropertiesReshape(newRefGrid);
        [newGrid]=EdgePropertiesReshape(newGrid);
        [newrestartsnake]=GenerateEdgeLoop(newRefGrid,boundstr,true);
    end
    
    [snaxel,snakposition,snakSave,loop,restartsnake,outTemp]=...
        ExecuteSnakes_Optim(newRefGrid,newrestartsnake,...
        newGrid,connectstructinfo,paramsnake,paramspline,outinfo,nIter,ii,nPop);
    population.location=outTemp.dirprofile;
    population.additional.snaxelVolRes=snakSave(end).currentConvVolume;
    population.additional.snaxelVelResV=snakSave(end).currentConvVelocity;
    
    supportstruct.loop=loop;
end

function population=EnforceConstraintViolation(population,defaultVal)
    % constraint is now something that will go from 1 (no violation) to 0
    % full violation numbers in between will be used as a ratio of the
    % default value that must be added.
    
    isConstraint=1-[population(:).constraint];
    for ii=1:length(population)
        if ~isempty(population(ii).objective)
            
            [population(ii).objective]=population(ii).objective+...
                isConstraint(ii)*(defaultVal);
            
        else
            [population(ii).objective]=isConstraint(ii)*(defaultVal);
        end
        
    end
    
end

function [iterstruct,paramoptim]...
        =GenerateNewPop(paramoptim,iterstruct,nIter,iterStart)
    procStr=['Generate New Population'];
    [tStart]=PrintStart(procStr,2);
    
    varExtract={'nPop','iterGap','optimMethod'};
    [nPop,iterGap,optimMethod]=ExtractVariables(varExtract,paramoptim);
    
    [isGradient]=CheckIfGradient(optimMethod);
    
    [newPop,iterstruct(nIter).population,paramoptim,deltas]=OptimisationMethod(paramoptim,...
        iterstruct(nIter).population,...
        iterstruct(max([nIter-iterGap,iterStart])).population);
    varExtract={'nPop'};
    [nPop]=ExtractVariables(varExtract,paramoptim);
    [iterstruct(nIter+1).population]=GeneratePopulationStruct(paramoptim);
    
    
    if ~isGradient
        for ii=1:nPop
            iterstruct(nIter+1).population(ii).fill=newPop(ii,:);
        end
    else
        for ii=1:nPop
            iterstruct(nIter+1).population(ii).fill=newPop(ii,:);
            iterstruct(nIter+1).population(ii).optimdat.var=deltas{ii}(1,:);
            iterstruct(nIter+1).population(ii).optimdat.value=deltas{ii}(2,:);
            
        end
    end
    
    iterstruct(nIter+1).population=ApplySymmetry(paramoptim,iterstruct(nIter+1).population);
    [~]=PrintEnd(procStr,2,tStart);
end


%% Parametrisation Interface

function [unstructured,unstructReshape,gridrefined,connectstructinfo,...
        unstructuredrefined,loop]=GridInitAndRefine(param)
    % Executes the Grid Initialisation process
    procStr='INITIAL GRID OPERATIONS';
    [tStart]=PrintStart(procStr,2);
    
    
    [unstructured,~,unstructReshape]=...
        GridInitialisationV2(param);
    [gridrefined,connectstructinfo,unstructuredrefined,loop]=...
        GridRefinement(unstructReshape,param);
    
    
    [~]=PrintEnd(procStr,2,tStart);
    
end

function [newGrid,newRefGrid,newRestart]=ReFillGrids(baseGrid,refinedGrid,...
        restartsnake,connectstructinfo,newFill)
    
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
    if sum(abs(newFill))==0
        newFill(round(numel(newFill)/2))=1e-3;
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

function [paramoptim]=OptimisationParametersModif(paramoptim,baseGrid)
    
    varExtract={'symType'};
    [symType]=ExtractVariables(varExtract,paramoptim);
    varExtract={'cellLevels','corneractive'};
    [cellLevels,corneractive]=ExtractVariables(varExtract,paramoptim.parametrisation);
    
    nDesVar=sum([baseGrid.cell(:).isactive]);
    paramoptim.general.nDesVar=nDesVar;
    
    paramoptim.general.symDesVarList...
        =BuildSymmetryLists(symType,cellLevels,corneractive);
    
    paramoptim.general.notDesInd...
        =BuildExclusionList(paramoptim.general.symDesVarList);
    
    [paramoptim]=CheckiterGap(paramoptim);
end

function [paramoptim]=CheckiterGap(paramoptim)
    varExtract={'optimMethod'};
    [optimMethod]=ExtractVariables(varExtract,paramoptim);
    
    switch optimMethod
        case 'conjgradls'
            paramoptim.general.iterGap=2;
        case 'conjgrad'
            paramoptim.general.iterGap=2;
        otherwise
            
            paramoptim.general.iterGap=1;
    end
end

function [symDesVarList]=BuildSymmetryLists(symType,cellLevels,corneractive)
    
    switch symType
        case 'none'
            symDesVarList=zeros([2,0]);
        case 'horz'
            warning('Symmetry assignement not robust when taken with the rest of the program')
            nRows=cellLevels(2);
            rowMatch=zeros([2,floor(nRows/2)]);
            for ii=1:floor(nRows/2)
                rowMatch(:,ii)=[ii;(nRows-ii+1)];
            end
            
            % Section to deal with inactive corners
            rowLength=ones([1,nRows])*cellLevels(1);
            if ~corneractive
                rowLength(1)=rowLength(1)-2;
                rowLength(end)=rowLength(end)-2;
            end
            indStart=cumsum([0,rowLength]);
            % Section to build index matching lists
            cellMatch=zeros([2,indStart(floor(nRows/2)+1)]);
            for ii=1:floor(nRows/2)
                for jj=1:rowLength(ii)
                    rootInd=jj+indStart(ii);
                    mirrorInd=jj+indStart(cellLevels(2)-(ii-1));
                    
                    cellMatch(:,rootInd)=[rootInd,mirrorInd];
                end
            end
            
            symDesVarList=cellMatch;
            
        case 'vert'
            error('Not coded yet')
            
    end
    
end

function notDesInd=BuildExclusionList(symDesVarList)
    
    notDesInd=[];
    notDesInd=[notDesInd,symDesVarList(2,:)];
    
    
    notDesInd=sort(RemoveIdenticalEntries(notDesInd));
    
end

function [desvarconnec,cellGrid,vertexGrid]=ExtractVolumeCellConnectivity(unstructReshape)
    % This relies on there being no empty indices for speed
    
    [cellGrid]=CellCentredGrid(unstructReshape);
    [vertexGrid]=VertexCentredGrid(unstructReshape);
    
    nCells=length(unstructReshape.cell);
    nEdges=length(unstructReshape.edge);
    nVerts=length(vertexGrid);
    cellList=[unstructReshape.cell(:).index];
    
    desvarconnec=struct('index',[],'neighbours',[],'corners',[]);
    desvarconnec=repmat(desvarconnec,[1,nCells]);
    for ii=1:nCells
        desvarconnec(ii).index=unstructReshape.cell(ii).index;
    end
    for ii=1:nEdges % Find edge connected cells
        cellInd=unstructReshape.edge(ii).cellindex;
        cellInd(cellInd==0)=[];
        cellSub=FindObjNum([],cellInd,cellList);
        for jj=1:length(cellInd)
            desvarconnec(cellSub(jj)).neighbours=...
                [desvarconnec(cellSub(jj)).neighbours,cellInd([1:jj-1,jj+1:end])];
        end
    end
    for ii=1:nVerts % Find vertex connected cells
        cellInd=[vertexGrid(ii).cell(:).index];
        cellInd(cellInd==0)=[];
        cellSub=FindObjNum([],cellInd,cellList);
        for jj=1:length(cellInd)
            desvarconnec(cellSub(jj)).corners=...
                [desvarconnec(cellSub(jj)).corners,cellInd([1:jj-1,jj+1:end])];
        end
    end
    for ii=1:nCells % Remove duplicates sort and build corners
        desvarconnec(ii).neighbours=unique(desvarconnec(ii).neighbours);
        desvarconnec(ii).corners=unique(desvarconnec(ii).corners);
        [~,~,IB] = intersect(desvarconnec(ii).neighbours,desvarconnec(ii).corners);
        desvarconnec(ii).corners(IB)=[];
    end
    
    
end

function [desvarconnec]=ExtractDesignVariableConnectivity(baseGrid,desvarconnec)
    
    actCell=[baseGrid.cell(:).isactive];
    indCell=[baseGrid.cell(:).index];
    desCellInd=cumsum(actCell).*actCell;
    desIndC{length(actCell)}=[];
    
    actSub=find(actCell);
    
    for ii=actSub
        desIndC{ii}=desCellInd(ii);
    end
    for ii=actSub
        desvarconnec(ii).index=desIndC{ii};
        desvarconnec(ii).neighbours=[desIndC{FindObjNum([],desvarconnec(ii).neighbours,indCell)}];
        desvarconnec(ii).corners=[desIndC{FindObjNum([],desvarconnec(ii).corners,indCell)}];
    end
    desvarconnec=desvarconnec(actSub);
end

%% Optimisation Specific Operations

function [iterstruct,paroptim]=InitialisePopulation(paroptim)
    
    varExtract={'nDesVar','nPop','startPop','desVarConstr','desVarVal',...
        'optimMethod','desvarconnec','specificFillName','initInterp'};
    [nDesVar,nPop,startPop,desVarConstr,desVarVal,optimMethod,desvarconnec,specificFillName,initInterp]...
        =ExtractVariables(varExtract,paroptim);
    varExtract={'cellLevels','corneractive'};
    [cellLevels,corneractive]=ExtractVariables(varExtract,paroptim.parametrisation);
    
    switch startPop
        case 'rand'
            origPop=rand([nPop,nDesVar]);
        case 'randuniform'
            
            origPop=rand([nPop,1])*ones([1 nDesVar]);
        case 'specificfill'
            
            [origPop]=StartFromFill(nDesVar,nPop,specificFillName);
            
        case 'randuniformsharp'
            
            LEind=1+(cellLevels(1)-2)*[0:(cellLevels(2)-1)];
            TEind=cellLevels(1)-2+(cellLevels(1)-2)*[0:(cellLevels(2)-1)];
            
            origPop=rand([nPop,1])*ones([1 nDesVar]);
            
            origPop(:,LEind)=origPop(:,LEind)/2;
            origPop(:,TEind)=origPop(:,TEind)/2;
            
        case 'halfuniformsharp'
            if ~corneractive
                LEind=1+(cellLevels(1)-2)*[0:(cellLevels(2)-1)];
                TEind=cellLevels(1)-2+(cellLevels(1)-2)*[0:(cellLevels(2)-1)];
            else
                LEind=1+(cellLevels(1))*[0:(cellLevels(2)-1)];
                TEind=cellLevels(1)+(cellLevels(1))*[0:(cellLevels(2)-1)];
            end
            
            origPop=ones([nPop nDesVar])*0.5;
            
            origPop(:,LEind)=origPop(:,LEind)/2;
            origPop(:,TEind)=origPop(:,TEind)/2;
            
        case 'outerbound'
            origPop=ones([nPop,nDesVar]);
            multiplier=zeros(size(origPop));
            for ii=1:length(desvarconnec)
                multiplier(:,ii)=numel(desvarconnec(ii).neighbours);
            end
            multiplier=2.^multiplier/2.^max(max(multiplier));
            origPop=multiplier.*origPop;
            
        case 'innerbound'
            origPop=zeros([nPop,nDesVar]);
            ii=1;
            try
                while isempty(regexp(desVarConstr{ii},'LocalVolFrac', 'once'))
                    ii=ii+1;
                end
            catch
                error('Invalid constraint - initialisation pair')
            end
            origPop(1,desVarVal{ii}{1})=desVarVal{ii}{2};
        case 'horzstrip'
            nStrips=cellLevels(2);
            origPop=zeros([nPop,nDesVar]);
            for ii=1:nPop
                pop=zeros(cellLevels);
                while sum(sum(pop))==0
                    nAct=randi(nStrips);
                    stripAct=randperm(nStrips,nAct);
                    
                    for jj=stripAct
                        
                        pop(:,jj)=rand;
                        
                    end
                end
                origPop(ii,1:nDesVar)=reshape(pop,[1,nDesVar]);
            end
        case 'initbusemann'
            [origPop]=InitialisePopBuseman(cellLevels,nPop,nDesVar,desVarConstr,...
                desVarVal);
            
        case 'initaeroshell'
            [origPop]=InitialiseAeroshell(cellLevels,nPop,nDesVar,desVarConstr,...
                desVarVal);
            
        case 'initinterp'
            [rootFill]=InitialiseFromFunction(cellLevels,corneractive,initInterp,nDesVar);
            origPop=ones([nPop 1])*rootFill;
    end
    
    
    [isGradient]=CheckIfGradient(optimMethod);
    if isGradient
        [origPop,nPop,deltas]=InitialiseGradientBased(origPop(1,:),paroptim);
        paroptim.general.nPop=nPop;
        [iterstruct]=InitialiseIterationStruct(paroptim);
        
        for ii=1:nPop
            iterstruct(1).population(ii).fill=origPop(ii,:);
            iterstruct(1).population(ii).optimdat.var=deltas{ii}(1,:);
            iterstruct(1).population(ii).optimdat.value=deltas{ii}(2,:);
        end
        
    else
        paroptim.general.nPop=nPop;
        [iterstruct]=InitialiseIterationStruct(paroptim);
        [origPop]=OverflowHandling(paroptim,origPop);
        for ii=1:nPop
            iterstruct(1).population(ii).fill=origPop(ii,:);
        end
    end
end

function [origPop,nPop,deltas]=InitialiseGradientBased(rootPop,paroptim)
    
    
    varExtract={'notDesInd','varActive','diffStepSize','desVarRange','desvarconnec','borderActivation'};
    [notDesInd,varActive,diffStepSize,desVarRange,desvarconnec,borderActivation]...
        =ExtractVariables(varExtract,paroptim);
    
    [inactiveVar]=SelectInactiveVariables(rootPop,varActive,desvarconnec,borderActivation);
    [desVarList]=ExtractActiveVariable(length(rootPop),notDesInd,inactiveVar);
    nPop=length(desVarList)*length(diffStepSize)+1;
    
    origPop=ones(nPop,1)*rootPop;
    deltas{nPop}=[];
    deltas{1}=[1;0];
    for jj=1:length(diffStepSize)
        origPop(2+((jj-1)*length(desVarList)):1+((jj)*length(desVarList)),desVarList)...
            =origPop(2+((jj-1)*length(desVarList)):1+((jj)*length(desVarList)),desVarList)...
            +eye(length(desVarList))*diffStepSize(jj);
        
        for ii=1:length(desVarList)
            deltas{1+ii+((jj-1)*length(desVarList))}=[desVarList(ii);diffStepSize(jj)];
        end
    end
    overFlowDiff=false(size(origPop));
    overFlowDiff(:,desVarList)=origPop(:,desVarList)>max(desVarRange);
    origPop(overFlowDiff)=max(desVarRange);
    
end

function [origPop]=InitialisePopBuseman(cellLevels,nPop,nDesVar,desVarConstr,...
        desVarVal)
    % Initialises a random number of aerodynamic looking strips in the
    % domain
    
    minTargFill=0;
    for ii=1:length(desVarConstr)
        if strcmp(desVarConstr{ii},'MinSumVolFrac')
            minTargFill=desVarVal{ii};
            
        end
    end
    
    
    nStrips=cellLevels(2);
    origPop=zeros([nPop,nDesVar]);
    for ii=1:nPop
        pop=zeros(cellLevels);
        while sum(sum(pop))==0
            nAct=randi(ceil(nStrips/4));
            stripAct=randperm(ceil(nStrips/2),ceil(nAct/2));
            stripAct=[stripAct,stripAct+1];
            stripAct(stripAct>nStrips)=nStrips;
            posPeak=randi(cellLevels(1)-1);
            hPeak=rand([length(stripAct),1]);
            ratio=minTargFill/(2*sum(hPeak));
            if ratio>1
                hPeak=ratio*hPeak;
            end
            
            ll=1;
            for jj=stripAct
                
                currPeak=hPeak(ll);
                ll=ll+1;
                totFrac=currPeak*cellLevels(1)/2;
                volLine=zeros([cellLevels(1)-2+1,1]);
                
                for kk=2:length(volLine)-1
                    if kk<=posPeak+1
                        volLine(kk)=currPeak*(kk-1)/posPeak;
                    else
                        volLine(kk)=currPeak-currPeak*((kk-1)-posPeak)...
                            /(length(volLine)-1-posPeak);
                    end
                end
                volFrac=zeros([cellLevels(1)-2,1]);
                for kk=1:length(volFrac)
                    volFrac(kk)=mean(volLine(kk:kk+1));
                end
                volFrac=[1e-3;volFrac;1e-3];
                volFrac(volFrac>1)=1;
                pop(:,jj)=volFrac;
                
            end
        end
        origPop(ii,1:nDesVar)=reshape(pop,[1,nDesVar]);
    end
end

function [origPop]=InitialiseAeroshell(cellLevels,nPop,nDesVar,desVarConstr,...
        desVarVal)
    % Initialises a random number of aerodynamic looking strips in the
    % domain
    
    minTargFill=0;
    for ii=1:length(desVarConstr)
        if strcmp(desVarConstr{ii},'MinSumVolFrac')
            minTargFill=desVarVal{ii};
            
        end
    end
    
    
    nStrips=cellLevels(2);
    origPop=zeros([nPop,nDesVar]);
    for ii=1:nPop
        pop=zeros(cellLevels);
        while sum(sum(pop))==0
            % Generate Active strips
            nAct=randi(ceil(nStrips));
            stripAct=randperm((nStrips),ceil(nAct/2));
            stripAct=[stripAct,stripAct+1];
            stripAct(stripAct>nStrips)=nStrips;
            stripAct=sort(RemoveIdenticalEntries(stripAct));
            % Find coefficients for active strips
            nAct=numel(stripAct);
            loStrip=[0,stripAct(1:end-1)];
            hiStrip=[stripAct(2:end),nStrips];
            coeffStrip=(hiStrip-loStrip)-1;
            % Generate peaks
            posPeak=randi(cellLevels(1)-1);
            hPeak=rand([length(stripAct),1]).*coeffStrip';
            ratio=minTargFill/(2*sum(hPeak));
            if ratio>1
                hPeak=ratio*hPeak;
            end
            
            ll=1;
            for jj=stripAct
                
                currPeak=hPeak(ll);
                ll=ll+1;
                totFrac=currPeak*cellLevels(1)/2;
                volLine=zeros([cellLevels(1)-2+1,1]);
                
                for kk=2:length(volLine)-1
                    if kk<=posPeak+1
                        volLine(kk)=currPeak*(kk-1)/posPeak;
                    else
                        volLine(kk)=currPeak-currPeak*((kk-1)-posPeak)...
                            /(length(volLine)-1-posPeak);
                    end
                end
                volFrac=zeros([cellLevels(1)-2,1]);
                for kk=1:length(volFrac)
                    volFrac(kk)=mean(volLine(kk:kk+1));
                end
                volFrac=[1e-3;volFrac;1e-3];
                
                pop(:,jj)=volFrac;
                
            end
        end
        origPop(ii,1:nDesVar)=reshape(pop,[1,nDesVar]);
    end
end

function [rootFill]=InitialiseFromFunction(cellLevels,corneractive,initInterp,nDesVar)
    
    load(initInterp{1}); % loads file with variable interpFunc with the same #rows as cellLevels
    
    
    [nFunc,nDes]=size(interpFunc);
    indFunc=linspace(0,1,nDes);
    indFill=linspace(0,1,cellLevels(1));
    rootFill=zeros([1,nDesVar]);
    
    if cellLevels(2)>nFunc
        warning('insufficient number of functions some will be repeated')
        
    end
    nCellLine=ones([1,cellLevels(2)])*cellLevels(1);
    if ~corneractive
        nCellLine(1)=cellLevels(1)-2;
        nCellLine(end)=cellLevels(1)-2;
    end
    startInd=cumsum([0,nCellLine(1:end-1)])+1;
    endInd=cumsum(nCellLine);
    for ii=1:cellLevels(2)
        
        interpFill=interp1(indFunc,interpFunc(ii,:),indFill);
        
        if ~corneractive && (ii==1 || ii==cellLevels(2))
            interpFill([1,end])=[];
        end
        
        rootFill(startInd(ii):endInd(ii))=interpFill;
        
    end
    
end


function [origPop]=StartFromFill(nDesVar,nPop,fillName)
    
    switch fillName
        case '24DVaverage'
            rootFill=[   0.072190973000000
                0.197280875000000
                0.307048435000000
                0.394557265000000
                0.458644380000000
                0.500609320000000
                0.523538210000000
                0.530283520000000
                0.523530135000000
                0.499683995000000
                0.452365100000000
                0.340267790000000
                0.072190973000000
                0.197280875000000
                0.307048435000000
                0.394557265000000
                0.458644380000000
                0.500609320000000
                0.523538210000000
                0.530283520000000
                0.523530135000000
                0.499683995000000
                0.452365100000000
                0.340267790000000]';
            
        otherwise
            
            error('invalid fill')
            
    end
    
    if nDesVar~=length(rootFill)
        error('invalid fill and parametrisation combination at initialisation')
    end
    
    origPop=ones([nPop,1])*rootFill;
    
end

function [iterstruct]=InitialiseIterationStruct(paramoptim)
    
    varExtract={'maxIter'};
    [maxIter]=ExtractVariables(varExtract,paramoptim);
    
    [popstruct]=GeneratePopulationStruct(paramoptim);
    
    [iterstruct(1:maxIter).population]=deal(popstruct);
    
end

function [popstruct]=GeneratePopulationStruct(paroptim)
    varExtract={'nPop','nDesVar','objectiveName'};
    [nPop,nDesVar,objectiveName]=ExtractVariables(varExtract,paroptim);
    
    [valFill{1:nPop}]=deal(zeros([1,nDesVar]));
    switch objectiveName
        case 'CutCellFlow'
            addstruct=struct('iter',[],'res',[],'cl',[],'cm',[],'cd',[],...
                'cx',[],'cy',[],'A',[],'L',[],'t',[],'c',[],'tc',[],'snaxelVolRes',[],'snaxelVelResV',[]);
            
        case 'LengthArea'
            addstruct=struct('A',[],'L',[],'t',[],'c',[],'tc',[],'snaxelVolRes',[],'snaxelVelResV',[]);
            
    end
    optimdatstruct=struct('var',[],'value',[]);
    popstruct=struct('fill',valFill,'location','','objective',[],'constraint'...
        ,true,'optimdat',optimdatstruct,'additional',addstruct,'exception','');
    
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

function [objValue,additional]=EvaluateObjective(objectiveName,paramoptim,member,loop)
    
    procStr=['Calculate Objective - ',objectiveName];
    [tStart]=PrintStart(procStr,2);
    
    objValue=[];
    [objValue,additional]=eval([objectiveName,'(paramoptim,member,loop);']);
    
    [tElapsed]=PrintEnd(procStr,2,tStart);
end

function [objValue,additional]=LengthArea(paramoptim,member,loop)
    for ii=1:length(loop)
        points=loop(ii).snaxel.coord(1:end-1,:);
        [A(ii)]=abs(CalculatePolyArea(points));
        vec=points([end,1:end-1],:)-points;
        L(ii)=sum(sqrt(sum(vec.^2,2)));
        t(ii)=max(points(:,2))-min(points(:,2));
        xMin(ii)=min(points(:,1));
        xMax(ii)=max(points(:,1));
        
    end
    objValue=sum(A)/sum(L);
    
    additional.A=sum(A);
    additional.L=sum(L);
    additional.t=sum(t);
    additional.c=max(xMax)-min(xMin);
    additional.tc=additional.t/additional.c;
end

function [objValue,additional]=CutCellFlow(paramoptim,member,loop)
    boundaryLoc=member.location;
    
    [obj]=CutCellFlow_Handler(paramoptim,boundaryLoc);
    [~,areaAdd]=LengthArea(paramoptim,member,loop);
    objValue=obj.cd;
    
    additional=obj;
    additional.A=areaAdd.A;
    additional.L=areaAdd.L;
    additional.t=areaAdd.t;
    additional.c=areaAdd.c;
    additional.tc=areaAdd.tc;
    
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

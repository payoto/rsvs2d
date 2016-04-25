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

function [iterstruct]=ExecuteOptimisation(caseStr,restartFromPop)
    close all
    clc
    procStr2=['OPTIMISATION - ',caseStr];
    [tStartOpt]=PrintStart(procStr2,0);
    %clusterObj=parcluster('OptimSnakes');
    [paramoptim,outinfo,iterstruct,unstrGrid,baseGrid,gridrefined,connectstructinfo,unstrRef,restartsnake]...
        =InitialiseOptimisation(caseStr);
    varExtract={'maxIter','restartSource'};
    [maxIter,restartSource]=ExtractVariables(varExtract,paramoptim);
    startIter=1;
    
    % Restart
    in2Flag=nargin==2;
    if in2Flag || ~isempty(restartSource)
        if in2Flag
            restartSource=restartFromPop;
        end
        load(restartSource)
        startIter=length(optimstruct);
        maxIter=startIter+maxIter;
        iterstruct=[optimstruct,iterstruct];
        [iterstruct]=GenerateNewPop(paramoptim,iterstruct,startIter);
       paramoptim.general.restartSource=restartSource;
       startIter=startIter+1;
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
    include_EdgeInformation
    include_SnakeParam
    include_EdgeInformation
    include_Utilities
    include_PostProcessing
    include_Mex_Wrapper
    
    diaryFile=[cd,'\Result_Template\Latest_Diary.log'];
    diaryFile=MakePathCompliant(diaryFile);
    fidDiary=fopen(diaryFile,'w');
    fclose(fidDiary);
    diary(diaryFile);
    
    
    % Initialise Optimisation
    % Get Parametrisation parameters
    paramoptim=StructOptimParam(caseStr);
    [outinfo]=OptimisationOutput('init',paramoptim);
    % Initialise Grid
    [unstrGrid,baseGrid,gridrefined,connectstructinfo,unstrRef,loop]...
        =GridInitAndRefine(paramoptim.parametrisation);
    
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
    [iterstruct]=InitialiseIterationStruct(paramoptim,paramoptim.general.nDesVar);
    [iterstruct]=InitialisePopulation(paramoptim,iterstruct);
    iterstruct(1).population=ApplySymmetry(paramoptim,iterstruct(1).population);
    [~,~,~,~,restartsnake]=ExecuteSnakes_Optim(gridrefined,loop,...
        baseGrid,connectstructinfo,paramoptim.initparam,paramoptim.spline,outinfo,0,0);
    
    [~]=PrintEnd(procStr,1,tStart);
end

function [population]=PerformIteration(paramoptim,outinfo,nIter,population,gridrefined,restartsnake,...
        baseGrid,connectstructinfo)
    
    
    varExtract={'nPop','objectiveName','defaultVal'};
    [nPop,objectiveName,defaultVal]=ExtractVariables(varExtract,paramoptim);
    
    paramsnake=paramoptim.parametrisation;
    paramspline=paramoptim.spline;
    [population]=ConstraintMethod('DesVar',paramoptim,population);
    
    
    [captureErrors{1:nPop}]=deal('');
    
    parfor ii=1:nPop
    %for ii=1:nPop
        
        currentMember=population(ii).fill;
        [newGrid,newRefGrid,newrestartsnake]=ReFillGrids(baseGrid,gridrefined,restartsnake,connectstructinfo,currentMember);
        try
            % Normal Execution
            population(ii)=NormalExecutionIteration(population(ii),newRefGrid,newrestartsnake,...
        newGrid,connectstructinfo,paramsnake,paramspline,outinfo,nIter,ii,objectiveName,paramoptim)
            
            
        catch MEexception
            % Error Capture
            population(ii).constraint=false;
            population(ii).exception=['error: ',MEexception.identifier];
            captureErrors{ii}=MEexception.getReport;
        end
    end
    
    [population]=ConstraintMethod('Res',paramoptim,population);
    population=EnforceConstraintViolation(population,defaultVal);
    [outinfo]=OptimisationOutput('iteration',paramoptim,nIter,outinfo,population,captureErrors);
    
end

function population=NormalExecutionIteration(population,newRefGrid,newrestartsnake,...
        newGrid,connectstructinfo,paramsnake,paramspline,outinfo,nIter,ii,objectiveName,paramoptim)
    
    varExtract={'restart','boundstr'};
    [isRestart]=ExtractVariables(varExtract,paramsnake);
    
    if ~isRestart
        [newrestartsnake]=GenerateEdgeLoop(newRefGrid,boundstr,true);
    end
    
    [~,~,snakSave,loop,~,outTemp]=ExecuteSnakes_Optim(newRefGrid,newrestartsnake,...
        newGrid,connectstructinfo,paramsnake,paramspline,outinfo,nIter,ii);
    population.location=outTemp.dirprofile;
    [population.objective,population.additional]=EvaluateObjective(objectiveName,paramoptim,population,loop);
    population.additional.snaxelVolRes=snakSave(end).currentConvVolume;
    population.additional.snaxelVelResV=snakSave(end).currentConvVelocity;
    
    
end

function population=EnforceConstraintViolation(population,defaultVal)
    
    isConstraint=[population(:).constraint];
    
    [population(~isConstraint).objective]=deal(defaultVal);
    
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
    iterstruct(nIter+1).population=ApplySymmetry(paramoptim,iterstruct(nIter+1).population);
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

function population=ApplySymmetry(paramoptim,population)
    
    varExtract={'symDesVarList'};
    [symDesVarList]=ExtractVariables(varExtract,paramoptim);
    
    for ii=1:length(population)
        population(ii).fill(symDesVarList(2,:))=...
            population(ii).fill(symDesVarList(1,:));
    end
    
    
end

function [paramoptim]=OptimisationParametersModif(paramoptim,baseGrid)
    
    varExtract={'symType'};
    [symType]=ExtractVariables(varExtract,paramoptim);
    varExtract={'cellLevels','corneractive'};
    [cellLevels,corneractive]=ExtractVariables(varExtract,paramoptim.parametrisation);
    
    nDesVar=sum([baseGrid.cell(:).isactive]);
    paramoptim.general.nDesVar=nDesVar;
    
    switch symType
        case 'none'
            paramoptim.general.symDesVarList=zeros([2,0]);
        case 'horz'
            
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
            
            paramoptim.general.symDesVarList=cellMatch;
            
        case 'vert'
            error('Not coded yet')
            
    end
            
    
end

%% Optimisation Specific Operations

function [iterstruct]=InitialisePopulation(paroptim,iterstruct)
    
    varExtract={'nDesVar','nPop','startPop'};
    [nDesVar,nPop,startPop]=ExtractVariables(varExtract,paroptim);
    varExtract={'cellLevels'};
    [cellLevels]=ExtractVariables(varExtract,paroptim.parametrisation);
    
    switch startPop
        case 'rand'
           origPop=rand([nPop,nDesVar]);
        case 'randuniform'
            
           origPop=rand([nPop,1])*ones([1 nDesVar]);
           
        case 'horzstrip'
            nStrips=cellLevels(2);
            origPop=zeros([nPop,nDesVar]);
            for ii=1:nPop
                pop=zeros(cellLevels);
                
                nAct=randi(nStrips);
                stripAct=randperm(nStrips,nAct);
                
                for jj=stripAct
                    
                    pop(:,jj)=rand;
                    
                end
                origPop(ii,1:nDesVar)=reshape(pop,[1,nDesVar]);
            end
    end
    
    
    for ii=1:nPop
        iterstruct(1).population(ii).fill=origPop(ii,:);
    end
end

function [iterstruct]=InitialiseIterationStruct(paroptim,nDesVar)
    
    varExtract={'nPop','maxIter'};
    [nPop,nIter]=ExtractVariables(varExtract,paroptim);
    
    [valFill{1:nPop}]=deal(zeros([1,nDesVar]));
    addstruct=struct('iter',[],'res',[],'cl',[],'cm',[],'cd',[],'cx',[],'cy',[],'A',[],'L',[],'snaxelVolRes',[],'snaxelVelResV',[]);
    population=struct('fill',valFill,'location','','objective',[],'constraint'...
        ,true,'additional',addstruct,'exception','none');
    
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
    
        
    end
    objValue=sum(A)/sum(L);

    additional.A=sum(A);
    additional.L=sum(L);
end

function [objValue,additional]=CutCellFlow(paramoptim,member,loop)
    boundaryLoc=member.location;
    
    [obj]=CutCellFlow_Handler(paramoptim,boundaryLoc);
    [~,areaAdd]=LengthArea(paramoptim,member,loop);
    objValue=obj.cd;
    additional=obj;
    additional.A=areaAdd.A;
    additional.L=areaAdd.L;
    
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

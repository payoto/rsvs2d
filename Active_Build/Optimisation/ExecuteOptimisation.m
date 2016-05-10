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
        [iterstruct,paramoptim]=GenerateNewPop(paramoptim,iterstruct,nIter);
        % create new population
        [~]=PrintEnd(procStr,1,tStart);
    end
    %% Finish Optimisation
    iterstruct(end)=[];
    [~]=PrintEnd(procStr2,0,tStartOpt);
    OptimisationOutput('final',paramoptim,outinfo,iterstruct);
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
    [iterstruct,paramoptim]=InitialisePopulation(paramoptim);
    
    iterstruct(1).population=ApplySymmetry(paramoptim,iterstruct(1).population);
    [~,~,~,~,restartsnake]=ExecuteSnakes_Optim(gridrefined,loop,...
        baseGrid,connectstructinfo,paramoptim.initparam,paramoptim.spline,outinfo,0,0,0);
    [outinfo]=OptimisationOutput('iteration',paramoptim,0,outinfo,iterstruct(1),{});
    
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
        newGrid,connectstructinfo,paramsnake,paramspline,outinfo,nIter,ii,objectiveName,paramoptim);
            
            
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
        newGrid,connectstructinfo,paramsnake,paramspline,outinfo,nIter,ii,...
        objectiveName,paramoptim)
    
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
    [population.objective,population.additional]=EvaluateObjective(objectiveName,paramoptim,population,loop);
    population.additional.snaxelVolRes=snakSave(end).currentConvVolume;
    population.additional.snaxelVelResV=snakSave(end).currentConvVelocity;
    
    
end

function population=EnforceConstraintViolation(population,defaultVal)
    
    isConstraint=[population(:).constraint];
    
    [population(~isConstraint).objective]=deal(defaultVal);
    
end

function [iterstruct,paramoptim]=GenerateNewPop(paramoptim,iterstruct,nIter)
    procStr=['Generate New Population'];
    [tStart]=PrintStart(procStr,2);
    
    varExtract={'nPop','iterGap'};
    [nPop,iterGap]=ExtractVariables(varExtract,paramoptim);
    
    [newPop,iterstruct(nIter).population,paramoptim]=OptimisationMethod(paramoptim,...
        iterstruct(nIter).population,...
        iterstruct(max([nIter-iterGap,1])).population);
   varExtract={'nPop'};
    [nPop]=ExtractVariables(varExtract,paramoptim);
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
    
    paramoptim.general.symDesVarList...
        =BuildSymmetryLists(symType,cellLevels,corneractive);
            
    paroptim.general.notDesInd...
        =BuildExclusionList(paramoptim.general.symDesVarList);
    
    [paramoptim]=CheckiterGap(paramoptim);
end

function [paramoptim]=CheckiterGap(paramoptim)
    varExtract={'optimMethod'};
    [optimMethod]=ExtractVariables(varExtract,paramoptim);
    
    switch optimMethod
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

%% Optimisation Specific Operations

function [iterstruct,paroptim]=InitialisePopulation(paroptim)
    
    varExtract={'nDesVar','nPop','startPop','desVarConstr','desVarVal','optimMethod'};
    [nDesVar,nPop,startPop,desVarConstr,desVarVal,optimMethod]...
        =ExtractVariables(varExtract,paroptim);
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
    end
    
    
    [isGradient]=CheckIfGradient(optimMethod);
    if isGradient
        [origPop,nPop]=InitialiseGradientBased(origPop(1,:),paroptim);
    end
    paroptim.general.nPop=nPop;
    
    [iterstruct]=InitialiseIterationStruct(paroptim,nDesVar);
    for ii=1:nPop
        iterstruct(1).population(ii).fill=origPop(ii,:);
    end
end

function [isGradient]=CheckIfGradient(optimMethod)
    
    switch optimMethod
        
        case 'DE'
            isGradient=false;
        case 'DEtan'
            isGradient=false;
        case 'conjgrad'
            isGradient=true;
        otherwise
            isGradient=false;
            warning('Optimisation method is not known as gradient based or otherwise, no gradient is assumed')
            
    end
end

function [origPop,nPop]=InitialiseGradientBased(rootPop,paroptim)
    varExtract={'notDesInd','nPop','diffStepSize','desVarRange'};
    [notDesInd,nPop,diffStepSize,desVarRange]=ExtractVariables(varExtract,paroptim);
    inactiveVar=[];
    [desVarList]=ExtractActiveVariable(length(rootPop),notDesInd,inactiveVar);
    nPop=length(desVarList)+1;
    
    origPop=ones(nPop,1)*rootPop;
    origPop(2:end,desVarList)=origPop(2:end,desVarList)+eye(nPop-1)*diffStepSize;
    
    overFlowDiff=false(size(origPop));
    overFlowDiff(:,desVarList)=origPop(:,desVarList)>max(desVarRange);
    origPop(overFlowDiff)=max(desVarRange);
    
end

function [origPop]=InitialisePopBuseman(cellLevels,nPop,nDesVar,desVarConstr,...
        desVarVal)
    
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

function [iterstruct]=InitialiseIterationStruct(paroptim,nDesVar)
    
    varExtract={'nPop','maxIter','objectiveName'};
    [nPop,nIter,objectiveName]=ExtractVariables(varExtract,paroptim);
    
    [valFill{1:nPop}]=deal(zeros([1,nDesVar]));
    switch objectiveName
        case 'CutCellFlow'
            addstruct=struct('iter',[],'res',[],'cl',[],'cm',[],'cd',[],...
                'cx',[],'cy',[],'A',[],'L',[],'t',[],'c',[],'tc',[],'snaxelVolRes',[],'snaxelVelResV',[]);
    
        case 'LengthArea'
            addstruct=struct('A',[],'L',[],'t',[],'c',[],'tc',[],'snaxelVolRes',[],'snaxelVelResV',[]);
    
    end
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

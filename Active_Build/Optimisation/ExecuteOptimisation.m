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

%{
function [] = ExecuteOptimisation()
    %FUNCTIONLIST allows local functions to be used globally once it has
    %been used.
    
    funcHandles=localfunctions;
    funcDir=[cd,'\Automated_Function_Directory_ExecOptim'];
    HeaderActivation(funcHandles,funcDir)
    
end
%}

function [iterstruct,outinfo]=ExecuteOptimisation(caseStr,restartFromPop,...
        debugArgIn)
    %close all
    clc
    
    if nargin==3
        OptimisationDebug(caseStr,debugArgIn)
        return
    end
    flagOut=true;
    procStr2=['OPTIMISATION - ',caseStr];
    [tStartOpt]=PrintStart(procStr2,0);
    %clusterObj=parcluster('OptimSnakes');
    [paramoptim,outinfo,iterstruct,~,baseGrid,gridrefined,...
        connectstructinfo,~,restartsnake]...
        =InitialiseOptimisation(caseStr);
    
    varExtract={'maxIter','restartSource','refineSteps','refineIter','restartIterNum'};
    [maxIter,restartSource,refineSteps,refineIter,restartIterNum]=ExtractVariables(varExtract...
        ,paramoptim);
    startIter=1;
    firstValidIter=1;
    refStart=1;
    % Restart
    inNFlag=nargin;
    if inNFlag==2 || ~isempty(restartSource{1})
        if inNFlag==2
            restartSource=restartFromPop;
        end
        [paramoptim,outinfo,iterstruct,baseGrid,gridrefined,...
            connectstructinfo,restartsnake,startIter,maxIter]=...
            HandleVariableRestart(restartSource{1},paramoptim,outinfo,...
            iterstruct,baseGrid,gridrefined,connectstructinfo,restartsnake,...
            startIter,maxIter,restartIterNum);
        if refineSteps>0
            refStart=numel(outinfo);
        end
        
        [iterstruct,paramoptim,firstValidIter]=GenerateRestartPop(paramoptim,...
            iterstruct,startIter,restartSource{2},baseGrid);
        
        paramoptim=SetVariables({'restartSource'},{restartSource},paramoptim);
        startIter=startIter+1;
    end
    
    % Introduce debug lines of code
    % debugScrip
    
    % Specify starting population
    for refStage=refStart:refStart+refineSteps
        
        % Start optimisation Loop
        for nIter=startIter:maxIter
            % Assign design variables to grid
            procStr=['ITERATION ',int2str(nIter)];
            [tStart]=PrintStart(procStr,1);
            % Compute Shape using snakes
            [iterstruct(nIter).population,restartsnake]=...
                PerformIteration(paramoptim,outinfo(end),nIter,...
                iterstruct(nIter).population,gridrefined,restartsnake,...
                baseGrid,connectstructinfo,iterstruct(1:nIter-1));
            OptimisationOutput('optstruct',paramoptim,outinfo(end),...
                iterstruct);
            % Evaluate Objective Function
            [iterstruct,paramoptim]=GenerateNewPop(paramoptim,iterstruct,...
                nIter,firstValidIter,baseGrid);
            % create new population
            DiskSpaceManagement(paramoptim,iterstruct)
            [~]=PrintEnd(procStr,1,tStart);
            % Convergence tests
            if ConvergenceTest_sloperefine(paramoptim,iterstruct,nIter,startIter) ...
                    && (mod(nIter,2)==0) && (refStage~=(refineSteps+refStart))
                fprintf('\n Optimisation Stopped By Slope convergence condition \n');
                break
            end
            if ConvergenceTest_static(paramoptim,iterstruct,nIter,startIter) ...
                    && (mod(nIter,2)==0)
                fprintf('\n Optimisation Stopped By convergence condition \n');
                break
            end
        end
        
        % Finish Optimisation
        iterstruct(end)=[];
        [~]=PrintEnd(procStr2,0,tStartOpt);
        pause(0.01)
        
        if isempty(nIter)
            nIter=numel(iterstruct);
            while isempty([iterstruct(nIter).population.objective]) && nIter>1
                nIter=nIter-1;
            end
            startIter=1;
            while numel(iterstruct(nIter).population(1).fill)>...
                    numel(iterstruct(startIter).population(1).fill)
                startIter=startIter+1;
            end
            flagOutSingle=false;
        end
        %try
        %save(['PreOptimOutFinal',int2str(refStage),'.mat'])
        try
            if (exist('flagOut','var') && ~flagOut) || ...
                    (exist('flagOutSingle','var') && ~flagOutSingle)
                disp('Output Skipped')
                flagOutSingle=true;
            else
                [paramoptim]=FindKnownOptim(paramoptim,iterstruct(1:nIter),...
                    baseGrid,gridrefined,restartsnake,connectstructinfo,...
                    outinfo(end));
                OptimisationOutput('final',paramoptim,outinfo,iterstruct(1:nIter));
            end
        catch ME;
            disp(ME.getReport),
        end
        
        %debugScript2
        
        save([outinfo(end).rootDir,filesep,'DvpAnisRefineMat.mat'])
        if refStage<(refineSteps+refStart)
            %save('DvpAnisRefineMat.mat')
            [paramoptim,outinfo(end+1),iterstruct2,~,baseGrid,gridrefined,...
                connectstructinfo,~,restartsnake]=...
                HandleRefinement(paramoptim,iterstruct(1:nIter),...
                outinfo(end),baseGrid,gridrefined,connectstructinfo,...
                refStage,nIter,startIter);
            
            if ~isempty(refineIter)
                startIter=nIter+1;
                maxIter=nIter+refineIter(min(refStage,numel(refineIter)));
            else
                maxIter=nIter+maxIter-(startIter-1);
                startIter=nIter+1;
            end
            iterstruct=iterstruct2;
        end
        
    end
    diary off
end

%% Handle Restarting of previous runs


function [iterstruct,paroptim,firstValidIter]=GenerateRestartPop(paroptim,...
        iterstruct,startIter,precOptMeth,baseGrid)
    
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
                iterstruct(max([startIter-iterGap,1])).population,baseGrid);
            firstValidIter=startIter;
            nLoc=(regexprep(iterstruct(firstValidIter).population(1).location,...
                'iteration_.*$',''));
            while firstValidIter>1 && strcmp(nLoc,regexprep(iterstruct(...
                    firstValidIter-1).population(1).location,'iteration_.*$',''))
                firstValidIter=firstValidIter-1;
            end
        else
            
            paroptim.optim.CG.lineSearch=true;
            paroptim.general.isRestart=true;
            [origPop,~,paroptim,deltas]...
                =OptimisationMethod(paroptim,iterstruct(startIter).population,...
                iterstruct(max([startIter-iterGap,startIter])).population,...
                baseGrid);
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
        length(iterstruct(startIter).population)
        if nPop==length(iterstruct(startIter).population)
            firstValidIter=1;
        end
        [iterstruct]=GenerateNewPop(paroptim,iterstruct,startIter,startIter,...
            baseGrid);
        
    end
    
end

function [paramoptim,outinfo,iterstruct,baseGrid,gridrefined,...
        connectstructinfo,restartsnake,startIter,maxIter]=...
        HandleVariableRestart(restartPath,paramoptim,outinfo,iterstruct,baseGrid,gridrefined,...
        connectstructinfo,restartsnake,startIter,maxIter,restartIterNum)
    % handles variability of restart inputs
    [optionalin]=LoadRestart(restartPath);
    fieldsIn=fieldnames(optionalin);
    maxIterInit=maxIter;
    for ii=1:numel(fieldsIn)
        repeatflag=true; % allows to rerun a field if they appear in wrong order
        while repeatflag
            switch fieldsIn{ii}
                case 'optimstruct'
                    startIter=min(length(optionalin.(fieldsIn{ii})),restartIterNum);
                    maxIter=startIter+maxIterInit;
                    if ~isfield(optionalin.(fieldsIn{ii})(1).population,...
                            'nonfillvar')
                        for kk=1:numel(optionalin.(fieldsIn{ii}))
                            [optionalin.(fieldsIn{ii})(kk).population(:)...
                                .nonfillvar]=deal([]);
                        end
                    end
                    iterstruct=[optionalin.(fieldsIn{ii})(1:startIter),iterstruct]; %#ok<AGROW>
                    isOptimStruct=true;
                    repeatflag=false;
                case 'grid'
                    baseGrid=optionalin.(fieldsIn{ii}).base;
                    gridrefined=optionalin.(fieldsIn{ii}).refined;
                    connectstructinfo=optionalin.(fieldsIn{ii}).connec;
                    connectstructinfo.edge=struct([]);
                    fieldsConnec=fieldnames(connectstructinfo.cell);
                    if ~strcmp(fieldsConnec{1},'old')
                        newFields={'old','new','targetfill'};
                        connectstructinfo.cell=cell2struct(struct2cell...
                            (connectstructinfo.cell),newFields(...
                            1:numel(fieldsConnec)),1);
                    end
                    [startVol,validVol]=ExtractVariables({'startVol',...
                        'validVol'},paramoptim);
                    paramoptim=SetVariables({'startVol'},{validVol},paramoptim);
                    [paramoptim,iterstruct]=InitialiseParamForGrid(baseGrid,...
                        paramoptim);
                    paramoptim=SetVariables({'startVol'},{startVol},paramoptim);
                    
                    isGrid=true;
                    needNewRestartSnake=true;
                    if exist('isOptimStruct','var');
                        repeatflag=isOptimStruct;
                        fieldsIn{ii}='optimstruct';
                    else
                        repeatflag=false;
                    end
                case 'paramoptim'
                    paramoptim=ImposePresets(paramoptim,...
                        optionalin.(fieldsIn{ii}));
                    paramoptim.parametrisation=ImposePresets(paramoptim...
                        .parametrisation,optionalin.(fieldsIn{ii}).parametrisation);
                    paramoptim.initparam=ImposePresets(paramoptim.initparam,...
                        optionalin.(fieldsIn{ii}).initparam);
                    
                    if exist('isSupport','var');
                        repeatflag=isSupport;
                        fieldsIn{ii}='supportOptim';
                    elseif exist('isGrid','var');
                        repeatflag=isGrid;
                        fieldsIn{ii}='grid';
                    else
                        repeatflag=false;
                    end
                case 'supportOptim'
                    paramoptim.optim.supportOptim=optionalin.(fieldsIn{ii});
                    repeatflag=false;
                    isSupport=true;
                case 'restartsnak'
                    restartsnake=optionalin.(fieldsIn{ii});
                    repeatflag=false;
                    isNewRestartSnake=true;
                case 'outinfo'
                    
                    outinfonew=outinfo;
                    outinfo=optionalin.(fieldsIn{ii});
                    outinfo(end+1)=outinfonew;
                    if ~isempty(regexp(outinfo(end-1).rootDir,'_[0-9]{1,2}$', 'once'))
                        numFinStage=regexprep(regexp(outinfo(end-1).rootDir,...
                            '_[0-9]{1,2}$','match'),'_','');
                        for jj=1:str2double(numFinStage{1})
                            outinfo(end).rootDir=[outinfo(end).rootDir,'_',int2str(jj)];
                            outinfo(end).marker=[outinfo(end).marker,'_',int2str(jj)];
                        end
                        copyfile(outinfonew.rootDir,outinfo(end).rootDir)
                    end
                    repeatflag=false;
                otherwise
                    warning([fieldsIn{ii},' is an unhandled restart input'])
                    repeatflag=false;
            end
        end
    end
    
    if exist('needNewRestartSnake','var')
        if exist('isNewRestartSnake','var')
            needNewRestartSnake=~isNewRestartSnake;
        end
        if needNewRestartSnake
%             varExtract={'boundstr'};
%             [boundstr]=ExtractVariables(varExtract,paramoptim.parametrisation);
%             [loop]=GenerateSnakStartLoop(gridrefined,boundstr);
%             [~,~,~,~,restartsnake]=ExecuteSnakes_Optim('snak',gridrefined,loop,...
%                 baseGrid,connectstructinfo,paramoptim.initparam,...
%                 paramoptim.spline,outinfo(end),0,0,0);
%             
            [restartsnake]=ReInitSnake(paramoptim,gridrefined,baseGrid,connectstructinfo,...
                    iterstruct(end),outinfo);
        end
    end
    if isOptimStruct
        if ~CheckIfGradient(ExtractVariables({'optimMethod'},paramoptim))
            nPop=numel(iterstruct(1).population);
            paramoptim=SetVariables({'nPop'},{nPop},paramoptim);
        end
    end
    
end


function param=ImposePresets(param,parampreset)
    
    for ii=1:length(parampreset.structdat.vars)
        param=SetVariables({parampreset.structdat.vars(ii).name},...
            {ExtractVariables({parampreset.structdat.vars(ii).name},parampreset)},param);
    end
    
end

function [optionalin]=LoadRestart(restartPath)
    
    load(restartPath)
    clear restartPath
    c=whos;
    for ii=1:length(c)
        optionalin.(c(ii).name)=eval(c(ii).name);
    end
end

%%  Optimisation Initialisation Blocks

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
    include_GridCheck
    
    
    % Initialise Optimisation
    % Get Parametrisation parameters
    paramoptim=StructOptimParam(caseStr);
    [~,paramoptim]=ConstraintMethod('init',paramoptim,[]); % Constraint setup early
    [paramoptim]=CheckConsistantParamoptim(paramoptim);
    
    [outinfo]=OptimisationOutput('init',paramoptim);
    diaryFile=[outinfo.rootDir,'\Latest_Diary.log'];
    diaryFile=MakePathCompliant(diaryFile);
    fidDiary=fopen(diaryFile,'w');
    fclose(fidDiary);
    diary(diaryFile);
    
    
    % Initialise Grid
    [unstrGrid,baseGrid,gridrefined,connectstructinfo,unstrRef,loop]...
        =GridInitAndRefine(paramoptim.parametrisation);
    % Update Parameters using grid information
    [paramoptim,iterstruct]=InitialiseParamForGrid(baseGrid,paramoptim);
    
    iterstruct(1).population=ApplySymmetry(paramoptim,iterstruct(1).population);
    % Start Parallel Pool
    [nWorker,nodeShareNum,node1ReserveNum]=ExtractVariables({'worker','nodeShareNum','node1ReserveNum'},paramoptim);
    workerList=StartParallelPool(nWorker,10);
    [machineList]=WriteMachineFiles(nWorker*nodeShareNum,outinfo.rootDir,node1ReserveNum);
    paramoptim=SetVariables({'workerList','machineList'},{workerList,machineList},paramoptim);
    
    if ExtractVariables({'useSnake'},paramoptim)
        
        [~,~,~,~,restartsnake]=ExecuteSnakes_Optim('snak',gridrefined,loop,...
            baseGrid,connectstructinfo,paramoptim.initparam,...
            paramoptim.spline,outinfo,0,0,0);
        
    else
        restartsnake=struct([]);
    end
    
    [outinfo]=OptimisationOutput('iteration',...
        paramoptim,0,outinfo,iterstruct(1),{});
    
    [~,paramoptim]=ConstraintMethod('init',paramoptim,iterstruct,baseGrid); % constraint set up post grid
    [~]=PrintEnd(procStr,1,tStart);
end

function [restartsnake]=ReInitSnake(paramoptim,gridrefined,baseGrid,connectstructinfo,...
        iterstruct,outinfo)
    
    rmdir([outinfo(end).rootDir,filesep,'iteration_0'], 's')
    [pltFile]=FindDir(outinfo(end).rootDir,'Tec360PLT',0);
    delete(pltFile{1})
    varExtract={'boundstr'};
    [boundstr]=ExtractVariables(varExtract,paramoptim.parametrisation);
    varExtract={'nPop'};
    [nPop]=ExtractVariables(varExtract,paramoptim);
    
    [gridrefined]=EdgePropertiesReshape(gridrefined);
    [baseGrid]=EdgePropertiesReshape(baseGrid);
    [newrestartsnake]=GenerateEdgeLoop(gridrefined,boundstr,true);
    
    [~,~,~,~,restartsnake]=ExecuteSnakes_Optim('snak',gridrefined,newrestartsnake,...
        baseGrid,connectstructinfo,paramoptim.initparam,...
        paramoptim.spline,outinfo(end),0,0,0);
    
    
    
    [~]=OptimisationOutput('iteration',...
        paramoptim,0,outinfo(end),iterstruct(1),{});
    
end


function [paramoptim]=CheckConsistantParamoptim(paramoptim)
    varExtract={'cellLevels','passDomBounds','corneractive','regulariseSnakeBounds'};
    [cellLevels,passDomBounds,corneractive]=ExtractVariables(varExtract(1:3),...
        paramoptim.parametrisation);
    [regulariseSnakeBounds]=ExtractVariables(varExtract(4),paramoptim);
    
    switch regulariseSnakeBounds
        case ''
        case '2c_orig'
            [passDomBounds]=MakeCartesianGridBounds(cellLevels);
        otherwise
            warning('Unrecognised regulariseSnakeBounds case')
    end
    paramoptim.parametrisation=SetVariables(varExtract(1:3),...
        {cellLevels,passDomBounds,corneractive},paramoptim.parametrisation);
    
end

function [paramoptim,iterstruct]=InitialiseParamForGrid(baseGrid,paramoptim)
    [desvarconnec,~,~]=ExtractVolumeCellConnectivity(baseGrid);
    [desvarconnec]=ExtractDesignVariableConnectivity(baseGrid,desvarconnec);
    paramoptim=SetVariables({'desvarconnec'},{desvarconnec},paramoptim);
    [paramoptim]=OptimisationParametersModif(paramoptim,baseGrid);
    [iterstruct,paramoptim]=InitialisePopulation(paramoptim,baseGrid);
    
    knownOptim=ExtractVariables({'knownOptim'},paramoptim);
    paramoptim=SetVariables({'knownOptimStart'},{knownOptim},paramoptim);
    
end

function [workerList]=StartParallelPool(nWorker,nTry)
    kk=0;
    while numel(gcp('nocreate'))==0 && kk<nTry
        comStr=computer;
        if strcmp(comStr(1:2),'PC')
            poolName=parallel.importProfile('ExportOptimSnakes.settings');
        else
            poolName=parallel.importProfile('ExportOptimSnakesLinux.settings');
        end
        clusterObj=parcluster(poolName);
        clusterObj.NumWorkers=nWorker;
        saveProfile(clusterObj);
        
        try
            p=parpool(poolName);
            p.IdleTimeout=Inf;
            
        catch ME
            
        end
        kk=kk+1;
    end
    
    if numel(gcp('nocreate'))~=0
        disp(['Parallel Pool succesfully created after ',int2str(kk),' attempt(s)'])
    else
        %error('Parrallel pool failed to start')
        throw(ME)
    end
    ll=1;
    workerList=[];
    while numel(workerList)<nWorker && ll<10
        thisworker=zeros([1 ll*nWorker]);
        parfor ii=1:ll*nWorker
            tempworker = getCurrentWorker;
            thisworker(ii)=tempworker.ProcessId;
        end
        workerList=unique(thisworker);
        ll=ll+1;
    end
    if numel(workerList)<nWorker && ~strcmp(computer,'PCWIN64')
        error('Failed to recover a list of workers')
    end
    pctRunOnAll ExecInclude
    pathStr=cd;
    parfor ii=1:ll*nWorker
        cd(pathStr);
    end
end

function [paramoptim]=InitialiseObjective(paramoptim)
    varExtract={'objectiveName'};
    [objectiveName]=ExtractVariables(varExtract,paramoptim);
    
    switch objectiveName
        case 'CutCellFlow'
            
        otherwise
            disp('Objective requires no initilisation')
            
    end
    
    
end

%% Iteration and convergence

function [population,restartsnake]=PerformIteration(paramoptim,outinfo,nIter,population,...
        gridrefined,restartsnake,baseGrid,connectstructinfo,iterstruct)
    
    
    varExtract={'nPop','objectiveName','defaultVal','lineSearch','useSnake'};
    [nPop,objectiveName,defaultVal,lineSearch,useSnake]=ExtractVariables(varExtract,paramoptim);
    
    isSensiv=~((~CheckSnakeSensitivityAlgorithm(paramoptim)) || lineSearch);
    
    if ~useSnake
        [population,supportstruct,captureErrors]=IterateNoSnake(paramoptim,population,baseGrid);
        
    elseif ~isSensiv
        [population,supportstruct,captureErrors]=IterateNoSensitivity(paramoptim,outinfo,nIter,population,...
            gridrefined,restartsnake,baseGrid,connectstructinfo);
        
    else
        %[population,supportstruct,captureErrors,restartsnake]=IterateSensitivity(paramoptim,outinfo,nIter,population,...
        [population,supportstruct,captureErrors,~]=IterateSensitivity(paramoptim,outinfo,nIter,population,...
            gridrefined,restartsnake,baseGrid,connectstructinfo);
    end
    nPop=numel(population);
    [paramoptim]=SetVariables({'nPop'},{nPop},paramoptim);
    
    %save('TestParamMeshMot.mat')
    
    [supportstruct,captureErrors,population]=ParamMeshMotion(paramoptim,population,...
        supportstruct,isSensiv,captureErrors,iterstruct);
    
    [population,captureErrors]=ParallelObjectiveCalc...
        (objectiveName,paramoptim,population,supportstruct,baseGrid,captureErrors);
    
    [population]=ConstraintMethod('Res',paramoptim,population);
    population=EnforceConstraintViolation(population,defaultVal);
    
    [outinfo]=OptimisationOutput('iteration',paramoptim,nIter,outinfo,population,captureErrors);
    
end

function [population,captureErrors]=ParallelObjectiveCalc...
        (objectiveName,paramoptim,population,supportstruct,baseGrid,captureErrors)
    
    nPop=numel(population);
    
    parfor ii=1:nPop %
        %for ii=1:nPop
        try
            [population(ii).objective,additional]=...
                EvaluateObjective(objectiveName,paramoptim,population(ii),...
                supportstruct(ii),baseGrid);
            % Overwrite fill if filldelta exists
            fieldsAdd=fieldnames(additional);
            filldeltaInd=find(~cellfun(@isempty,regexp(fieldsAdd,'filldelta')));
            if ~isempty(filldeltaInd)
                population(ii).fill=population(ii).fill+additional.filldelta;
            end
            % Additional Data to be saved
            fieldsAdd(filldeltaInd)=[];
            for jj=1:numel(fieldsAdd)
                population(ii).additional.(fieldsAdd{jj})=additional.(fieldsAdd{jj});
            end
        catch MEexception
            population(ii).constraint=false;
            population(ii).exception=[population(ii).exception,'error: ',MEexception.identifier ,' - '];
            captureErrors{ii}=[captureErrors{ii},MEexception.getReport];
        end
    end
    
end


function [population,captureErrors]=SerialisedParallelObjectiveCalc...
        (objectiveName,paramoptim,population,supportstruct,captureErrors)
    
    nPop=numel(population);
    nRunning=0;
    objCalcFinished=false;
    statusstruct=repmat(struct('state','ready','PID',[],'stage','init',...
        'repeat',1,'dat',struct([])),[1 nPop]);
    
    while ~objCalcFinished
        for ii=1:nWorker-nRunning %
            
            try
                [population(ii).objective,additional,statusstruct(ii)]=...
                    EvaluateObjectiveSerial(objectiveName,paramoptim,population(ii),...
                    supportstruct(ii),statusstruct(ii));
                if strcmp(statusstruct(ii).state,'finished')
                    fieldsAdd=fieldnames(additional);
                    for jj=1:numel(fieldsAdd)
                        population(ii).additional.(fieldsAdd{jj})=additional.(fieldsAdd{jj});
                    end
                end
            catch MEexception
                population(ii).constraint=false;
                statusstruct(ii).state='error';
                population(ii).exception=[population(ii).exception,'error: ',MEexception.identifier ,' - '];
                captureErrors{ii}=[captureErrors{ii},MEexception.getReport];
            end
        end
    end
end

function [isConv]=ConvergenceTest_static(paramoptim,iterstruct,nIter,startIter)
    comEps=@(f1,f2,eps) all(abs(f1-f2)<eps);
    varExtract={'optimMethod','iterGap'};
    [optimMethod,iterGap]=ExtractVariables(varExtract,paramoptim);
    isConv=false;
    eps=1e-10;
    if CheckIfGradient(optimMethod) && ((nIter-startIter)>(4*iterGap+1))
        
        if numel(iterstruct(nIter).population(1).fill)==numel(iterstruct(nIter-iterGap).population(1).fill) ...
                && numel(iterstruct(nIter).population(1).fill)==numel(iterstruct(nIter-2*iterGap).population(1).fill)...
                && numel(iterstruct(nIter).population(1).fill)==numel(iterstruct(nIter-3*iterGap).population(1).fill)...
                && numel(iterstruct(nIter).population(1).fill)==numel(iterstruct(nIter-4*iterGap).population(1).fill)
            isConv=comEps(iterstruct(nIter).population(1).fill,iterstruct(nIter-iterGap).population(1).fill,eps) ...
                && comEps(iterstruct(nIter).population(1).fill,iterstruct(nIter-2*iterGap).population(1).fill,eps)...
                && comEps(iterstruct(nIter).population(1).fill,iterstruct(nIter-3*iterGap).population(1).fill,eps)...
                && comEps(iterstruct(nIter).population(1).fill,iterstruct(nIter-4*iterGap).population(1).fill,eps);
        end
        
    else
        
        isConv=false; % Need to find a convergence condition for non grad optim
        
    end
    
end

function [isConv]=ConvergenceTest_sloperefine(paramoptim,iterstruct,nIter,startIter)
    comEps=@(f1,f2,eps) all(abs(f1-f2)<eps);
    varExtract={'optimMethod','iterGap','knownOptimStart','slopeConv'};
    [optimMethod,iterGap,knownOptimStart,slopeConv]=ExtractVariables(varExtract,paramoptim);
    mmaSpan=3;
    
    isConv=false;
    if CheckIfGradient(optimMethod) && ((nIter-startIter)>=(mmaSpan*iterGap))
        
        %objDat=zeros([1,nIter-startIter+1]);
        
        refMin=[];
        for ii=flip(1:startIter)
            refMin=min([iterstruct(ii).population(:).objective,refMin]);
        end
        
        kk=1;
        for ii=flip(nIter:-iterGap:startIter)
            objDat(kk)=min([iterstruct(ii).population(:).objective]);kk=kk+1;
        end
        objDat=min(objDat,refMin);
        if knownOptimStart~=0
            mma=MovingAverage(objDat,mmaSpan);
            mmaSlope=abs(mma(2:end)-mma(1:end-1));
            isConv=(max(mmaSlope(end-(mmaSpan-1):end))/max(mmaSlope))<slopeConv;
        else
            mma=MovingAverage(log10(objDat),mmaSpan);
            mmaSlope=abs(mma(2:end)-mma(1:end-1));
            isConv=(max(mmaSlope(end-(mmaSpan-1):end))/max(mmaSlope))<slopeConv;
        end
    else
        
        isConv=false; % Need to find a convergence condition for non grad optim
        
    end
    
end

function [population,supportstruct,captureErrors]=IterateNoSensitivity(paramoptim,outinfo,nIter,population,...
        gridrefined,restartsnake,baseGrid,connectstructinfo)
    
    varExtract={'nPop'};
    [nPop]=ExtractVariables(varExtract,paramoptim);
    
    paramsnake=paramoptim.parametrisation;
    paramspline=paramoptim.spline;
    [population]=ConstraintMethod('DesVar',paramoptim,population,baseGrid);
    
    [captureErrors{1:nPop}]=deal('');
    supportstruct=repmat(struct('loop',[]),[1,nPop]);
    
    parfor ii=1:nPop
        %for ii=1:nPop
        %for ii=flip(1:nPop)
        worker = getCurrentWorker
        currentMember=population(ii).fill;
        [newGrid,newRefGrid,newrestartsnake]=ReFillGrids(baseGrid,gridrefined,...
            restartsnake,connectstructinfo,currentMember);
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
    [supportstruct(:).parentMesh]=deal('');
end

function [population,supportstruct,captureErrors]=IterateNoSnake(paramoptim,population,baseGrid)
    
    varExtract={'nPop'};
    [nPop]=ExtractVariables(varExtract,paramoptim);
    
    [population]=ConstraintMethod('DesVar',paramoptim,population,baseGrid);
    
    [captureErrors{1:nPop}]=deal('');
    supportstruct=repmat(struct('loop',[],'parentMesh',''),[1,nPop]);
    
    
end

%% Support for  mesh motion

function [supportstruct,captureErrors,population]=ParamMeshMotion(paramoptim,population,...
        supportstruct,isSensiv,captureErrors,iterstruct)
    % DVP Find folders and paths required for mesh motion
    % Also need separate function to Output the surface and displacements
    % single load operation needed to get the root loop that generated the
    % mesh
    
    
    varExtract={'meshDefSens','meshDefNorm','objectiveName','useSnake','rootMesh'};
    [meshDefSens,meshDefNorm,objectiveName,useSnake,rootMesh]=...
        ExtractVariables(varExtract,paramoptim);
    
    if (isSensiv && meshDefSens) || (~isSensiv && meshDefNorm) && useSnake
        % Execute the method for selection
        [supportstruct]=SelectPrevMesh(paramoptim,supportstruct,population,iterstruct);
        % Generate the profiles
        [captureErrors,supportstruct,population]=GenerateSurfaceMotionFiles...
            (paramoptim,population,supportstruct,captureErrors);
        
    end
    
end

function [supportstruct]=SelectPrevMesh(paramoptim,supportstruct,population,iterstruct)
    
    varExtract={'optimMethod','rootMesh','lineSearch','direction'};
    [optimMethod,rootMesh,lineSearch,direction]=ExtractVariables(varExtract,paramoptim);
    
    switch rootMesh{1}
        case 'path'
            [supportstruct(:).parentMesh]=deal(rootMesh{2});
            
        case 'previter'
            [prevMesh]=FindPrevIterMesh(optimMethod,direction,lineSearch,iterstruct,rootMesh);
            [supportstruct(:).parentMesh]=deal(prevMesh);
        case 'none'
            [supportstruct(:).parentMesh]=deal('');
        otherwise
            error('Invalid mesh motion specified')
    end
    
end

function [captureErrors,supportstruct,population]=GenerateSurfaceMotionFiles...
        (paramoptim,population,supportstruct,captureErrors)
    % DVP Generates the displacements and the surface.xyz and
    % displacements.xyz files for the files concerned.
    %
    
    varExtract={'typeLoop'};
    [typeLoop]=...
        ExtractVariables(varExtract,paramoptim.parametrisation);
    
    for ii=1:numel(supportstruct)
        try
            if ~isempty(supportstruct(ii).parentMesh)
                if ii==1 || ~strcmp(supportstruct(ii).parentMesh,...
                        supportstruct(max(ii-1,1)).parentMesh)
                    [rootProfile]=FindDir(supportstruct(ii).parentMesh,'restart',false);
                    indat=load(rootProfile{1},'loop');
                end
                looproot=indat.loop;
                loop=supportstruct(ii).loop;
                if numel(loop)~=numel(looproot)
                    fprintf('Profile %i : Different Topologies, re-meshing required\n',ii);
                    supportstruct(ii).parentMesh='';
                else
                    [catloop]=MatchLoops(looproot,loop,typeLoop);
                    fidSurf=fopen([population(ii).location,filesep,'surface.xyz'],'w');
                    fidDisp=fopen([population(ii).location,filesep,'displacements.xyz'],'w');
                    
                    DisplacementsOutput(catloop,fidSurf,fidDisp,typeLoop,false);
                    fclose(fidSurf);
                    fclose(fidDisp);
                end
                
            end
        catch MEexception
            %population(ii).constraint=false;
            population(ii).exception=[population(ii).exception,'Warning: ',MEexception.identifier ,' - '];
            captureErrors{ii}=[captureErrors{ii},MEexception.getReport];
            supportstruct(ii).parentMesh='';
        end
    end
end

function [catloop]=MatchLoops(loop1,loop2,typeLoop)
    % loop1 is the reference loop
    % loop2 is the one that must be transformed in displacements
    
    [nel1]=(cellfun(@numel,{loop1(:).(typeLoop)}));
    nel1(2,:)=[loop1(:).isccw];
    %nel1(3,:)=[loop1(:).isinternal];
    [nel2]=(cellfun(@numel,{loop2(:).(typeLoop)}));
    nel2(2,:)=[loop2(:).isccw];
    %nel2(3,:)=[loop2(:).isinternal];
    
    for ii=1:2
        [~,orderL1]=sort(nel1(ii,:));
        nel1=nel1(:,orderL1);
        [~,orderL2]=sort(nel2(ii,:));
        nel2=nel2(:,orderL2);
        loop1=loop1(orderL1);
        loop2=loop2(orderL2);
    end
    
    catloop=repmat(struct(typeLoop,[],'isccw',false),[numel(loop1),3]);
    if size(nel1,2)==size(nel2,2)
        for ii=1:numel(loop1)
            catloop(ii,1).(typeLoop)=loop1(ii).(typeLoop);
            catloop(ii,2).(typeLoop)=loop2(ii).(typeLoop);
            catloop(ii,1).isccw=loop1(ii).isccw;
            catloop(ii,2).isccw=loop2(ii).isccw;
        end
        for ii=1:size(nel1,2)
            if all(nel1(1,:)==nel2(1,:)) ...
                    && ~xor(catloop(ii,1).isccw,catloop(ii,2).isccw)
                % if they are all the same assume the loops combien
                % directly
                
                catloop(ii,3).(typeLoop)=catloop(ii,2).(typeLoop)-catloop(ii,1).(typeLoop);
                catloop(ii,3).isccw=catloop(ii,1).isccw;
            elseif all(nel1(1,:)==nel2(1,:)) ...
                    && xor(catloop(ii,1).isccw,catloop(ii,2).isccw)
                % if they are all the same assume the loops combine
                % directly
                
                catloop(ii,3).(typeLoop)=flip(catloop(ii,2).(typeLoop))-catloop(ii,1).(typeLoop);
                catloop(ii,3).isccw=catloop(ii,1).isccw;
                warning('support for misalligned loops not included yet')
            else
                % The loops do not all match
                error('Non matching loops provided, robust process must be coded')
            end
        end
        
    else
        disp('Different topologies, differences cannot be assigned')
    end
    
    
end

function [prevMesh]=FindPrevIterMesh(objectiveName,direction,lineSearch,iterstruct,rootMesh)
    if numel(iterstruct)>0
        isGradient=CheckIfGradient(objectiveName);
        if isGradient
            if lineSearch
                prevMesh=iterstruct(end).population(1).location;
            else
                obj=[iterstruct(end).population(:).objective];
                isObjEmpty=find(~cellfun(@isempty,{iterstruct(end).population(:).objective}));
                switch direction
                    case 'min'
                        [~,minIterPos]=min(obj);
                    case 'max'
                        [~,minIterPos]=max(obj);
                        
                end
                prevMesh=iterstruct(end).population(isObjEmpty(minIterPos)).location;
            end
        else
            error('Global optimisation cannot be selected for previous iteration')
        end
        
    else
        prevMesh=rootMesh{2};
        
    end
end
%% Gradient Iteration
function [population,supportstruct,captureErrors,restartsnake]=IterateSensitivity(paramoptim,outinfo,nIter,population,...
        gridrefined,restartsnake,baseGrid,connectstructinfo)
    
    
    paramsnake=paramoptim.parametrisation;
    paramspline=paramoptim.spline;
    
    [population]=ConstraintMethod('DesVar',paramoptim,population,baseGrid);
    
    [population,supportstruct,restartsnake,paramsnake,paramoptim,captureErrors]=ComputeRootSensitivityPopFill...
        (paramsnake,paramspline,paramoptim,population,baseGrid,gridrefined,...
        restartsnake,connectstructinfo,outinfo,nIter);
    
    population=ApplySymmetry(paramoptim,population);
    [population]=ConstraintMethod('DesVar',paramoptim,population,baseGrid);
    
    varExtract={'nPop','sensCalc'};
    [nPop,sensCalc]=ExtractVariables(varExtract,paramoptim);
    [captureErrors{2:nPop}]=deal('');
    
    supportstruct=[supportstruct,repmat(struct('loop',[]),[1,nPop-1])];
    
    supportstructsens=repmat(struct('loop',[],'rootMesh','','loopsens',[],'volumefraction',[]),[1,nPop-1]);
    
    if strcmp('analytical',sensCalc)
        
        [paramoptim,paramsnake]=HandleNonFillVar(population(1),paramoptim,paramsnake);
        [supportstructsens,population]=ComputeRootSensitivityPopProfiles(paramsnake,paramoptim,...
            population,baseGrid,gridrefined,restartsnake,connectstructinfo,supportstructsens);
        
    end
    
    rootPop=population(1);
    
    parfor ii=1:nPop-1
        %for ii=1:nPop-1
        currentMember=population(ii+1).fill;
        [newGrid,newRefGrid,newrestartsnake]=ReFillGrids(baseGrid,gridrefined,...
            restartsnake,connectstructinfo,currentMember);
        try
            % Normal Execution
            % DVP check for output type and if sens need to make sure the
            % comparison profile is passed for reference.
            % Potentially develop a code that if two profiles have
            % different number of points it resamples using an intersection
            % if there are no intersections it resamples both profiles
            % starting from the closest points.
            % The main question is do I go and collect the data from the
            % grid source or pass it around.
            switch sensCalc
                case'snake'
                    [population(ii+1),supportstruct(ii+1)]=NormalExecutionIteration(...
                        population(ii+1),newRefGrid,newrestartsnake,newGrid,...
                        connectstructinfo,paramsnake,paramspline,outinfo,nIter,ii+1,paramoptim);
                case 'analytical'
                    [population(ii+1),supportstruct(ii+1)]=SensitivityExecutionIteration(...
                        population(ii+1),newRefGrid,supportstructsens(ii),newGrid,...
                        connectstructinfo,paramsnake,paramspline,outinfo,...
                        nIter,ii+1,paramoptim,rootPop);
            end
            
        catch MEexception
            
            population(ii+1).constraint=false;
            population(ii+1).exception=['error: ',MEexception.identifier, ' - '];
            captureErrors{ii+1}=MEexception.getReport;
        end
    end
    [supportstruct(:).parentMesh]=deal('');
    
end

function [newpopulation,supportstruct,restartsnake,paramsnake,paramoptim,captureErrors]=...
        ComputeRootSensitivityPopFill(paramsnake,paramspline,paramoptim,...
        population,baseGrid,gridrefined,restartsnake,connectstructinfo,outinfo,nIter)
    
    captureErrors{1}='';
    nDesVar=ExtractVariables({'nDesVar'},paramoptim);
    try
        % Compute root member profile
        currentMember=population(1).fill;
        [newGrid,newRefGrid,newrestartsnake]=ReFillGrids(baseGrid,gridrefined,restartsnake,connectstructinfo,currentMember);
        
        [population(1),supportstruct,restartsnake]=NormalExecutionIteration(population(1),newRefGrid,...
            newrestartsnake,newGrid,connectstructinfo,paramsnake,paramspline...
            ,outinfo,nIter,1,paramoptim);
        
        % Compute sensitivity
        
        [paramoptim,paramsnake]=HandleNonFillVar(population(1),paramoptim,paramsnake);
        newfillstruct=SnakesSensitivity('fill',newRefGrid,restartsnake,...
            newGrid,paramsnake,paramoptim,currentMember);
        
        keepPop=false(numel(population));
        for ii=1:numel(population)
            keepPop(ii)=any(population(ii).optimdat.var>nDesVar);
        end
        keepPop=find(keepPop);
        keepPop(keepPop==1)=[];
        nPop=numel(newfillstruct)+1+numel(keepPop);
        [paramoptim]=SetVariables({'nPop'},{nPop},paramoptim);
        [newpopulation]=GeneratePopulationStruct(paramoptim);
        newpopulation(1)=population(1);
        for ii=2:numel(newfillstruct)+1
            newpopulation(ii).fill=newfillstruct(ii-1).fill;
            newpopulation(ii).nonfillvar=newpopulation(1).nonfillvar;
            newpopulation(ii).optimdat.var=newfillstruct(ii-1).optimdat.var;
            newpopulation(ii).optimdat.value=newfillstruct(ii-1).optimdat.val;
        end
        nPrc=numel(newfillstruct)+1;
        for ii=1:numel(keepPop)
            newpopulation(nPrc+ii).fill=newpopulation(1).fill;
            newpopulation(nPrc+ii).nonfillvar=population(keepPop(ii)).nonfillvar;
            newpopulation(nPrc+ii).optimdat.var=population(keepPop(ii)).optimdat.var;
            newpopulation(nPrc+ii).optimdat.value=population(keepPop(ii)).optimdat.value;
        end
        
        varExtract={'restart'};
        [paramsnake]=SetVariables(varExtract,{true},paramsnake);
    catch MEexception
        population(1).constraint=false;
        population(1).exception=['error: ',MEexception.identifier];
        captureErrors{1}=MEexception.getReport;
        disp(captureErrors{1})
        newpopulation=population;
        supportstruct.loop=[];
        disp(MEexception.getReport)
        warning('Sensitivity Extraction has failed, basis will not be smoothed.')
    end
    
    
end


function [supportstruct,population]=...
        ComputeRootSensitivityPopProfiles(paramsnake,paramoptim,...
        population,baseGrid,gridrefined,restartsnake,connectstructinfo,supportstruct)
    procStr=['Compute new profiles from sensitivity'];
    [tStart]=PrintStart(procStr,2);
    
    
    % Compute root member profile
    currentMember=population(1).fill;
    [newGrid,newRefGrid,~]=ReFillGrids(baseGrid,gridrefined,restartsnake,connectstructinfo,currentMember);
    
    % Compute sensitivity
    newProfileLoops=SnakesSensitivity('profile',newRefGrid,restartsnake,...
        newGrid,paramsnake,paramoptim,currentMember,population(2:end));
    
    
    [supportstruct(:).loopsens]=deal(newProfileLoops(:).loopsens);
    [supportstruct(:).volumefraction]=deal(newProfileLoops(:).volumefraction);
    varExtract={'modeSmoothType'};
    [modeSmoothType]=ExtractVariables(varExtract,paramsnake);
    if strcmp(modeSmoothType,'optimlinprog')
        [population(2:numel(newProfileLoops)+1).fill]=deal(newProfileLoops(:).fill);
        [population(2:numel(newProfileLoops)+1).optimdat]=deal(newProfileLoops(:).optimdat);
    end
    
    
    [~]=PrintEnd(procStr,2,tStart);
end

function [population,supportstruct,restartsnake]=NormalExecutionIteration(population,newRefGrid,...
        newrestartsnake,newGrid,connectstructinfo,paramsnake,paramspline...
        ,outinfo,nIter,ii,paramoptim)
    
    [paramoptim,paramsnake]=HandleNonFillVar(population,paramoptim,paramsnake);
    
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
        ExecuteSnakes_Optim('snak',newRefGrid,newrestartsnake,...
        newGrid,connectstructinfo,paramsnake,paramspline,outinfo,nIter,ii,nPop);
    population.location=outTemp.dirprofile;
    population.additional.snaxelVolRes=snakSave(end).currentConvVolume;
    population.additional.snaxelVelResV=snakSave(end).currentConvVelocity;
    
    supportstruct.loop=loop;
end


function [population,supportstruct,restartsnake]=PostExecutionIteration(population,newRefGrid,...
        newrestartsnake,newGrid,connectstructinfo,paramsnake,paramspline...
        ,outinfo,nIter,ii,paramoptim)
    
    [paramoptim,paramsnake]=HandleNonFillVar(population,paramoptim,paramsnake);
    
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
        ExecuteSnakes_Optim('post',newRefGrid,newrestartsnake,...
        newGrid,connectstructinfo,paramsnake,paramspline,outinfo,nIter,ii,nPop);
    population.location=outTemp.dirprofile;
    population.additional.snaxelVolRes=snakSave(end).currentConvVolume;
    population.additional.snaxelVelResV=snakSave(end).currentConvVelocity;
    
    supportstruct.loop=loop;
    supportstruct.parentMesh='';
end

function [population,supportstruct,restartsnake]=SensitivityExecutionIteration(population,newRefGrid,...
        supportstructsens,newGrid,connectstructinfo,paramsnake,paramspline...
        ,outinfo,nIter,ii,paramoptim,rootmember)
    
    varExtract={'nPop'};
    [nPop]=ExtractVariables(varExtract,paramoptim);
    varExtract={'axisRatio'};
    [axisRatio]=ExtractVariables(varExtract,paramsnake);
    [paramoptim,paramsnake]=HandleNonFillVar(population,paramoptim,paramsnake);
    varExtract={'axisRatio'};
    [axisRatioNew]=ExtractVariables(varExtract,paramsnake);
    axisRatioNew=axisRatioNew/axisRatio;
    [paramsnake]=SetVariables(varExtract,{axisRatioNew},paramsnake);
    
    
    [~,~,~,loop,restartsnake,outTemp]=...
        ExecuteSnakes_Optim('sens',newRefGrid,supportstructsens,...
        newGrid,connectstructinfo,paramsnake,paramspline,outinfo,nIter,ii,nPop);
    
    population.location=outTemp.dirprofile;
    population.additional.snaxelVolRes=rootmember.additional.snaxelVolRes;
    population.additional.snaxelVelResV=rootmember.additional.snaxelVelResV;
    
    supportstruct.loop=loop;
end

function population=EnforceConstraintViolation(population,defaultVal)
    % constraint is now something that will go from 1 (no violation) to 0
    % full violation numbers in between will be used as a ratio of the
    % default value that must be added.
    
    
    for ii=1:length(population)
        isConstraint=1-[population(ii).constraint];
        if isempty(isConstraint)
            isConstraint=1;
        end
        if ~isempty(population(ii).objective)
            
            [population(ii).objective]=population(ii).objective+...
                isConstraint*(defaultVal);
            
        else
            [population(ii).objective]=isConstraint*(defaultVal);
        end
        
    end
    
end

function [iterstruct,paramoptim]...
        =GenerateNewPop(paramoptim,iterstruct,nIter,iterStart,baseGrid)
    procStr=['Generate New Population'];
    [tStart]=PrintStart(procStr,2);
    
    varExtract={'nPop','iterGap','optimMethod','nDesVar'};
    [nPop,iterGap,optimMethod,nDesVar]=ExtractVariables(varExtract,paramoptim);
    
    [isGradient]=CheckIfGradient(optimMethod);
    iterCurr=RescaleDesVarNoFill('tofill',paramoptim,iterstruct(nIter).population);
    iterM1=RescaleDesVarNoFill('tofill',paramoptim,iterstruct(max([nIter-iterGap,...
        iterStart])).population);
    [newPop,iterstruct(nIter).population,paramoptim,deltas]=OptimisationMethod(paramoptim,...
        iterCurr,iterM1,baseGrid);
    [iterstruct(nIter).population]=RescaleDesVarNoFill('tovar',...
        paramoptim,iterstruct(nIter).population);
    
    varExtract={'nPop'};
    [nPop]=ExtractVariables(varExtract,paramoptim);
    [iterstruct(nIter+1).population]=GeneratePopulationStruct(paramoptim);
    
    nFill=nDesVar;
    if ~isGradient
        for ii=1:nPop
            iterstruct(nIter+1).population(ii).fill=newPop(ii,1:nFill);
            iterstruct(nIter+1).population(ii).nonfillvar=newPop(ii,nFill+1:end);
        end
    else
        for ii=1:nPop
            iterstruct(nIter+1).population(ii).fill=newPop(ii,1:nFill);
            iterstruct(nIter+1).population(ii).nonfillvar=newPop(ii,nFill+1:end);
            iterstruct(nIter+1).population(ii).optimdat.var=deltas{ii}(1,:);
            iterstruct(nIter+1).population(ii).optimdat.value=deltas{ii}(2,:);
            
        end
    end
    [iterstruct(nIter+1).population]=RescaleDesVarNoFill('tovar',...
        paramoptim,iterstruct(nIter+1).population);
    iterstruct(nIter+1).population=ApplySymmetry(paramoptim,iterstruct(nIter+1).population);
    [~]=PrintEnd(procStr,2,tStart);
end



%% Parametrisation Interface

function [unstructured,unstructReshape,gridrefined,connectstructinfo,...
        unstructuredrefined,loop]=GridInitAndRefine(param,unstructReshape)
    % Executes the Grid Initialisation process
    procStr='INITIAL GRID OPERATIONS';
    [tStart]=PrintStart(procStr,2);
    
    if nargin<2
        [unstructured,~,unstructReshape]=GridInitialisationV2(param);
    else
        [unstructured]=ModifReshape(unstructReshape);
    end
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
    
    varExtract={'symType','useSnake','startVol','validVol','diffStepSize',...
        'maxDiffStep','minDiffStep','gradScaleType','nDesVar','numNonFillVar'};
    [symType,useSnake,startVol,validVol,diffStepSize,maxDiffStep,minDiffStep...
        ,gradScaleType,nDesVar,numNonFillVar]=ExtractVariables(varExtract,paramoptim);
    varExtract={'cellLevels','corneractive'};
    [cellLevels,corneractive]=ExtractVariables(varExtract,paramoptim.parametrisation);
    if useSnake
        nDesVar=sum([baseGrid.cell(:).isactive]);
        paramoptim=SetVariables({'nDesVar'},{nDesVar},paramoptim);
        paramoptim=SetVariables({'gradScale'},{[BuildCellRatios(baseGrid,gradScaleType),...
            ones([1,sum(numNonFillVar)])]},paramoptim);
    else
        paramoptim=SetVariables({'gradScale'},{ones(1,nDesVar+sum(numNonFillVar))}...
            ,paramoptim);
    end
    if ~isempty(startVol) && startVol~=0
        validVol=max(validVol,startVol);
    end
    if isempty(maxDiffStep) || maxDiffStep==0
        maxDiffStep=max([abs(diffStepSize),minDiffStep]);
    end
    diffStepSize=diffStepSize/max(abs(diffStepSize))*maxDiffStep;
    paramoptim=SetVariables({'maxDiffStep','validVol','diffStepSize'},...
        {maxDiffStep,validVol,diffStepSize},paramoptim);
    
    symDesVarList=BuildSymmetryLists(symType,cellLevels,corneractive,baseGrid);
    notDesInd=BuildExclusionList(symDesVarList);
    paramoptim=SetVariables({'notDesInd','symDesVarList'},...
        {notDesInd,symDesVarList},paramoptim);
    [paramoptim]=CheckiterGap(paramoptim);
    
    
    % Set resampleSnak to match iin both with precedence in the
    % optimisation option
    paramoptim.parametrisation=SetVariables({'resampleSnak'},...
        {ExtractVariables({'resampleSnak'},paramoptim)},paramoptim.parametrisation);
    
    
end

function [gradScale]=BuildCellRatios(baseGrid,gradScaleType)
    
    cellCentredInf=CellCentreGridInformation((baseGrid));
    isAct=logical([cellCentredInf(:).isactive]);
    volFill=[cellCentredInf(isAct).volume];
    gradScale=ones(size(volFill));
    switch gradScaleType
        case 'volume'
            
            % Check for negative Volume
            if any(sign(volFill)==-1)
                warning('Negative Volumes in grid data')
                volFill=abs(volFill);
            end
            
            gradScale=1./volFill;
            gradScale=gradScale/mean(gradScale);
            
        otherwise
            disp('Gradient scaling inactive')
            
    end
    
    
end

function [paramoptim]=CheckiterGap(paramoptim)
    varExtract={'optimMethod'};
    [optimMethod]=ExtractVariables(varExtract,paramoptim);
    
    switch optimMethod
        case 'conjgradls'
            paramoptim=SetVariables({'iterGap'},{2},paramoptim);
        case 'conjgrad'
            paramoptim=SetVariables({'iterGap'},{2},paramoptim);
        otherwise
            paramoptim=SetVariables({'iterGap'},{1},paramoptim);
    end
end

function [symDesVarList]=BuildSymmetryLists(symType,cellLevels,corneractive,baseGrid)
    
    switch symType
        case 'none'
            symDesVarList=zeros([2,0]);
            
        case 'horz'
            [symDesVarList]=RobustSymmetryAssignement(baseGrid,2);
        case 'horzold'
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
            
            [symDesVarList]=RobustSymmetryAssignement(baseGrid,1);
            
    end
    
end

function [symList]=RobustSymmetryAssignement(baseGrid,dim)
    % finds symmetric sites by matching coordinates and volumes.
    
    cellBaseGrid=CellCentreGridInformation(baseGrid);
    %cellInd=[cellBaseGrid(:).index];
    for ii=1:numel(cellBaseGrid)
        
        cellBaseGrid(ii).coord=mean(vertcat(cellBaseGrid(ii).vertex(:).coord));
        
    end
    cellCoords=vertcat(cellBaseGrid(:).coord);
    symLines=mean(cellCoords);
    
    symCompCoord=cellCoords;
    symCompCoord(:,dim)=abs(symCompCoord(:,dim)-symLines(dim));
    symCompCoord=[symCompCoord,vertcat(cellBaseGrid(:).volume)];
    compMat1=repmat(reshape(symCompCoord,[1 size(symCompCoord)]),[numel(cellBaseGrid),1]);
    compMat2=repmat(reshape(symCompCoord,[size(symCompCoord,1) 1 size(symCompCoord,2)]),[1 numel(cellBaseGrid)]);
    
    [sub1,sub2]=find(triu(all((abs(compMat1-compMat2)<1e-14),3),1));
    
    
    subFillPos=find(logical([baseGrid.cell(:).isactive]));
    
    subFill1=FindObjNum([],sub1,subFillPos);
    subFill2=FindObjNum([],sub2,subFillPos);
    
    symList=[subFill1,subFill2];
    symList=symList(find(all(symList~=0,2)),:);
    if isempty(symList)
        warning('Symmetry list is empty (because the grid is not symmetric)')
        symList=zeros([2,0]);
    else
        symList=RemoveIdenticalVectors(symList)';
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

function [iterstruct,paroptim]=InitialisePopulation(paroptim,baseGrid)
    
    varExtract={'nDesVar','nPop','startPop','desVarConstr','desVarVal',...
        'optimMethod','desvarconnec','specificFillName','initInterp','numNonFillVar'};
    [nDesVar,nPop,startPop,desVarConstr,desVarVal,optimMethod,desvarconnec,...
        specificFillName,initInterp,numNonFillVar]=ExtractVariables(varExtract,paroptim);
    varExtract={'cellLevels','corneractive','defaultCorner','typDat'};
    [cellLevels,corneractive,defaultCorner,typDat]=ExtractVariables(varExtract,paroptim.parametrisation);
    
    switch startPop
        case 'image'
            fill=Image2Fill(typDat);
            origPop=repmat(fill(:)',[nPop 1]);
        case 'rand'
            origPop=rand([nPop,nDesVar]);
        case 'Rosenrand'
            origPop=rand([nPop,nDesVar])*4-2;
        case 'Rosen'
            origPop=ones([nPop 1])*[-0.527550124697614,0.621201147399720,...
                0.814520320829864, 0.133167225580308,-1.030580011564433,...
                1.199348401514573,-1.143489960105117,0.541368333874143];
        case 'halfuniform'
            origPop=ones([nPop nDesVar])*0.5;
        case 'halfuniformthin'
            origPop=ones([nPop nDesVar])*0.2;
            
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
                origPop=ones([nPop nDesVar])*0.5;
                
                origPop(:,LEind)=origPop(:,LEind)/2;
                origPop(:,TEind)=origPop(:,TEind)/2;
            else
                LEind=1+(cellLevels(1))*[0:(cellLevels(2)-1)];
                TEind=cellLevels(1)+(cellLevels(1))*[0:(cellLevels(2)-1)];
                origPop=ones([nPop nDesVar])*0.5;
                
                origPop(:,LEind)=origPop(:,LEind)/5;
                origPop(:,TEind)=origPop(:,TEind)/5;
                
                %                 origPop(:,LEind+1)=origPop(:,LEind+1)/2;
                %                 origPop(:,TEind-1)=origPop(:,TEind-1)/2;
                %                 origPop(:,LEind)=defaultCorner;
                %                 origPop(:,TEind)=defaultCorner;
            end
            
            
            
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
        case 'NACA0012'
            [rootFill]=NacaOuterLimit0012(baseGrid,paroptim);
            origPop=ones([nPop 1])*rootFill{2};
        case 'NACAmulti'
            
            for ii=1:numel(initInterp)
                [rootFill]=NacaOuterLimit4d(baseGrid,paroptim,initInterp{ii});
                origPop(ii,1:numel(rootFill{2}))=rootFill{2};
            end
            nPop=numel(initInterp);
        case 'AEROmulti'
            for ii=1:numel(initInterp)
                [rootFill]=OuterLimitInverse(baseGrid,paroptim,initInterp{ii});
                origPop(ii,1:numel(rootFill{2}))=rootFill{2};
            end
            nPop=numel(initInterp);
        case 'loadshape'
            specificFillName=MakePathCompliant(specificFillName);
            [rootFill]=MatchVoltoShape(baseGrid,paroptim,specificFillName);
            origPop=ones([nPop 1])*rootFill{2};
        case 'initfreefempop'
            
            [origPop]=InitFreefemPopulation(cellLevels,nPop,nDesVar,desVarConstr,...
                desVarVal);
    end
    
    
    [nonFillPop]=InitialiseNonFillPop(paroptim,size(origPop,1));
    origPop=[origPop,nonFillPop];
    [isGradient]=CheckIfGradient(optimMethod);
    if isGradient
        [origPop,nPop,deltas]=InitialiseGradientBased(origPop(1,:),paroptim);
        paroptim.general.nPop=nPop;
        [iterstruct]=InitialiseIterationStruct(paroptim);
        
        for ii=1:nPop
            iterstruct(1).population(ii).fill=origPop(ii,1:nDesVar);
            iterstruct(1).population(ii).nonfillvar=origPop(ii,nDesVar+1:end);
            iterstruct(1).population(ii).optimdat.var=deltas{ii}(1,:);
            iterstruct(1).population(ii).optimdat.value=deltas{ii}(2,:);
        end
        
    else
        paroptim.general.nPop=nPop;
        [iterstruct]=InitialiseIterationStruct(paroptim);
        [origPop]=OverflowHandling(paroptim,origPop);
        for ii=1:nPop
            iterstruct(1).population(ii).fill=origPop(ii,1:nDesVar);
            iterstruct(1).population(ii).nonfillvar=origPop(ii,nDesVar+1:end);
        end
        if strcmp('NACAmulti',startPop) || strcmp('AEROmulti',startPop)
            for ii=1:nPop
                iterstruct(1).population(ii).additional.Airfoil=initInterp{ii};
            end
        end
    end
    
    iterstruct(1).population=RescaleDesVarNoFill('tovar',paroptim,iterstruct(1).population);
end

function [nonfillPop]=InitialiseNonFillPop(paramoptim,nPop)
    varExtract={'desVarRange','desVarRangeNoFill','nonFillVar','numNonFillVar','startPopNonFill'};
    
    [desVarRange,desVarRangeNoFill,nonFillVar,numNonFillVar,startPopNonFill]=...
        ExtractVariables(varExtract,paramoptim);
    nonfillPop=zeros([nPop,sum(numNonFillVar)]);
    
    for jj=1:numel(nonFillVar)
        nS=1+sum(numNonFillVar(1:(jj-1)));
        nE=sum(numNonFillVar(1:jj));
        
        if isnumeric(startPopNonFill{jj})
            nonfillPop(:,nS:nE)=repmat(startPopNonFill{jj},[nPop,1]);
        else
            switch startPopNonFill{jj}
                case 'rand'
                    nonfillPop(:,nS:nE)=rand([nPop,numNonFillVar(jj)]);
                case 'mid'
                    nonfillPop(:,nS:nE)=ones([nPop,numNonFillVar(jj)])*0.5;
                case 'low'
                    nonfillPop(:,nS:nE)=zeros([nPop,numNonFillVar(jj)]);
                case 'high'
                    nonfillPop(:,nS:nE)=ones([nPop,numNonFillVar(jj)]);
                    
                otherwise
                    
                    error('Unknown starting non fill population')
                    
            end
        end
    end
end

function [origPop,nPop,deltas]=InitialiseGradientBased(rootPop,paroptim)
    
    
    varExtract={'notDesInd','varActive','diffStepSize','desVarRange',...
        'desvarconnec','borderActivation','nDesVar'};
    [notDesInd,varActive,diffStepSize,desVarRange,desvarconnec,borderActivation,nDesVar]...
        =ExtractVariables(varExtract,paroptim);
    
    [inactiveVar]=SelectInactiveVariables(rootPop(1:nDesVar),varActive,desvarconnec,borderActivation);
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
                volLine=currPeak*ones([cellLevels(1)-2+1,1]);
                if numel(volLine)>2
                    volLine([1,end])=0;
                end
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
                deltaLE=5e-4+rand*1e-3;
                volFrac=[deltaLE;volFrac;deltaLE];
                volFrac(volFrac>1.5)=1.5;
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
            stripAct=randperm((nStrips)+1,ceil(nAct/2));
            stripAct=[stripAct,stripAct+1]-1;
            stripAct(stripAct>nStrips)=nStrips;
            stripAct(stripAct<1)=1;
            stripAct=unique(stripAct);
            % Find coefficients for active strips
            nAct=numel(stripAct);
            loStrip=[0,stripAct(1:end-1)];
            hiStrip=[stripAct(2:end),nStrips+1];
            coeffStrip=(hiStrip-loStrip)-1;
            % Generate peaks
            %posPeak=randi(cellLevels(1)-1);
            posPeak = randi(cellLevels(1)-1,length(stripAct),1);
            hPeak=rand([length(stripAct),1]).*coeffStrip';
            ratio=minTargFill/(2*sum(hPeak));
            if ratio>1
                hPeak=ratio*hPeak;
            end
            
            ll=1;
            for jj=stripAct
                
                currPeak=hPeak(ll);
                
                totFrac=currPeak*cellLevels(1)/2;
                volLine=zeros([cellLevels(1)-2+1,1]);
                
                for kk=2:length(volLine)-1
                    if kk<=posPeak(ll)
                        volLine(kk)=currPeak*(kk)/posPeak(ll);
                    else
                        volLine(kk)=currPeak-currPeak*((kk-1)-posPeak(ll))...
                            /(length(volLine)-1-posPeak(ll));
                    end
                end
                volFrac=zeros([cellLevels(1)-2,1]);
                for kk=1:length(volFrac)
                    volFrac(kk)=mean(volLine(kk:kk+1));
                end
                volFrac=[1e-3;volFrac;1e-3];
                
                pop(:,jj)=volFrac;
                ll=ll+1;
            end
        end
        origPop(ii,1:nDesVar)=reshape(pop,[1,nDesVar])*(rand+0.11)/1.1;
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
        case 'testgrad_busemann'
            % Dir_2016-10-14T195656_AreaM2sweep_Nc_1
            % iter 62 pop 5
            rootFill=[0.017758864412533
                0.040091115054392
                0.088013447073493
                0.158576855436767
                0.244766013034044
                0.341655165814136
                0.437886860400727
                0.526237687126564
                0.596263255527965
                0.647600512408900
                0.675462927779081
                0.677062059567565
                0.655705870095454
                0.621064073965241
                0.577486441976678
                0.518179829527757
                0.459524434442214
                0.432528992901254
                0.414613928163020
                0.326783439010478
                0.017758864412461
                0.040091115054298
                0.088013447073474
                0.158576855436773
                0.244766013033964
                0.341655165814072
                0.437886860400853
                0.526237687126764
                0.596263255528053
                0.647600512408853
                0.675462927778959
                0.677062059567526
                0.655705870095593
                0.621064073965397
                0.577486441976701
                0.518179829527685
                0.459524434442118
                0.432528992901215
                0.414613928163050
                0.326783439010473]';
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
    varExtract={'nPop','nDesVar','objectiveName','numNonFillVar'};
    [nPop,nDesVar,objectiveName,numNonFillVar]=ExtractVariables(varExtract,paroptim);
    
    [valFill{1:nPop}]=deal(zeros([1,nDesVar]));
    [valFill2{1:nPop}]=deal(zeros([1,sum(numNonFillVar)]));
    switch objectiveName
        case 'CutCellFlow'
            addstruct=struct('iter',[],'res',[],'cl',[],'cm',[],'cd',[],...
                'cx',[],'cy',[],'A',[],'L',[],'t',[],'c',[],'tc',[],'snaxelVolRes',[],'snaxelVelResV',[]);
            
        case 'LengthArea'
            addstruct=struct('A',[],'L',[],'t',[],'c',[],'tc',[],'snaxelVolRes',[],'snaxelVelResV',[]);
            
        case 'InverseDesign'
            addstruct=struct('sum',[],'mean',[],'std',[],'max',[],'min',[],'A',...
                [],'L',[],'t',[],'c',[],'tc',[],'snaxelVolRes',[],'snaxelVelResV',[]);
            
        case 'InverseDesignTopo'
            addstruct=struct('sum',[],'mean',[],'std',[],'max',[],'min',[],'A',...
                [],'L',[],'t',[],'c',[],'tc',[],'snaxelVolRes',[],'snaxelVelResV',[]);
            
        case 'InverseDesignBulk'
            addstruct=struct('sum',[],'mean',[],'std',[],'max',[],'min',[],'A',...
                [],'L',[],'t',[],'c',[],'tc',[],'snaxelVolRes',[],'snaxelVelResV',[]...
                ,'Airfoil','');
        case 'ASOFlow'
            addstruct=struct('CL',[],'CD',[],'CM',[],'CD0',[],'directResidual',[],...
                'directNIter',[],'adjointResidual',[],'adjointNIter'...
                ,[],'majorIt',[],'minorIt',[],'A',[],'L',[],'t',[],'c',[],...
                'tc',[],'snaxelVolRes',[],'snaxelVelResV',[]);
            
        case 'FreeFemPPTest'
            addstruct=struct('A',...
                [],'L',[],'t',[],'c',[],'tc',[],'snaxelVolRes',[],'snaxelVelResV',[]...
                ,'Airfoil','');
            
        case 'FreeFemPP'
            addstruct=struct('A',...
                [],'L',[],'t',[],'c',[],'tc',[],'snaxelVolRes',[],'snaxelVelResV',[]...
                ,'Airfoil','');
            
            
        otherwise
            addstruct=struct('y',[]);
            
    end
    addstruct.objTime=[];
    
    optimdatstruct=struct('var',[],'value',[],'mutVecPopInd',[],'mutVecFillInd',[]);
    popstruct=struct('fill',valFill,'nonfillvar',valFill2,'location','','objective',[],'constraint'...
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
            disp(['    Time Elapsed:',int2str(floor(tElapsed*24)),datestr(tElapsed,':MM:SS:FFF')]);
            disp(datestr(tStart,0))
            disp('--------------------------------------------------------------------------------------------')
            disp('********************************************************************************************')
            disp('--------------------------------------------------------------------------------------------')
            disp('  ')
        case 1
            disp('  ')
            disp('--------------------------------------------------------------------------------------------')
            disp(['    Time Elapsed:',int2str(floor(tElapsed*24)),datestr(tElapsed,':MM:SS:FFF')]);
            disp(procStart)
            disp('--------------------------------------------------------------------------------------------')
            disp('--------------------------------------------------------------------------------------------')
            disp('  ')
            
        case 2
            
            
            disp(['    Time Elapsed:',int2str(floor(tElapsed*24)),datestr(tElapsed,':MM:SS:FFF')]);
            disp(procStart)
            disp('-----------------------')
            disp('  ')
        case 3
            disp('----------')
            
    end
    
end

%% Objective Function

function [objValue,additional]=EvaluateObjective(objectiveName,paramoptim,...
        member,supportstruct,baseGrid)
    
    procStr=['Calculate Objective - ',objectiveName];
    [tStart]=PrintStart(procStr,2);
    
    objInput=ExtractVariables({'objInput'},paramoptim);
    
    % DVP : replace loop by full support struct add the mesh motion to
    % paramoptim
    loop=supportstruct.loop;
    paramoptim=SetVariables({'parentMesh'},{supportstruct.parentMesh},paramoptim);
    objValue=[];
    
    
    [objValue,additional]=eval([objectiveName,'(paramoptim,member,',objInput,');']);
    
    [tElapsed]=PrintEnd(procStr,2,tStart);
    additional.objTime=datestr(tElapsed,'dd-HH:MM:SS');
end

function [objValue,additional]=LengthArea(paramoptim,member,loop)
    for ii=1:length(loop)
        
        [xMin(ii),xMax(ii),t(ii),L(ii),A(ii)]=...
            ClosedLoopProperties(loop(ii).snaxel.coord(1:end-1,:));
        
    end
    objValue=sum(A)/sum(L);
    
    additional.A=sum(A);
    additional.L=sum(L);
    additional.t=sum(t);
    additional.c=max(xMax)-min(xMin);
    additional.tc=additional.t/additional.c;
end

% Aerodynamics

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

% Inverse design

function [objValue,additional]=InverseDesign(paramoptim,member,loop)
    
    [obj,h]=InverseDesign_Error(paramoptim,loop);
    hgsave(h,[member.location,filesep,'prof.fig']);
    close(h)
    [~,areaAdd]=LengthArea(paramoptim,member,loop);
    objValue=obj.sum;
    
    additional=obj;
    additional.A=areaAdd.A;
    additional.L=areaAdd.L;
    additional.t=areaAdd.t;
    additional.c=areaAdd.c;
    additional.tc=areaAdd.tc;
end

function [objValue,additional]=InverseDesignBulk(paramoptim,member,loop)
    
    paramoptim=SetVariables({'aeroName'},{member.additional.Airfoil},paramoptim);
    
    [obj,h]=InverseDesign_Error(paramoptim,loop);
    hgsave(h,[member.location,filesep,'prof.fig']);
    close(h)
    [~,areaAdd]=LengthArea(paramoptim,member,loop);
    objValue=obj.max;
    
    additional=obj;
    additional.A=areaAdd.A;
    additional.L=areaAdd.L;
    additional.t=areaAdd.t;
    additional.c=areaAdd.c;
    additional.tc=areaAdd.tc;
end

function [objValue,additional]=InverseDesignTopo(paramoptim,member,loop)
    
    [obj,h]=InverseDesign_ErrorTopo(paramoptim,loop);
    hgsave(h,[member.location,filesep,'prof.fig']);
    close(h)
    [~,areaAdd]=LengthArea(paramoptim,member,loop);
    objValue=obj.sum;
    
    additional=obj;
    additional.A=areaAdd.A;
    additional.L=areaAdd.L;
    additional.t=areaAdd.t;
    additional.c=areaAdd.c;
    additional.tc=areaAdd.tc;
end

% Analytical test

function [objValue,additional]=Rosenbrock(paramoptim,member,loop)
    
    [objValue] = RosenbrockFunction(member.fill);
    
    
    additional.y=objValue;
end

% Structures

function [objValue,additional]=FreeFemPP(paramoptim,member,loop)
    
    varExtract={'typeLoop'};
    [typeLoop]=ExtractVariables(varExtract,paramoptim.parametrisation);
    [isInternal]=mod(FindInternalLoop(loop),2);
    coords{numel(loop),11}=[];
    for ii=1:numel(loop)
        coords{ii,1}=loop(ii).(typeLoop)(:,1);
        coords{ii,2}=loop(ii).(typeLoop)(:,2);
        coords{ii,3}=loop(ii).isccw;
        coords{ii,4}=isInternal(ii);
    end
    
    [objValue,outstruct]=CallFreeFemPP(coords,paramoptim);
    
    
    [~,areaAdd]=LengthArea(paramoptim,member,loop);
    
    
    additional=outstruct;
    additional.A=areaAdd.A;
    additional.L=areaAdd.L;
    additional.t=areaAdd.t;
    additional.c=areaAdd.c;
    additional.tc=areaAdd.tc;
end

function [objValue,additional]=FreeFemPPTest(paramoptim,member,loop)
    
    varExtract={'typeLoop'};
    [typeLoop]=ExtractVariables(varExtract,paramoptim.parametrisation);
    [isInternal]=mod(FindInternalLoop(loop),2);
    coords{numel(loop),11}=[];
    for ii=1:numel(loop)
        coords{ii,1}=loop(ii).(typeLoop)(:,1);
        coords{ii,2}=loop(ii).(typeLoop)(:,2);
        coords{ii,3}=loop(ii).isccw;
        coords{ii,4}=isInternal(ii);
    end
    
    objValue=1;
    
    [~,areaAdd]=LengthArea(paramoptim,member,loop);
    
    
    additional.A=areaAdd.A;
    additional.L=areaAdd.L;
    additional.t=areaAdd.t;
    additional.c=areaAdd.c;
    additional.tc=areaAdd.tc;
end

%% Find known Optim

function [paramoptim]=FindKnownOptim(paramoptim,iterstruct,baseGrid,gridrefined,...
        restartsnake,connectstructinfo,outinfo)
    
    varExtract={'objectiveName','knownOptim'};
    [objectiveName,knownOptim]=ExtractVariables(varExtract,paramoptim);
    try
        switch objectiveName
            case 'InverseDesign'
                [paramoptim]=FindKnownOptimInvDesign(paramoptim,baseGrid,gridrefined,...
                    restartsnake,connectstructinfo,outinfo);
            case 'InverseDesignTopo'
                [paramoptim]=FindKnownOptimInvDesign(paramoptim,baseGrid,gridrefined,...
                    restartsnake,connectstructinfo,outinfo);
            case 'CutCellFlow'
                [paramoptim]=FindKnownOptimCutCell(paramoptim,outinfo,iterstruct);
        end
        
    catch ME
        disp(ME.getReport)
        paramoptim=SetVariables({'knownOptim'},{knownOptim},paramoptim);
    end
    
    
end


function [paramoptim]=FindKnownOptimInvDesign(paramoptim,baseGrid,gridrefined,...
        restartsnake,connectstructinfo,outinfo)
    
    varExtract={'aeroClass','aeroName','objectiveName'};
    [aeroClass,aeroName,objectiveName]=ExtractVariables(varExtract,paramoptim);
    
    switch aeroClass
        case 'lib'
            [rootFill]=OuterLimitInverse(baseGrid,paramoptim,aeroName);
        case 'NACA'
            if numel(aeroName)==4
                [rootFill]=NacaOuterLimit4d(baseGrid,paramoptim,aeroName);
            else
                [rootFill]=NacaOuterLimit4d(baseGrid,paramoptim,'0012');
                
            end
    end
    [iterstruct]=InitialiseIterationStruct(paramoptim);
    population=iterstruct(1).population(1);
    population.fill=rootFill{2};
    [~,paramoptim]=ConstraintMethod('init',paramoptim,iterstruct,baseGrid);
    [population]=ConstraintMethod('DesVar',paramoptim,population,baseGrid);
    
    currentMember=population.fill;
    [newGrid,newRefGrid,newrestartsnake]=ReFillGrids(baseGrid,gridrefined,...
        restartsnake,connectstructinfo,currentMember);
    
    paramsnake=SetVariables({'restart','snakData'},{false,'all'},paramoptim.parametrisation);
    paramspline=paramoptim.spline;
    
    % Normal Execution
    [population,supportstruct]=PostExecutionIteration(population,newRefGrid,newrestartsnake,...
        newGrid,connectstructinfo,paramsnake,paramspline,outinfo,0,1,paramoptim);
    [population.objective,population.additional]=...
        EvaluateObjective(objectiveName,paramoptim,population,...
        supportstruct,baseGrid);
    knownOptim=population.objective;
    paramoptim=SetVariables({'knownOptim'},{knownOptim},paramoptim);
end

function [paramoptim]=FindKnownOptimCutCell(paramoptim,outinfo,iterstruct)
    
    varExtract={'direction'};
    [direction]=ExtractVariables(varExtract,paramoptim);
    [dat]=GenerateOptimalSolDat(direction,iterstruct);
    
    [knownOptim]=SupersonicOptimLinRes(paramoptim,outinfo(end).rootDir,...
        dat.xMin,dat.xMax,dat.A,dat.nPoints);
    
    paramoptim=SetVariables({'knownOptim'},{knownOptim},paramoptim);
end

function [dat]=GenerateOptimalSolDat(optimDirection,optimstruct)
    
    [~,posOpt]=eval([optimDirection,'([optimstruct(end).population(:).objective])']);
    optimsolution=optimstruct(end).population(posOpt);
    
    profileDir=optimsolution.location;
    
    c=dir(profileDir);
    isFileName=false;
    ii=0;
    while(~isFileName)
        ii=ii+1;
        isFileName=~isempty(regexp(c(ii).name,'restart', 'once'));
    end
    outload=load([profileDir,filesep,c(ii).name]);
    
    dat.A=optimstruct(end).population(posOpt).additional.A;
    coord=vertcat(outload.loop(:).subdivision);
    dat.xMin=min(coord(:,1));
    dat.xMax=max(coord(:,1));
    dat.t=max(coord(:,2))-min(coord(:,2));
    dat.nPoints=length(coord(:,1));
    
end

%% Refined Optimisation


function [paramoptim,outinfo,iterstruct,unstrGrid,baseGrid,gridrefined,...
        connectstructinfo,unstrRef,restartsnake]=...
        HandleRefinement(paramoptim,iterstruct,outinfo,oldBase,gridrefined,...
        connectstructinfo,refStep,nIter,firstValidIter)
    % This must only recieve a portion of iterstruct
    
    oldGrid.base=oldBase;
    oldGrid.refined=gridrefined;
    oldGrid.connec=connectstructinfo;
    oldGrid.cellrefined=CellCentreGridInformation(gridrefined);
    
    % grid refinement
    [paramoptim,outinfo,iterstruct,unstrGrid,baseGrid,gridrefined,...
        connectstructinfo,unstrRef,restartsnake]=...
        InitialiseRefinement(paramoptim,iterstruct,outinfo,oldGrid,refStep,firstValidIter);
    
    paramoptim=SetVariables({'isRestart'},{true},paramoptim);
    % This line causes problems and exceeding boundmax
    %paramoptim.parametrisation.general.refineSteps=paramoptim.parametrisation.general.refineSteps-1;
    [iterstruct,paramoptim]=GenerateNewPop(paramoptim,iterstruct,nIter,firstValidIter,baseGrid);
    
end

%

function [paramoptim,outinfo,iterstruct,unstrGrid,baseGrid,gridrefined,...
        connectstructinfo,unstrRef,restartsnake]=...
        InitialiseRefinement(paramoptim,iterstruct,outinfoOld,oldGrid,refStep,firstValidIter)
    
    procStr='REFINE OPTIMISATION PROCESS';
    [tStart]=PrintStart(procStr,1);
    supportOptim=paramoptim.optim.supportOptim;
    % Initialise Optimisation
    % Get Parametrisation parameters
    varNames={'optimCase','desvarconnec'};
    [optimCase]=ExtractVariables(varNames(1),paramoptim);
    paramoptim=SetVariables(varNames,{[optimCase,'_',int2str(refStep)],[]},paramoptim);
    varNames={'boundstr','corneractive','defaultCorner'};
    [boundstr,corneractive,defaultCorner]=ExtractVariables(varNames,paramoptim.parametrisation);
    
    
    [outinfo]=OptimisationOutput('init',paramoptim);
    
    diaryFile=[outinfo.rootDir,'\Latest_Diary.log'];
    diaryFile=MakePathCompliant(diaryFile);
    fidDiary=fopen(diaryFile,'w');
    fclose(fidDiary);
    diary off
    diary(diaryFile);
    
    warning(['[~,paramoptim]=ConstraintMethod(''init'',paramoptim,[]);',...
        ' partially supported no refinement in constrained cells']);
    
    % Refine Grid
    varNames={'refineOptim'};
    refineCellLvl=ExtractVariables(varNames,paramoptim);
    %     refparamsnake=SetVariables({'refineGrid','typeRefine'},{refineCellLvl(refStep,:),'automatic'},...
    %         paramoptim.parametrisation);
    % Need to add here the additional refinement options
    refinementStruct.oldgrid=oldGrid;
    refinementStruct.pop=iterstruct(end).population;
    refinementStruct.param=paramoptim;
    oldGrid=SelectRefinementCells(iterstruct(end).population,oldGrid,paramoptim);
    
    [gridrefined,connectstructinfo,oldGrid,refCellLevels]=AnisotropicRefinement...
        (oldGrid.base,oldGrid,paramoptim,iterstruct,refStep);
    refinementStruct.connecdesvar=connectstructinfo;
    %     [~,baseGrid,gridrefined,connectstructinfo,~,~]...
    %         =GridInitAndRefine(refparamsnake,oldGrid.base);
    
    
    defaultFillRefMat=BuildDefaultFillRef(oldGrid,gridrefined,connectstructinfo);
    
    newgrid.base=oldGrid.base;
    newgrid.refined=gridrefined;
    newgrid.connec=connectstructinfo;
    newgrid.cellrefined=CellCentreGridInformation(gridrefined);
    
    [unstrGrid,baseGrid,gridrefined,connectstructinfo,unstrRef,loop]...
        =GridInitAndRefine(paramoptim.parametrisation,gridrefined);
    
    % Update some parameters
    
    [desvarconnec,~,~]=ExtractVolumeCellConnectivity(baseGrid);
    [desvarconnec]=ExtractDesignVariableConnectivity(baseGrid,desvarconnec);
    paramoptim=SetVariables({'desvarconnec'},{desvarconnec},paramoptim);
    [paramoptim]=OptimisationParametersModif(paramoptim,baseGrid);
    
    %iterstruct(1).population=ApplySymmetry(paramoptim,iterstruct(1).population);
    
    varNames={'lineSearch'};
    lineSearch=ExtractVariables(varNames,paramoptim);
    paramoptim=SetVariables(varNames,{~lineSearch},...
        paramoptim);
    
    % Update Fill Information to match snake
    varNames={'refineGrid'};
    refineGrid=ExtractVariables(varNames,paramoptim.parametrisation);
    if numel(refineGrid)==1; refineGrid=ones(1,2)*refineGrid;end
    [gridmatch,~]=GridMatching(oldGrid,newgrid,refineGrid,refCellLevels);
    
    % Fill follows the order of activeCells
    %outinfoOld.rootDir='C:\Users\ap1949\Local Documents\PhD\Development Work\Snakes Volume Parametrisation\results\Optimisation\HPC\Debug_17-02-11\Dir_2017-02-11T150949_volsweeplocal_1_000000e-01_cv_'
    [profloops]=ExtractVolInfo(outinfoOld.rootDir);
    
    actFill=logical([baseGrid.cell(:).isactive]);
    [profloops,transformstruct]=ExtractVolumeFraction(gridmatch,profloops,...
        firstValidIter,defaultFillRefMat,iterstruct,actFill);
    
    
    % Generate new snake restarts.
    [newFrac]=BuildNewRestartFrac(iterstruct,profloops,baseGrid);
    [baseGrid,gridrefined]=ReFracGrids(baseGrid,gridrefined,...
        connectstructinfo,newFrac);
    
    if false % ~corneractive
        [baseGrid,gridrefined]=ReduceCornerFrac(baseGrid,gridrefined,...
            connectstructinfo,defaultCorner);
    end
    
    [gridrefined]=EdgePropertiesReshape(gridrefined);
    [loop]=GenerateSnakStartLoop(gridrefined,boundstr);
    
    %     paramoptim.optim.supportOptim=UpdateStepDir(supportOptim,...
    %         profloops,transformstruct,oldGrid.connec,iterstruct,oldGrid.base,baseGrid);
    [paramoptim]=UpdateSupportOptim(paramoptim,profloops,transformstruct,oldGrid.connec,iterstruct,oldGrid.base,baseGrid);
    iterstruct=RewriteHistory(iterstruct,profloops,baseGrid,firstValidIter,defaultFillRefMat);
    
    [~,~,~,~,restartsnake]=ExecuteSnakes_Optim('snak',gridrefined,loop,...
        baseGrid,connectstructinfo,paramoptim.initparam,...
        paramoptim.spline,outinfo,0,0,0);
    [~,paramoptim]=ConstraintMethod('init',paramoptim,iterstruct,baseGrid); % constraint set up post grid
    [outinfo]=OptimisationOutput('iteration',...
        paramoptim,0,outinfo,iterstruct(1),{});
    save([outinfo.rootDir,filesep,'iteration_0',filesep,'refinementinfo.mat'],'refinementStruct')
    [~]=PrintEnd(procStr,1,tStart);
end

function [paramoptim]=UpdateSupportOptim(paramoptim,profloops,transformstruct,connectstructinfo,...
        iterstruct,oldBaseGrid,newBaseGrid)
    
    
    varNames={'optimMethod','stepAlgo','numNonFillVar'};
    [optimMethod,stepAlgo,numNonFillVar]=ExtractVariables(varNames,paramoptim);
    
    switch optimMethod
        case 'conjgrad'
            supportOptim=paramoptim.optim.supportOptim;
            switch stepAlgo
                
                case 'conjgrad'
                    [newPrevStep]=UpdateStepDir(...
                        supportOptim.curr.prevStep(1:end-(sum(numNonFillVar))),profloops,transformstruct,...
                        connectstructinfo,iterstruct,oldBaseGrid,newBaseGrid);
                    supportOptim.curr.prevStep=[newPrevStep,...
                        supportOptim.curr.prevStep(end-(sum(numNonFillVar)-1):end)];
                case 'BFGS'
                    [newPrevStep]=UpdateStepDir(...
                        supportOptim.curr.prevDir(1:end-(sum(numNonFillVar))),profloops,transformstruct,...
                        connectstructinfo,iterstruct,oldBaseGrid,newBaseGrid);
                    supportOptim.curr.prevDir=[newPrevStep,...
                        supportOptim.curr.prevDir(end-(sum(numNonFillVar)-1):end)];
                    supportOptim.curr.Bk=eye(numel(supportOptim.curr.prevDir));
                    supportOptim.curr.Bkinv=eye(numel(supportOptim.curr.prevDir));
                    supportOptim.curr.iter=1;
                otherwise
                    error('Unspecified behaviour in refinement for a conjgrad step algorithm')
            end
            
            paramoptim.optim.supportOptim=supportOptim;
            
        case 'DE'
            disp('Support Optim not used')
        otherwise
            warning('Unspecified Refinement behaviour for the support structure to the optimisation')
            
    end
    
    
end

function [profloops]=ExtractVolInfo(optimRootDir)
    
    [iterationPath,iterName]=FindDirRegex(optimRootDir,'iteration_[0-9]+',true);
    iterNum=regexp(iterName,'iteration_','split');
    
    profloops=repmat(struct('iter',[],'prof',[],'refinevolfrac',[],'isfilled',false),[1 0]);
    
    for ii=1:length(iterationPath)
        
        [profloops]=[profloops,FindProfileLoops(iterationPath{ii},...
            str2double(iterNum{ii}{2}))];
    end
    
end

function [returnPath,returnName]=FindDir(rootDir,strDir,isTargDir)
    returnPath={};
    returnName={};
    subDir=dir(rootDir);
    subDir(1:2)=[];
    nameCell={subDir(:).name};
    isprofileDirCell=strfind(nameCell,strDir);
    for ii=1:length(subDir)
        subDir(ii).isProfile=(~isempty(isprofileDirCell{ii})) && ...
            ~xor(subDir(ii).isdir,isTargDir);
    end
    
    returnSub=find([subDir(:).isProfile]);
    
    
    if isempty(returnSub)
        fprintf('FindDir Could not find requested item %s in:\n%s \n',strDir,rootDir)
    end
    for ii=1:length(returnSub)
        returnPath{ii}=[rootDir,filesep,subDir(returnSub(ii)).name];
        returnName{ii}=subDir(returnSub(ii)).name;
        
    end
    
end


function [profloops]=FindProfileLoops(rootDir,iterNum)
    
    [profilePath,profileName]=FindDir(rootDir,'profile',true);
    profNum=regexp(profileName,'profile_','split');
    %profNum=profNum(:,2);
    kk=1;
    for ii=1:length(profilePath)
        [restartPath,restartName]=FindDir(profilePath{ii},'restart',false);
        load(restartPath{1},'snakSave')
        
        profloops(kk).iter=iterNum;
        profloops(kk).prof=str2num(profNum{ii}{2});
        
        if isfield(snakSave(end),'volumefraction') %#ok<COLND>
            profloops(kk).refinevolfrac=snakSave(end).volumefraction.refinedInfo; %#ok<COLND>
            profloops(kk).isfilled=true;
        else
            profloops(kk).refinevolfrac=[];
            profloops(kk).isfilled=false;
        end
        kk=kk+1;
    end
    isFilledSub=find(~[profloops(:).isfilled]);
    iterInd=[profloops(:).iter];
    profInd=[profloops(:).prof];
    for ii=isFilledSub
        refProf=find(iterInd==iterInd(ii) & profInd==1);
        profloops(ii).refinevolfrac=profloops(refProf(1)).refinevolfrac;
    end
end

function [profloops,transformstruct]=ExtractVolumeFraction(gridmatch,...
        profloops,firstValidIter,defaultFillRefMat,iterstruct,actFill)
    
    [transformstruct,~]=BuildMatrix(gridmatch);
    
    [profloops]=ConvertProfToFill(profloops,transformstruct,firstValidIter,...
        defaultFillRefMat,iterstruct,actFill);
end

function [transformstruct,coeffMat]=BuildMatrix(gridmatch)
    
    newGridIndsMulti=cell2mat(cellfun(@(new,old)old*ones([1,numel(new)]),...
        {gridmatch.matchstruct(:).oldGridInd},...
        {gridmatch.matchstruct(:).newGridInd},'UniformOutput',false));
    oldGridInds=[gridmatch.matchstruct(:).oldGridInd];
    oldGridCoeff=[gridmatch.matchstruct(:).coeff];
    newGridInds=[gridmatch.matchstruct(:).newGridInd];
    
    oldGridUniq=RemoveIdenticalEntries(sort(oldGridInds));
    oldGridSub=FindObjNum([],oldGridInds,oldGridUniq);
    newGridMultiSub=FindObjNum([],newGridIndsMulti,newGridInds);
    
    coeffMat=zeros([numel(newGridInds),numel(oldGridUniq)]);
    coeffMat(sub2ind(size(coeffMat),newGridMultiSub,oldGridSub))=oldGridCoeff;
    transformstruct.coeff=coeffMat;
    transformstruct.indNew=newGridInds;
    transformstruct.indOld=oldGridUniq;
    transformstruct.volumeNew=[gridmatch.matchstruct(:).newvolume]';
end

function [newFrac]=BuildNewRestartFrac(iterstruct,profloops,baseGrid)
    
    iterProf=[profloops(:).iter];
    profProf=[profloops(:).prof];
    iterProfSubs=FindObjNum([],max(iterProf),iterProf);
    profSub=find(profProf(iterProfSubs)==1);
    newFrac=profloops(iterProfSubs(profSub)).newFracs;
    
    newFrac(newFrac~=0 & logical([baseGrid.cell(:).isactive]'))=0.5;
    
end

function [profloops]=ConvertProfToFill(profloops,transformstruct,...
        firstValidIter,defaultFillRefMat,iterstruct,actFill)
    
    iterProf=[profloops(:).iter];
    
    for ii=find((iterProf>=firstValidIter) | (0==iterProf))
        volSubs=FindObjNum([],transformstruct.indOld,profloops(ii).refinevolfrac.index);
        profloops(ii).newFracs=(transformstruct.coeff*...
            profloops(ii).refinevolfrac.fractionvol(volSubs)')./ transformstruct.volumeNew;
    end
    
    isFilledSub=find(~[profloops(:).isfilled]);
    
    for ii=isFilledSub
        
        profloops(ii).newFracs(actFill)=profloops(ii).newFracs(actFill)+(defaultFillRefMat*...
            (iterstruct(profloops(ii).iter).population(profloops(ii).prof).fill...
            -iterstruct(profloops(ii).iter).population(1).fill)');
    end
    
end

function [stepModif]=UpdateStepDir(stepModif,profloops,transformstruct,connectstructinfo,...
        iterstruct,oldBaseGrid,newBaseGrid)
    
    profProf=[profloops(:).prof];
    iterProf=[profloops(:).iter];
    indFindRoot=find(iterProf==max(iterProf) & profProf==1);
    newInd=[connectstructinfo.cell(:).new];
    matOld2New=zeros(numel(connectstructinfo),numel(newInd));
    
    
    for ii=1:numel(connectstructinfo.cell)
        matOld2New(ii,FindObjNum([],connectstructinfo.cell(ii).new,...
            profloops(indFindRoot).refinevolfrac.index))=1;
    end
    
    profFineOldVol=profloops(indFindRoot).refinevolfrac.fractionvol';
    iterFillOldVol=vertcat(oldBaseGrid.cell(:).fill);
    iterFillOldVol(logical([oldBaseGrid.cell(:).isactive]))=iterstruct(end).population(1).fill;
    
    ratiosFill2Prof=profFineOldVol./(matOld2New'*iterFillOldVol);
    ratiosFill2Prof(isnan(ratiosFill2Prof))=0;
    ratiosFill2Prof(~isfinite(ratiosFill2Prof))=0;
    matOld2New=matOld2New*diag(ratiosFill2Prof);
    
    stepOldVol=zeros([numel(oldBaseGrid.cell),1]);
    stepOldVol(logical([oldBaseGrid.cell(:).isactive]))=stepModif;%supportOptim.curr.prevStep;
    
    volSubs=FindObjNum([],transformstruct.indOld,profloops(indFindRoot).refinevolfrac.index);
    profStepFill=(matOld2New'*stepOldVol);
    newStep=(transformstruct.coeff*profStepFill(volSubs))./ transformstruct.volumeNew;
    
    actFill=logical([newBaseGrid.cell(:).isactive]);
    
    stepModif=newStep(actFill)'; %supportOptim.curr.prevStep=...
end

function [newGrid,newRefGrid]=ReFracGrids(baseGrid,refinedGrid,...
        connectstructinfo,newFracs)
    
    activeCell=true(size(baseGrid.cell));
    activeInd=[baseGrid.cell((activeCell)).index];
    
    connecInd=[connectstructinfo.cell(:).old];
    activConnecSub=FindObjNum([],activeInd,connecInd);
    activeCellSub=find(activeCell);
    refCellInd=[refinedGrid.cell(:).index];
    
    if numel(newFracs)~=numel(activeCellSub)
        error('Fill and Active Set do not match in size')
    end
    if sum(abs(newFracs))==0
        newFracs(round(numel(newFracs)/2))=1e-3;
    end
    
    newGrid=baseGrid;
    newRefGrid=refinedGrid;
    
    for ii=1:length(activeCellSub)
        newGrid.cell(activeCellSub(ii)).fill=newFracs(ii);
        newSub=FindObjNum([],[connectstructinfo.cell(activConnecSub(ii)).new],refCellInd);
        [newRefGrid.cell(newSub).fill]=deal(newFracs(ii));
    end
    
end

function [newGrid,newRefGrid]=ReduceCornerFrac(baseGrid,refinedGrid,...
        connectstructinfo,cornerFill)
    
    activeCell=~logical([baseGrid.cell(:).isactive]);
    activeFill=logical([baseGrid.cell(:).fill]~=0);
    activeCell=activeCell & activeFill;
    activeInd=[baseGrid.cell((activeCell)).index];
    
    connecInd=[connectstructinfo.cell(:).old];
    activConnecSub=FindObjNum([],activeInd,connecInd);
    activeCellSub=find(activeCell);
    refCellInd=[refinedGrid.cell(:).index];
    
    newGrid=baseGrid;
    newRefGrid=refinedGrid;
    
    for ii=1:length(activeCellSub)
        newGrid.cell(activeCellSub(ii)).fill=cornerFill;
        newSub=FindObjNum([],[connectstructinfo.cell(activConnecSub(ii)).new],refCellInd);
        [newRefGrid.cell(newSub).fill]=deal(cornerFill);
    end
    
end

function [loop]=GenerateSnakStartLoop(gridrefined2,boundstr)
    
    isEdge=[gridrefined2.edge(:).(boundstr{1})];
    cond=boundstr{3};
    %[loop]=OrderSurfaceVertexReshape(gridrefined2,isEdge,cond);
    [loop]=GenerateEdgeLoop(gridrefined2,boundstr,1);
end

function [iterstruct]=RewriteHistory(iterstruct,profloops,baseGrid,firstValidIter,...
        defaultFillRefMat)
    
    actFill=logical([baseGrid.cell(:).isactive]);
    iterProf=[profloops(:).iter];
    profProf=[profloops(:).prof];
    for ii=firstValidIter:numel(iterstruct) % enforce
        for jj=1:numel(iterstruct(ii).population)
            iterstruct(ii).population(jj).fill=(defaultFillRefMat*iterstruct(ii).population(jj).fill')';
        end
    end
    for ii=find((iterProf>=firstValidIter))
        
        profSub=find(sort(profProf(FindObjNum([],profloops(ii).iter,iterProf)))==profloops(ii).prof);
        
        iterstruct(profloops(ii).iter).population(profSub).fill=...
            profloops(ii).newFracs(actFill)';
        
    end
    
    % rewrite optimdat
    iterstruct(1).population(1).optimdat.value=0;
    iterstruct(1).population(1).optimdat.var=1;
    for ii=(firstValidIter+2):2:length(iterstruct)
        iterstruct(ii).population(1).optimdat.value=[iterstruct(ii).population(1).fill...
            -iterstruct(ii-1).population(1).fill,iterstruct(ii).population(1).nonfillvar...
            -iterstruct(ii-1).population(1).nonfillvar];
        
        iterstruct(ii).population(1).optimdat.var=1:(numel(iterstruct(ii-1).population(1).fill)+...
            numel(iterstruct(ii-1).population(1).nonfillvar));
    end
    for ii=firstValidIter:length(iterstruct)
        for jj=2:length(iterstruct(ii).population)
            iterstruct(ii).population(jj).optimdat.value=...
                [iterstruct(ii).population(jj).fill-iterstruct(ii).population(1).fill,...
                iterstruct(ii).population(jj).nonfillvar-iterstruct(ii).population(1).nonfillvar];
            iterstruct(ii).population(jj).optimdat.var=find(...
                iterstruct(ii).population(jj).optimdat.value~=0);
            iterstruct(ii).population(jj).optimdat.value=...
                iterstruct(ii).population(jj).optimdat.value(...
                iterstruct(ii).population(jj).optimdat.var);
        end
        iterstruct(ii).population(1).optimdat.value=...
            iterstruct(ii).population(jj).optimdat.value;
        iterstruct(ii).population(1).optimdat.var=...
            iterstruct(ii).population(jj).optimdat.var;
    end
    
    % supportOptim.curr.prevStep
end

function [defaultFillRef]=BuildDefaultFillRef(grid,newGrid,connectToNew)
    
    actSubOld=find([grid.base.cell(:).isactive]);
    actSubNew=logical([newGrid.cell(:).isactive]);
    newFillInd=[connectToNew.cell(actSubOld).new];
    newFillIndOrd=[newGrid.cell(actSubNew).index];
    
    
    defaultFillRef=zeros([numel(newFillInd),numel(actSubOld)]);
    for ii=1:numel(actSubOld)
        nRef=numel(connectToNew.cell(actSubOld(ii)).new);
        newFillSub=FindObjNum([],connectToNew.cell(actSubOld(ii)).new,newFillIndOrd);
        defaultFillRef(newFillSub,ii)=1/nRef;
    end
    
end
%% Refinement Cells Selection

function oldGrid=SelectRefinementCells(lastpopulation,oldGrid,paramoptim)
    % function
    
    varNames={'refineOptimType','refineOptimRatio','refineOptimPopRatio',...
        'symDesVarList','direction','rankType'};
    [refineOptimType,refineOptimRatio,refineOptimPopRatio,symDesVarList,...
        direction,rankType]=ExtractVariables(varNames,paramoptim);
    
    % for all gradient based and non gradient use the entire last
    % population has a guide. ie evaluate the condition for each member.
    
    isRefine=false(size(oldGrid.base.cell));
    [indexPop]=SelectRefinementForPop(lastpopulation,refineOptimPopRatio,direction);
    
    switch refineOptimType
        case 'all'
            % refines all cells
            isRefine=true(size(isRefine));
        case 'notempty'
            activeCell=logical([oldGrid.base(:).isactive]);
            
            activeCellSub=find(activeCell);
            
            for ii=indexPop;
                isRefine(activeCellSub)= isRefine(activeCellSub) | (lastpopulation(ii).fill>0);
            end
        case 'contour'
            % refines all cells with snaxels
            for ii=indexPop;
                isRefine=isRefine | REFINE_contour(lastpopulation(ii));
            end
        case 'desvargrad'
            % refines cells based on the gradient between design variables
            actInd=find([oldGrid.base.cell(:).isactive]);
            cellInd=[oldGrid.base.cell(:).index];
            edgeCellSub=FindObjNum([],[oldGrid.base.edge(:).cellindex],cellInd);
            
            for ii=indexPop;
                cellRank=reshape(REFINE_desvargrad(lastpopulation(ii),...
                    oldGrid.base,actInd,cellInd,edgeCellSub,refineOptimRatio),size(isRefine));
                isRefine=isRefine | RefineSortMethod(cellRank,refineOptimRatio,rankType);
            end
        case 'desvargradadvanced'
            % refines cells based on the gradient between design variables
            actInd=find([oldGrid.base.cell(:).isactive]);
            cellInd=[oldGrid.base.cell(:).index];
            edgeCellSub=FindObjNum([],[oldGrid.base.edge(:).cellindex],cellInd);
            
            oldIndsNewOrd=cell2mat(cellfun(@(new,old)old*ones([1,numel(new)]),...
                {oldGrid.connec.cell(:).new},...
                {oldGrid.connec.cell(:).old},'UniformOutput',false));
            for ii=indexPop;
                cellRank=reshape(REFINE_desvargradadvanced(lastpopulation(ii),...
                    oldGrid.base,actInd,cellInd,edgeCellSub,refineOptimRatio,...
                    oldIndsNewOrd),size(isRefine));
                isRefine=isRefine | RefineSortMethod(cellRank,refineOptimRatio,rankType);
            end
        case 'contlength'
            % refines cells based on the L/sqrt(A)
            
            % refines cells based on the gradient between design variables
            actInd=find([oldGrid.base.cell(:).isactive]);
            cellInd=[oldGrid.base.cell(:).index];
            
            oldIndsNewOrd=cell2mat(cellfun(@(new,old)old*ones([1,numel(new)]),...
                {oldGrid.connec.cell(:).new},...
                {oldGrid.connec.cell(:).old},'UniformOutput',false));
            
            cellCentredCoarse=CellCentreGridInformation(oldGrid.base);
            newIndsCell=[oldGrid.connec.cell(:).new];
            for ii=indexPop;
                cellRank=reshape(REFINE_contlength(lastpopulation(ii),...
                    oldGrid,actInd,cellInd,refineOptimRatio,...
                    oldIndsNewOrd,cellCentredCoarse,newIndsCell),size(isRefine));
                isRefine=isRefine | RefineSortMethod(cellRank,refineOptimRatio,rankType);
            end
            
        case 'contlengthnorm'
            % refines cells based on the Length normalised to the cell
            % size inside the length calculation
            
            % refines cells based on the gradient between design variables
            actInd=find([oldGrid.base.cell(:).isactive]);
            cellInd=[oldGrid.base.cell(:).index];
            
            oldIndsNewOrd=cell2mat(cellfun(@(new,old)old*ones([1,numel(new)]),...
                {oldGrid.connec.cell(:).new},...
                {oldGrid.connec.cell(:).old},'UniformOutput',false));
            
            cellCentredCoarse=CellCentreGridInformation(oldGrid.base);
            newIndsCell=[oldGrid.connec.cell(:).new];
            for ii=indexPop;
                cellRank=reshape(REFINE_contlengthnorm(lastpopulation(ii),...
                    oldGrid,actInd,cellInd,refineOptimRatio,...
                    oldIndsNewOrd,cellCentredCoarse,newIndsCell),size(isRefine));
                isRefine=isRefine | RefineSortMethod(cellRank,refineOptimRatio,rankType);
            end
        case 'contcurvevol'
            % refines cells based on the sqrt(curvature)*A
            % refines cells based on the gradient between design variables
            actInd=find([oldGrid.base.cell(:).isactive]);
            cellInd=[oldGrid.base.cell(:).index];
            
            oldIndsNewOrd=cell2mat(cellfun(@(new,old)old*ones([1,numel(new)]),...
                {oldGrid.connec.cell(:).new},...
                {oldGrid.connec.cell(:).old},'UniformOutput',false));
            
            cellCentredCoarse=CellCentreGridInformation(oldGrid.base);
            newIndsCell=[oldGrid.connec.cell(:).new];
            for ii=indexPop;
                cellRank=reshape(REFINE_contcurvevol(lastpopulation(ii),...
                    oldGrid,actInd,cellInd,refineOptimRatio,...
                    oldIndsNewOrd,cellCentredCoarse,newIndsCell),size(isRefine));
                isRefine=isRefine | RefineSortMethod(cellRank,refineOptimRatio,rankType);
            end
        case 'curvelength'
            % refines cells based on the sqrt(curvature)*A
            % refines cells based on the gradient between design variables
            actInd=find([oldGrid.base.cell(:).isactive]);
            cellInd=[oldGrid.base.cell(:).index];
            
            oldIndsNewOrd=cell2mat(cellfun(@(new,old)old*ones([1,numel(new)]),...
                {oldGrid.connec.cell(:).new},...
                {oldGrid.connec.cell(:).old},'UniformOutput',false));
            
            cellCentredCoarse=CellCentreGridInformation(oldGrid.base);
            newIndsCell=[oldGrid.connec.cell(:).new];
            for ii=indexPop;
                cellRank=reshape(REFINE_curvelength(lastpopulation(ii),...
                    oldGrid,actInd,cellInd,refineOptimRatio,...
                    oldIndsNewOrd,cellCentredCoarse,newIndsCell),size(isRefine));
                isRefine=isRefine | RefineSortMethod(cellRank,refineOptimRatio,rankType);
            end
        case 'contcurvenoedge'
            % refines cells based on the sqrt(curvature)*A
            % refines cells based on the gradient between design variables
            actInd=find([oldGrid.base.cell(:).isactive]);
            cellInd=[oldGrid.base.cell(:).index];
            
            oldIndsNewOrd=cell2mat(cellfun(@(new,old)old*ones([1,numel(new)]),...
                {oldGrid.connec.cell(:).new},...
                {oldGrid.connec.cell(:).old},'UniformOutput',false));
            
            cellCentredCoarse=CellCentreGridInformation(oldGrid.base);
            newIndsCell=[oldGrid.connec.cell(:).new];
            for ii=indexPop;
                cellRank=reshape(REFINE_contcurvevolnoedge(lastpopulation(ii),...
                    oldGrid,actInd,cellInd,refineOptimRatio,...
                    oldIndsNewOrd,cellCentredCoarse,newIndsCell),size(isRefine));
                isRefine=isRefine | RefineSortMethod(cellRank,refineOptimRatio,rankType);
            end
        case 'contcurve'
            % refines cells based on the sqrt(curvature)*A
            % refines cells based on the gradient between design variables
            actInd=find([oldGrid.base.cell(:).isactive]);
            cellInd=[oldGrid.base.cell(:).index];
            
            oldIndsNewOrd=cell2mat(cellfun(@(new,old)old*ones([1,numel(new)]),...
                {oldGrid.connec.cell(:).new},...
                {oldGrid.connec.cell(:).old},'UniformOutput',false));
            
            cellCentredCoarse=CellCentreGridInformation(oldGrid.base);
            newIndsCell=[oldGrid.connec.cell(:).new];
            for ii=indexPop;
                cellRank=reshape(REFINE_contcurve(lastpopulation(ii),...
                    oldGrid,actInd,cellInd,refineOptimRatio,...
                    oldIndsNewOrd,cellCentredCoarse,newIndsCell),size(isRefine));
                isRefine=isRefine | RefineSortMethod(cellRank,refineOptimRatio,rankType);
            end
        case 'contcurvescale'
            % refines cells based on the sqrt(curvature)*A
            % refines cells based on the gradient between design variables
            actInd=find([oldGrid.base.cell(:).isactive]);
            cellInd=[oldGrid.base.cell(:).index];
            
            oldIndsNewOrd=cell2mat(cellfun(@(new,old)old*ones([1,numel(new)]),...
                {oldGrid.connec.cell(:).new},...
                {oldGrid.connec.cell(:).old},'UniformOutput',false));
            
            cellCentredCoarse=CellCentreGridInformation(oldGrid.base);
            newIndsCell=[oldGrid.connec.cell(:).new];
            for ii=indexPop;
                cellRank=reshape(REFINE_contcurvescale(lastpopulation(ii),...
                    oldGrid,actInd,cellInd,refineOptimRatio,...
                    oldIndsNewOrd,cellCentredCoarse,newIndsCell),size(isRefine));
                isRefine=isRefine | RefineSortMethod(cellRank,refineOptimRatio,rankType);
            end
            
        otherwise
            error('Unsupported Cell Selection refinement')
    end
    
    if all(~isRefine)
        error('isRefine not set properly, no refinement has taken place')
    end
    
    % Match Symmetry
    if isempty(symDesVarList)
        symDesVarList=zeros([2,0]);
    end
    isActSub=find([oldGrid.base.cell(:).isactive]);
    isRefine(isActSub(symDesVarList(1,:)))=isRefine(isActSub(symDesVarList(1,:))) ...
        | isRefine(isActSub(symDesVarList(2,:)));
    isRefine(isActSub(symDesVarList(2,:)))=isRefine(isActSub(symDesVarList(1,:))) ...
        | isRefine(isActSub(symDesVarList(2,:)));
    
    for ii=1:numel(oldGrid.base.cell)
        oldGrid.base.cell(ii).isrefine=isRefine(ii);
    end
end

function [indexPop]=SelectRefinementForPop(population,ratio,optimDir)
    
    objPop=[population(:).objective];
    
    switch optimDir
        case 'min'
            [~,indexPop]=sort(objPop,'ascend');
        case 'max'
            
            [~,indexPop]=sort(objPop,'descend');
    end
    indexPop=indexPop(1:max(round(numel(objPop)*ratio),1));
end

function [isRefine]=RefineSortMethod(cellRank,refineOptimRatio,rankType)
    
    switch rankType
        case 'value'
            cellRank=cellRank/max(cellRank);
            isRefine=cellRank>=(1-refineOptimRatio);
        case 'rank'
            numRefine=ceil(sum(cellRank~=0)*refineOptimRatio);
            [~,cellOrd]=sort(cellRank);
            isRefine=false(size(cellRank));
            isRefine((cellOrd((numel(cellRank)+1-(numRefine)):end)))=true;
        case 'rankvalue'
            numRefine=ceil(sum(cellRank~=0)*refineOptimRatio);
            [cellRank,cellOrd]=sort(cellRank);
            isRefine=false(size(cellRank));
            isRefine((cellOrd((numel(cellRank)+1-(numRefine)):end)))=true;
            isRefine=isRefine | (cellRank>(cellRank((numel(cellRank)+1-(numRefine)))*0.95));
        case 'number'
            numRefine=min([numel(cellRank),refineOptimRatio]);
            [~,cellOrd]=sort(cellRank);
            isRefine=false(size(cellRank));
            isRefine((cellOrd((numel(cellRank)+1-(numRefine)):end)))=true;
        otherwise
            error('Unknown ranking type')
    end
    
end

function [isSnax]=REFINE_contour(population)
    % refines all cells with snaxels
    
    [restartPath,~]=FindDir(population.location,'restart',false);
    
    load(restartPath{1},'snakSave')
    snakSave=snakSave;
    isSnax=snakSave(end).volumefraction.isSnax;
    
end

function [cellRank,isRefine]=REFINE_desvargrad(population,gridBase,actInd,cellInd,...
        edgeCellSub,refineOptimRatio)
    % refines cells based on the gradient between design variables
    % this system is edge based It evaluates the gradients between every
    % edges
    % Then rejects edges where any cell is inactive or/and without snaxel
    
    [restartPath,~]=FindDir(population.location,'restart',false);
    
    load(restartPath{1},'snakSave')
    snakSave=snakSave;
    isSnax=snakSave(end).volumefraction.isSnax;
    
    for ii=1:numel(actInd)
        gridBase.cell(actInd(ii)).fill=population.fill(ii);
    end
    
    
    isAct=false([1,numel(gridBase.cell)]);
    isAct(actInd)=true;
    
    cellFillExt=[0,[gridBase.cell(:).fill]];
    isSnaxExt=[false,isSnax];
    isActExt=[false,isAct];
    
    edgeFill=cellFillExt(edgeCellSub+1);
    edgeIsSnax=isSnaxExt(edgeCellSub+1);
    edgeIsAct=isActExt(edgeCellSub+1);
    
    edgeSnax=edgeIsSnax(1:2:end) & edgeIsSnax(2:2:end);
    edgeAct=edgeIsAct(1:2:end) & edgeIsAct(2:2:end);
    edgeGrad=abs(edgeFill(1:2:end)-edgeFill(2:2:end));
    
    cond=~edgeSnax;
    %cond=~edgeSnax | ~edgeAct;
    
    cond=cond | edgeGrad==0;
    edgeSub=1:numel(edgeGrad);
    edgeGrad(cond)=[];
    edgeSub(cond)=[];
    
    
    %     [~,~,edgeRank(edgeSub)]=unique(edgeGrad);
    %     edgeRank=edgeRank/max(edgeRank);
    
    
    cellGrad=zeros([1,numel(gridBase.cell)+1]);
    for jj=1:numel(edgeGrad)
        ii=edgeSub(jj);
        cellGrad(edgeCellSub((ii-1)*2+(1:2))+1)=max([edgeGrad(jj),edgeGrad(jj);
            cellGrad(edgeCellSub((ii-1)*2+(1:2))+1)],[],1);
    end
    %[~,~,cellRank]=unique(cellGrad(2:end));
    cellRank=cellGrad(2:end);
    cellRank(~isAct)=0;
    cellRank=cellRank/max(cellRank);
    isRefine=cellRank>=(1-refineOptimRatio);
    
    %     cellSub=FindObjNum([],[gridBase.edge(isEdgeRefine).cellindex],cellInd);
    %     isRefine=false(size(gridBase.cell));
    %    isRefine(cellSub)=true;
end


function [cellRank,isRefine]=REFINE_desvargradadvanced(population,gridBase,actInd,cellInd,...
        edgeCellSub,refineOptimRatio,oldIndsNewOrd)
    % refines cells based on the gradient between design variables
    % this system is edge based It evaluates the gradients between every
    % edges
    % Then rejects edges where any cell is inactive or/and without snaxel
    
    [restartPath,~]=FindDir(population.location,'restart',false);
    
    load(restartPath{1},'snakSave','restartsnak')
    snakSave=snakSave;
    isSnax=snakSave(end).volumefraction.isSnax;
    
    for ii=1:numel(actInd)
        gridBase.cell(actInd(ii)).fill=population.fill(ii);
    end
    
    isAct=false([1,numel(gridBase.cell)]);
    isAct(actInd)=true;
    
    cellFillExt=[0,[gridBase.cell(:).fill]];
    
    [cellGrad]=IdentifiedCrossedEdge(restartsnak.volfracconnec,restartsnak.snaxel...
        ,cellInd,cellFillExt,oldIndsNewOrd);
    cellRank=cellGrad(2:end);
    cellRank(~isAct)=0;
    cellRank=cellRank/max(cellRank);
    isRefine=cellRank>=(1-refineOptimRatio);
    isRefine(~isAct)=false;
    %     cellSub=FindObjNum([],[gridBase.edge(isEdgeRefine).cellindex],cellInd);
    %     isRefine=false(size(gridBase.cell));
    %    isRefine(cellSub)=true;
end

function [cellGrad]=IdentifiedCrossedEdge(connec,snaxel,cellInd,fillExt,oldIndsNewOrd)
    % THis is going to return all the edges along which the gradient should
    % be calculated
    
    connecEdgeNewInd=[connec.edge(:).newedge];
    connecCellNew=[connec.cell(:).newCellInd];
    snaxEdge=[snaxel(:).edge];
    snaxNewCell=[connec.edge(FindObjNum([],snaxEdge,connecEdgeNewInd)).newcell];
    snaxCell=oldIndsNewOrd(FindObjNum([],snaxNewCell,connecCellNew));
    snaxCellSub=FindObjNum([],snaxCell,cellInd);
    
    snaxFill=fillExt(snaxCellSub+1);
    snaxGrad=abs(snaxFill(1:2:end)-snaxFill(2:2:end));
    
    cellGrad=zeros([1,numel(cellInd)+1]);
    for jj=1:numel(snaxGrad)
        ii=(jj);
        cellGrad(snaxCellSub((ii-1)*2+(1:2))+1)=max([snaxGrad(jj),snaxGrad(jj);
            cellGrad(snaxCellSub((ii-1)*2+(1:2))+1)],[],1);
    end
    
end

function [cellRank,isRefine]=REFINE_contlength(population,grid,actInd,cellInd,...
        refineOptimRatio,oldIndsNewOrd,cellCentredCoarse,newIndsCell)
    % refines cells based on the gradient between design variables
    % this system is edge based It evaluates the gradients between every
    % edges
    % Then rejects edges where any cell is inactive or/and without snaxel
    
    [restartPath,~]=FindDir(population.location,'restart',false);
    
    load(restartPath{1},'snakSave','restartsnak')
    snakSave=snakSave;
    isSnax=snakSave(end).volumefraction.isSnax;
    
    for ii=1:numel(actInd)
        grid.base.cell(actInd(ii)).fill=population.fill(ii);
    end
    
    isAct=false([1,numel(grid.base.cell)]);
    isAct(actInd)=true;
    
    [cellCentredCoarse]=CoarseCellCentred(restartsnak.snaxel,grid.refined,...
        grid.cellrefined,cellCentredCoarse,oldIndsNewOrd,newIndsCell);
    
    cellRank=[cellCentredCoarse(:).lSnax]./sqrt([cellCentredCoarse(:).volume]);
    cellRank(~isfinite(cellRank))=0;
    cellRank(~isAct)=0;
    cellRank=cellRank/max(cellRank);
    isRefine=cellRank>=(1-refineOptimRatio);
    isRefine(~isAct)=false;
    
end


function [cellRank,isRefine]=REFINE_contlengthnorm(population,grid,actInd,cellInd,...
        refineOptimRatio,oldIndsNewOrd,cellCentredCoarse,newIndsCell)
    % refines cells based on the gradient between design variables
    % this system is edge based It evaluates the gradients between every
    % edges
    % Then rejects edges where any cell is inactive or/and without snaxel
    
    [restartPath,~]=FindDir(population.location,'restart',false);
    
    load(restartPath{1},'snakSave','restartsnak')
    snakSave=snakSave;
    isSnax=snakSave(end).volumefraction.isSnax;
    
    for ii=1:numel(actInd)
        grid.base.cell(actInd(ii)).fill=population.fill(ii);
    end
    
    isAct=false([1,numel(grid.base.cell)]);
    isAct(actInd)=true;
    
    [cellCentredCoarse]=CoarseCellCentred(restartsnak.snaxel,grid.refined,...
        grid.cellrefined,cellCentredCoarse,oldIndsNewOrd,newIndsCell);
    
    cellRank=[cellCentredCoarse(:).lSnaxNorm];
    cellRank(~isfinite(cellRank))=0;
    cellRank(~isAct)=0;
    cellRank=cellRank/max(cellRank);
    isRefine=cellRank>=(1-refineOptimRatio);
    isRefine(~isAct)=false;
    
end

function [cellRank,isRefine]=REFINE_contcurve(population,grid,actInd,cellInd,...
        refineOptimRatio,oldIndsNewOrd,cellCentredCoarse,newIndsCell)
    % refines cells based on the gradient between design variables
    % this system is edge based It evaluates the gradients between every
    % edges
    % Then rejects edges where any cell is inactive or/and without snaxel
    
    [restartPath,~]=FindDir(population.location,'restart',false);
    
    load(restartPath{1},'snakSave','restartsnak')
    snakSave=snakSave;
    isSnax=snakSave(end).volumefraction.isSnax;
    
    for ii=1:numel(actInd)
        grid.base.cell(actInd(ii)).fill=population.fill(ii);
    end
    
    isAct=false([1,numel(grid.base.cell)]);
    isAct(actInd)=true;
    
    [cellCentredCoarse]=CoarseCellCentred(restartsnak.snaxel,grid.refined,...
        grid.cellrefined,cellCentredCoarse,oldIndsNewOrd,newIndsCell);
    
    % This was shown to leave no effect of cell size on the refinement
    % criterion "TestCurvChangeArea"
    cellRank=([cellCentredCoarse(:).curvSnax]).*([cellCentredCoarse(:).volume]);
    
    cellRank(~isfinite(cellRank))=0;
    cellRank(~isAct)=0;
    cellRank=cellRank/max(cellRank);
    isRefine=cellRank>=(1-refineOptimRatio);
    isRefine(~isAct)=false;
    
end

function [cellRank,isRefine]=REFINE_contcurvevol(population,grid,actInd,cellInd,...
        refineOptimRatio,oldIndsNewOrd,cellCentredCoarse,newIndsCell)
    % refines cells based on the gradient between design variables
    % this system is edge based It evaluates the gradients between every
    % edges
    % Then rejects edges where any cell is inactive or/and without snaxel
    
    [restartPath,~]=FindDir(population.location,'restart',false);
    
    load(restartPath{1},'snakSave','restartsnak')
    snakSave=snakSave;
    isSnax=snakSave(end).volumefraction.isSnax;
    
    for ii=1:numel(actInd)
        grid.base.cell(actInd(ii)).fill=population.fill(ii);
    end
    
    isAct=false([1,numel(grid.base.cell)]);
    isAct(actInd)=true;
    
    [cellCentredCoarse]=CoarseCellCentred(restartsnak.snaxel,grid.refined,...
        grid.cellrefined,cellCentredCoarse,oldIndsNewOrd,newIndsCell);
    
    % This was shown to leave no effect of cell size on the refinement
    % criterion "TestCurvChangeArea"
    cellRank=sqrt(([cellCentredCoarse(:).curvSnax])).*([cellCentredCoarse(:).volume]);
    
    cellRank(~isfinite(cellRank))=0;
    cellRank(~isAct)=0;
    cellRank=cellRank/max(cellRank);
    isRefine=cellRank>=(1-refineOptimRatio);
    isRefine(~isAct)=false;
    
end


function [cellRank,isRefine]=REFINE_curvelength(population,grid,actInd,cellInd,...
        refineOptimRatio,oldIndsNewOrd,cellCentredCoarse,newIndsCell)
    % refines cells based on the gradient between design variables
    % this system is edge based It evaluates the gradients between every
    % edges
    % Then rejects edges where any cell is inactive or/and without snaxel
    
    [restartPath,~]=FindDir(population.location,'restart',false);
    
    load(restartPath{1},'snakSave','restartsnak')
    snakSave=snakSave;
    isSnax=snakSave(end).volumefraction.isSnax;
    
    for ii=1:numel(actInd)
        grid.base.cell(actInd(ii)).fill=population.fill(ii);
    end
    
    isAct=false([1,numel(grid.base.cell)]);
    isAct(actInd)=true;
    
    [cellCentredCoarse]=CoarseCellCentred(restartsnak.snaxel,grid.refined,...
        grid.cellrefined,cellCentredCoarse,oldIndsNewOrd,newIndsCell);
    
    % This was shown to leave no effect of cell size on the refinement
    % criterion "TestCurvChangeArea"
    cellRank=sqrt(([cellCentredCoarse(:).curvSnax])).*([cellCentredCoarse(:).lSnax]).*([cellCentredCoarse(:).volume]);
    
    cellRank(~isfinite(cellRank))=0;
    cellRank(~isAct)=0;
    cellRank=cellRank/max(cellRank);
    isRefine=cellRank>=(1-refineOptimRatio);
    isRefine(~isAct)=false;
    
end


function [cellRank,isRefine]=REFINE_contcurvevolnoedge(population,grid,actInd,cellInd,...
        refineOptimRatio,oldIndsNewOrd,cellCentredCoarse,newIndsCell)
    % refines cells based on the gradient between design variables
    % this system is edge based It evaluates the gradients between every
    % edges
    % Then rejects edges where any cell is inactive or/and without snaxel
    
    [restartPath,~]=FindDir(population.location,'restart',false);
    
    load(restartPath{1},'snakSave','restartsnak')
    snakSave=snakSave;
    isSnax=snakSave(end).volumefraction.isSnax;
    
    for ii=1:numel(actInd)
        grid.base.cell(actInd(ii)).fill=population.fill(ii);
    end
    
    isAct=false([1,numel(grid.base.cell)]);
    isAct(actInd)=true;
    
    [cellCentredCoarse]=CoarseCellCentred(restartsnak.snaxel,grid.refined,...
        grid.cellrefined,cellCentredCoarse,oldIndsNewOrd,newIndsCell);
    
    % This was shown to leave no effect of cell size on the refinement
    % criterion "TestCurvChangeArea"
    cellRank=sqrt(([cellCentredCoarse(:).curvNoEdge])).*([cellCentredCoarse(:).volume]);
    
    cellRank(~isfinite(cellRank))=0;
    cellRank(~isAct)=0;
    cellRank=cellRank/max(cellRank);
    isRefine=cellRank>=(1-refineOptimRatio);
    isRefine(~isAct)=false;
    
end


function [cellRank,isRefine]=REFINE_contcurvescale(population,grid,actInd,cellInd,...
        refineOptimRatio,oldIndsNewOrd,cellCentredCoarse,newIndsCell)
    % refines cells based on the gradient between design variables
    % this system is edge based It evaluates the gradients between every
    % edges
    % Then rejects edges where any cell is inactive or/and without snaxel
    
    [restartPath,~]=FindDir(population.location,'restart',false);
    
    load(restartPath{1},'snakSave','restartsnak')
    snakSave=snakSave;
    isSnax=snakSave(end).volumefraction.isSnax;
    
    for ii=1:numel(actInd)
        grid.base.cell(actInd(ii)).fill=population.fill(ii);
    end
    
    isAct=false([1,numel(grid.base.cell)]);
    isAct(actInd)=true;
    cellFill=[grid.base.cell(:).fill];
    [cellCentredCoarse]=CoarseCellCentred(restartsnak.snaxel,grid.refined,...
        grid.cellrefined,cellCentredCoarse,oldIndsNewOrd,newIndsCell);
    
    % This was shown to leave no effect of cell size on the refinement
    % criterion "TestCurvChangeArea"
    cellRank=([cellCentredCoarse(:).curvSnax]).*([cellCentredCoarse(:).volume])...
        .*min(cellFill,1-cellFill)*0.5;
    
    cellRank(~isfinite(cellRank))=0;
    cellRank(~isAct)=0;
    cellRank=cellRank/max(cellRank);
    isRefine=cellRank>=(1-refineOptimRatio);
    isRefine(~isAct)=false;
    
end

% Build a cellCentredGrid with Snaxels

function [cellCentredCoarse]=CoarseCellCentred(snaxel,refinedGrid,...
        cellCentredFine,cellCentredCoarse,oldIndsNewOrd,newIndsCell)
    
    
    [cellCentredFine]=IdentifyCellSnaxel(snaxel,refinedGrid,cellCentredFine);
    newEdgeInd=[refinedGrid.edge(:).index];
    % this line matches each Fine cell to its coarse cell in the
    % cellCentredGrid
    newToOldCell=FindObjNum([],oldIndsNewOrd(FindObjNum([],[cellCentredFine(:).index],...
        newIndsCell)),[cellCentredCoarse(:).index]);
    cellCentredCoarse(1).snaxel=struct([]);
    for ii=1:numel(newToOldCell)
        cellCentredCoarse(newToOldCell(ii)).snaxel=[cellCentredCoarse(newToOldCell(ii)).snaxel, ...
            cellCentredFine(ii).snaxel];
    end
    for ii=1:numel(cellCentredCoarse)
        cellCentredCoarse(ii).lSnax=0;
        cellCentredCoarse(ii).curvSnax=0;
        cellCentredCoarse(ii).curvNoEdge=0;
        cellCentredCoarse(ii).lSnaxNorm=0;
        coord=vertcat(cellCentredCoarse(ii).vertex(:).coord);
        
        cellCentredCoarse(ii).cellLength=max(coord)-min(coord);
        if ~isempty(cellCentredCoarse(ii).snaxel)
            [cellCentredCoarse(ii).snaxel]...
                =CompressSnaxelChain(cellCentredCoarse(ii).snaxel);
            % mark snaxel which are not on the border of the cell (ie who's edge has only one cell)
            for jj=1:numel(cellCentredCoarse(ii).snaxel)
                snaxCellNew=refinedGrid.edge((FindObjNum([],cellCentredCoarse(ii).snaxel(jj).edge,newEdgeInd))).cellindex;
                snaxCellOld=unique(oldIndsNewOrd(FindObjNum([],snaxCellNew,newIndsCell)));
                cellCentredCoarse(ii).snaxel(jj).isborder=numel(snaxCellOld)>1;
            end
        end
    end
    for ii=1:numel(cellCentredCoarse)
        if ~isempty(cellCentredCoarse(ii).snaxel)
            [cellCentredCoarse(ii).lSnax,cellCentredCoarse(ii).curvSnax,...
                cellCentredCoarse(ii).lSnaxNorm,cellCentredCoarse(ii).curvNoEdge]...
                =ExploreSnaxelChain(cellCentredCoarse(ii).snaxel,cellCentredCoarse(ii).cellLength);
        end
    end
    
end

function [cellCentredCoarse]=CellCentredSnaxelInfo(snaxel,refinedGrid,...
        cellCentredFine,cellCentredCoarse,connecstruct)
    
    oldIndsNewOrd=cell2mat(cellfun(@(new,old)old*ones([1,numel(new)]),...
        {connecstruct.cell(:).newCellInd},...
        {connecstruct.cell(:).oldCellInd},'UniformOutput',false));
    newIndsCell=[connecstruct.cell(:).newCellInd];
    newEdgeInd=[refinedGrid.edge(:).index];
    [cellCentredFine]=IdentifyCellSnaxel(snaxel,refinedGrid,cellCentredFine);
    
    % this line matches each Fine cell to its coarse cell in the
    % cellCentredGrid
    newToOldCell=FindObjNum([],oldIndsNewOrd(FindObjNum([],[cellCentredFine(:).index],...
        newIndsCell)),[cellCentredCoarse(:).index]);
    cellCentredCoarse(1).snaxel=struct([]);
    for ii=1:numel(newToOldCell)
        cellCentredCoarse(newToOldCell(ii)).snaxel=[cellCentredCoarse(newToOldCell(ii)).snaxel, ...
            cellCentredFine(ii).snaxel];
    end
    for ii=1:numel(cellCentredCoarse)
        cellCentredCoarse(ii).lSnax=0;
        cellCentredCoarse(ii).curvSnax=0;
        cellCentredCoarse(ii).curvNoEdge=0;
        cellCentredCoarse(ii).lSnaxNorm=0;
        coord=vertcat(cellCentredCoarse(ii).vertex(:).coord);
        
        cellCentredCoarse(ii).cellLength=max(coord)-min(coord);
        if ~isempty(cellCentredCoarse(ii).snaxel)
            [cellCentredCoarse(ii).snaxel]...
                =CompressSnaxelChain(cellCentredCoarse(ii).snaxel);
            % mark snaxel which are not on the border of the cell (ie who's edge has only one cell)
            for jj=1:numel(cellCentredCoarse(ii).snaxel)
                snaxCellNew=refinedGrid.edge((FindObjNum([],cellCentredCoarse(ii).snaxel(jj).edge,newEdgeInd))).cellindex;
                snaxCellOld=unique(oldIndsNewOrd(FindObjNum([],snaxCellNew,newIndsCell)));
                cellCentredCoarse(ii).snaxel(jj).isborder=numel(snaxCellOld)>1;
            end
        end
        
        
        
    end
    for ii=1:numel(cellCentredCoarse)
        if ~isempty(cellCentredCoarse(ii).snaxel)
            [cellCentredCoarse(ii).lSnax,cellCentredCoarse(ii).curvSnax,...
                cellCentredCoarse(ii).lSnaxNorm,cellCentredCoarse(ii).curvNoEdge]...
                =ExploreSnaxelChain(cellCentredCoarse(ii).snaxel,cellCentredCoarse(ii).cellLength);
        end
    end
    
end

% These functions are now in include_SnakeParam
%{
function [cellCentredGrid]=IdentifyCellSnaxel(snaxel,refinedGrid,cellCentredGrid)
    % Extracts the snaxel data and matches it to the cells
    
    cellCentredGrid(1).snaxel=struct([]);
    [snakposition]=PositionSnakesStruct(snaxel,refinedGrid);
    [snaxel]=CalculateSnaxelCurvature(snaxel,snakposition);
    %     [snakPosInd]=ReferenceCompArray(snakposition,inf,'inf','index');
    %     [edgeInd]=ReferenceCompArray(refinedGrid.edge,inf,'inf','index');
    %     [cellInd]=ReferenceCompArray(cellCentredGrid,inf,'inf','index');
    %     [vertIndex]=ReferenceCompArray(refinedGrid.vertex,inf,'inf','index');
    %     [vertCoord]=ReferenceCompArrayVertCat(refinedGrid.vertex,inf,'inf','coord');
    
    snakPosInd=[snakposition(:).index];
    edgeInd=[refinedGrid.edge(:).index];
    cellInd=[cellCentredGrid(:).index];
    vertIndex=[refinedGrid.vertex(:).index];
    vertCoord=vertcat(refinedGrid.vertex(:).coord);
    
    
    
    %     snakPosFields='coord,vector,vectorprec,vectornext,normvector,edgelength';
    %     snakPosFields=VerticalStringArray(snakPosFields,',');
    snakPosFields={'coord','vector','vectorprec','vectornext','normvector','edgelength'};
    for ii=1:length(snaxel)
        snaxEdge=snaxel(ii).edge;
        snaxEdgeSub=FindObjNum(refinedGrid.edge,snaxEdge,edgeInd);
        snaxCells=refinedGrid.edge(snaxEdgeSub).cellindex;
        snaxCells=snaxCells(snaxCells~=0);
        snaxCellsSub=FindObjNum(cellCentredGrid,snaxCells,cellInd);
        snaxPosSub=FindObjNum(snakposition,snaxel(ii).index,snakPosInd);
        
        for jj=1:length(snaxCellsSub)
            iToVert=vertCoord(find(vertIndex==snaxel(ii).tovertex),:); %#ok<FNDSB> % extract vertex coordinates
            iFromVert=vertCoord(find(vertIndex==snaxel(ii).fromvertex),:); %#ok<FNDSB>
            
            snaxelCell=snaxel(ii);
            for kk=1:length(snakPosFields(:,1))
                fieldNam=deblank(snakPosFields{kk});
                snaxelCell.(fieldNam)=0;
            end
            
            snaxelCell.coord=snakposition(snaxPosSub).coord;
            snaxelCell.vector=snakposition(snaxPosSub).vector;
%             snaxelCell.vectorprec=snakposition(snaxPosSub).vectorprec;
%             snaxelCell.vectornext=snakposition(snaxPosSub).vectornext;
%             snaxelCell.normvector=snakposition(snaxPosSub).normvector;
            snaxelCell.edgelength=norm(iToVert-iFromVert);
            
            cellCentredGrid(snaxCellsSub(jj)).snaxel=...
                [cellCentredGrid(snaxCellsSub(jj)).snaxel,snaxelCell];
        end
    end
    
end

function [snakposition]=PositionSnakesStruct(snaxel,unstructured)
    % Returns an array with Snaxel coordinates preceded by snaxel indices
    vertIndex=[unstructured.vertex(:).index];
    vertCoord=vertcat(unstructured.vertex(:).coord);
    fromVertex=[snaxel(:).fromvertex];
    toVertex=[snaxel(:).tovertex];
    
    nSnaxel=length(snaxel);
    
    for ii=nSnaxel:-1:1
        iToVert=vertCoord(find(vertIndex==toVertex(ii)),:); %#ok<FNDSB> % extract vertex coordinates
        iFromVert=vertCoord(find(vertIndex==fromVertex(ii)),:); %#ok<FNDSB>
        
        snakposition(ii).index=snaxel(ii).index;
        snakposition(ii).coord=iFromVert+(iToVert-iFromVert)*snaxel(ii).d;
        snakposition(ii).vectornotnorm=(iToVert-iFromVert);
        snakposition(ii).vertInit=iFromVert;
        snakposition(ii).vector=(iToVert-iFromVert)/norm(iToVert-iFromVert);
    end
    
end

function [lSnax,curvSnax,lSnaxNorm]=ExploreSnaxelChain(snaxelpart,cellLength)
    
    [snaxInd]=[snaxelpart(:).index];
    [snaxPrec]=[snaxelpart(:).snaxprec];
    [snaxNext]=[snaxelpart(:).snaxnext];
    snaxCon=[snaxPrec;snaxNext];
    isVisited=false(size(snaxInd));
    lSnax=0;
    curvSnax=sum(sqrt(sum(vertcat(snaxelpart(:).curv).^2,2)));
%     plotPoints= @(points) plot(points([1:end],1),points([1:end],2));
%     plotPoints(vertcat(snaxelpart(:).coord));
    while ~all(isVisited)
        
        startSub=find(~isVisited,1,'first');
        
        for isfwd=0:1
            currSub=startSub;
            flag=true;
            while flag
                isVisited(currSub)=true;
                nextSub=FindObjNum([],snaxCon(isfwd+1,currSub),snaxInd);
                if nextSub==0
                    break
                end
                coord1=snaxelpart(currSub).coord;
                coord2=snaxelpart(nextSub).coord;
                lSnax=lSnax+sqrt(sum((coord2-coord1).^2));
                lSnaxNorm=lSnax+sqrt(sum(((coord2-coord1)./cellLength).^2));
                currSub=nextSub;
                flag=~isVisited(currSub);
            end
        end
        
    end
    
end

function [snaxel]=CompressSnaxelChain(snaxel)
    
    snaxInd=[snaxel(:).index];
    
    [~,uniqSub]=unique(snaxInd);
    
    snaxel=snaxel(uniqSub);
    
end

function [snaxel]=CalculateSnaxelCurvature(snaxel,snakposition)
    
    curvFunc=@(pi,pip1,pim1,s1,s2)(-pi*(s1+s2)+pip1*s2+pim1*s1)/(s1^2*s2+s2^2*s1);
    
    
    snaxInd=[snakposition(:).index];
    snaxPrecSub=FindObjNum([],[snaxel(:).snaxprec],snaxInd);
    snaxNextSub=FindObjNum([],[snaxel(:).snaxnext],snaxInd);
%     figure,hold on
    for ii=1:numel(snaxel)
        pi1=snakposition(ii).coord;
        pip1=snakposition(snaxNextSub(ii)).coord;
        pim1=snakposition(snaxPrecSub(ii)).coord;
        snaxel(ii).curv=curvFunc(pi1,pip1,pim1,norm(pi1-pip1),norm(pi1-pim1));
%         quiver(pi1(1),pi1(2),snaxel(ii).curv(1)/50,snaxel(ii).curv(2)/50)
%         plot(pi1(1),pi1(2),'r*')
    end
    
    
end
%}

%% Cell Refinement Orientation

function [gridrefined,fullconnect,oldGrid,refCellLevels]=AnisotropicRefinement(baseGrid,oldGrid,...
        paramoptim,iterstruct,refStep)
    
    varNames={'refinePattern','refineOptim','refineOptimPopRatio','direction','desVarVal','desVarConstr'};
    [refinePattern,refineOptim,refineOptimPopRatio,direction,desVarVal,...
        desVarConstr]=ExtractVariables(varNames,paramoptim);
    
    oldIndsNewOrd=cell2mat(cellfun(@(new,old)old*ones([1,numel(new)]),...
        {oldGrid.connec.cell(:).new},...
        {oldGrid.connec.cell(:).old},'UniformOutput',false));
    
    
    % Do not refine any cell which is subject to a local constraint
    actSub=find(logical([baseGrid.cell(:).isactive]));
    listLocConstr=find(~cellfun(@isempty,regexp(desVarConstr,'LocalVolFrac_')));
    for ii=listLocConstr
        [baseGrid.cell(actSub(desVarVal{ii}{1})).isrefine]=deal(false);
    end
    
    % implement different refinement patterns
    switch refinePattern
        case 'preset'
            refineList=[[baseGrid.cell(:).index];[baseGrid.cell(:).isrefine]];
            refCellLevels=refineOptim(min(refStep,size(refineOptim,1)),1:2);
            refSubSteps=1;
        case 'edgecross'
            [refineList,refCellLevels,refSubSteps]=Refinement_edgecross(baseGrid,...
                oldGrid.refined,iterstruct(end).population,refineOptimPopRatio,...
                direction,oldIndsNewOrd);
        case 'curvature'
            
        otherwise
            error('Refinement Pattern does not exist')
    end
    
    
    % Need to add here the additional refinement options
    
    gridrefined=baseGrid;
    for ii=1:refSubSteps
        refparamsnake=SetVariables({'refineGrid','typeRefine'},...
            {refCellLevels(ii,:),'automatic'},...
            paramoptim.parametrisation);
        cellSub=FindObjNum([],refineList(1,:),[gridrefined.cell(:).index]);
        [gridrefined.cell(:).isrefine]=deal(false);
        for jj=1:numel(cellSub)
            gridrefined.cell(cellSub(jj)).isrefine=refineList(2,jj)==ii;
        end
        [~,newbaseGrid{ii},gridrefinedcell{ii},connectstructinfo{ii},~,~]...
            =GridInitAndRefine(refparamsnake,gridrefined); %#ok<AGROW,ASGLU>
        gridrefined=gridrefinedcell{ii};
    end
    for ii=1:numel(oldGrid.base.cell)
        oldGrid.base.cell(ii).isrefine=refineList(2,ii);
    end
    fullconnect=connectstructinfo{1};
    for ii=2:refSubSteps
        % Build connection information
        cellSub1=FindObjNum([],refineList(1,find(refineList(2,:)==ii)),...
            [fullconnect.cell(:).old]);
        cellSub2=FindObjNum([],refineList(1,find(refineList(2,:)==ii)),...
            [connectstructinfo{ii}.cell(:).old]);
        [fullconnect.cell(cellSub1).new]=...
            deal(connectstructinfo{ii}.cell(cellSub2).new);
    end
    for ii=1:refSubSteps
        % Build the correct refinement vectors
        currActiveSub=find(FindObjNum([],[gridrefinedcell{ii}.cell(:).index],...
            [newbaseGrid{ii}.cell(:).index])==0)';
        currActiveSub=sort([currActiveSub,find(refineList(2,:)==ii)]);
        [gridrefined.cell(currActiveSub).refineVec]=...
            deal(gridrefinedcell{ii}.cell(currActiveSub).refineVec);
    end
end

function [refineList,refCellLevels,refSubSteps]=Refinement_edgecross(baseGrid,...
        refinedGrid,population,refineOptimPopRatio,direction,oldIndsNewOrd)
    
    refineList=[[baseGrid.cell(:).index];[baseGrid.cell(:).isrefine]];
    cellGrad=false([2,size(refineList,2)]);
    [indexPop]=SelectRefinementForPop(population,refineOptimPopRatio,direction);
    
    for ii=indexPop
        [restartPath,~]=FindDir(population(ii).location,'restart',false);
        load(restartPath{1},'snakSave','restartsnak')
        [cellGrad]=cellGrad | CrossedEdgeRefinementOrientation(restartsnak.volfracconnec,restartsnak.snaxel,...
            refinedGrid,refineList(1,:),oldIndsNewOrd);
    end
    
    cellGrad=(cellGrad'*[1;2])';
    
    refineList(2,find(refineList(2,:)))=cellGrad(find(refineList(2,:)));
    refSubSteps=3;
    refCellLevels=[1 2;2 1;2 2];
    actLevels=unique(refineList(2,:));
    actLevels(actLevels==0)=[];
    refineList(2,:)=reshape(FindObjNum([],refineList(2,:),actLevels),...
        [1 size(refineList,2)]);
    refSubSteps=numel(actLevels);
    refCellLevels=refCellLevels(actLevels,:);
    
end

function [cellEdgeLog]=CrossedEdgeRefinementOrientation(connec,snaxel,...
        refineGrid,cellInd,oldIndsNewOrd)
    % Returns the orientation of the active edges of each cell
    
    connecEdgeNewInd=[connec.edge(:).newedge];
    connecCellNew=[connec.cell(:).newCellInd];
    snaxEdge=[snaxel(:).edge];
    snaxNewCell=[connec.edge(FindObjNum([],snaxEdge,connecEdgeNewInd)).newcell];
    snaxCell=oldIndsNewOrd(FindObjNum([],snaxNewCell,connecCellNew));
    snaxCellSub=FindObjNum([],snaxCell,cellInd);
    snaxOrient=[refineGrid.edge(FindObjNum([],snaxEdge,[refineGrid.edge(:).index])).orientation]+1;
    
    cellEdgeOrient=zeros([2,numel(cellInd)+1]);
    for jj=1:numel(snaxOrient)
        ii=(jj);
        curSub=snaxCellSub((ii-1)*2+(1:2))+1;
        if curSub(1)~=curSub(2)
            cellEdgeOrient(snaxOrient(ii),curSub)=cellEdgeOrient(snaxOrient(ii),curSub)+1;
        end
    end
    cellEdgeOrient(:,1)=[];
    cellEdgeLog=cellEdgeOrient>0;
    
    for ii=1:size(cellEdgeOrient,2)
        cellEdgeLog(:,ii)=cellEdgeLog(:,ii) | any(cellEdgeOrient(:,ii)>3);
    end
end
%% Test Function

function []=OptimisationDebug(caseStr,debugArgIn)
    
    
    [paramoptim,outinfo,iterstruct]=OptimisationDebugStart(caseStr);
    if numel(debugArgIn)>5
        paramoptim=debugArgIn{6};
    end
    [grid,loop,restartsnake,snakSave,newFill]=RestartSnakeFill(debugArgIn{1:5});
    popuDebug=iterstruct(1).population(1);
    
    if iscell(newFill)
        popuDebug=repmat(popuDebug,[1,numel(newFill)]);
        for ii=1:numel(newFill)
            popuDebug(ii).fill=newFill{ii};
        end
    elseif size(newFill,1)>1
        popuDebug=repmat(popuDebug,[1,size(newFill,1)]);
        for ii=1:size(newFill,1)
            popuDebug(ii).fill=newFill(ii,:);
            newFill2{ii}=newFill(ii,:);
        end
        newFill=newFill2;
    else
        popuDebug.fill=newFill;
    end
    paramoptim.parametrisation.snakes.step.snakesSteps=100;
    paramoptim.parametrisation.snakes.step.snakData='all';
    paramoptim.parametrisation.snakes.step.snakesConsole=true;
    
    paramsnake=paramoptim.parametrisation;
    paramspline=paramoptim.spline;
    connectstructinfo=grid.connec;
    oldField = 'oldCellInd';
    newField = 'old';
    [connectstructinfo.cell.(newField)] = connectstructinfo.cell.(oldField);...
        connectstructinfo.cell = rmfield(connectstructinfo.cell,oldField);
    oldField = 'newCellInd';
    newField = 'new';
    [connectstructinfo.cell.(newField)] = connectstructinfo.cell.(oldField);...
        connectstructinfo.cell = rmfield(connectstructinfo.cell,oldField);
    
    baseGrid=grid.base;
    gridrefined=grid.refined;
    
    if numel(popuDebug)>1
        supportstruct=cell([1 numel(popuDebug)]);
        parfor ii=1:numel(popuDebug)
            [newGrid,newRefGrid,newrestartsnake]=ReFillGrids(baseGrid,gridrefined,...
                restartsnake,connectstructinfo,newFill{ii});
            [popuDebug(ii),supportstruct{ii}]=NormalExecutionIteration(...
                popuDebug(ii),newRefGrid,newrestartsnake,newGrid,...
                connectstructinfo,paramsnake,paramspline,outinfo,debugArgIn{2}(1),debugArgIn{3}(1)+ii-1,paramoptim);
        end
    else
        
        [newGrid,newRefGrid,newrestartsnake]=ReFillGrids(baseGrid,gridrefined,...
            restartsnake,connectstructinfo,newFill);
        [popuDebug,supportstruct]=NormalExecutionIteration(...
            popuDebug,newRefGrid,newrestartsnake,newGrid,...
            connectstructinfo,paramsnake,paramspline,outinfo,debugArgIn{2},debugArgIn{3},paramoptim);
    end
    
    
    
    %     [paramoptim]=FindKnownOptimInvDesign(paramoptim,baseGrid,gridrefined,...
    %         restartsnake,connectstructinfo,outinfo);
end

function [paramoptim,outinfo,iterstruct]=OptimisationDebugStart(caseStr)
    
    procStr='INITIALISE DEBUG PROCESS';
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
    
    [unstrGrid,baseGrid,gridrefined,connectstructinfo,unstrRef,loop]...
        =GridInitAndRefine(paramoptim.parametrisation);
    [desvarconnec,~,~]=ExtractVolumeCellConnectivity(baseGrid);
    [desvarconnec]=ExtractDesignVariableConnectivity(baseGrid,desvarconnec);
    paramoptim=SetVariables({'desvarconnec'},{desvarconnec},paramoptim);
    
    [paramoptim]=OptimisationParametersModif(paramoptim,baseGrid);
    [iterstruct,paramoptim]=InitialisePopulation(paramoptim,baseGrid);
    
    %iterstruct(1).population=ApplySymmetry(paramoptim,iterstruct(1).population);
    
end


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
function [iterstruct,outinfo]=ExecuteOptimisation(caseStr,restartFromPop,debugArgIn)
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
    
    varExtract={'maxIter','restartSource','refineSteps','refineIter'};
    [maxIter,restartSource,refineSteps,refineIter]=ExtractVariables(varExtract,paramoptim);
    startIter=1;
    firstValidIter=1;
    % Restart
    inNFlag=nargin;
    if inNFlag==2 || ~isempty(restartSource{1})
        if inNFlag==2
            restartSource=restartFromPop;
        end
        [iterstruct,startIter,maxIter,paramoptim,firstValidIter]=RestartOptions(paramoptim,inNFlag,...
            restartSource,maxIter,iterstruct,baseGrid);
    end
    
    % Introduce debug lines of code
    %debugScrip
    
    % Specify starting population
    for refStage=1:refineSteps+1
        
        % Start optimisation Loop
        for nIter=startIter:maxIter
            % Assign design variables to grid
            procStr=['ITERATION ',int2str(nIter)];
            [tStart]=PrintStart(procStr,1);
            % Compute Shape using snakes
            [iterstruct(nIter).population,restartsnake]=PerformIteration(paramoptim,outinfo(refStage),nIter,iterstruct(nIter).population,gridrefined,restartsnake,...
                baseGrid,connectstructinfo);
            OptimisationOutput('optstruct',paramoptim,outinfo(refStage),iterstruct);
            % Evaluate Objective Function
            [iterstruct,paramoptim]=GenerateNewPop(paramoptim,iterstruct,nIter,firstValidIter,baseGrid);
            % create new population
            
            [~]=PrintEnd(procStr,1,tStart);
            % Convergence tests
            if ConvergenceTest_sloperefine(paramoptim,iterstruct,nIter,startIter) ...
                    && (mod(nIter,2)==0) && (refStage~=refineSteps+1)
                fprintf('\n Optimisation Stopped By Slope convergence condition \n');
                break
            end
            if ConvergenceTest_static(paramoptim,iterstruct,nIter,startIter) && (mod(nIter,2)==0)
                fprintf('\n Optimisation Stopped By convergence condition \n');
                break
            end
        end
        
        % Finish Optimisation
        iterstruct(end)=[];
        [~]=PrintEnd(procStr2,0,tStartOpt);
        pause(0.01)
        diary off
        if isempty(nIter)
            nIter=numel(iterstruct);
            startIter=1;
        end
        %try
        %save(['PreOptimOutFinal',int2str(refStage),'.mat'])
         try
            if exist('flagOut','var') && ~flagOut
                disp('Output Skipped')
            else
                OptimisationOutput('final',paramoptim,outinfo,iterstruct(1:nIter));
            end
         catch ME;
             disp(ME.getReport),
         end

        %debugScript2
        
        save('DvpAnisRefineMat.mat')
        if refStage<(refineSteps+1)
            save('DvpAnisRefineMat.mat')
            [paramoptim,outinfo(refStage+1),iterstruct2,~,baseGrid,gridrefined,...
                connectstructinfo,~,restartsnake]=...
                HandleRefinement(paramoptim,iterstruct(1:nIter),outinfo(refStage),baseGrid,gridrefined,...
                connectstructinfo,refStage,nIter,startIter);
            
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
    
end

function [iterstruct,startIter,maxIter,paramoptim,firstValidIter]=RestartOptions(paramoptim,inNFlag,...
        restartSource,maxIter,iterstruct,baseGrid)
    
    
    
    load(restartSource{1})
    if exist('supportOptim','var');paramoptim.optim.supportOptim=supportOptim;end
    startIter=length(optimstruct);
    maxIter=startIter+maxIter;
    iterstruct=[optimstruct,iterstruct];
    
    [iterstruct,paramoptim,firstValidIter]=GenerateRestartPop(paramoptim,...
        iterstruct,startIter,restartSource{2},baseGrid);
    
    paramoptim.general.restartSource=restartSource;
    startIter=startIter+1;
    
    
    
    
end

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
        [iterstruct]=GenerateNewPop(paroptim,iterstruct,startIter,startIter,baseGrid);
        
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
    include_GridCheck
    
    
    
    % Initialise Optimisation
    % Get Parametrisation parameters
    paramoptim=StructOptimParam(caseStr);
    [outinfo]=OptimisationOutput('init',paramoptim);
    diaryFile=[outinfo.rootDir,'\Latest_Diary.log'];
    diaryFile=MakePathCompliant(diaryFile);
    fidDiary=fopen(diaryFile,'w');
    fclose(fidDiary);
    diary(diaryFile);
    
    [~,paramoptim]=ConstraintMethod('init',paramoptim,[]);
    % Initialise Grid
    [unstrGrid,baseGrid,gridrefined,connectstructinfo,unstrRef,loop]...
        =GridInitAndRefine(paramoptim.parametrisation);
    [desvarconnec,~,~]=ExtractVolumeCellConnectivity(baseGrid);
    [paramoptim.general.desvarconnec]=...
        ExtractDesignVariableConnectivity(baseGrid,desvarconnec);
    % Start Parallel Pool
    StartParallelPool(ExtractVariables({'worker'},paramoptim),10);
    
    [paramoptim]=OptimisationParametersModif(paramoptim,baseGrid);
    [iterstruct,paramoptim]=InitialisePopulation(paramoptim,baseGrid);
    
    iterstruct(1).population=ApplySymmetry(paramoptim,iterstruct(1).population);
    
    if ExtractVariables({'useSnake'},paramoptim)
        [~,~,~,~,restartsnake]=ExecuteSnakes_Optim('snak',gridrefined,loop,...
            baseGrid,connectstructinfo,paramoptim.initparam,...
            paramoptim.spline,outinfo,0,0,0);
    else
        restartsnake=struct([]);
    end
    
    [outinfo]=OptimisationOutput('iteration',...
        paramoptim,0,outinfo,iterstruct(1),{});
    
    [~]=PrintEnd(procStr,1,tStart);
end

function []=StartParallelPool(nWorker,nTry)
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
            parpool(poolName)
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
    
end

function [population,restartsnake]=PerformIteration(paramoptim,outinfo,nIter,population,...
        gridrefined,restartsnake,baseGrid,connectstructinfo)
    
    
    varExtract={'nPop','objectiveName','defaultVal','lineSearch','useSnake'};
    [nPop,objectiveName,defaultVal,lineSearch,useSnake]=ExtractVariables(varExtract,paramoptim);
    
    if ~useSnake
        [population,supportstruct,captureErrors]=IterateNoSnake(paramoptim,population,baseGrid);
        
    elseif (~CheckSnakeSensitivityAlgorithm(paramoptim)) || lineSearch
        [population,supportstruct,captureErrors]=IterateNoSensitivity(paramoptim,outinfo,nIter,population,...
            gridrefined,restartsnake,baseGrid,connectstructinfo);
        
    else
        %[population,supportstruct,captureErrors,restartsnake]=IterateSensitivity(paramoptim,outinfo,nIter,population,...
        [population,supportstruct,captureErrors,~]=IterateSensitivity(paramoptim,outinfo,nIter,population,...
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
    %for ii=1:nPop
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
    varExtract={'optimMethod','iterGap','knownOptim','slopeConv'};
    [optimMethod,iterGap,knownOptim,slopeConv]=ExtractVariables(varExtract,paramoptim);
    mmaSpan=3;
    
    isConv=false;
    if CheckIfGradient(optimMethod) && ((nIter-startIter)>=(mmaSpan*iterGap))
        
        objDat=zeros([1,nIter-startIter+1]);
        for ii=startIter:nIter
            objDat(ii-startIter+1)=min([iterstruct(ii).population(:).objective]);
        end
        if knownOptim~=0
            mma=MovingAverage(objDat(1:iterGap:end),mmaSpan);
            mmaSlope=abs(mma(2:end)-mma(1:end-1));
            isConv=(max(mmaSlope(end-(mmaSpan-1):end))/max(mmaSlope))<slopeConv;
        else
            mma=MovingAverage(log10(objDat(1:iterGap:end)),mmaSpan);
            mmaSlope=abs(mma(2:end)-mma(1:end-1));
            isConv=(max(mmaSlope(end-(mmaSpan-1):end))/max(mmaSlope))<slopeConv;
        end
    else
        
        isConv=false; % Need to find a convergence condition for non grad optim
        
    end
    
end
%% Normal Iteration

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
    
end


function [population,supportstruct,captureErrors]=IterateNoSnake(paramoptim,population,baseGrid)
    
    varExtract={'nPop'};
    [nPop]=ExtractVariables(varExtract,paramoptim);
    
    [population]=ConstraintMethod('DesVar',paramoptim,population,baseGrid);
    
    [captureErrors{1:nPop}]=deal('');
    supportstruct=repmat(struct('loop',[]),[1,nPop]);
    
    
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
    supportstructsens=repmat(struct('loop',[],'loopsens',[],'volumefraction',[]),[1,nPop-1]);
    
    if strcmp('analytical',sensCalc)
        [supportstructsens]=ComputeRootSensitivityPopProfiles(paramsnake,paramoptim,...
            population,baseGrid,gridrefined,restartsnake,connectstructinfo,supportstructsens);
        
    end
    
    rootPop=population(1);
    
    parfor ii=1:nPop-1
    %for ii=1:nPop-1
        currentMember=population(ii+1).fill;
        [newGrid,newRefGrid,newrestartsnake]=ReFillGrids(baseGrid,gridrefined,restartsnake,connectstructinfo,currentMember);
        try
            % Normal Execution
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
            % Error Capture
            population(ii+1).constraint=false;
            population(ii+1).exception=['error: ',MEexception.identifier];
            captureErrors{ii+1}=MEexception.getReport;
        end
    end
    
end

function [newpopulation,supportstruct,restartsnake,paramsnake,paramoptim,captureErrors]=...
        ComputeRootSensitivityPopFill(paramsnake,paramspline,paramoptim,...
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
        newfillstruct=SnakesSensitivity('fill',newRefGrid,restartsnake,...
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


function [supportstruct]=...
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
    
    
    
    [~]=PrintEnd(procStr,2,tStart);
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
        ExecuteSnakes_Optim('snak',newRefGrid,newrestartsnake,...
        newGrid,connectstructinfo,paramsnake,paramspline,outinfo,nIter,ii,nPop);
    population.location=outTemp.dirprofile;
    population.additional.snaxelVolRes=snakSave(end).currentConvVolume;
    population.additional.snaxelVelResV=snakSave(end).currentConvVelocity;
    
    supportstruct.loop=loop;
end

function [population,supportstruct,restartsnake]=SensitivityExecutionIteration(population,newRefGrid,...
        supportstructsens,newGrid,connectstructinfo,paramsnake,paramspline...
        ,outinfo,nIter,ii,paramoptim,rootmember)
    
    varExtract={'nPop'};
    [nPop]=ExtractVariables(varExtract,paramoptim);
    
    
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
        =GenerateNewPop(paramoptim,iterstruct,nIter,iterStart,baseGrid)
    procStr=['Generate New Population'];
    [tStart]=PrintStart(procStr,2);
    
    varExtract={'nPop','iterGap','optimMethod'};
    [nPop,iterGap,optimMethod]=ExtractVariables(varExtract,paramoptim);
    
    [isGradient]=CheckIfGradient(optimMethod);
    
    [newPop,iterstruct(nIter).population,paramoptim,deltas]=OptimisationMethod(paramoptim,...
        iterstruct(nIter).population,...
        iterstruct(max([nIter-iterGap,iterStart])).population,baseGrid);
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
        'maxDiffStep','minDiffStep','gradScaleType','nDesVar'};
    [symType,useSnake,startVol,validVol,diffStepSize,maxDiffStep,minDiffStep...
        ,gradScaleType,nDesVar]=ExtractVariables(varExtract,paramoptim);
    varExtract={'cellLevels','corneractive'};
    [cellLevels,corneractive]=ExtractVariables(varExtract,paramoptim.parametrisation);
    if useSnake
        nDesVar=sum([baseGrid.cell(:).isactive]);
        paramoptim.general.nDesVar=nDesVar;
        paramoptim=SetVariables({'gradScale'},{BuildCellRatios(baseGrid,gradScaleType)}...
            ,paramoptim);
    else
        paramoptim=SetVariables({'gradScale'},{ones(1,nDesVar)}...
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
    
    paramoptim.general.symDesVarList...
        =BuildSymmetryLists(symType,cellLevels,corneractive,baseGrid);
    paramoptim.general.notDesInd...
        =BuildExclusionList(paramoptim.general.symDesVarList);
    
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
            paramoptim.general.iterGap=2;
        case 'conjgrad'
            paramoptim.general.iterGap=2;
        otherwise
            
            paramoptim.general.iterGap=1;
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
    
    symList=RemoveIdenticalVectors(symList(find(all(symList~=0,2)),:))';
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
        'optimMethod','desvarconnec','specificFillName','initInterp'};
    [nDesVar,nPop,startPop,desVarConstr,desVarVal,optimMethod,desvarconnec,...
        specificFillName,initInterp]=ExtractVariables(varExtract,paroptim);
    varExtract={'cellLevels','corneractive'};
    [cellLevels,corneractive]=ExtractVariables(varExtract,paroptim.parametrisation);
    
    switch startPop
        case 'rand'
            origPop=rand([nPop,nDesVar]);
        case 'Rosenrand'
            origPop=rand([nPop,nDesVar])*4-2;
        case 'Rosen'
            origPop=ones([nPop 1])*[-0.527550124697614,0.621201147399720,0.814520320829864,...
                0.133167225580308,-1.030580011564433,1.199348401514573,-1.143489960105117,0.541368333874143];
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
                
                origPop(:,LEind)=origPop(:,LEind)/4;
                origPop(:,TEind)=origPop(:,TEind)/4;
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
        if strcmp('NACAmulti',startPop) || strcmp('AEROmulti',startPop)
            for ii=1:nPop
                iterstruct(1).population(ii).additional.Airfoil=initInterp{ii};
            end
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
    varExtract={'nPop','nDesVar','objectiveName'};
    [nPop,nDesVar,objectiveName]=ExtractVariables(varExtract,paroptim);
    
    [valFill{1:nPop}]=deal(zeros([1,nDesVar]));
    switch objectiveName
        case 'CutCellFlow'
            addstruct=struct('iter',[],'res',[],'cl',[],'cm',[],'cd',[],...
                'cx',[],'cy',[],'A',[],'L',[],'t',[],'c',[],'tc',[],'snaxelVolRes',[],'snaxelVelResV',[]);
            
        case 'LengthArea'
            addstruct=struct('A',[],'L',[],'t',[],'c',[],'tc',[],'snaxelVolRes',[],'snaxelVelResV',[]);
            
        case 'InverseDesign'
            addstruct=struct('sum',[],'mean',[],'std',[],'max',[],'min',[],'A',...
                [],'L',[],'t',[],'c',[],'tc',[],'snaxelVolRes',[],'snaxelVelResV',[]);
            
        case 'InverseDesignBulk'
            addstruct=struct('sum',[],'mean',[],'std',[],'max',[],'min',[],'A',...
                [],'L',[],'t',[],'c',[],'tc',[],'snaxelVolRes',[],'snaxelVelResV',[]...
                ,'Airfoil','');
            
            
        otherwise
            addstruct=struct('y',[]);
            
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
    [~,areaAdd]=LengthArea(paramoptim,member,loop);
    objValue=obj.max;
    
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
    
    paramoptim.general.isRestart=true;
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
    diary(diaryFile);
    warning ('[~,paramoptim]=ConstraintMethod(''init'',paramoptim,[]); Not supported');
    
    % Refine Grid
    varNames={'refineOptim'};
    refineCellLvl=ExtractVariables(varNames,paramoptim);
%     refparamsnake=SetVariables({'refineGrid','typeRefine'},{refineCellLvl(refStep,:),'automatic'},...
%         paramoptim.parametrisation);
    % Need to add here the additional refinement options
    refinementStruct.oldgrid=oldGrid;
    refinementStruct.pop=iterstruct(end).population;
    refinementStruct.param=paramoptim; %#ok<STRNU>
    oldGrid=SelectRefinementCells(iterstruct(end).population,oldGrid,paramoptim);
    
    [gridrefined,connectstructinfo,oldGrid,refCellLevels]=AnisotropicRefinement...
        (oldGrid.base,oldGrid,paramoptim,iterstruct,refStep);
    
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
    [paramoptim.general.desvarconnec]=...
        ExtractDesignVariableConnectivity(baseGrid,desvarconnec);
    
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
    
    [profloops,transformstruct]=ExtractVolumeFraction(gridmatch,profloops,firstValidIter);
    
    
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
    
    [outinfo]=OptimisationOutput('iteration',...
        paramoptim,0,outinfo,iterstruct(1),{});
    save([outinfo.rootDir,filesep,'iteration_0',filesep,'refinementinfo.mat'],'refinementStruct')
    [~]=PrintEnd(procStr,1,tStart);
end

function [paramoptim]=UpdateSupportOptim(paramoptim,profloops,transformstruct,connectstructinfo,...
        iterstruct,oldBaseGrid,newBaseGrid)
    
    
    varNames={'optimMethod','stepAlgo'};
    [optimMethod,stepAlgo]=ExtractVariables(varNames,paramoptim);
    
    switch optimMethod
        case 'conjgrad'
            supportOptim=paramoptim.optim.supportOptim;
            switch stepAlgo
                
                case 'conjgrad'
                    [supportOptim.curr.prevStep]=UpdateStepDir(...
                        supportOptim.curr.prevStep,profloops,transformstruct,...
                        connectstructinfo,iterstruct,oldBaseGrid,newBaseGrid);
                case 'BFGS'
                    [supportOptim.curr.prevDir]=UpdateStepDir(...
                        supportOptim.curr.prevDir,profloops,transformstruct,...
                        connectstructinfo,iterstruct,oldBaseGrid,newBaseGrid);
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
    
    [iterationPath,iterName]=FindDir(optimRootDir,'iteration',true);
    iterNum=regexp(iterName,'iteration_','split');
    
    profloops=repmat(struct('iter',[],'prof',[],'refinevolfrac',[]),[1 0]);
    
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
        disp('FindDir Could not find requested item')
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
        profloops(kk).refinevolfrac=snakSave(end).volumefraction.refinedInfo; %#ok<COLND>
        kk=kk+1;
    end
    
    
    
end

function [profloops,transformstruct]=ExtractVolumeFraction(gridmatch,profloops,firstValidIter)
    
    [transformstruct,~]=BuildMatrix(gridmatch);
    
    [profloops]=ConvertProfToFill(profloops,transformstruct,firstValidIter);
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

function [profloops]=ConvertProfToFill(profloops,transformstruct,firstValidIter)
    
    iterProf=[profloops(:).iter];
    
    for ii=find((iterProf>=firstValidIter) | (0==iterProf))
        volSubs=FindObjNum([],transformstruct.indOld,profloops(ii).refinevolfrac.index);
        profloops(ii).newFracs=(transformstruct.coeff*...
            profloops(ii).refinevolfrac.fractionvol(volSubs)')./ transformstruct.volumeNew;
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
    [loop]=OrderSurfaceVertexReshape(gridrefined2,isEdge,cond);
    
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
        iterstruct(ii).population(1).optimdat.value=iterstruct(ii).population(1).fill...
            -iterstruct(ii-1).population(1).fill;
        
        iterstruct(ii).population(1).optimdat.var=1:numel(iterstruct(ii-1).population(1).fill);
    end
    for ii=firstValidIter:length(iterstruct)
        for jj=2:length(iterstruct(ii).population)
            iterstruct(ii).population(jj).optimdat.value=...
                iterstruct(ii).population(jj).fill-iterstruct(ii).population(1).fill;
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
    
    varNames={'refineOptimType','refineOptimRatio','refineOptimPopRatio','symDesVarList','direction'};
    [refineOptimType,refineOptimRatio,refineOptimPopRatio,symDesVarList,...
        direction]=ExtractVariables(varNames,paramoptim);
    
    % for all gradient based and non gradient use the entire last
    % population has a guide. ie evaluate the condition for each member.
    
    isRefine=false(size(oldGrid.base.cell));
    [indexPop]=SelectRefinementForPop(lastpopulation,refineOptimPopRatio,direction);
    
    switch refineOptimType
        case 'all'
            % refines all cells
            isRefine=true(size(isRefine));
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
                isRefine=isRefine | reshape(REFINE_desvargrad(lastpopulation(ii),...
                    oldGrid.base,actInd,cellInd,edgeCellSub,refineOptimRatio),size(isRefine));
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
                isRefine=isRefine | reshape(REFINE_desvargradadvanced(lastpopulation(ii),...
                    oldGrid.base,actInd,cellInd,edgeCellSub,refineOptimRatio,...
                    oldIndsNewOrd),size(isRefine));
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
                isRefine=isRefine | reshape(REFINE_contlength(lastpopulation(ii),...
                    oldGrid,actInd,cellInd,refineOptimRatio,...
                    oldIndsNewOrd,cellCentredCoarse,newIndsCell),size(isRefine));
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
                isRefine=isRefine | reshape(REFINE_contlengthnorm(lastpopulation(ii),...
                    oldGrid,actInd,cellInd,refineOptimRatio,...
                    oldIndsNewOrd,cellCentredCoarse,newIndsCell),size(isRefine));
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
                isRefine=isRefine | reshape(REFINE_contcurve(lastpopulation(ii),...
                    oldGrid,actInd,cellInd,refineOptimRatio,...
                    oldIndsNewOrd,cellCentredCoarse,newIndsCell),size(isRefine));
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
                isRefine=isRefine | reshape(REFINE_contcurvescale(lastpopulation(ii),...
                    oldGrid,actInd,cellInd,refineOptimRatio,...
                    oldIndsNewOrd,cellCentredCoarse,newIndsCell),size(isRefine));
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

function [isSnax]=REFINE_contour(population)
            % refines all cells with snaxels
    
    [restartPath,~]=FindDir(population.location,'restart',false);
    
    load(restartPath{1},'snakSave')
    snakSave=snakSave;
    isSnax=snakSave(end).volumefraction.isSnax;
    
end

function [isRefine]=REFINE_desvargrad(population,gridBase,actInd,cellInd,...
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


function [isRefine]=REFINE_desvargradadvanced(population,gridBase,actInd,cellInd,...
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

function [isRefine]=REFINE_contlength(population,grid,actInd,cellInd,...
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


function [isRefine]=REFINE_contlengthnorm(population,grid,actInd,cellInd,...
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

function [isRefine]=REFINE_contcurve(population,grid,actInd,cellInd,...
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

function [isRefine]=REFINE_contcurvescale(population,grid,actInd,cellInd,...
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
        cellCentredCoarse(ii).lSnaxNorm=0;
        coord=vertcat(cellCentredCoarse(ii).vertex(:).coord);
        
        cellCentredCoarse(ii).cellLength=max(coord)-min(coord);
        if ~isempty(cellCentredCoarse(ii).snaxel)
            [cellCentredCoarse(ii).snaxel]...
                =CompressSnaxelChain(cellCentredCoarse(ii).snaxel);
        end
    end
    for ii=1:numel(cellCentredCoarse)
        if ~isempty(cellCentredCoarse(ii).snaxel)
            [cellCentredCoarse(ii).lSnax,cellCentredCoarse(ii).curvSnax,cellCentredCoarse(ii).lSnaxNorm]...
                =ExploreSnaxelChain(cellCentredCoarse(ii).snaxel,cellCentredCoarse(ii).cellLength);
        end
    end
    
end

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


%% Cell Refinement Orientation

function [gridrefined,fullconnect,oldGrid,refCellLevels]=AnisotropicRefinement(baseGrid,oldGrid,...
        paramoptim,iterstruct,refStep)
    
    varNames={'refinePattern','refineOptim','refineOptimPopRatio','direction'};
    [refinePattern,refineOptim,refineOptimPopRatio,direction]=ExtractVariables(varNames,paramoptim);
    
    oldIndsNewOrd=cell2mat(cellfun(@(new,old)old*ones([1,numel(new)]),...
                {oldGrid.connec.cell(:).new},...
                {oldGrid.connec.cell(:).old},'UniformOutput',false));
            
    switch refinePattern
        case 'preset'
            refineList=[[baseGrid.cell(:).index];[baseGrid.cell(:).isrefine]];
            refCellLevels=refineOptim(refStep,1:2);
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

function [cellEdgeOrient]=CrossedEdgeRefinementOrientation(connec,snaxel,...
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
            cellEdgeOrient(snaxOrient(ii),curSub)=true;
        end
    end
    cellEdgeOrient(:,1)=[];
end
%% Test Function 

function []=OptimisationDebug(caseStr,debugArgIn)
    
    
    [paramoptim,outinfo,iterstruct]=OptimisationDebugStart(caseStr);
    [grid,loop,restartsnake,snakSave,newFill]=RestartSnakeFill(debugArgIn{:});
    popuDebug=iterstruct(1).population(1);
    
    popuDebug.fill=newFill;
    
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
    
    [newGrid,newRefGrid,newrestartsnake]=ReFillGrids(baseGrid,gridrefined,...
        restartsnake,connectstructinfo,newFill);
        
    [popuDebug,supportstruct]=NormalExecutionIteration(...
        popuDebug,newRefGrid,newrestartsnake,newGrid,...
        connectstructinfo,paramsnake,paramspline,outinfo,debugArgIn{2},debugArgIn{3},paramoptim);
    
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
    [paramoptim.general.desvarconnec]=...
        ExtractDesignVariableConnectivity(baseGrid,desvarconnec);
    [paramoptim]=OptimisationParametersModif(paramoptim,baseGrid);
    [iterstruct,paramoptim]=InitialisePopulation(paramoptim,baseGrid);
    
    %iterstruct(1).population=ApplySymmetry(paramoptim,iterstruct(1).population);
    
end


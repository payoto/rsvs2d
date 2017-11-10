% clear
% load('C:\Users\ap1949\Local Documents\PhD\Development Work\Snakes Volume Parametrisation\results\Optimisation\Archive_2017_07\Day_2017-07-24\Dir_2017-07-24T143924_TestProfConstr_1\DvpAnisRefineMat.mat')
if refStage<(refineSteps+1)
    %save('DvpAnisRefineMat.mat')
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
refStage=refStage+1;

save('DebugAtStartIter.mat')

for nIter=startIter:maxIter
    % Assign design variables to grid
    procStr=['ITERATION ',int2str(nIter)];
    [tStart]=PrintStart(procStr,1);
    % Compute Shape using snakes
    [iterstruct(nIter).population,restartsnake]=...
        PerformIteration(paramoptim,outinfo(refStage),nIter,...
        iterstruct(nIter).population,gridrefined,restartsnake,...
        baseGrid,connectstructinfo,iterstruct(1:nIter-1));
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
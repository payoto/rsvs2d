%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%              University of Bristol
%      Department of Aerospace Engineering
%                     2016
%
%            Optimisation Using
%          Parametric Snakes for
%         for Aerodynamic shape
%         parametrisation
%
%             Alexandre Payot
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [newPop,iterCurr,paramoptim,deltas]=OptimisationMethod(paramoptim,varargin)
    % Function distributing the optimisation to various methods
    
    deltas={};
    varExtract={'optimMethod'};
    [optimMethod]=ExtractVariables(varExtract,paramoptim);
    
    switch optimMethod
        case 'DE'
            proj='cut';
            [newPop,iterCurr,paramoptim]=DifferentialEvolution(paramoptim,proj,varargin{1},varargin{2});
        case 'DEtan'
            proj='tan';
            [newPop,iterCurr,paramoptim]=DifferentialEvolution(paramoptim,proj,varargin{1},varargin{2});
            
        case 'conjgrad'
            [newPop,iterCurr,paramoptim,deltas]=ConjugateGradient2(paramoptim,varargin{1},varargin{2});
        case 'conjgradls'
            [newPop,iterCurr,paramoptim,deltas]=ConjugateGradient2(paramoptim,varargin{1},varargin{2});
        case 'DEStrip'
            
        case 'GA'
            
        case 'GSA'
            
        case 'SQP'
            
    end
    
end


%% Optimisation Methods

function [newPop,iterCurr,paramoptim]=DifferentialEvolution(paramoptim,proj,iterCurr,iterm1)
    switch proj
        case 'cut'
            projFunc=@(x) x;%(atan(x)/pi+1/2);
            projInv=@(x) x;%(tan((x-1/2)*pi));
        case 'tan'
            projFunc=@(x) (atan(x)/pi+1/2);
            projInv=@(x) (tan((x-1/2)*pi));
    end
    
    nPop=length(iterCurr);
    DEstruct=paramoptim.optim.DE;
    diffAmp=DEstruct.diffAmplification;
    CR=DEstruct.xOverRatio;
    varExtract={'desVarRange','direction','geneType'};
    [desVarRange,direction,geneType]=ExtractVariables(varExtract,paramoptim);
    
    varExtract={'cellLevels'};
    [cellLevels]=ExtractVariables(varExtract,paramoptim.parametrisation);
    nFill=length(iterCurr(1).fill);
    
    newPop=zeros([nPop,nFill]);
    
    % selection
    switch direction
        case 'min'
            noImprovLog=[iterm1(:).objective]<[iterCurr(:).objective];
        case 'max'
            noImprovLog=[iterm1(:).objective]>[iterCurr(:).objective];
    end
    
    [iterCurr(noImprovLog)]=deal(iterm1(noImprovLog));
    
    switch geneType
        case 'single'
            nDes=nFill;
            for ii=nFill:-1:1
                indexmap{ii}=ii;
            end
            
        case 'horz'
            nDes=cellLevels(2);
            indices=reshape(1:nFill,cellLevels)';
            for ii=nDes:-1:1
                indexmap{ii}=indices(ii,:);
            end
            
        case 'vert'
            nDes=cellLevels(1);
            indices=reshape(1:nFill,cellLevels);
            for ii=nDes:-1:1
                indexmap{ii}=indices(ii,:);
            end
    end
    
    
    for ii=1:nPop
        % Mutation
        rInd=randperm(nPop-1,3);
        rInd(rInd>=ii)=rInd(rInd>=ii)+1;
        mutVec=projFunc(projInv(iterCurr(rInd(1)).fill)+diffAmp*...
            (projInv(iterCurr(rInd(2)).fill)-projInv(iterCurr(rInd(3)).fill)));
        %mutVec=mutVec*(max(desVarRange)-min(desVarRange))+min(desVarRange);
        mutVec(mutVec>max(desVarRange))=max(desVarRange);
        mutVec(mutVec<min(desVarRange))=min(desVarRange);
        % Crossover
        crossVec=-ones([1,nFill]);
        
        [fromMutVecLog]=ExtractDEIndices(nDes,nFill,CR,indexmap);
        
        crossVec(fromMutVecLog)=mutVec(fromMutVecLog);
        crossVec(~fromMutVecLog)=iterCurr(ii).fill(~fromMutVecLog);
        
        newPop(ii,:)=crossVec;
        
    end
    
end

function [fromMutVecLog]=ExtractDEIndices(nDes,nFill,CR,indexmap)
    
    fixInd=randi(nDes,1);
    fromMutVecLogDes=(rand([1,nDes])<=CR);
    fromMutVecLogDes(fixInd)=true;
    fromMutVecLog=false([1,nFill]);
    
    fromMutVecLog([indexmap{fromMutVecLogDes}])=true;
end

% Conjugate Gradient Method

function [newPop,iterCurr,paramoptim]=ConjugateGradient(paramoptim,iterCurr,iterm1)
    % To build more functionality need a universal design variable to fill
    % process construction
    
    varExtract={'diffStepSize','direction','notDesInd','desVarRange',...
        'lineSearch','worker','nPop'};
    [diffStepSize,direction,notDesInd,desVarRange,lineSearch,worker,nPop]...
        =ExtractVariables(varExtract,paramoptim);
    [desVarList]=ExtractActiveVariable(length(iterCurr(1).fill),notDesInd,[]);
    
    
    [gradFt]=ExtractFiniteGradient(iterCurr,desVarList);
    [gradFm1]=ExtractFiniteGradient(iterm1,desVarList);
    Dxtm1=iterCurr(1).fill(desVarList)-iterm1(1).fill(desVarList);
    switch direction
        case 'min'
            DxT=-gradFt+(sum(gradFt.^2)/sum(gradFm1.^2))*Dxtm1;
        case 'max'
            DxT=gradFt+(sum(gradFt.^2)/sum(gradFm1.^2))*Dxtm1;
    end
    
    
    newRootPop=iterCurr(1).fill(desVarList)+DxT;
    newRootFill=iterCurr(1).fill;
    newRootFill(desVarList)=newRootPop;
    [newRootFill,desVarList]=OverflowHandling(paramoptim,desVarList,newRootFill);
    
    newDiffFill=(ones(length(desVarList),1)*newRootFill);
    negMov=find(newRootFill(desVarList)>=max(desVarRange));
    signMat=eye(length(desVarList));
    signMat(:,negMov)=signMat(:,negMov)*-1;
    newDiffFill(:,desVarList)=newDiffFill(:,desVarList)+signMat...
        *diffStepSize;
    overFlowDiff=false(size(newDiffFill));
    overFlowDiff(:,desVarList)=newDiffFill(:,desVarList)>max(desVarRange);
    newDiffFill(overFlowDiff)=max(desVarRange);
    
    newPop=[newRootFill;newDiffFill];
    
    [nPop,nFill]=size(newPop);
    paramoptim.general.nPop=nPop;
    paramoptim.optim.CG.lineSearch=~lineSearch;
end

function [newPop,iterCurr,paramoptim]=ConjugateGradientLS(paramoptim,iterCurr,iterm1)
    % To build more functionality need a universal design variable to fill
    % process construction
    
    varExtract={'diffStepSize','direction','notDesInd','desVarRange',...
        'lineSearch','worker','nPop'};
    [diffStepSize,direction,notDesInd,desVarRange,lineSearch,worker,nPop]...
        =ExtractVariables(varExtract,paramoptim);
    [desVarList]=ExtractActiveVariable(length(iterCurr(1).fill),notDesInd,[]);
    if ~lineSearch
        
        [gradFt]=ExtractFiniteGradient(iterCurr,desVarList);
        [gradFm1]=ExtractFiniteGradient(iterm1,desVarList);
        Dxtm1=iterCurr(1).fill(desVarList)-iterm1(1).fill(desVarList);
        switch direction
            case 'min'
                DxT=-gradFt+(sum(gradFt.^2)/sum(gradFm1.^2))*Dxtm1;
            case 'max'
                DxT=gradFt+(sum(gradFt.^2)/sum(gradFm1.^2))*Dxtm1;
        end
    else
        [DxT]=FindStepLength(iterCurr,direction,nPop);
    end
    
    
    if lineSearch
        
        newRootPop=iterCurr(1).fill(desVarList)+DxT(desVarList);
        newRootFill=iterCurr(1).fill;
        newRootFill(desVarList)=newRootPop;
        [newRootFill,desVarList]=OverflowHandling(paramoptim,desVarList,newRootFill);
        
        newDiffFill=(ones(length(desVarList),1)*newRootFill);
        negMov=find(newRootFill(desVarList)>=max(desVarRange));
        signMat=eye(length(desVarList));
        signMat(:,negMov)=signMat(:,negMov)*-1;
        newDiffFill(:,desVarList)=newDiffFill(:,desVarList)+signMat...
            *diffStepSize;
        overFlowDiff=false(size(newDiffFill));
        overFlowDiff(:,desVarList)=newDiffFill(:,desVarList)>max(desVarRange);
        newDiffFill(overFlowDiff)=max(desVarRange);
        
        newPop=[newRootFill;newDiffFill];
    else
        newRootFill=iterCurr(1).fill(desVarList);
        minChange=min(desVarRange)-newRootFill-0.1;
        maxChange=max(desVarRange)-newRootFill+0.1;
        minChange=minChange./DxT;
        maxChange=maxChange./DxT;
        minMin=min(minChange(minChange>0));
        minMax=min(maxChange(maxChange>0));
        
        changeMax=min([minMin,minMax]);
        
        changeLenght=linspace(0,changeMax,worker);
        
        newPop=ones(worker,1)*iterCurr(1).fill;
        newPop(:,desVarList)=newPop(:,desVarList)+changeLenght'*DxT;
        
        for ii=1:worker
            newPop(ii,:)=OverflowHandling(paramoptim,desVarList,newPop(ii,:));
        end
    end
    [nPop,nFill]=size(newPop);
    paramoptim.general.nPop=nPop;
    paramoptim.optim.CG.lineSearch=~lineSearch;
end


function [gradF]=ExtractFiniteGradient(popstruct,desVarList)
    
    nPop=length(popstruct);
    nFill=length(desVarList);
    
    Dx=vertcat(popstruct(2:end).fill)-(ones((nPop-1),1)*popstruct(1).fill);
    Dx=Dx(:,desVarList);
    isDeriv=Dx~=0;
    flag1=false;
    flag2=false;
    % Checks
    if ~sum(sum(isDeriv,1)~=1)==0
        if ~sum(sum(isDeriv,1)~=1)==nFill
            warning('Conjugate gradient checksums not true')
            flag1=true;
        end
        flag2=true;
    end
    if sum(sum(isDeriv,2)~=1)~=0
        if sum(sum(isDeriv,2)~=2)~=(nPop-1)
            warning('Conjugate gradient checksums not true')
            flag1=true;
        end
        flag2=true;
    end
    if ~flag2 && ~flag1
        [rowPos,colPos]=find(isDeriv);
    else
        rowPos=1:(nFill);
        colPos=1:(nFill);
    end
    %
    
    Df=vertcat(popstruct(2:end).objective)-popstruct(1).objective;
    
    
    gradFMat=Df*ones(1,nFill)./Dx;
    gradF=zeros(1,nFill);
    
    if ~flag1
        for ii=1:length(colPos)
            gradF(colPos(ii))=gradFMat(rowPos(ii),colPos(ii));
        end
    else
        warning('Backup not coded yet')
        for ii=1:length(colPos)
            gradF(colPos(ii))=gradFMat(rowPos(ii),colPos(ii));
        end
    end
    
end

function [DxT]=FindStepLength(popstruct,direction,nPop)
    
    Dx=vertcat(popstruct(2:nPop).fill)-(ones((nPop-1),1)*popstruct(1).fill);
    
    f=vertcat(popstruct(1:end).objective);
    pp=spline(1:nPop,f);
    
    iTest=linspace(1,nPop,1000);
    
    switch direction
        case 'min'
            [targObj,indexLoc]=min(ppval(pp,iTest));
        case 'max'
            [targObj,indexLoc]=max(ppval(pp,iTest));
    end
    
    iNew=iTest(indexLoc);
    prec=floor(iNew);
    next=ceil(iNew);
    vec=popstruct(next).fill-popstruct(prec).fill;
    DxT=vec*(iNew-prec)+popstruct(prec).fill;
    
end

function [newRootFill,desVarList]=OverflowHandling(paramoptim,desVarList,newRootFill)
    % Function which handles steps overflowing the fill constraint
    
    varExtract={'desVarRange','varOverflow'};
    [desVarRange,varOverflow]=ExtractVariables(varExtract,paramoptim);
    
    switch varOverflow
        case 'truncate'
            minD=min(desVarRange);
            maxD=max(desVarRange);
            newRootFill(newRootFill<minD)=minD;
            newRootFill(newRootFill>maxD)=maxD;
            
        otherwise
            error('No valid design variable overflow mechanism')
            
    end
    
    
    
end

%%

function [newPop,iterCurr,paramoptim,deltas]=ConjugateGradient2(paramoptim,iterCurr,iterm1)
    
    varExtract={'diffStepSize','direction','notDesInd','desVarRange',...
        'lineSearch','worker','nPop','validVol','varActive'};
    
    [diffStepSize,direction,notDesInd,desVarRange,lineSearch,worker,...
        nPop,validVol,varActive]=ExtractVariables(varExtract,paramoptim);
    
    
    % Extract previous iteration information
    [gradstruct_curr]=GetIterationInformation(iterCurr);
    [gradstruct_m1]=GetIterationInformation(iterm1);
    
    rootPop=iterCurr(1).fill;
    
    % Case dependant statements
    if lineSearch
        [stepVector]=FindOptimalStepVector(iterCurr,worker,direction);
        [newRoot,deltaRoot]=GenerateNewRootFill(rootPop,stepVector,paramoptim);
        % Need to build function for activation and deactivation of variables
        [inactiveVar]=SelectInactiveVariables(newRoot,varActive);
        [desVarList]=ExtractActiveVariable(length(iterCurr(1).fill),notDesInd,inactiveVar);
        [newGradPop,deltaGrad]=GenerateNewGradientPop(newRoot,desVarRange,diffStepSize,desVarList);
        newPop=[newRoot;newGradPop];
        deltas=[deltaRoot,deltaGrad];
        
        paramoptim.optim.CG.lineSearch=false;
    else % Direction Search
        
        % Get component change
        % The assumption is the the modes have been selected sensibly and
        % will be close to orthogonal
        [modestruct]=ExtractModes(gradstruct_curr,gradstruct_m1);
        [gradF_curr,gradF_m1]...
            =BuildGradientVectors(gradstruct_curr,gradstruct_m1,modestruct);
        
        % Get Corresponding design vector direction
        [stepVector]=NewStepDirection...
            (gradF_curr,gradF_m1,modestruct,iterCurr,direction);
        % Generate Linesearch Distances
        rootPop=iterCurr(1).fill;
        [stepLengths]=StepLengthsForLS(rootPop,stepVector,worker,validVol,desVarRange);
        [newPop,deltas]=...
            GenerateNewLineSearchPop(rootPop,stepVector,stepLengths,paramoptim);
        % Population trimming for invalid values
        
        % Declare linesearch
        paramoptim.optim.CG.lineSearch=true;
    end
    
    % Securing the end
    [nPop,~]=size(newPop);
    paramoptim.general.nPop=nPop;
    
end

function [gradientopt]=GetIterationInformation(population)
    
    nPop=length(population);
    nGrad=nPop-1;
    nFill=length(population(1).fill);
    
    rootPop=population(1).fill;
    rootObj=population(1).objective;
    
    gradientopt=struct('design',zeros(1,nFill),'objective',0,'fill',zeros(1,nFill));
    gradientopt=repmat(gradientopt,[1,nGrad]);
    
    for ii=1:nGrad
        
        gradientopt(ii).design(population(ii+1).optimdat.var)...
            =population(ii+1).optimdat.value;
        gradientopt(ii).objective=population(ii+1).objective-rootObj;
        gradientopt(ii).fill=population(ii+1).fill-rootPop;
        
    end
    
    
    
end

function [modestruct]=ExtractModes(gradstruct_curr,gradstruct_m1)
    
    nCurr=length(gradstruct_curr);
    desModes_curr=vertcat(gradstruct_curr(:).design);
    nM1=length(gradstruct_m1);
    desModes_m1=vertcat(gradstruct_m1(:).design);
    
    allModes=[desModes_curr;desModes_m1];
    
    modeSimilarity=FindIdenticalVectorOrd(allModes);
    
    nModes=length(modeSimilarity);
    
    modestruct=struct('mode',zeros(size(gradstruct_curr(1).design)),...
        'curr',[],'m1',[]);
    
    modestruct=repmat(modestruct,[1,nModes]);
    
    for ii=1:nModes
        modestruct(ii).mode=allModes(modeSimilarity{ii}(1),:);
        for jj=1:length(modeSimilarity{ii})
            if modeSimilarity{ii}(jj)<=nCurr
                modestruct(ii).curr=[modestruct(ii).curr,modeSimilarity{ii}(jj)];
            elseif modeSimilarity{ii}(jj)<=(nCurr+nM1)
                modestruct(ii).m1=[modestruct(ii).m1,modeSimilarity{ii}(jj)-nCurr];
            else
                error('Invalid Index was produced by the mode identification')
            end
        end
    end
    
    
end

function [gradF_curr,gradF_m1]=...
        BuildGradientVectors(gradstruct_curr,gradstruct_m1,modestruct)
    
    
    nModes=length(modestruct);
    
    gradF_curr=zeros([nModes,1]);
    gradF_m1=zeros([nModes,1]);
    
    for ii=1:nModes
        
        currInd=modestruct(ii).curr;
        m1Ind=modestruct(ii).curr;
        if ~isempty(currInd)
            gradF_curr(ii)=gradstruct_curr(currInd(1)).objective;
        end
        if ~isempty(m1Ind)
            gradF_m1(ii)=gradstruct_m1(m1Ind(1)).objective;
        end
    end
    
    
    
end

function [stepVector]=NewStepDirection(gradF_curr,gradF_m1,modestruct,iterCurr,direction)
    
    normVec=@(v) sqrt(sum(v.^2,2));
    
    prevStep=zeros(size(iterCurr(1).fill));
    prevStep(iterCurr(1).optimdat.var)=iterCurr(1).optimdat.value;
    
    modes=vertcat(modestruct(:).mode)';
    
    gradDes_curr=(modes*gradF_curr)';
    gradDes_m1=(modes*gradF_m1)';
    
    switch direction
        case 'min'
            signD=-1;
        case 'max'
            signD=1;
    end
    scale=(normVec(gradDes_curr)/normVec(prevStep))^2;
    %scale=dot(gradDes_curr,gradDes_curr-gradDes_m1)/dot(gradDes_m1,gradDes_m1);
    if ~isfinite(scale)
        scale=1;
    end
    
    stepVector=signD*gradDes_curr...
        +scale*prevStep;
    
    
end

function [stepLengths]=StepLengthsForLS(rootPop,stepVector,worker,validVol,desVarRange)
    
    
    
    % Limit Imposed by design variable range
    minChange=min(desVarRange)-rootPop;
    maxChange=max(desVarRange)-rootPop;
    minMin=FindMinPosMov(minChange,stepVector);
    minMax=FindMinPosMov(maxChange,stepVector);
    % Limit imposed by the maximum step size
    stepToLim=validVol*ones(size(stepVector))./abs(stepVector);
    stepToLim=min(stepToLim(isfinite(stepToLim)));
    
    %
    stepLow=0;
    stepHigh=stepToLim;
    
    % Step Length will be bisecting backtrack
    stepLengths=zeros([1,worker]);
    stepLengths(end)=stepHigh;
    
    for ii=worker-1:-1:2
        stepLengths(ii)=stepLengths(ii+1)/2;
    end
    
end

function [minChange]=FindMinPosMov(change,step)
    
    change=change./step;
    change=change(isfinite(change));
    change=change(change>0);
    minChange=min(change);
end

function [newPop,deltas]=...
        GenerateNewLineSearchPop(rootPop,stepVector,stepLengths,paramoptim)
    
    nNew=length(stepLengths);
    nFill=length(rootPop);
    deltaDes=(stepLengths'*stepVector);
    newPop=(ones([nNew,1])*rootPop)+deltaDes;
    
    for ii=nNew:-1:1
        actVar=find(deltaDes(ii,:));
        deltas{ii}=[actVar;deltaDes(ii,actVar)];
        [newPop(ii,:)]=OverflowHandling2(paramoptim,newPop(ii,:));
    end
    
end

function [newRootFill]=OverflowHandling2(paramoptim,newRootFill)
    % Function which handles steps overflowing the fill constraint
    
    varExtract={'desVarRange','varOverflow'};
    [desVarRange,varOverflow]=ExtractVariables(varExtract,paramoptim);
    
    switch varOverflow
        case 'truncate'
            minD=min(desVarRange);
            maxD=max(desVarRange);
            newRootFill(newRootFill<minD)=minD;
            newRootFill(newRootFill>maxD)=maxD;
            
        otherwise
            error('No valid design variable overflow mechanism')
            
    end
    
    
    
end

function [stepVector]=FindOptimalStepVector(iterstruct,worker,direction)
    
    f=[iterstruct(:).objective];
    g=[iterstruct(:).constraint];
    vec=zeros(size(iterstruct(1).fill));
    vec(iterstruct(end).optimdat.var)=iterstruct(end).optimdat.value;
    
    invalidPoints=g<0.9;
    
    [~,bestPoint]=min(f);
    stepLengths=1./2.^[inf,(worker-2):-1:0];
    
    if bestPoint==1
        bestPoint=2;
    end
    if bestPoint==worker
        bestPoint=worker-1;
    end
    prevPoint=bestPoint-1;
    nextPoint=bestPoint+1;
    
    indSearch=[prevPoint,bestPoint,nextPoint];
    
    
    %pp=spline(stepLengths(indSearch),f(indSearch));
    [pp]=ParabolicFit(stepLengths(indSearch),f(indSearch));
    iTest=linspace(stepLengths(indSearch(1)),stepLengths(indSearch(3)),1000);
    switch direction
        case 'min'
            [targObj,indexLoc]=min(ParabolicVal(pp,iTest));
        case 'max'
            [targObj,indexLoc]=max(ParabolicVal(pp,iTest));
    end
    stepLength=iTest(indexLoc);
    
    if stepLength==0
        warning('Step Length is stagnant this iteration')
    end
    
    stepVector=vec*stepLength;
    
end

function [coeff]=ParabolicFit(xI,yI)
    
    parabola=@(x) [x.^2, x ,ones(size(x))];
    
    xI=reshape(xI,[numel(xI),1]);
    yI=reshape(yI,[numel(yI),1]);
    R=parabola(xI);
    
    if numel(xI)==3
        coeff=R\yI;
    else
        warning('Least Square fit used')
        coeff=(R'*R)\R'*yI;
        
    end
    
end

function [yy]=ParabolicVal(coeff,xx)
    parabola=@(x) [x.^2, x ,ones(size(x))];
    
    xx=reshape(xx,[numel(xx),1]);
    
    yy=parabola(xx)*coeff;
end

function [newRoot,deltas]=GenerateNewRootFill(rootFill,stepVector,paramoptim)
    
    newRoot=rootFill+stepVector;
    [newRoot]=OverflowHandling2(paramoptim,newRoot);
    
    realStep=newRoot-rootFill;
    
    actVar=find(realStep~=0);
    deltas{1}=[actVar;realStep(actVar)];
    
end

function [newGradPop,deltas]=GenerateNewGradientPop(rootFill,desVarRange,stepSize,desVarList)
    
    nActVar=length(desVarList);
    
    newGradPop=(ones(nActVar,1)*rootFill);
    
    negMov=find(rootFill(desVarList)<=(min(desVarRange)+stepSize));
    signMat=-eye(nActVar);
    signMat(:,negMov)=signMat(:,negMov)*-1;
    
    partialSteps=signMat*stepSize;
    newGradPop(:,desVarList)=newGradPop(:,desVarList)+partialSteps;
    
    for ii=nActVar:-1:1
        actVar=find(partialSteps(ii,:)~=0);
        deltas{ii}=[desVarList(actVar);partialSteps(ii,actVar)];
    end
    
end



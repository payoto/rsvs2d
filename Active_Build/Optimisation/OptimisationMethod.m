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
            [newPop,iterCurr,paramoptim,deltas]=ConjugateGradient(paramoptim,varargin{1},varargin{2},varargin{3});
        
        case 'DEStrip'
            
        case 'GA'
            
        case 'GSA'
            
        case 'SQP'
            
        otherwise
            warning('unrecognised optimisation')
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
    warning('Use of cellLevels here, needs to be deprecated')
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
%         mutVec(mutVec>max(desVarRange))=max(desVarRange);
%         mutVec(mutVec<min(desVarRange))=min(desVarRange);
        % Crossover
        crossVec=-ones([1,nFill]);
        
        [fromMutVecLog]=ExtractDEIndices(nDes,nFill,CR,indexmap);
        
        crossVec(fromMutVecLog)=mutVec(fromMutVecLog);
        crossVec(~fromMutVecLog)=iterCurr(ii).fill(~fromMutVecLog);
        
        [crossVec]=OverflowHandling(paramoptim,crossVec);
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

%% Conjugate gradient

function [newPop,iterCurr,paramoptim,deltas]=ConjugateGradient(paramoptim,iterCurr,iterm1,baseGrid)
    
    varExtract={'diffStepSize','direction','notDesInd','desVarRange',...
        'lineSearch','nLineSearch','nPop','validVol','varActive','desvarconnec',...
        'isRestart','borderActivation','lineSearchType','minDiffStep','stepAlgo'};
    
    [diffStepSize,direction,notDesInd,desVarRange,lineSearch,nLineSearch,...
        nPop,validVol,varActive,desvarconnec,isRestart,borderActivation,...
        lineSearchType,minDiffStep,stepAlgo]...
        =ExtractVariables(varExtract,paramoptim);
    supportOptim=paramoptim.optim.supportOptim;
    iterOrig=iterCurr;
    % Extract previous iteration information
    [iterCurr,validCurr]=ExtractValidIter(iterCurr);
    [iterm1,validM1]=ExtractValidIter(iterm1);
    [gradstruct_curr]=GetIterationInformation(iterCurr);
    [gradstruct_m1]=GetIterationInformation(iterm1);
    
    rootPop=iterCurr(1).fill;
    
    % Case dependant statements
    if lineSearch
        if ~isRestart
            [stepVector,validVol,diffStepSize]=FindOptimalStepVector(iterOrig,...
                UnitStepLength(nLineSearch,lineSearchType),direction,validVol,diffStepSize,minDiffStep);
            [newRoot,deltaRoot]=GenerateNewRootFill(rootPop,stepVector,paramoptim,baseGrid);
        else
            [newRoot,deltaRoot]=FindOptimalRestartPop(iterCurr,direction);
            paramoptim.general.isRestart=false;
        end
        
        % Need to build function for activation and deactivation of variables
        [inactiveVar]=SelectInactiveVariables(newRoot,varActive,desvarconnec,borderActivation);
        [desVarList]=ExtractActiveVariable(length(iterCurr(1).fill),notDesInd,inactiveVar);
        [newGradPop,deltaGrad]=GenerateNewGradientPop(newRoot,desVarRange,diffStepSize,desVarList);
        newPop=[newRoot;vertcat(newGradPop{:})];
        deltas=[deltaRoot,deltaGrad{:}];
        paramoptim.optim.CG.validVol=validVol;
        paramoptim.optim.CG.diffStepSize=diffStepSize;
        paramoptim.optim.CG.lineSearch=false;
    else % Direction Search
        
        % Get component change
        % The assumption is the the modes have been selected sensibly and
        % will be close to orthogonal
        [modestruct]=ExtractModes(gradstruct_curr,gradstruct_m1);
        [modestruct]=RemoveFailedModes(modestruct,gradstruct_curr,iterCurr(1).fill...
            ,desVarRange,direction);
        [gradF_curr,gradF_m1]...
            =BuildGradientVectors(gradstruct_curr,gradstruct_m1,modestruct);
        
        % Get Corresponding design vector direction
        [stepVector,supportOptim]=NewStepDirection(stepAlgo,gradF_curr,...
        gradF_m1,modestruct,iterCurr,direction,supportOptim);
        
        % Generate Linesearch Distances
        rootPop=iterCurr(1).fill;
        [unitSteps]=UnitStepLength(nLineSearch,lineSearchType);
        
        [stepLengths]=StepLengthsForLS(rootPop,stepVector,unitSteps,validVol,desVarRange);
            
        [newPop,deltas]=...
            GenerateNewLineSearchPop(rootPop,stepVector,stepLengths,paramoptim);
        % Population trimming for invalid values
        
        % Declare linesearch
        paramoptim.optim.CG.lineSearch=true;
    end
    
    % Securing the end
    [nPop,~]=size(newPop);
    paramoptim.general.nPop=nPop;
    paramoptim.optim.supportOptim=supportOptim;
end

function [population,validIter]=ExtractValidIter(population)
    
    validIter=[1,find([population(2:end).constraint]>0.8)+1];
    population=population(validIter);
    
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
    
    normVec=@(v) sqrt(sum(v.^2,2));
    
    nCurr=length(gradstruct_curr);
    desModes_curr=vertcat(gradstruct_curr(:).design);
    nM1=length(gradstruct_m1);
    desModes_m1=vertcat(gradstruct_m1(:).design);
    
    allModes=[desModes_curr;desModes_m1];
    modeCoeff=zeros([(nCurr+nM1),1]);
    for ii=1:length(allModes(:,1))
        modeCoeff(ii)=allModes(ii,find(allModes(ii,:)~=0,1,'first'));
    end
    modeMultiplier=modeCoeff*ones(size(allModes(1,:)));
    allModes=allModes./modeMultiplier;
    modeCoeffNorm=normVec(allModes);
    modeMultiplier=modeCoeffNorm*ones(size(allModes(1,:)));
    allModes=allModes./modeMultiplier;
    modeCoeff=modeCoeff.*modeCoeffNorm;
    modeSimilarity=FindIdenticalVectorOrd(allModes);
    
    nModes=length(modeSimilarity);
    
    modestruct=struct('mode',zeros(size(gradstruct_curr(1).design)),...
        'curr',struct('ind',[],'coeff',[]),'m1',struct('ind',[],'coeff',[]));
    
    modestruct=repmat(modestruct,[1,nModes]);
    
    for ii=1:nModes
        modestruct(ii).mode=allModes(modeSimilarity{ii}(1),:);
        for jj=1:length(modeSimilarity{ii})
            if modeSimilarity{ii}(jj)<=nCurr
                modestruct(ii).curr.ind=[modestruct(ii).curr.ind,modeSimilarity{ii}(jj)];
                modestruct(ii).curr.coeff(end+1)=modeCoeff(modeSimilarity{ii}(jj));
            elseif modeSimilarity{ii}(jj)<=(nCurr+nM1)
                modestruct(ii).m1.ind=[modestruct(ii).m1.ind,modeSimilarity{ii}(jj)-nCurr];
                modestruct(ii).m1.coeff(end+1)=modeCoeff(modeSimilarity{ii}(jj));
            else
                error('Invalid Index was produced by the mode identification')
            end
        end
    end
    
    
    
    
end

function [modestruct]=RemoveFailedModes(modestruct,gradstruct,rootfill,desVarRange,direction)
    
    [newModeInd]=FindNewModes(modestruct,rootfill,desVarRange);
    
    [indModeFailing]=FindFailingMode(modestruct(newModeInd),gradstruct,direction);
    
    modestruct(newModeInd(indModeFailing))=[];
    
    
end

function [newModeInd]=FindNewModes(modestruct,rootfill,desVarRange)
    % New Modes are caracterised by being 0 or 1 in the root member
    
    activeVar=find((rootfill>min(desVarRange)) & (rootfill<max(desVarRange)));
    isNewMode=false(size(modestruct));
    for ii=1:length(modestruct)
        modeVar=find(modestruct(ii).mode);
        
        isNewMode(ii)=sum((FindObjNum([],modeVar,activeVar)==0))>=ceil(numel(modeVar)/2);
        
    end
    newModeInd=find(isNewMode);
end

function [indModeFailing]=FindFailingMode(modestruct,gradstruct,direction)
    
    switch direction
        case 'min'
            multiplier=1;
        case 'max'
            multiplier=-1;
    end
    isFailing=false(size(modestruct));
    for ii=1:length(modestruct)
        
        obj=[gradstruct(modestruct(ii).curr.ind).objective];
        coeffs=modestruct(ii).curr.coeff;
        if numel(obj)>=2
            if numel(obj)>2
                obj=obj(1:2);
                coeffs=coeffs(1:2);
            end
            slope=multiplier*((obj(1)-obj(2))/(coeffs(1)-coeffs(2)))>=0;
        else
            slope=true;
        end
        val=sum(multiplier*obj>0)>0;
        isFailing(ii)=slope && val;
    end
    indModeFailing=find(isFailing);
end

function [gradF_curr,gradF_m1]=...
        BuildGradientVectors(gradstruct_curr,gradstruct_m1,modestruct)
    
    FDO2= @(fi,fm,fp,dxm,dxp) (fp.*dxm.^2-fm.*dxp.^2+fi.*(dxp.^2-dxm.^2))...
        ./((dxm*dxp.^2)+(dxm.^2*dxp)); % assumes dxm and dxp are on different sides and absolute
    FDO2V2= @(Dfm,Dfp,dxm,dxp) (Dfp.*dxm.^2-Dfm.*dxp.^2)...
        ./((dxp*dxm.^2)-(dxp.^2*dxm));
    
    nModes=length(modestruct);
    
    gradF_curr=zeros([nModes,1]);
    gradF_m1=zeros([nModes,1]);
    
    for ii=1:nModes
        
        currInd=modestruct(ii).curr.ind;
        currCoeff=modestruct(ii).curr.coeff;
        [gradF_curr(ii)]=GenerateGradientEntry(currInd,currCoeff,...
            gradstruct_curr,FDO2V2);
        
        m1Ind=modestruct(ii).m1.ind;
        m1Coeff=modestruct(ii).m1.coeff;
        [gradF_m1(ii)]=GenerateGradientEntry(m1Ind,m1Coeff,...
            gradstruct_m1,FDO2V2);
    end
    
end

function [gradF]=GenerateGradientEntry(ind,coeff,gradstruct,FDO2V2)
    
    gradF=0;
    if ~isempty(ind)
        if numel(ind)==1
            gradF=gradstruct(ind(1)).objective/coeff(1);
        else % Additional functionality would be to build recursive process
            % to find right indices of currInd
            Dfm=gradstruct(ind(1)).objective;
            Dfp=gradstruct(ind(2)).objective;
            dxm=coeff(1);
            dxp=coeff(2);
            gradF=FDO2V2(Dfm,Dfp,dxm,dxp);
        end
    end
    
end

function [stepVector,supportOptim]=NewStepDirection(stepAlgo,gradF_curr,...
        gradF_m1,modestruct,iterCurr,direction,supportOptim)
    
    switch stepAlgo
        case 'conjgrad'

             [stepVector,supportOptim]=NewStepDirectionConjGrad(gradF_curr,...
                 gradF_m1,modestruct,iterCurr,direction,supportOptim);
    
        case 'BFGS'
            [stepVector,supportOptim]=NewStepDirectionBFGS(gradF_curr,...
                gradF_m1,modestruct,iterCurr,direction,supportOptim);
    
            
    end
    
end


function [stepVector,supportOptim]=NewStepDirectionConjGrad(gradF_curr,gradF_m1,modestruct,iterCurr,direction,supportOptim)
    
    normVec=@(v) sqrt(sum(v.^2,2));
    if isempty(supportOptim)
    prevStep=zeros(size(iterCurr(1).fill));
    prevStep(iterCurr(1).optimdat.var)=iterCurr(1).optimdat.value;
    else
        prevStep=supportOptim.prevStep;
    end
    
    modes=vertcat(modestruct(:).mode)';
    gradF_curr(isnan(gradF_curr))=0;
    gradF_m1(isnan(gradF_m1))=0;
    gradDes_curr=(modes*gradF_curr)';
    gradDes_m1=(modes*gradF_m1)';
    
    switch direction
        case 'min'
            signD=-1;
        case 'max'
            signD=1;
    end
    %scale=(normVec(gradDes_curr)/normVec(prevStep))^2;
    scale=dot(gradDes_curr,gradDes_curr-gradDes_m1)/dot(gradDes_m1,gradDes_m1);
    if ~isfinite(scale)
        scale=1;
    end
    
    stepVector=signD*gradDes_curr+scale*prevStep;
    supportOptim.prevStep=stepVector;
    
end

function [stepVector,supportOptim]=NewStepDirectionBFGS(gradF_curr,gradF_m1,modestruct,iterCurr,direction,supportOptim)
    
    normVec=@(v) sqrt(sum(v.^2,2));
    
    
    modes=vertcat(modestruct(:).mode)';
    gradF_curr(isnan(gradF_curr))=0;
    gradF_m1(isnan(gradF_m1))=0;
    gradDes_curr=(modes*gradF_curr)';
    gradDes_m1=(modes*gradF_m1)';
    if isempty(supportOptim)
        
        Bk=eye(numel(gradDes_curr));
        Bkinv=eye(numel(gradDes_curr));
        gradDes_m1=gradDes_m1*0;
        prevDir=ones([1,numel(gradDes_curr)]);
    else
        Bk=supportOptim.Bk;
        Bkinv=supportOptim.Bkinv;
        
        prevDir=supportOptim.prevDir;
        prevDir=ones(size(iterCurr(1).fill));
        prevDir(iterCurr(1).optimdat.var)=iterCurr(1).optimdat.value;
    end
    
    switch direction
        case 'min'
            signD=-1;
        case 'max'
            signD=1;
    end
    %scale=(normVec(gradDes_curr)/normVec(prevStep))^2;
    yk=(gradDes_curr-gradDes_m1)';
    sk=prevDir';
    Bk=Bk+(yk*yk')/(yk'*sk)-(Bk*sk*sk'*Bk)/(sk'*Bk*sk);
    Bkinv=Bkinv+(sk'*yk+yk'*Bkinv*yk)*(sk*sk')/(sk'*yk)^2 ...
        -(Bkinv*yk*sk'+sk*yk'*Bkinv)/(sk'*yk);
    
    stepVector=(Bkinv*signD*gradDes_curr')';
    
    supportOptim.prevDir=stepVector;
    supportOptim.Bk=Bk;
    supportOptim.Bkinv=Bkinv;
    
end

function [stepLengths]=StepLengthsForLS(rootPop,stepVector,unitStep,validVol,desVarRange)
    
    
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
    stepLengths=unitStep*stepHigh;
    
end

function [stepLengths]=UnitStepLength(worker,lineSearchType)
    
    switch lineSearchType
        case 'backbisection'
            stepLengths=zeros([1,worker]);
            stepLengths(end)=1;
            
            for ii=worker-1:-1:2
                stepLengths(ii)=stepLengths(ii+1)/2;
            end
        case 'polydistrib'
            coeffs=@(vol,v0,g1,g0) [g1+g0-2*vol,3*vol-2*g0-g1,g0,v0]';
            Yfunc=@(x,coeffs) [x.^3 x.^2 x.^1 x.^0]*coeffs;
            stepLengths=Yfunc(linspace(0,1,worker)',coeffs(1,0,1/3,0))';
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
        [newPop(ii,:)]=OverflowHandling(paramoptim,newPop(ii,:));
    end
    deltas{1}=deltas{end};
end

function [stepVector,validVol,diffStepSize]=FindOptimalStepVector(...
        iterstruct,unitSteps,direction,validVol,diffStepSize,minDiffStep)
    
    f=[iterstruct(:).objective];
    g=[iterstruct(:).constraint];
    vec=zeros(size(iterstruct(1).fill));
    vec(iterstruct(1).optimdat.var)=iterstruct(1).optimdat.value;
    
    invalidPoints=g<0.9;
    
    switch direction
        case 'min'
            [~,bestPoint]=min(f);
        case 'max'
            [~,bestPoint]=max(f);
    end
    
    stepLengths=unitSteps;
    stepLength=stepLengths(bestPoint);
    %{
    if bestPoint==1
        bestPoint=2;
    end
    if bestPoint==numel(f)
        bestPoint=numel(f)-1;
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
    %stepLength=iTest(indexLoc);
    %}
    
    
    if stepLength==0
        warning('Step Length is stagnant this iteration')
        validVol=validVol*stepLengths(end-2);
    end
    if (stepLength==0) || (validVol*stepLength<10*max(abs(diffStepSize)))
        
        diffStepSize=sign(diffStepSize).*max(abs(diffStepSize)/4,minDiffStep);
        
    end
    volMulti=1+round((bestPoint)/numel(unitSteps)*2-1)/2;
    
    validVol=validVol*volMulti;
    stepVector=vec*stepLength;
    
end

function [newRoot,deltaRoot]=FindOptimalRestartPop(iterstruct,direction)
    
    f=[iterstruct(:).objective];
    
    switch direction
        case 'min'
            [targObj,indexLoc]=min(f);
        case 'max'
            [targObj,indexLoc]=max(f);
    end
    
    
    newRoot=iterstruct(indexLoc).fill;
    deltaRoot={[1;0]};
    
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

function [newRoot,deltas]=GenerateNewRootFill(rootFill,stepVector,paramoptim,baseGrid)
    
    newRoot=rootFill+stepVector;
    popVec.fill=newRoot;
    popVec=ApplySymmetry(paramoptim,popVec);
    [popVec]=ConstraintMethod('DesVar',paramoptim,popVec,baseGrid);
    newRoot=popVec.fill;
    
    [newRoot]=OverflowHandling(paramoptim,newRoot);
    
    realStep=newRoot-rootFill;
    
    actVar=find(realStep~=0);
    deltas{1}=[actVar;realStep(actVar)];
    
end

function [newGradPop,deltas]=GenerateNewGradientPop(rootFill,desVarRange,stepSize,desVarList)
    
    nActVar=length(desVarList);
    
    rootFillMaxUp=(1-rootFill)/2;
    rootFillMaxDown=-(rootFill)/2;
    
    for jj=1:length(stepSize)
        
        ratioUp=rootFillMaxUp/stepSize(jj);
        ratioUp(ratioUp<0)=1;
        ratioUp(ratioUp==0)=-0.1;
        ratioDown=rootFillMaxDown/stepSize(jj);
        ratioDown(ratioDown==0)=-0.1;
        if stepSize(jj)>=0;
            ratioStep=min([ratioUp;ones(size(ratioUp))]);
        else
            ratioStep=min([ratioDown;ones(size(ratioUp))]);
        end
        
        newGradPop{jj}=(ones(nActVar,1)*rootFill);
        
        %negMov=find(rootFill(desVarList)<=(min(desVarRange)+stepSize(jj)));
        signMat=eye(nActVar);
        for ii=1:nActVar
            signMat(ii,ii)=ratioStep(desVarList(ii));
        end

        partialSteps=signMat*stepSize(jj);
        newGradPop{jj}(:,desVarList)=newGradPop{jj}(:,desVarList)+partialSteps;
        
        
        
        for ii=nActVar:-1:1
            actVar=find(partialSteps(ii,:)~=0);
            deltas{jj}{ii}=[desVarList(actVar);partialSteps(ii,actVar)];
        end
    end
    
end

%% various

%{
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

%}

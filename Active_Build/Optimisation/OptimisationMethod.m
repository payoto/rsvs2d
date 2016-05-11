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

function [newPop,iterCurr,paramoptim]=OptimisationMethod(paramoptim,varargin)
    % Function distributing the optimisation to various methods
    
    
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
            [newPop,iterCurr,paramoptim]=ConjugateGradient(paramoptim,varargin{1},varargin{2});
        case 'conjgradls'
            [newPop,iterCurr,paramoptim]=ConjugateGradientLS(paramoptim,varargin{1},varargin{2});
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

function [newRootFill,desVarList]=OverflowHandling(paramoptim,desVarList,newRootFill)
    % Function which handles steps overflowing the fill constraint
    
    varExtract={'desVarRange','varOverflow','notDesInd'};
    [desVarRange,varOverflow,notDesInd]=ExtractVariables(varExtract,paramoptim);
    
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

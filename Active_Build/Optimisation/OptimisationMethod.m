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
        case 'none'
            iterCurr=varargin{1};
            newPop=vertcat(iterCurr.fill);
            
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
    nNonFill=length(iterCurr(1).nonfillvar);
    warning('Use of cellLevels here, needs to be deprecated')
    newPop=zeros([nPop,nFill+nNonFill]);
    
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
            geneTypeList={'single'};
        case 'horz'
            geneTypeList={'horz'};
        case 'vert'
            geneTypeList={'vert'};
        case 'horzvert'
            geneTypeList={'horz','vert'};
        case 'all'
            geneTypeList={'horz','vert','single'};

    end
    
    for jj=1:numel(geneTypeList)
        switch geneTypeList{jj}
            case 'single'
                nDes=nFill;
                for ii=nFill:-1:1
                    indexmap{jj}{ii}=ii;
                end
                
            case 'horz'
                nDes=cellLevels(2);
                indices=reshape(1:nFill,cellLevels)';
                for ii=nDes:-1:1
                    indexmap{jj}{ii}=indices(ii,:);
                end
                
            case 'vert'
                nDes=cellLevels(1);
                indices=reshape(1:nFill,cellLevels);
                for ii=nDes:-1:1
                    indexmap{jj}{ii}=indices(ii,:);
                end
        end
    end
    
    if nPop>=4
        for ii=1:nPop
            % Mutation
            rInd=randperm(nPop-1,3);
            rInd(rInd>=ii)=rInd(rInd>=ii)+1;
            currRandVecs=[vertcat(iterCurr(rInd).fill),vertcat(iterCurr(rInd).nonfillvar)];
            currVec=[(iterCurr(ii).fill),(iterCurr(ii).nonfillvar)];
            mutVec=projFunc(projInv(currRandVecs(1,:))+diffAmp*...
                (projInv(currRandVecs(2,:))-projInv(currRandVecs(3,:))));
            %mutVec=mutVec*(max(desVarRange)-min(desVarRange))+min(desVarRange);
            %         mutVec(mutVec>max(desVarRange))=max(desVarRange);
            %         mutVec(mutVec<min(desVarRange))=min(desVarRange);
            % Crossover
            crossVec=-ones(size(currVec));
            
            indMapChoice=randi(numel(indexmap),1);
            [fromMutVecLog]=ExtractDEIndices(numel(indexmap{indMapChoice}),...
                nFill,nNonFill,CR,indexmap{indMapChoice});
            
            crossVec(fromMutVecLog)=mutVec(fromMutVecLog);
            crossVec(~fromMutVecLog)=currVec(~fromMutVecLog);
            
            [crossVec]=OverflowHandling(paramoptim,crossVec);
            newPop(ii,:)=crossVec;
            
        end
    else
        warning('Population size is too small for DE - new Pop will be same as old')
        newPop=[vertcat(iterCurr(:).fill),vertcat(iterCurr(:).nonfillvar)];
    end
    
end

function [fromMutVecLog]=ExtractDEIndices(nDes,nFill,nNonFill,CR,indexmap)
    
    fixInd=randi(nDes,1);
    fromMutVecLogDes=(rand([1,nDes])<=CR);
    fromMutVecLogDes(fixInd)=true;
    fromMutVecLogNonFill=(rand([1,nNonFill])<=CR);
    fromMutVecLog=false([1,nFill]);
    
    fromMutVecLog([indexmap{fromMutVecLogDes}])=true;
    fromMutVecLog=[fromMutVecLog,fromMutVecLogNonFill];
end

%% Conjugate gradient

function [newPop,iterOrig,paramoptim,deltas]=ConjugateGradient(paramoptim,iterCurr,iterm1,baseGrid)
    
    varExtract={'diffStepSize','direction','notDesInd','desVarRange',...
        'lineSearch','nLineSearch','nPop','validVol','varActive','desvarconnec',...
        'isRestart','borderActivation','lineSearchType','minDiffStep',...
        'stepAlgo','minVol','gradScale'};
    
    [diffStepSize,direction,notDesInd,desVarRange,lineSearch,nLineSearch,...
        nPop,validVol,varActive,desvarconnec,isRestart,borderActivation,...
        lineSearchType,minDiffStep,stepAlgo,minVol,gradScale]...
        =ExtractVariables(varExtract,paramoptim);
    
    
    supportOptim=paramoptim.optim.supportOptim;
    iterOrig=iterCurr;
    % Extract previous iteration information
    [iterCurr,validCurr]=ExtractValidIter(iterCurr);
    [constrCurr]=ConstraintDistanceExtraction(paramoptim,iterCurr,baseGrid);
    [iterm1,validM1]=ExtractValidIter(iterm1);
    [constrm1]=ConstraintDistanceExtraction(paramoptim,iterm1,baseGrid);
    [gradstruct_curr]=GetIterationInformation(iterCurr,constrCurr);
    [gradstruct_m1]=GetIterationInformation(iterm1,constrm1);
    
    rootPop=[iterCurr(1).fill,iterCurr(1).nonfillvar];
    prevStep=[iterCurr(1).fill,iterCurr(1).nonfillvar]-[iterm1(1).fill,iterm1(1).nonfillvar];
    obj_curr=iterCurr(1).objective;
    obj_m1=iterm1(1).objective;
    lll=0;
    precRoot=[];
    while isempty(precRoot) && lll<numel(iterCurr)
        lll=lll+1;precRoot=iterCurr(lll).location;
    end
    % Case dependant statements
    if lineSearch
        if ~isRestart
            [stepVector,validVol,diffStepSize]=FindOptimalStepVector(iterOrig,...
                UnitStepLength(nLineSearch,lineSearchType),direction,...
                validVol,diffStepSize,minDiffStep,minVol);
            [newRoot,deltaRoot]=GenerateNewRootFill(rootPop,stepVector,paramoptim,baseGrid,precRoot);
        else
            [newRoot,deltaRoot]=FindOptimalRestartPop(iterCurr,direction);
            paramoptim=SetVariables({'isRestart'},{false},paramoptim);
        end
        
        % Need to build function for activation and deactivation of variables
        [inactiveVar]=SelectInactiveVariables(newRoot,varActive,desvarconnec,borderActivation);
        [desVarList]=ExtractActiveVariable(length([iterCurr(1).fill,iterCurr(1).nonfillvar]),notDesInd,inactiveVar);
        [newGradPop,deltaGrad]=GenerateNewGradientPop(newRoot,desVarRange,diffStepSize,desVarList);
        newPop=[newRoot;vertcat(newGradPop{:})];
        deltas=[deltaRoot,deltaGrad{:}];
        validVol
        
        paramoptim=SetVariables({'validVol','diffStepSize','lineSearch'},...
            {validVol,diffStepSize,false},paramoptim);
    else % Direction Search
        
        % Get component change
        % The assumption is the the modes have been selected sensibly and
        % will be close to orthogonal
        [modestruct]=ExtractModes(gradstruct_curr,gradstruct_m1);
        [modestruct]=RemoveFailedModes(modestruct,gradstruct_curr,[iterCurr(1).fill,iterCurr(1).nonfillvar]...
            ,desVarRange,direction);
        [gradF_curr,gradF_m1,gradAdd_curr,gradAdd_m1]...
            =BuildGradientVectors(gradstruct_curr,gradstruct_m1,modestruct,supportOptim);
        
        gradScale=gradScale/sum(abs(diffStepSize)); % Scale the gradients by the size of the steps taken.
        [gradDes_curr,gradDes_m1]=GradFtoGradDes(gradF_curr,gradF_m1,modestruct,gradScale);
        [gradDes_curr]=ConstraintScaleGradient(gradDes_curr,gradAdd_curr,constrCurr,direction,validVol);
        [gradDes_m1]=ConstraintScaleGradient(gradDes_m1,gradAdd_m1,constrm1,direction,validVol);
        % Get Corresponding design vector direction
        
        [stepVector,supportOptim]=NewStepDirection(stepAlgo,gradDes_curr,...
            gradDes_m1,prevStep,direction,supportOptim,obj_curr,obj_m1);
        % Generate Linesearch Distances
        rootPop=[iterCurr(1).fill,iterCurr(1).nonfillvar];
        [unitSteps]=UnitStepLength(nLineSearch,lineSearchType);
        
        [stepLengths]=StepLengthsForLS(rootPop,stepVector,unitSteps,validVol,desVarRange);
        
        [newPop,deltas]=...
            GenerateNewLineSearchPop(rootPop,stepVector,stepLengths,paramoptim,precRoot);
        % Population trimming for invalid values
        
        % Declare linesearch
        paramoptim=SetVariables({'lineSearch'},{true},paramoptim);
    end
    
    % Securing the end
    [nPop,~]=size(newPop);
    paramoptim=SetVariables({'nPop'},{nPop},paramoptim);
    paramoptim.optim.supportOptim=supportOptim;
end

function [gradF]=ConstraintScaleGradient(gradF,gradAdd,constr,direction,validVol)
    s=@(x,eps) sqrt(1-((max(0,min(x,eps))-eps)/eps).^2); % elliptical equation governing the scale
    
    switch direction
        case 'min'
            dir=-1;
        case 'max'
            dir=1;
    end
    
    dist=(constr(1).desvarconstr);
    distGrad=gradAdd.desvarconstr;
    
    scaleGrad=min(((((dir*ones(size(distGrad,1),1)*sign(gradF)).*sign(distGrad))>=0)...
        +(ones(size(distGrad,1),1)*s(dist,validVol))),1);
    scaleGrad=min(scaleGrad,[],1);
    gradF=gradF.*scaleGrad;
end

function [constrRoot]=ConstraintDistanceExtraction(paramoptim,population,baseGrid)
    
    constrRoot=repmat(struct('desvarconstr',...
        ones([1,numel([population(1).fill,population(1).nonfillvar])])),[1,numel(population)]);
    [~,constrDistance]=ConstraintMethod('DesVar',paramoptim,population,baseGrid);
    
    if ~isempty(constrDistance{1})
        constrRoot=constrDistance{1};
        for ii=2:length(constrDistance)
            for jj=1:numel(population)
                constrRoot(jj).desvarconstr=min(constrRoot(jj).desvarconstr,constrDistance{1}(jj).desvarconstr);
            end
        end
        
    end
end

function [population,validIter]=ExtractValidIter(population)
    
    validIter=[1,find([population(2:end).constraint]>0.8)+1];
    population=population(validIter);
    
end

function [gradientopt]=GetIterationInformation(population,constr)
    
    
    
    nPop=length(population);
    nGrad=nPop-1;
    nFill=length([population(1).fill,population(1).nonfillvar]);
    
    rootPop=[population(1).fill,population(1).nonfillvar];
    rootObj=population(1).objective;
    rootAdd=catstruct(population(1).additional,constr(1));
    rootAdd.desvar=rootPop;
    
    gradientopt=struct('design',zeros(1,nFill),'objective',0,'fill',...
        zeros(1,nFill),'additional',rootAdd);
    gradientopt=repmat(gradientopt,[1,nGrad]);
    
    fieldsAdd=fieldnames(population(1).additional);
    fieldsConstr=fieldnames(constr(1));
    for ii=1:nGrad
        
        gradientopt(ii).design(population(ii+1).optimdat.var)...
            =population(ii+1).optimdat.value;
        gradientopt(ii).objective=population(ii+1).objective-rootObj;
        gradientopt(ii).fill=[population(ii+1).fill,population(ii+1).nonfillvar]-rootPop;
        for jj=1:numel(fieldsAdd)
            gradientopt(ii).additional.(fieldsAdd{jj})=...
                population(ii+1).additional.(fieldsAdd{jj})-rootAdd.(fieldsAdd{jj});
        end
        for jj=1:numel(fieldsConstr)
            gradientopt(ii).additional.(fieldsConstr{jj})=...
                constr(ii+1).(fieldsConstr{jj})-rootAdd.(fieldsConstr{jj});
        end
        gradientopt(ii).additional.desvar=[population(ii+1).fill,population(ii+1).nonfillvar]-rootPop;
    end
    
    % Remove modes which don't have any fill change
    % This can happen because of the approximate refinement of designs when
    % the snake is calculated through sensitivity.
    gradientopt=gradientopt(any(vertcat(gradientopt(:).design),2));
    
    
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
        modeCoeff(ii)=sign(allModes(ii,find(allModes(ii,:)~=0,1,'first')));
    end
    modeMultiplier=modeCoeff*ones(size(allModes(1,:)));
    allModes=allModes./modeMultiplier;
    %     modeCoeffNorm=normVec(allModes);
    %     modeMultiplier=modeCoeffNorm*ones(size(allModes(1,:)));
    %     allModes=allModes./modeMultiplier;
    %     modeCoeff=modeCoeff.*modeCoeffNorm;
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

function [modestruct]=ExtractModes_OLD(gradstruct_curr,gradstruct_m1)
    
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
    % Removes modes which are "failing" from the gradients. This is to
    % avoid one discontinuity dominating the gradients and causing the
    % optimisation to stall. This is used to avoid "The opening of new
    % cells" to be detrimental
    
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

function [gradF_curr,gradF_m1,gradAdd_curr,gradAdd_m1]=...
        BuildGradientVectors(gradstruct_curr,gradstruct_m1,modestruct,supportOptim)
    % This build includes a gradient rejection criteria
    
    
    FDO2= @(fi,fm,fp,dxm,dxp) (fp.*dxm.^2-fm.*dxp.^2+fi.*(dxp.^2-dxm.^2))...
        ./((dxm*dxp.^2)+(dxm.^2*dxp)); % assumes dxm and dxp are on different sides and absolute
    FDO2V2= @(Dfm,Dfp,dxm,dxp) (Dfp.*dxm.^2-Dfm.*dxp.^2)...
        ./((dxp*dxm.^2)-(dxp.^2*dxm));
    normVec=@(v) sqrt(sum(v.^2,2));
    
    
    
    nModes=length(modestruct);
    if ~isempty(supportOptim);
        
        prevGradNorm=(normVec(supportOptim.hist(end).gradfk));
    else
        prevGradNorm=1;
    end
    
    gradF_curr=zeros([nModes,1]);
    gradF_m1=zeros([nModes,1]);
    gradAdd_curr=repmat(opstruct('*',gradstruct_curr(1).additional,0),[nModes,1]);
    gradAdd_m1=gradAdd_curr;
    
    for ii=1:nModes
        
        currInd=modestruct(ii).curr.ind;
        currCoeff=modestruct(ii).curr.coeff;
        [gradF_curr(ii),gradAdd_curr(ii)]=GenerateGradientEntry(currInd,currCoeff,...
            gradstruct_curr,FDO2V2,prevGradNorm,gradAdd_curr(ii));
        
        m1Ind=modestruct(ii).m1.ind;
        m1Coeff=modestruct(ii).m1.coeff;
        [gradF_m1(ii),gradAdd_m1(ii)]=GenerateGradientEntry(m1Ind,m1Coeff,...
            gradstruct_m1,FDO2V2,prevGradNorm,gradAdd_m1(ii));
    end
    gradAdd_curr=opstruct('vertcat',gradAdd_curr);
    gradAdd_m1=opstruct('vertcat',gradAdd_m1);
end


function [gradF,gradAdd]=GenerateGradientEntry(ind,coeff,gradstruct,FDO2V2,prevGradNorm,gradAdd)
    
    isOutLier=@(gradVec,prevGradNorm) log10(abs(gradVec))> (max(mean(log10(abs(gradVec))),log10(prevGradNorm))+1);   % Takes in absolute values
    
    gradF=0;
    if ~isempty(ind)
        gradTest=zeros(size(ind));
        for ii=1:numel(ind)
            gradTest(ii)=gradstruct(ind(ii)).objective/coeff(ii);
        end
        validStep=~isOutLier(gradTest,prevGradNorm);
        ind=ind(validStep);
        coeff=coeff(validStep);
        
        if numel(ind)==1
            gradF=gradstruct(ind(1)).objective/coeff(1);
            gradAdd=opstruct('/',gradstruct(ind(1)).additional,coeff(1));
        else % Additional functionality would be to build recursive process
            % to find right indices of currInd
            
            
            Dfm=gradstruct(ind(1)).objective;
            Dfp=gradstruct(ind(2)).objective;
            dxm=coeff(1);
            dxp=coeff(2);
            gradF=FDO2V2(Dfm,Dfp,dxm,dxp);
            gradAdd=opstruct(FDO2V2,gradstruct(ind(1)).additional,gradstruct(ind(2)).additional,dxm,dxp);
        end
    else
        %disp('Warning: Empty Mode')
    end
    
end

function [gradDes_curr,gradDes_m1]=GradFtoGradDes(gradF_curr,gradF_m1,modestruct,gradScale)
    
    modes=vertcat(modestruct(:).mode)';
    gradF_curr(isnan(gradF_curr))=0;
    gradF_m1(isnan(gradF_m1))=0;
    gradDes_curr=((modes*gradF_curr)'.*gradScale);
    gradDes_m1=((modes*gradF_m1)'.*gradScale);
end

function [stepVector,supportOptim]=NewStepDirection(stepAlgo,gradDes_curr,...
        gradDes_m1,prevStep,direction,supportOptim,obj_curr,obj_m1)
    
    
    
    switch stepAlgo
        case 'conjgrad'
            
            [stepVector,supportOptim]=NewStepDirectionConjGrad(gradDes_curr,...
                gradDes_m1,prevStep,direction,supportOptim,obj_curr,obj_m1);
            
        case 'BFGS'
            [stepVector,supportOptim]=NewStepDirectionBFGS(gradDes_curr,...
                gradDes_m1,prevStep,direction,supportOptim,obj_curr,obj_m1);
            
            
    end
    
end


function [stepVector,supportOptim]=NewStepDirectionConjGrad(gradDes_curr,gradDes_m1,...
        prevStep,direction,supportOptim,obj_curr,obj_m1)
    
    normVec=@(v) sqrt(sum(v.^2,2));
    
    
    isEmptSupport=isempty(supportOptim);
    if isEmptSupport
        prevDir=zeros(size(gradDes_m1));
    else
        prevDir=supportOptim.curr.prevStep;
    end
    
    if isEmptSupport
        ii=1;
    else
        ii=numel(supportOptim.hist)+1;
    end
    supportOptim.hist(ii).gradfk=gradDes_curr;
    supportOptim.hist(ii).gradfkm1=gradDes_m1;
    supportOptim.hist(ii).prevDir =prevDir;
    supportOptim.hist(ii).prevStep=prevStep;
    supportOptim.hist(ii).obj_curr=obj_curr;
    supportOptim.hist(ii).obj_m1=obj_m1;
    [supportOptim.hist(ii).isWolfe,supportOptim.hist(ii).isStrongWolfe]=...
        WolfeCondition(1e-4,0.9,prevDir,prevStep,...
        gradDes_m1,obj_m1,gradDes_curr,obj_curr);
    switch direction
        case 'min'
            signD=-1;
        case 'max'
            signD=1;
    end
    %scale=(normVec(gradDes_curr)/normVec(prevStep))^2;
    scale=dot(gradDes_curr,gradDes_curr-gradDes_m1)/dot(gradDes_m1,gradDes_m1); % Polak Ribiere
    if ~isfinite(scale) || all(prevStep==0)
        scale=0;
    end
    supportOptim.hist(ii).scale=scale;
    stepVector=signD*gradDes_curr+scale*prevDir;
    supportOptim.curr.prevStep=stepVector;
    
end

function [stepVector,supportOptim]=NewStepDirectionBFGS(gradDes_curr,gradDes_m1,...
        prevStep,direction,supportOptim,obj_curr,obj_m1)
    
    normVec=@(v) sqrt(sum(v.^2,2));
    
    isEmptSupport=isempty(supportOptim);
    if isEmptSupport
        
        Bk=eye(numel(gradDes_curr));
        Bkinv=eye(numel(gradDes_curr));
        gradDes_m1=gradDes_m1*0;
        prevDir=ones([1,numel(gradDes_curr)]);
        iter=0;
    else
        Bk=supportOptim.curr.Bk;
        Bkinv=supportOptim.curr.Bkinv;
        iter=supportOptim.curr.iter;
        prevDir=supportOptim.curr.prevDir;
        %prevDir=ones(size([iterCurr(1).fill,iterCurr(1).nonfillvar]));
        %prevDir(iterCurr(1).optimdat.var)=iterCurr(1).optimdat.value;
    end
    
    if isEmptSupport
        ii=1;
    else
        ii=numel(supportOptim.hist)+1;
    end
    supportOptim.hist(ii).Bk=Bk;
    supportOptim.hist(ii).Bkinv=Bkinv;
    supportOptim.hist(ii).gradfk=gradDes_curr;
    supportOptim.hist(ii).gradfkm1=gradDes_m1;
    supportOptim.hist(ii).prevDir =prevDir;
    supportOptim.hist(ii).prevStep=prevStep;
    supportOptim.hist(ii).obj_curr=obj_curr;
    supportOptim.hist(ii).obj_m1=obj_m1;
    supportOptim.hist(ii).iter=iter;
    [supportOptim.hist(ii).isWolfe,supportOptim.hist(ii).isStrongWolfe]=...
        WolfeCondition(1e-4,0.9,prevDir,prevStep,...
        gradDes_m1,obj_m1,gradDes_curr,obj_curr);
    switch direction
        case 'min'
            signD=-1;
        case 'max'
            signD=1;
    end
    %scale=(normVec(gradDes_curr)/normVec(prevStep))^2;
    yk=(gradDes_curr-gradDes_m1)';
    sk=prevStep';
    
    
    
    if isEmptSupport
        stepVector=(Bkinv*(signD*gradDes_curr)')';
    else
        Bk=Bk+(yk*yk')/(yk'*sk)-(Bk*(sk*sk')*Bk)/(sk'*Bk*sk);
        Bkinv=Bkinv+(sk'*yk+yk'*Bkinv*yk)*(sk*sk')/(sk'*yk)^2 ...
            -(Bkinv*yk*sk'+sk*yk'*Bkinv)/(sk'*yk);
        
        if  ~supportOptim.hist(ii).isStrongWolfe || iter>8
            %         Bkinv=eye(numel(gradDes_curr));
            %         Bk=eye(numel(gradDes_curr));
            Bkinv=diag(diag(Bkinv));
            Bk=diag(diag(Bk));
            iter=0;
            
        end
        
        if any(isnan(Bkinv(:)))
            Bkinv=eye(numel(gradDes_curr));
            Bk=eye(numel(gradDes_curr));
            iter=0;
        end
        stepVector=(Bkinv*(signD*gradDes_curr)')';
    end
    if any(~isfinite(stepVector))
        warning('stepVector not finite')
    end
    
    supportOptim.curr.prevDir=stepVector;
    supportOptim.curr.Bk=Bk;
    supportOptim.curr.Bkinv=Bkinv;
    supportOptim.curr.iter=iter+1;
    
end


function [isWolfe,isStrongWolfe]=WolfeCondition(c1,c2,prevDir,prevStep,gradfk,fk,gradfkp1,fkp1)
    
    isWolfe1=fkp1<=(fk+c1*dot(prevStep,gradfk));
    isWolfe2=(dot(prevDir,gradfkp1))>=(c2*dot(prevDir,gradfk));
    isWolfe3=abs(dot(prevDir,gradfkp1))<=c2*abs(dot(prevDir,gradfk));
    
    isWolfe=isWolfe1 && isWolfe2;
    isStrongWolfe=isWolfe1 && isWolfe3;
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
        GenerateNewLineSearchPop(rootPop,stepVector,stepLengths,paramoptim,precRoot)
    
    nNew=length(stepLengths);
    nFill=length(rootPop);
    deltaDes=(stepLengths'*stepVector);
    newPop=(ones([nNew,1])*rootPop)+deltaDes;
    
    for ii=nNew:-1:1
        actVar=find(deltaDes(ii,:));
        deltas{ii}=[actVar;deltaDes(ii,actVar)];
        [newPop(ii,:),precRoot]=OverflowHandling(paramoptim,newPop(ii,:),precRoot);
    end
    deltas{1}=deltas{end};
end

function [stepVector,validVol,diffStepSize]=FindOptimalStepVector(...
        iterstruct,unitSteps,direction,validVol,diffStepSize,minDiffStep,minVol)
    
    f=[iterstruct(:).objective];
    g=[iterstruct(:).constraint];
    vec=zeros(size([iterstruct(1).fill,iterstruct(1).nonfillvar]));
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
    
    
    
    if (stepLength==0) || (validVol*stepLength<10*max(abs(diffStepSize)))
        
        diffStepSize=sign(diffStepSize).*max(abs(diffStepSize)/(10),minDiffStep);
        
    end
    if stepLength==0
        warning('Step Length is stagnant this iteration')
        %validVol=max(validVol*stepLengths(end-2),minVol);
        validVol=validVol/4;
    else
        volMulti=1+round((bestPoint)/numel(unitSteps)*2-1)/2;
        validVol=validVol*volMulti;
        validVol=max(validVol,minVol);
    end
    
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
    
    
    newRoot=[iterstruct(indexLoc).fill,iterstruct(indexLoc).nonfillvar];
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

function [newRoot,deltas]=GenerateNewRootFill(rootFill,stepVector,paramoptim,baseGrid,precRoot)
    
    newRoot=rootFill+stepVector;
    [newRoot]=OverflowHandling(paramoptim,newRoot,precRoot);
    popVec.fill=newRoot;
    popVec=ApplySymmetry(paramoptim,popVec);
    [popVec]=ConstraintMethod('DesVar',paramoptim,popVec,baseGrid);
    newRoot=popVec.fill;
    
    
    
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
    %paramoptim.general.nPop=nPop;
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
    %paramoptim.general.nPop=nPop;
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

<<<<<<< HEAD


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

function [inactiveVar]=SelectInactiveVariables(newFill,varActive)
    
    switch varActive
        case 'all'
            inactiveVar=[];
        case 'border'
            error('Not coded yet')
            
        case 'wideborder'
            error('Not coded yet')
            
        otherwise
            error('unrecognised variable activation criterion')
    end
    
    
end

function [newGradPop,deltas]=GenerateNewGradientPop(rootFill,desVarRange,stepSize,desVarList)
    
    nActVar=length(desVarList);
    
    newGradPop=(ones(nActVar,1)*rootFill);
    
    negMov=find(rootFill(desVarList)>=(max(desVarRange)-stepSize));
    signMat=eye(nActVar);
    signMat(:,negMov)=signMat(:,negMov)*-1;
    
    partialSteps=signMat*stepSize;
    newGradPop(:,desVarList)=newGradPop(:,desVarList)+partialSteps;
    
    for ii=nActVar:-1:1
        actVar=find(partialSteps(ii,:)~=0);
        deltas{ii}=[desVarList(actVar);partialSteps(ii,actVar)];
    end
    
end
>>>>>>> 2f6fdfae0c257288b7da037c663953265f6f756a
%}
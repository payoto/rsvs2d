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

function [newPop,iterCurr]=OptimisationMethod(paramoptim,varargin)
    % Function distributing the optimisation to various methods
    
    
    varExtract={'optimMethod'};
    [optimMethod]=ExtractVariables(varExtract,paramoptim);
    
    switch optimMethod
        case 'DE'
            proj='cut';
            [newPop,iterCurr]=DifferentialEvolution(paramoptim,proj,varargin{1},varargin{2});
        case 'DEtan'
            proj='tan';
            [newPop,iterCurr]=DifferentialEvolution(paramoptim,proj,varargin{1},varargin{2});
            
        case 'DEStrip'
            
        case 'GA'
            
        case 'GSA'
            
        case 'SQP'
            
    end
    
end


%% Optimisation Methods

function [newPop,iterCurr]=DifferentialEvolution(paramoptim,proj,iterCurr,iterm1)
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
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
    varExtract={'desVarRange','direction'};
    [desVarRange,direction]=ExtractVariables(varExtract,paramoptim);
    
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
        fixInd=randi(nFill,1);
        fromMutVecLog=(rand([1,nFill])<=CR);
        fromMutVecLog(fixInd)=true;
        
        crossVec(fromMutVecLog)=mutVec(fromMutVecLog);
        crossVec(~fromMutVecLog)=iterCurr(ii).fill(~fromMutVecLog);
        
        newPop(ii,:)=crossVec;
        
    end
    
end
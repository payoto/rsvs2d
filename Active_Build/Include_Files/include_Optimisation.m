function [] = include_Optimisation()
    %FUNCTIONLIST allows local functions to be used globally once it has
    %been used.
    
    funcHandles=localfunctions;
    funcDir=[cd,'\Automated_Function_Directory'];
    HeaderActivation(funcHandles,funcDir)
    
end


%%

function [desVarList]=ExtractActiveVariable(nFill,notDesInd,inactiveVar)
    
    desVarList=1:nFill;
    desVarList([notDesInd,inactiveVar])=[];
    
end

function [inactiveVar]=SelectInactiveVariables(newFill,varActive,derivtenscalc,cutoff)
    
    switch varActive
        case 'all'
            inactiveVar=[];
        case 'border'
            [inactiveVar]=InactiveVariables_border(newFill);
            
        case 'wideborder'
            [inactiveVar,activeVar]=InactiveVariables_wideborder(newFill,derivtenscalc,cutoff);
            
        case 'snaksensiv'
            
            [inactiveVar]=InactiveVariables_border(newFill);
        otherwise
            error('unrecognised variable activation criterion')
    end
    
    
end

function [inactiveVar,activeVar]=InactiveVariables_border(newFill)
    
    is0=newFill==0;
    is1=newFill==1;
    
    inactiveVar=find(is0);
    inactiveVar=[inactiveVar,find(is1)];
    
    isBord=~is0 & ~is1;
    activeVar=find(isBord);
    
    
end

function [inactiveVar,activeVar]=InactiveVariables_wideborder(newFill,derivtenscalc,cutoff)
    
    [inactiveVar,activeVar]=InactiveVariables_border(newFill);
    
    indexDeriv=[derivtenscalc(:).index];
    actVarSub=FindObjNum([],activeVar,indexDeriv);
    actFill=newFill(activeVar);
    cutActVar=actFill>cutoff;
    activeVar=unique([activeVar,[derivtenscalc(actVarSub(cutActVar)).neighbours]]);
    inactiveVar=1:length(newFill);
    inactiveVar(activeVar)=[];
end

function [isGradient]=CheckIfGradient(optimMethod)
    
    switch optimMethod
        
        case 'DE'
            isGradient=false;
        case 'DEtan'
            isGradient=false;
        case 'conjgrad'
            isGradient=true;
        case 'conjgradls'
            isGradient=true;
        otherwise
            isGradient=false;
            warning('Optimisation method is not known as gradient based or otherwise, no gradient is assumed')
            
    end
end

function [isSnakeSensitivity]=CheckSnakeSensitivityAlgorithm(paramoptim)
    
    varExtract={'varActive','optimMethod'};
    [varActive,optimMethod]=ExtractVariables(varExtract,paramoptim);
    
    isSnakeSensitivity=false;
    [isGradient]=CheckIfGradient(optimMethod);
    
    if isGradient
        isSnakeSensitivity=strcmp(varActive,'snaksensiv');
    end
end

function [newRootFill]=OverflowHandling(paramoptim,newRootFill)
    % Function which handles steps overflowing the fill constraint
    
    varExtract={'desVarRange','varOverflow'};
    [desVarRange,varOverflow]=ExtractVariables(varExtract,paramoptim);
    
    switch varOverflow
        case 'truncate'
            minD=min(desVarRange);
            maxD=max(desVarRange);
            newRootFill(newRootFill<minD)=minD;
            newRootFill(newRootFill>maxD)=maxD;
            
        case 'spill'
            newFill=zeros(size(newRootFill));
            for ii=1:size(newRootFill,1)
                [newFill(ii,:)]=SpillOverflow(paramoptim,newRootFill(ii,:));
            end
            newRootFill=newFill;
            minD=min(desVarRange);
            maxD=max(desVarRange);
            newRootFill(newRootFill<minD)=minD;
            newRootFill(newRootFill>maxD)=maxD;
        otherwise
            error('No valid design variable overflow mechanism')
            
    end
    
    
    
end

function [newFill]=SpillOverflow(paramoptim,newRootFill)
    
    varExtract={'desVarRange','desvarconnec'};
    [desVarRange,desvarconnec]=ExtractVariables(varExtract,paramoptim);
    
    normRootFill=(newRootFill-min(desVarRange))/(max(desVarRange)-min(desVarRange));
    
    overVar=find(normRootFill>1);
    underVar=find(normRootFill<0);
    exitFlag=false;
    kk=0;
    while (~isempty(overVar) || ~isempty(underVar)) && ~exitFlag
        
        
        [normRootFill]=SpillOverflowVarHandling(normRootFill,desvarconnec,overVar);
        
        normRootFill=1-normRootFill;
        [normRootFill]=SpillOverflowVarHandling(normRootFill,desvarconnec,underVar);
        normRootFill=1-normRootFill;
        
        
        
        overVar=find(normRootFill>1);
        underVar=find(normRootFill<0);
        meanFill=mean(normRootFill);
        stdFill=std(normRootFill);
        kk=kk+1;
        if ((meanFill>=1 || meanFill<=0)  && (abs(stdFill)<1e-5)) || kk>100
            exitFlag=true;
        end
        
    end
    
    newFill=normRootFill*(max(desVarRange)-min(desVarRange))+min(desVarRange);
end

function [newRootFill]=SpillOverflowVarHandling(newRootFill,desvarconnec,flowVar)
    
    desVarIndList=[desvarconnec(:).index];
    
    overFlowMat=zeros([length(flowVar),length(newRootFill)]);
    
    for ii=1:length(flowVar)
        
        currSub=FindObjNum([],flowVar(ii),desVarIndList);
        
        neighInd=desvarconnec(currSub).neighbours;
        cornInd=desvarconnec(currSub).corners;
        
        neighSub=FindObjNum([],neighInd,desVarIndList);
        cornSub=FindObjNum([],cornInd,desVarIndList);
        
        currVol=newRootFill(currSub);
        neighVol=newRootFill(neighSub);
        cornVol=newRootFill(cornSub);
        
        neighEmpt=neighSub(neighVol<=0);
        cornEmpt=cornSub(cornVol<=0);
        neighGreyCell={[]};
        for jj=1:length(cornEmpt)
            neighGreyCell{jj}=FindObjNum([],desvarconnec(cornEmpt(jj)).neighbours,neighInd)';
        end
        neighAll=[neighGreyCell{:}];
        neighSubAll=neighSub(RemoveIdenticalEntries(neighAll(neighAll~=0)));
        neighSubAll=RemoveIdenticalEntries([neighEmpt;neighSubAll]);
        neighSubAll=neighSubAll(newRootFill(neighSubAll)<1);
        
        nCorn=numel(cornEmpt);
        nNeigh=sum(1-newRootFill(neighSubAll));
        isNormProb=true;
        
        if nCorn>0
            baseRate=(-nNeigh+sqrt(nNeigh^2+4*nCorn))/(2*nCorn);
            
        elseif nNeigh>0
            baseRate=1/nNeigh;
            
        elseif numel(neighSub(newRootFill(neighSub)<1))>0 ...
                && numel(cornSub(newRootFill(cornSub)<1))>0
            
            neighSubAll=neighSub(newRootFill(neighSub)<1);
            cornEmpt=cornSub(newRootFill(cornSub)<1);
            nCorn=numel(cornEmpt);
            nNeigh=sum(1-newRootFill(neighSubAll));
            
            baseRate=(-nNeigh+sqrt(nNeigh^2+4*nCorn))/(2*nCorn);
            
        elseif numel(neighSub(newRootFill(neighSub)<1))>0
            neighSubAll=neighSub(newRootFill(neighSub)<1);
            nNeigh=sum(1-newRootFill(neighSubAll));
            baseRate=1/nNeigh;
            
        elseif numel(cornSub(newRootFill(cornSub)<1))>0
            neighSubAll=cornSub(newRootFill(cornSub)<1);
            nNeigh=sum(1-newRootFill(neighSubAll));
            baseRate=1/nNeigh;
            
        elseif numel(neighSub(newRootFill(neighSub)<currVol))>0 || ...
                numel(cornSub(newRootFill(cornSub)<currVol))>0
            neighSubAll=[neighSub(newRootFill(neighSub)<currVol);cornSub(newRootFill(cornSub)<currVol)];
            nNeigh=sum(currVol-newRootFill(neighSubAll));
            baseRate=1/nNeigh;
            
            isNormProb=false;
        else
            neighSubAll=[];
            nNeigh=sum(-newRootFill(neighSubAll));
            
            baseRate=0;
            
            isNormProb=false;
        end
        
        if isNormProb
            overFlowVol=zeros(size(newRootFill));
            overFlowVol(neighSubAll)=(1-newRootFill(neighSubAll))*baseRate;
            overFlowVol(cornEmpt)=baseRate^2;
            overFlowVol(currSub)=-1;
            overFlowMat(ii,:)=overFlowVol*(currVol-1);
            %newRootFill=newRootFill+overFlowVol*(currVol-1);
        else
            overFlowVol=zeros(size(newRootFill));
            overFlowVol(neighSubAll)=baseRate*(currVol-newRootFill(neighSubAll));
            overFlowVol(cornEmpt)=baseRate^2;
            overFlowVol(currSub)=-sum(overFlowVol);
            overFlowMat(ii,:)=overFlowVol*(currVol-max([min(newRootFill(neighSubAll)),1]))*2/3;
            %newRootFill=newRootFill+overFlowVol*(currVol-max([min(newRootFill(neighSubAll)),1]))*2/3;
        end
        
    end
    
    newRootFill=newRootFill+sum(overFlowMat,1);
    
end


function population=ApplySymmetry(paramoptim,population)
    
    varExtract={'symDesVarList'};
    [symDesVarList]=ExtractVariables(varExtract,paramoptim);
    
    for ii=1:length(population)
        population(ii).fill(symDesVarList(2,:))=...
            population(ii).fill(symDesVarList(1,:));
    end
    
    
end











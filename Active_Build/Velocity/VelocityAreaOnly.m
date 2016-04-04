%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%              University of Bristol
%      Department of Aerospace Engineering
%                     2015
%
%          Subdivision of Surfaces
%      for Aerodynamic shape parametrisation
%                - Snakes -
%        - Velocity calculation -
%               Using Area 
%             Alexandre Payot
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%#codegen

function [snaxel,snakposition,snaxelmodvel]=VelocityAreaOnly(snaxel,snakposition,volumefraction,coeffstructure,forceparam)
    
    [velAverage]=CalculateAverageVelocities(volumefraction,coeffstructure,forceparam);
    [snaxelvel]=DistributeVelocityToSnaxel(velAverage);
    
    [velstruct]=ConvertToVelStructure(snaxelvel);
    [newcelltoold]=StructureNewToOldcell(volumefraction);
    [velstruct]=CalculateDeviationVelocity(velstruct,coeffstructure,newcelltoold,volumefraction);
    [snaxel,snaxelvel2]=AssignVelocityToSnaxel(velstruct,snaxel);
    
    snaxelmodvel=snaxel;
    for ii=1:length(snaxel)
        snaxelmodvel(ii).v=snaxelvel2(ii).average;
    end
end

function [velstruct]=ConvertToVelStructure(snaxelvel)
    
    nVel=length([snaxelvel(:).cellindex]);
    velstruct(nVel).index=[];
    
    velSub=1;
    for ii=1:length(snaxelvel)
        for jj=1:length(snaxelvel(ii).cellindex)
            velstruct(velSub).index=velSub;
            velstruct(velSub).cellindex=snaxelvel(ii).cellindex(jj);
            velstruct(velSub).snaxelindex=snaxelvel(ii).index;
            %velstruct(velSub).tensvel=snaxelvel(ii).tensvel;
            velstruct(velSub).averagevel=snaxelvel(ii).averagevel(jj);
            velSub=velSub+1;
        end
    end
    
end

function [velAverage]=CalculateAverageVelocities2(volumefraction,coeffstructure)
    % Calculates the average velocities required to match the condition
    [coeffCellInd]=[coeffstructure(:).cellindex];
    [coeffSnaxInd]=[coeffstructure(:).snaxelindex];
    
    isNeg=zeros(length(volumefraction),1);
    is0=zeros(length(volumefraction),1);
    for ii=length(volumefraction):-1:1
        
        coeffSubs=FindObjNum([],volumefraction(ii).newCellInd,coeffCellInd);
        velAverage(ii).oldCellInd=volumefraction(ii).oldCellInd;
        velAverage(ii).deltaArea=volumefraction(ii).targetfill*volumefraction(ii).totalvolume...
            -volumefraction(ii).totalfraction;
        
        if sum(coeffSubs)==0
            
        else
            coeffSubs(coeffSubs==0)=[];
            velAverage(ii).sumcoeff=sum([coeffstructure(coeffSubs).value]);
            
            [velAverage(ii).velocity,isNeg(ii),is0(ii),velAverage(ii).sumcoeff]...
                =LogicBlockForAverageVel(velAverage(ii).sumcoeff,...
                velAverage(ii).deltaArea);
            
            snaxIndList=coeffSnaxInd(coeffSubs);
            velAverage(ii).snaxelindexlist=RemoveIdenticalEntries(snaxIndList);
        end
    end
    if sum(is0)>0
        disp(['    The sum of coefficients was equal to zero ',num2str(sum(is0)),' Times'])
    end
    if sum(isNeg)>0
        disp(['    The sum of coefficients was negative ',num2str(sum(isNeg)),' Times'])
    end
end

function [velocity,isNeg,is0,sumCoeff]=LogicBlockForAverageVel...
        (sumCoeff,deltaArea)
    % Logic block for the calculation of the average vel required for the
    % counteracting of the tensile velocities
    isNeg=0;
    is0=0;
    if abs(sumCoeff)<10^-15
        sumCoeff=0;
    end
    if sumCoeff<0
        %warning('The sum of derivative coefficients is negative ignoring value (average velocity will not drive convergence)')
        isNeg=1;
        sumCoeff=1;
    end
    if sumCoeff==0
        is0=1;
        %warning('The sum of derivative coefficients is 0 ignoring value (average velocity will not drive convergence)')
        
        velocity=(deltaArea)/1;
    else
        velocity=(deltaArea)/sumCoeff;
    end
    
    
end

function [velstruct]=CalculateDeviationVelocity(velstruct,coeffstructure,newcelltoold,volumefraction)
    
    % [A;C]*v_dev+[B;0]*v_bar=0
    
    velEqualMat=VelocityEqualityConditions(velstruct);
    deviationCondMat=DeviationVelocityNoAreaCondition(velstruct,coeffstructure,newcelltoold);
    subRmvRows=find(sum(abs(deviationCondMat),2)==0);
    deviationCondMat(subRmvRows,:)=[];
    [fitConditionInd]=Force0DeviationVelOnEmpty(volumefraction,velstruct,velEqualMat);
    
    vBar=[velstruct(:).averagevel];
    condMatrix=[velEqualMat;deviationCondMat];
    [nCond,nUnknown]=size(condMatrix);
    nSet=nUnknown-nCond;
    [condMatrix2,rmvColumns]=RemoveUnknowns(0,condMatrix,fitConditionInd);
    
    condSatisfiedMat=condMatrix(:,rmvColumns);
    condMatrix=condMatrix2;
    
    ControlMat=[velEqualMat;zeros(size(deviationCondMat))];
    controlVec=ControlMat*vBar';
    
    vDev=pinv(condMatrix)*(-controlVec);
    %vDev=(condMatrix)\(-controlVec);
    ll=1;
    for ii=1:length(velstruct)
        if sum(ii==rmvColumns)>0
            vDevComplete(ii)=0;
        else
            vDevComplete(ii)=vDev(ll);
            ll=ll+1;
        end
    end
    for ii=1:length(velstruct)
        velstruct(ii).deviationvel=vDevComplete(ii);
    end
    
end

function [condMatrix,rmvColumns]=RemoveUnknowns(nSet,condMatrix,fitConditionInd)
    % Removes columns associated with possible to remove variables without
    % removing an equation
    
    condMatrixWorking=condMatrix;
    fitConditionIndWorking=fitConditionInd;
    rmvColumns=[];
    for ii=1:nSet
        condition=true;
        ll=1;
        forceExit=length(fitConditionIndWorking)>ll;
        while condition && forceExit
            
            indWork=fitConditionIndWorking(ll);
            logicalElmCond=(sum(...
                abs(condMatrixWorking(:,[1:indWork-1,indWork+1:end]))~=0,2)~=0)...
                ==(sum(abs(condMatrix)~=0,2)~=0);
            
            condition=prod(logicalElmCond)==0; % true when index not to be removed
            
            if ~condition
                rmvColumns(ii)=indWork;
                fitConditionIndWorking(ll)=[];
                condMatrixWorking(:,indWork)=0;
            end
            ll=ll+1;
            forceExit=length(fitConditionIndWorking)>ll;
        end
        if ~forceExit
            disp('No More column may be removed')
            break
        end
        
    end
    condMatrix(:,rmvColumns)=[];
    
end

function [velEqualMat]=VelocityEqualityConditions(velstruct)
    % Calculates the matrix governing which velocities need to be equal
    snaxelIndList=[velstruct(:).snaxelindex];
    nVel=length(velstruct);
    
    velEqualMat=zeros(nVel);
    
    
    for ii=1:nVel
        otherSideVel=FindObjNum([],velstruct(ii).snaxelindex,snaxelIndList);
        
        otherSideVel(otherSideVel<=ii)=[];
        if numel(otherSideVel)>0
            velEqualMat(ii,ii)=1;
            velEqualMat(ii,otherSideVel)=-1;
        end
        
        
    end
    velEqualMat(find(sum(abs(velEqualMat),2)==0),:)=[];
    
end

function [newcelltoold]=StructureNewToOldcell(volumefraction)
    
    newIndexCell=[volumefraction(:).newCellInd];
    newcelltoold(length(newIndexCell)).index=[];
    ll=1;
    
    for ii=1:length(volumefraction)
        for jj=1:numel(volumefraction(ii).newCellInd)
            [newcelltoold(ll).newCell]=volumefraction(ii).newCellInd(jj);
            [newcelltoold(ll).oldCell]=volumefraction(ii).oldCellInd;
            ll=ll+1;
        end
    end
    
end

function [deviationCondMat]=DeviationVelocityNoAreaCondition(velstruct,coeffstructure,newcelltoold)
    
    velCellInd=[velstruct(:).cellindex];
    coeffNewCellInd=[coeffstructure(:).cellindex];
    newToOldNewInd=[newcelltoold(:).newCell];
    oldIndexSub=FindObjNum([],coeffNewCellInd,newToOldNewInd);
    coeffOldCellInd=[newcelltoold(oldIndexSub).oldCell];
    
    velSnaxelInd=[velstruct(:).snaxelindex];
    coeffSnaxelInd=[coeffstructure(:).snaxelindex];
    cellIndUnique=RemoveIdenticalEntries(velCellInd);
    
    
    nCell=length(cellIndUnique);
    nVel=length(velstruct);
    
    
    hashKeyVel=[velCellInd',velSnaxelInd'];
    hashKeyCoeff=[coeffOldCellInd',coeffSnaxelInd'];
    
    [positionHash]=CompareHashKeys(hashKeyCoeff,hashKeyVel);
    
    for ii=1:nVel
        velstruct(ii).velcoeff=sum([coeffstructure(positionHash{ii}).value]);
    end
    deviationCondMat=zeros([nCell,nVel]);
    for ii=1:nVel
        cellSubs=FindObjNum([],velstruct(ii).cellindex,cellIndUnique);
        deviationCondMat(cellSubs,ii)=velstruct(ii).velcoeff;
    end
    
end

function [snaxelvel]=DistributeVelocityToSnaxel(velaverage)
    % Distributes velocities to a snaxel centred structure
    snaxelList=[velaverage(:).snaxelindexlist];
    snaxelList=RemoveIdenticalEntries(snaxelList);
    
    for ii=length(snaxelList):-1:1
        snaxelvel(ii).index=snaxelList(ii);
        snaxelvel(ii).averagevel=[];
        snaxelvel(ii).averagevelbase=[];
        snaxelvel(ii).cellindex=[];
    end
    
    for ii=1:length(velaverage)
        snaxelVelSub=FindObjNum([],velaverage(ii).snaxelindexlist,snaxelList);
        for jj=1:length(snaxelVelSub)
            [snaxelvel(snaxelVelSub(jj)).cellindex(end+1)]=velaverage(ii).oldCellInd;
            [snaxelvel(snaxelVelSub(jj)).averagevelbase(end+1)]=velaverage(ii).velocity;
            
            snaxelvel(snaxelVelSub(jj)).averagevel(end+1)=...
                snaxelvel(snaxelVelSub(jj)).averagevelbase(end);
        end
    end
    
end

function [velAverage]=CalculateAverageVelocities(volumefraction,coeffstructure,forceparam)
    % Calculates the average velocities required to match the condition
    [coeffCellInd]=[coeffstructure(:).cellindex];
    [coeffSnaxInd]=[coeffstructure(:).snaxelindex];
    isNeg=zeros(length(volumefraction),1);
    is0=zeros(length(volumefraction),1);
    
    for ii=length(volumefraction):-1:1
        
        coeffSubs=FindObjNum([],volumefraction(ii).newCellInd,coeffCellInd);
        velAverage(ii).oldCellInd=volumefraction(ii).oldCellInd;
        velAverage(ii).deltaArea=volumefraction(ii).targetfill*volumefraction(ii).totalvolume...
            -volumefraction(ii).totalfraction;
        
        if sum(coeffSubs)==0
            
        else
            coeffSubs(coeffSubs==0)=[];
            
            coeffsWorking=[coeffstructure(coeffSubs).value];
            snaxIndList=coeffSnaxInd(coeffSubs);
            
            sumcoeff=sum(coeffsWorking);
            
            [velAverage(ii).velocity,isNeg(ii),is0(ii),sumcoeff]...
                =LogicBlockForAverageVel(sumcoeff,...
                velAverage(ii).deltaArea);
            
            velAverage(ii).snaxelindexlist=RemoveIdenticalEntries(snaxIndList);
        end
    end
    if sum(is0)>0
        disp(['    The sum of coefficients was equal to zero ',num2str(sum(is0)),' Times'])
    end
    if sum(isNeg)>0
        disp(['    The sum of coefficients was negative ',num2str(sum(isNeg)),' Times'])
    end
end

function [fitConditionInd]=Force0DeviationVelOnEmpty(volumefraction,velstruct,velEqualMat)
    
    oldCellVolFrac=[volumefraction(:).oldCellInd];
    oldCellVelStruct=[velstruct(:).cellindex];
    cellVellStructInVolFrac=FindObjNum([],oldCellVelStruct,oldCellVolFrac);
    
    logicFill0=[volumefraction(cellVellStructInVolFrac).targetfill]==0;
    logicAveVelis0=[velstruct(:).averagevel]==0;
    logicBorderVel=sum(abs(velEqualMat))>0;
    
    fitConditionLogic= logicFill0 & logicAveVelis0 & logicBorderVel;
    fitConditionInd=[velstruct(fitConditionLogic).index];
    
    
end

function [snaxel,snaxelvel]=AssignVelocityToSnaxel(velstruct,snaxel)
    % Assigns velocities to the snaxel structure
    
    snaxInd=[velstruct(:).snaxelindex];
    
    for ii=1:length(snaxel)
        
        velSub=FindObjNum([],snaxel(ii).index,snaxInd);
        
        snaxelVel=[velstruct(velSub).averagevel]+[velstruct(velSub).deviationvel];
        if numel(snaxelVel)>1
            if abs(snaxelVel(1)-snaxelVel(2))>10^-12
                disp(['Difference in Vel is ',num2str(snaxelVel(1)-snaxelVel(2))])
                warning('Averaging difference in velocities')
            end
        end
        snaxelvel(ii).index=snaxel(ii).index;
        snaxelvel(ii).average=velstruct(velSub(1)).averagevel;
        snaxelvel(ii).deviation=velstruct(velSub(1)).deviationvel;
        %snaxelvel(ii).force=snaxeltensvel(velTensSub).forcevel;
        snaxel(ii).v=mean(snaxelVel);
        
        if abs(snaxel(ii).v)<10^-12
            snaxel(ii).v=0;
        end
    end
    
    
end

function [positionHash]=CompareHashKeys(hashList,hashTest)
    % Matches the hashkeys fron hashTest to hashList returning the position
    % on in hashList
    
    [nTest,nDim]=size(hashTest);
    positionHash{nTest}=[];
    for jj=1:nTest
        hashLogical=false(size(hashList));
        for ii=1:nDim
            hashLogical(:,ii)=hashList(:,ii)==hashTest(jj,ii);
        end
        positionHash{jj}=find(prod(hashLogical,2));
    end
    
end

%% Find Correct position single step


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
%      Using Area and geometry forcing
%             Alexandre Payot
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [snaxel,snakposition,snaxelmodvel]=VelocityForce_FE(snaxel,snakposition,volumefraction,coeffstructure,forceparam)
    
    [snaxeltensvel,snakposition]=GeometryForcingVelocity(snaxel,snakposition,forceparam);
    [velAverage]=CalculateAverageVelocities(volumefraction,coeffstructure,snaxeltensvel,forceparam);
    [snaxelvel]=DistributeVelocityToSnaxel(velAverage,snaxeltensvel);
    
    [velstruct]=ConvertToVelStructure(snaxelvel);
    [newcelltoold]=StructureNewToOldcell(volumefraction);
    [velstruct]=CalculateDeviationVelocity(velstruct,coeffstructure,newcelltoold,volumefraction);
    [snaxel,snaxelvel2]=AssignVelocityToSnaxel(velstruct,snaxel,snaxelvel);
    
    snaxelmodvel=snaxel;
    for ii=1:length(snaxel)
        snaxelmodvel(ii).v=snaxeltensvel(ii).bendforce;
        %snaxel(ii).v=snaxeltensvel(ii).forcevel;
        snaxel(ii).acc=snaxeltensvel(ii).forceacc;
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

function [snaxelvel]=DistributeVelocityToSnaxel(velaverage,snaxeltensvel)
    % Distributes velocities to a snaxel centred structure
    snaxelList=[velaverage(:).snaxelindexlist];
    snaxelList=RemoveIdenticalEntries(snaxelList);
    snaxTensList=[snaxeltensvel(:).index];
    
    for ii=length(snaxelList):-1:1
        snaxelvel(ii).index=snaxelList(ii);
        snaxelvel(ii).averagevel=[];
        snaxelvel(ii).averagevelbase=[];
        snaxelvel(ii).averagevelforce=[];
        snaxelvel(ii).forcevel=[];
        snaxelvel(ii).cellindex=[];
    end
    
    for ii=1:length(velaverage)
        snaxelVelSub=FindObjNum([],velaverage(ii).snaxelindexlist,snaxelList);
        snaxelTensVelSub=FindObjNum([],velaverage(ii).snaxelindexlist,snaxTensList);
        for jj=1:length(snaxelVelSub)
            [snaxelvel(snaxelVelSub(jj)).cellindex(end+1)]=velaverage(ii).oldCellInd;
            %% Watch modification
            [snaxelvel(snaxelVelSub(jj)).averagevelbase(end+1)]=velaverage(ii).velocity;
            %%
            [snaxelvel(snaxelVelSub(jj)).averagevelforce(end+1)]=velaverage(ii).averagevelforce;
            [snaxelvel(snaxelVelSub(jj)).forcevel]=snaxeltensvel(snaxelTensVelSub(jj)).forcevel*velaverage(ii).forcevelcoeff;
            
            snaxelvel(snaxelVelSub(jj)).averagevel(end+1)=...
                snaxelvel(snaxelVelSub(jj)).averagevelbase(end)...
                +snaxelvel(snaxelVelSub(jj)).averagevelforce(end);
        end
    end
    
end

function [velAverage]=CalculateAverageVelocities(volumefraction,coeffstructure,snaxelvel,forceparam)
    % Calculates the average velocities required to match the condition
    [coeffCellInd]=[coeffstructure(:).cellindex];
    [coeffSnaxInd]=[coeffstructure(:).snaxelindex];
    snaxVelInd=[snaxelvel(:).index];
    isNeg=zeros(length(volumefraction),1);
    is0=zeros(length(volumefraction),1);
    
    maxVelRatio=forceparam.maxVelRatio;
    
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
            snaxSubWorking=FindObjNum([],snaxIndList,snaxVelInd);
            tensVelWorking=snaxelvel(snaxSubWorking).forcevel;
            
            tensVelCoeffProd=sum(coeffsWorking.*tensVelWorking);
            sumcoeff=sum(coeffsWorking);
            
            [velAverage(ii).velocity,isNeg(ii),is0(ii),sumcoeff]...
                =LogicBlockForAverageVel(sumcoeff,...
                velAverage(ii).deltaArea);
            [velAverage(ii).sumcoeff,velAverage(ii).forcevelCoeffProd,...
                velAverage(ii).averagevelforce]=LogicBlockForAverageTensVel...
                (sumcoeff,tensVelCoeffProd,tensVelWorking);
            velAverage(ii).forcevelcoeff=1;
%             if abs(velAverage(ii).averagevelforce)>maxVelRatio;%abs(velAverage(ii).velocity)*
%                 newAverageTensVel=(sign(velAverage(ii).averagevelforce))*maxVelRatio; %velAverage(ii).velocity*
%                 velAverage(ii).forcevelcoeff=newAverageTensVel/velAverage(ii).averagevelforce;
%                 velAverage(ii).averagevelforce=newAverageTensVel;
%             end
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

function [sumcoeff,tensVelCoeffProd,averageveltens]=LogicBlockForAverageTensVel...
        (sumcoeff,tensVelCoeffProd,tensVelWorking)
    % Logic block for the calculation of the average vel required for the
    % counteracting of the tensile velocities
    
    if abs(sumcoeff)<10^-15
        sumcoeff=0;
    end
    if abs(tensVelCoeffProd)<10^-15
        tensVelCoeffProd=0;
    end
    if sumcoeff<0
        %warning('The sum of derivative coefficients is negative ignoring value (average velocity will not drive convergence)')
        sumcoeff=1;
    end
    
    if sumcoeff==0 && tensVelCoeffProd==0
        %warning('The sum of derivative coefficients is 0 ignoring value (average velocity will not drive convergence)')
        averageveltens=0;% -sum(tensVelWorking);
    elseif sumcoeff==0 && tensVelCoeffProd~=0
        averageveltens=0; %-tensVelCoeffProd;
        warning('This is weird')
    else
        averageveltens=-tensVelCoeffProd/...
            sumcoeff;
        
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

function [snaxel,snaxelvel]=AssignVelocityToSnaxel(velstruct,snaxel,snaxeltensvel)
    % Assigns velocities to the snaxel structure
    
    snaxInd=[velstruct(:).snaxelindex];
    snaxTensInd=[snaxeltensvel(:).index];
    
    for ii=1:length(snaxel)
        
        velSub=FindObjNum([],snaxel(ii).index,snaxInd);
        velTensSub=FindObjNum([],snaxel(ii).index,snaxTensInd);
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
        snaxelvel(ii).force=snaxeltensvel(velTensSub).forcevel;
        snaxel(ii).v=mean(snaxelVel)+snaxeltensvel(velTensSub).forcevel;
        
        
    end
    
    
end

%% Geometry feature forcing

function [snaxeltensvel,snakposition]=GeometryForcingVelocity(snaxel,snakposition,forceparam)
    % Calculate a tensile velocity to influence the solution
    snakPosIndex=[snakposition(:).index];
    snaxIndex=[snaxel(:).index];
    
    if sum(snakPosIndex~=snaxIndex)>0
        error('Indices do not match')
    end
    [snaxeltensvel,snakposition]=CalculateTensileVelocity(snaxel,snakposition,snakPosIndex);
    [snaxeltensvel,snakposition]=CalculateBendingVelocity(snaxel,snakposition,snakPosIndex,snaxeltensvel);
    [derivtenscalc]=CalculateTensileVelocity2(snaxel,snakposition,snakPosIndex);
    [implicitMatTens,implicitMatBend]=BuildImplicitMatrix(derivtenscalc);
    
    %[bendforce]=CalculateBendForce(snaxel,snakposition,snakPosIndex);
    [bendforce]=CalculateShearForce_fullderiv(snaxel,snakposition,snakPosIndex);
    [bendforce]=SnaxelToForce_curvature(snaxel,snakposition,bendforce,forceparam);
    
    bendingVelInfluence=forceparam.bendingVelInfluence;
    tensVelInfluence=forceparam.tensVelInfluence;    
    maxForceVel=forceparam.maxForceVel;
 
    forceVelScaling=maxForceVel/(bendingVelInfluence+tensVelInfluence);
    tensCoeff=forceVelScaling*tensVelInfluence;
    bendCoeff=forceVelScaling*bendingVelInfluence;
    
    %% This is not zorking
    [tensVelVec]=GeometryForcingVelCalc(derivtenscalc,implicitMatTens,implicitMatBend,tensCoeff,bendCoeff);
    %%
    for ii=1:length(snaxeltensvel)
        %snaxeltensvel(ii).forcevel= snaxeltensvel(ii).tensvel*tensCoeff+ snaxeltensvel(ii).bendvel*bendCoeff;
        snaxeltensvel(ii).bendforce=bendforce(ii).d2pds2; % watch the change of force
        snaxeltensvel(ii).forcevel=bendforce(ii).vel;
        snaxeltensvel(ii).forceacc=bendforce(ii).acc;
    end
    
end

function [snaxeltensvel,snakposition]=CalculateTensileVelocity(snaxel,snakposition,snakPosIndex)
    
    
    rotCW=[0 -1; 1 0];
    rotCCW=[0 1; -1 0];
    velRatio=1/2; % limits the velocity to 0.1
    for ii=1:length(snakposition)
        neighbourSub=FindObjNum([],[snaxel(ii).snaxprec,snaxel(ii).snaxnext],snakPosIndex);
        coordNeighbour=vertcat(snakposition(neighbourSub).coord);
        testVecs=coordNeighbour-([1;1]*snakposition(ii).coord);
%         testPrec=norm(testVecs(1,:))>0;
%         testNext=norm(testVecs(2,:))>0;
%         dirVecs(1,:)=rotCW*snakposition(ii).vectorprec'*testPrec;
%         dirVecs(2,:)=rotCCW*snakposition(ii).vectornext'*testNext;
%         dirVecNorm=[norm(dirVecs(1,:));norm(dirVecs(2,:))];
%         dirVecNorm(dirVecNorm==0)=1;
%         
%         unitDirVecs=dirVecs./(dirVecNorm*[1 1]);
        
        snakposition(ii).tensVector=sum(testVecs);
        snaxeltensvel(ii).index=snakposition(ii).index;
        snaxeltensvel(ii).tensvel=dot(snakposition(ii).tensVector,snakposition(ii).vector)/(2*velRatio);
    end
    
end

function [implicitMatTens,implicitMatBend]=BuildImplicitMatrix(derivtenscalc)
    numCoeff=length(derivtenscalc);
    
    implicitMatTens=zeros(numCoeff);
    implicitMatBend=implicitMatTens;
    implicitMatBendm1=implicitMatTens;
    implicitMatBendp1=implicitMatTens;
    % Shifted indices lists
%     shiftM2=[numCoeff-1:numCoeff,1:numCoeff-2];
%     shiftM1=[numCoeff,1:numCoeff-1];
%     shiftP2=[3:numCoeff,1:2];
%     shiftP1=[2:numCoeffoeff,1];
    
    snaxPrecList=[derivtenscalc(:).snaxprec];
    snaxNextList=[derivtenscalc(:).snaxnext];
    snaxIndList=[derivtenscalc(:).index];
    snaxPrecSub=FindObjNum([],snaxPrecList,snaxIndList);
    snaxPrec2Sub=snaxPrecSub(snaxPrecSub);
    snaxNextSub=FindObjNum([],snaxNextList,snaxIndList);
    snaxNext2Sub=snaxNextSub(snaxNextSub);
    
    for ii=1:length(derivtenscalc)
        implicitMatTens(ii,ii)=derivtenscalc(ii).velcoeff_i;
        implicitMatTens(ii,snaxPrecSub(ii))=derivtenscalc(ii).velcoeff_m;
        implicitMatTens(ii,snaxNextSub(ii))=derivtenscalc(ii).velcoeff_p;
    end
    
    for ii=1:length(derivtenscalc)
        % Base
        implicitMatBend(ii,ii)=derivtenscalc(ii).velcoeff_i;
        implicitMatBend(ii,snaxPrecSub(ii))=derivtenscalc(ii).velcoeff_m;
        implicitMatBend(ii,snaxNextSub(ii))=derivtenscalc(ii).velcoeff_p;
        %p1
        implicitMatBendp1(ii,snaxNextSub(ii))=derivtenscalc(snaxNextSub(ii)).velcoeff_i;
        implicitMatBendp1(ii,snaxPrecSub(snaxNextSub(ii)))=derivtenscalc(snaxNextSub(ii)).velcoeff_m;
        implicitMatBendp1(ii,snaxNextSub(snaxNextSub(ii)))=derivtenscalc(snaxNextSub(ii)).velcoeff_p;
        %m1
        implicitMatBendm1(ii,snaxPrecSub(ii))=derivtenscalc(snaxPrecSub(ii)).velcoeff_i;
        implicitMatBendm1(ii,snaxPrecSub(snaxPrecSub(ii)))=derivtenscalc(snaxPrecSub(ii)).velcoeff_m;
        implicitMatBendm1(ii,snaxNextSub(snaxPrecSub(ii)))=derivtenscalc(snaxPrecSub(ii)).velcoeff_p;
    end
    implicitMatBend=-2*implicitMatBend+implicitMatBendm1+implicitMatBendp1;
end

function [derivtenscalc]=CalculateTensileVelocity2(snaxel,snakposition,snakPosIndex)
    
    
    rotCW=[0 -1; 1 0];
    rotCCW=[0 1; -1 0];
    velRatio=1; % limits the velocity to 0.1
    for ii=length(snakposition):-1:1
        neighSub=FindObjNum([],[snaxel(ii).snaxprec,snaxel(ii).snaxnext],snakPosIndex);
        
        derivtenscalc(ii).index=snakposition(ii).index;
        derivtenscalc(ii).snaxnext=snaxel(ii).snaxnext;
        derivtenscalc(ii).snaxprec=snaxel(ii).snaxprec;
        % extracting data from preexisting arrays
        derivtenscalc(ii).dir_i=snakposition(ii).vectornotnorm;
        derivtenscalc(ii).dir_p=snakposition(neighSub(2)).vectornotnorm;
        derivtenscalc(ii).dir_m=snakposition(neighSub(1)).vectornotnorm;
        derivtenscalc(ii).dir_norm=snakposition(ii).vector;
        derivtenscalc(ii).pos_i=snakposition(ii).coord;
        derivtenscalc(ii).pos_p=snakposition(neighSub(2)).coord;
        derivtenscalc(ii).pos_m=snakposition(neighSub(1)).coord;
        derivtenscalc(ii).tension=dot(snakposition(ii).tensVector,snakposition(ii).vector);
        bendingVec=2*snakposition(ii).tensVector-(snakposition(neighSub(1)).tensVector+snakposition(neighSub(2)).tensVector);
        derivtenscalc(ii).bending=dot(bendingVec,snakposition(ii).vector);
        % calculating data
        derivtenscalc(ii).Dvec_p=derivtenscalc(ii).pos_p-derivtenscalc(ii).pos_i;
        derivtenscalc(ii).Dvec_m=derivtenscalc(ii).pos_m-derivtenscalc(ii).pos_i;
        derivtenscalc(ii).Dnorm_p=norm(derivtenscalc(ii).Dvec_p);
        derivtenscalc(ii).Dnorm_m=norm(derivtenscalc(ii).Dvec_m);
        
        [velcoeff_i,velcoeff_p,velcoeff_m]=TensileVelDerivCoeff(derivtenscalc(ii));
        
        derivtenscalc(ii).velcoeff_i=velcoeff_i;
        derivtenscalc(ii).velcoeff_p=velcoeff_p;
        derivtenscalc(ii).velcoeff_m=velcoeff_m;
        
    end
    
end

function [velcoeff_i,velcoeff_p,velcoeff_m]=TensileVelDerivCoeff(derivtenscalc)
    
    fields=fieldnames(derivtenscalc);
    
    for ii=1:length(fields)
        eval([fields{ii},'=derivtenscalc.',fields{ii},';']);
    end
    
    velcoeff_i=-TensileVelDerivCoeffBreakdown(dir_i,dir_norm,Dnorm_p,Dvec_p)...
        -TensileVelDerivCoeffBreakdown(dir_i,dir_norm,Dnorm_m,Dvec_m);
    
    velcoeff_p=TensileVelDerivCoeffBreakdown(dir_p,dir_norm,Dnorm_p,Dvec_p);
    
    velcoeff_m=TensileVelDerivCoeffBreakdown(dir_m,dir_norm,Dnorm_m,Dvec_m);
    
    
    if isnan(velcoeff_p)
        velcoeff_p=0;
    end
    if isnan(velcoeff_m)
        velcoeff_m=0;
    end
    if isnan(velcoeff_i)
        velcoeff_i=-velcoeff_p-velcoeff_m;
    end
end

function [velcoeff]=TensileVelDerivCoeffBreakdown(Dg_ipm,Dgi_unit,Dpnorm,Dp_mp)

    % Equation for normalised bending
%     velcoeff=Dpnorm*((dot(Dg_ipm,Dgi_unit)/Dpnorm)...
%         -((dot(Dp_mp,Dgi_unit)/(Dpnorm^3))...
%         *sum(Dp_mp.*Dg_ipm)));

% Equation for actual bending
velcoeff=dot(Dg_ipm,Dgi_unit);
    
    
end

function [snaxeltensvel,snakposition]=CalculateBendingVelocity(snaxel,snakposition,snakPosIndex,snaxeltensvel)
    
    
    
    velRatio=8; % limits the velocity to the same range as the 
    for ii=1:length(snakposition)
        neighbourSub=FindObjNum([],[snaxel(ii).snaxprec,snaxel(ii).snaxnext],snakPosIndex);
        tensVecNeighbour=vertcat(snakposition(neighbourSub).tensVector);

        snakposition(ii).bendVector=-sum(tensVecNeighbour)+2*snakposition(ii).tensVector;
        snaxeltensvel(ii).bendvel=dot(snakposition(ii).bendVector,snakposition(ii).vector)/(velRatio);
    end
    
end

function [tensVelVec]=GeometryForcingVelCalc(derivtenscalc,implicitMatTens,implicitMatBend,tensCoeff,bendCoeff)
    
    T=[derivtenscalc(:).tension]';
    B=[derivtenscalc(:).bending]';
    
    F=T*tensCoeff+B*bendCoeff;
    
    dFdt=implicitMatTens*tensCoeff+implicitMatBend*bendCoeff;
    
    if rcond(dFdt)>1e-10
        tensVelVec=(dFdt)\(-F);
    else
        tensVelVec=pinv(dFdt)*(-F);
    end
    %tensVelVec=tensVelVec/mean(abs(tensVelVec));
    
    
    
end


%% Bending - Tom's equation

function [bendforce]=CalculateBendForce(snaxel,snakposition,snakPosIndex)

    snaxIndList=[snaxel(:).index];
    for ii=length(snakposition):-1:1
        %Find closest neighbour with different cooerdinates
         [precInd]=FindNextNonIdenticalSnaxel(snaxIndList,snakposition(ii).index,...
             snaxel,snakposition,'prec');
         [nextInd]=FindNextNonIdenticalSnaxel(snaxIndList,snakposition(ii).index,...
             snaxel,snakposition,'next');
         simFlag=false;
        if precInd~=snaxel(ii).snaxprec || nextInd~=snaxel(ii).snaxnext
            simFlag=true;
        end
        [bendforce(ii)]=BendForceDataExtraction(ii,nextInd,precInd,snakposition,snakPosIndex,simFlag);
    end
    
    % Calculating second difference    
    for ii=length(bendforce):-1:1

        bendforce(ii).d2pds2=SecondDerivativeVarStencil...
            (bendforce(ii).pos_i,bendforce(ii).pos_p,bendforce(ii).pos_m,...
            bendforce(ii).dist_p,bendforce(ii).dist_m);
        %{
        % extracting data from preexisting arrays
        derivtenscalc(ii).dir_i=snakposition(ii).vectornotnorm;
        derivtenscalc(ii).dir_p=snakposition(neighSub(2)).vectornotnorm;
        derivtenscalc(ii).dir_m=snakposition(neighSub(1)).vectornotnorm;
        derivtenscalc(ii).dir_norm=snakposition(ii).vector;
        derivtenscalc(ii).pos_i=snakposition(ii).coord;
        derivtenscalc(ii).pos_p=snakposition(neighSub(2)).coord;
        derivtenscalc(ii).pos_m=snakposition(neighSub(1)).coord;
        derivtenscalc(ii).tension=dot(snakposition(ii).tensVector,snakposition(ii).vector);
        bendingVec=2*snakposition(ii).tensVector-(snakposition(neighSub(1)).tensVector+snakposition(neighSub(2)).tensVector);
        derivtenscalc(ii).bending=dot(bendingVec,snakposition(ii).vector);
        % calculating data
        derivtenscalc(ii).Dvec_p=derivtenscalc(ii).pos_p-derivtenscalc(ii).pos_i;
        derivtenscalc(ii).Dvec_m=derivtenscalc(ii).pos_m-derivtenscalc(ii).pos_i;
        derivtenscalc(ii).Dnorm_p=norm(derivtenscalc(ii).Dvec_p);
        derivtenscalc(ii).Dnorm_m=norm(derivtenscalc(ii).Dvec_m);
        
        [velcoeff_i,velcoeff_p,velcoeff_m]=TensileVelDerivCoeff(derivtenscalc(ii));
        
        derivtenscalc(ii).velcoeff_i=velcoeff_i;
        derivtenscalc(ii).velcoeff_p=velcoeff_p;
        derivtenscalc(ii).velcoeff_m=velcoeff_m;
        %}
    end
    % Calculate 3rd, edge centred difference
    for ii=length(bendforce):-1:1
        neighSub=FindObjNum([],bendforce(ii).snaxnext,snakPosIndex);
        
        bendforce(ii).edge_d3pds3=(bendforce(neighSub).d2pds2-...
            bendforce(ii).d2pds2)/bendforce(ii).dist_p;
    end
    % Calculate Shear force at node
    for ii=length(bendforce):-1:1
        neighSub=FindObjNum([],bendforce(ii).snaxprec,snakPosIndex);
        
        bendforce(ii).shearforce=(bendforce(neighSub).edge_d3pds3-...
            bendforce(ii).edge_d3pds3);
    end
    for ii=length(bendforce):-1:1
        neighSub=FindObjNum([],[bendforce(ii).snaxprec,bendforce(ii).snaxnext],snakPosIndex);
        
        bendforce(ii).spuriousloading=SecondDerivativeVarStencil...
            (bendforce(ii).d2pds2,bendforce(neighSub(2)).d2pds2,bendforce(neighSub(1)).d2pds2,...
            bendforce(ii).dist_p,bendforce(ii).dist_m);
    end
end

function [bendforce]=CalculateShearForce_fullderiv(snaxel,snakposition,snakPosIndex)

    snaxIndList=[snaxel(:).index];
    for ii=length(snakposition):-1:1
        %Find closest neighbour with different cooerdinates
         [precInd]=FindNextNonIdenticalSnaxel(snaxIndList,snakposition(ii).index,...
             snaxel,snakposition,'prec');
         [nextInd]=FindNextNonIdenticalSnaxel(snaxIndList,snakposition(ii).index,...
             snaxel,snakposition,'next');
         simFlag=false;
        if precInd~=snaxel(ii).snaxprec || nextInd~=snaxel(ii).snaxnext
            simFlag=true;
        end
        [bendforce(ii)]=ShearForceDataExtraction(ii,nextInd,precInd,snakposition,snakPosIndex,simFlag);
    end
    
    % Calculating second difference    
    for ii=length(bendforce):-1:1

        bendforce(ii).d2pds2=SecondDerivativeVarStencil...
            (bendforce(ii).pos_i,bendforce(ii).pos_p,bendforce(ii).pos_m,...
            bendforce(ii).dist_p,bendforce(ii).dist_m);
        %{
        % extracting data from preexisting arrays
        derivtenscalc(ii).dir_i=snakposition(ii).vectornotnorm;
        derivtenscalc(ii).dir_p=snakposition(neighSub(2)).vectornotnorm;
        derivtenscalc(ii).dir_m=snakposition(neighSub(1)).vectornotnorm;
        derivtenscalc(ii).dir_norm=snakposition(ii).vector;
        derivtenscalc(ii).pos_i=snakposition(ii).coord;
        derivtenscalc(ii).pos_p=snakposition(neighSub(2)).coord;
        derivtenscalc(ii).pos_m=snakposition(neighSub(1)).coord;
        derivtenscalc(ii).tension=dot(snakposition(ii).tensVector,snakposition(ii).vector);
        bendingVec=2*snakposition(ii).tensVector-(snakposition(neighSub(1)).tensVector+snakposition(neighSub(2)).tensVector);
        derivtenscalc(ii).bending=dot(bendingVec,snakposition(ii).vector);
        % calculating data
        derivtenscalc(ii).Dvec_p=derivtenscalc(ii).pos_p-derivtenscalc(ii).pos_i;
        derivtenscalc(ii).Dvec_m=derivtenscalc(ii).pos_m-derivtenscalc(ii).pos_i;
        derivtenscalc(ii).Dnorm_p=norm(derivtenscalc(ii).Dvec_p);
        derivtenscalc(ii).Dnorm_m=norm(derivtenscalc(ii).Dvec_m);
        
        [velcoeff_i,velcoeff_p,velcoeff_m]=TensileVelDerivCoeff(derivtenscalc(ii));
        
        derivtenscalc(ii).velcoeff_i=velcoeff_i;
        derivtenscalc(ii).velcoeff_p=velcoeff_p;
        derivtenscalc(ii).velcoeff_m=velcoeff_m;
        %}
    end
    % Calculate 3rd derivative
    for ii=length(bendforce):-1:1
        neighSub=FindObjNum([],[bendforce(ii).snaxprec,bendforce(ii).snaxnext],snakPosIndex);
        
        bendforce(ii).snax_pp=bendforce(neighSub(2)).snaxnext;
        bendforce(ii).snax_mm=bendforce(neighSub(1)).snaxprec;
        
        bendforce(ii).pos_pp=bendforce(neighSub(2)).pos_p;
        bendforce(ii).pos_mm=bendforce(neighSub(1)).pos_m;
        bendforce(ii).dist_pp=bendforce(neighSub(2)).dist_p;
        bendforce(ii).dist_mm=bendforce(neighSub(1)).dist_m;
        
        bendforce(ii).dist_ep=bendforce(ii).dist_p+bendforce(ii).dist_pp;
        bendforce(ii).dist_em=bendforce(ii).dist_m+bendforce(ii).dist_mm;
        
        bendforce(ii).shearforce=ThrirdDerivativeVarStencil(...
            bendforce(ii).pos_i,...
            bendforce(ii).pos_p,...
            bendforce(ii).pos_m,...
            bendforce(ii).pos_pp,...
            bendforce(ii).pos_mm,...
            bendforce(ii).dist_p,...
            bendforce(ii).dist_m,...
            bendforce(ii).dist_ep,...
            bendforce(ii).dist_em...
            );
        
        bendforce(ii).d4pds4=SecondDerivativeVarStencil...
            (bendforce(ii).d2pds2,bendforce(neighSub(1)).d2pds2,bendforce(neighSub(2)).d2pds2,...
            bendforce(ii).dist_p,bendforce(ii).dist_m);
    end

end

function [neighInd]=FindNextNonIdenticalSnaxel(snaxIndList,rootInd,snaxel,snakposition,dirSnax)
% Returns the closest non similar snaxel
% Where dirSnax is a string: either 'prec' or 'next'

snaxSub=FindObjNum([],rootInd,snaxIndList);
rootPos=snakposition(snaxSub).coord;

neighInd=snaxel(snaxSub).(['snax',dirSnax]);
neighSub=FindObjNum([],neighInd,snaxIndList);
neighPos=snakposition(neighSub).coord;

while sum(neighPos~=rootPos)==0
    neighInd=snaxel(neighSub).(['snax',dirSnax]);
    neighSub=FindObjNum([],neighInd,snaxIndList);
    neighPos=snakposition(neighSub).coord;
end

end

function [derivtenscalc]=BendForceDataExtraction(ii,nextInd,precInd,snakposition,snakPosIndex,simFlag)
    % Extracts data for BendForce
    
    derivtenscalc.index=snakposition(ii).index;
    derivtenscalc.snaxnext=nextInd;
    derivtenscalc.snaxprec=precInd;
    derivtenscalc.samepos=simFlag;
    neighSub=FindObjNum([],[precInd,nextInd],snakPosIndex);
    % Extracting positions
    derivtenscalc.pos_i=snakposition(ii).coord;
    derivtenscalc.pos_p=snakposition(neighSub(2)).coord; % p for plus: i+1
    derivtenscalc.pos_m=snakposition(neighSub(1)).coord; % m for minus: i-1

    % Calculating distances
    derivtenscalc.dist_p=sqrt(sum((derivtenscalc.pos_i-derivtenscalc.pos_p).^2));
    derivtenscalc.dist_m=sqrt(sum((derivtenscalc.pos_i-derivtenscalc.pos_m).^2));
end

function [derivtenscalc]=ShearForceDataExtraction(ii,nextInd,precInd,snakposition,snakPosIndex,simFlag)
    % Extracts data for BendForce
    
    derivtenscalc.index=snakposition(ii).index;
    derivtenscalc.snaxnext=nextInd;
    derivtenscalc.snaxprec=precInd;
    derivtenscalc.samepos=simFlag;
    neighSub=FindObjNum([],[precInd,nextInd],snakPosIndex);
    % Extracting positions
    derivtenscalc.pos_i=snakposition(ii).coord;
    derivtenscalc.pos_p=snakposition(neighSub(2)).coord; % p for plus: i+1
    derivtenscalc.pos_m=snakposition(neighSub(1)).coord; % m for minus: i-1

    % Calculating distances
    derivtenscalc.dist_p=sqrt(sum((derivtenscalc.pos_i-derivtenscalc.pos_p).^2));
    derivtenscalc.dist_m=sqrt(sum((derivtenscalc.pos_i-derivtenscalc.pos_m).^2));
end

function [d2pds2]=SecondDerivativeVarStencil(Pi,Pp,Pm,Dp,Dm)
    % Calculation of a second difference with variable stencil
    
    d2pds2=2*(-Pi*(Dp+Dm)+Dp*Pm+Dm*Pp)...
        /(Dp^2*Dm+Dm^2*Dp);
    
end

function [d3pds3]=ThrirdDerivativeVarStencil(Pi,Pp,Pm,Ppp,Pmm,Dp,Dm,Ep,Em)
    % Calculation of a second difference with variable stencil
    
    d3pds3=6*(((Pp*Dm^2-Pm*Dp^2)*(Ep*Em^2+Em*Ep^2))...
        -((Em^2*Ppp-Ep^2*Pmm)*(Dp*Dm^2+Dm*Dp^2))-...
        Pi*((Dm^2-Dp^2)*(Ep^2*Em+Ep*Em^2)-(Em^2-Ep^2)*(Dp^2*Dm+Dp*Dm^2)))...
        /...
        ((Dp^2*Dm^3+Dp^3*Dm^2)*(Ep^2*Em+Ep*Em^2)-(Dp*Dm^2+Dp^2*Dm)*(Ep^2*Em^3+Ep^3*Em^2));
    
end


%% Dynamic equations

function [bendforce]=SnaxelToForce_curvature(snaxel,snakposition,bendforce,forceparam)
    % need to pass parameter array
    
    % Unpack constants
    global maxDt
    Dt=maxDt;
    dampBase=forceparam.dampBase;
    dampSides=forceparam.dampSides;
    vectorMagAveraging=forceparam.vectorMagAveraging;
    % Extract structure data
    vt=[snaxel(:).v]';
    Fvec=vertcat(bendforce(:).d2pds2);
    Ft=[snaxel(:).acc]';
    dirSnak=vertcat(snakposition(:).vector);
    
    if vectorMagAveraging
        Fd=MagnitudeAveraging(Fvec,dirSnak,snakposition);
    else
        Fd=Fvec;
    end
    
    for ii=1:length(Fvec)
        Ftp1(ii)=dot(Fd(ii,:),dirSnak(ii,:));
    end
    simPos=find([bendforce(:).samepos]);
    for ii=simPos
        Ftp1(ii)=1/dot(Fvec(ii,:),dirSnak(ii,:));
    end
    Ftp1(isinf(Ftp1))=sign(Ftp1(isinf(Ftp1)))*1000;
    Ftp1=Ftp1';
    
    [implicitMatDamp]=BuildDampingMatrix(bendforce,snaxel,dampBase,dampSides);
    
    [vtp1,Ftp1]=ForceToVelocity(vt,Ftp1,Ft,Dt,implicitMatDamp);
    for ii=1:length(vtp1)
        bendforce(ii).vel=vtp1(ii);
        bendforce(ii).acc=Ftp1(ii);
    end
end

function FAveraged=MagnitudeAveraging(Fvec,dirSnak,snakposition)
    % Sets the mean signed magnitude to 0 for a group of vectors Fvec where
    % the sign is defined relative to dirSnak
    
    [m,n]=size(Fvec);
    FMag=zeros([1,m]);
    FMagSign=zeros([m,1]);
    FUnit=zeros([m n]);
    FAveraged=zeros([m n]);
    
    for ii=1:length(Fvec)
        FMag(ii)=sqrt(sum(Fvec(ii,:).^2)); % Get magnitude
    end
    flat=find(FMag==0);
    for ii=1:length(Fvec)
        FUnit(ii,1:2)=Fvec(ii,:)/sqrt(sum(Fvec(ii,:).^2));% Calculate unit vector
    end
    for ii=flat
        FUnit(ii,1:2)=snakposition(ii).normvector{1}; % Extract unit vector if [0 0]
    end
    for ii=1:length(Fvec)
        FMagSign(ii)=sign(dot(FUnit(ii,:),dirSnak(ii,:)))*FMag(ii); % % Signed magnitude
    end
    meanMag=mean(FMagSign);
    for ii=1:length(Fvec)
        FAveraged(ii,1:2)=Fvec(ii,:)... % Remove som eof th emagnitude along unit axis
            -sign(dot(FUnit(ii,:),dirSnak(ii,:)))*meanMag*FUnit(ii,:);
    end
end

function [vtp1,Ftp1]=ForceToVelocity(vt,Ftp1,Ft,Dt,dampMat)
    %rhonm=1 keeps energy bounded, and conserves in the linear case
    rhonm=1.0;
    gam=(3-rhonm)/(2*(1+rhonm));
    
    maxForce=1;
    forceScaling=maxForce/max([mean(abs(Ftp1)),maxForce]);
    Fd=dampMat*vt;
    %Ftp1=Ftp1*forceScaling+Fd;
    Ftp1(Ftp1>maxForce)=maxForce;
    Ftp1(Ftp1<-maxForce)=-maxForce;
    Ftp1=Ftp1+Fd;
    vtp1=vt+((1-gam)*Ft+gam*Ftp1)*Dt;
    
    %vtp1=(eye(length(Fd))-dampMat)\(Dt*Ft*forceScaling+vt);
end

function [implicitMatDamp]=BuildDampingMatrix(derivtenscalc,snaxel,dampBase,dampSides)
    numCoeff=length(derivtenscalc);
    
    implicitMatDamp=zeros(numCoeff);


    
    snaxPrecList=[snaxel(:).snaxprec];
    snaxNextList=[snaxel(:).snaxnext];
    snaxIndList=[derivtenscalc(:).index];
    snaxPrecSub=FindObjNum([],snaxPrecList,snaxIndList);
    
    snaxNextSub=FindObjNum([],snaxNextList,snaxIndList);
    
    
    for ii=1:length(derivtenscalc)
        implicitMatDamp(ii,ii)=-1*dampBase;
        implicitMatDamp(ii,snaxPrecSub(ii))=0.5*dampSides;
        implicitMatDamp(ii,snaxNextSub(ii))=0.5*dampSides;
    end
end

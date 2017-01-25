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
% %#codegen
function [snaxel,snakposition,snaxelmodvel,velcalcinfostruct,sensSnax,forceparam]...
        =VelocityLengthMinimisationSQP(snaxel,snakposition,volumefraction,coeffstructure,forceparam)
    
    [snaxeltensvel,snakposition,velcalcinfostruct,sensSnax,forceparam]=GeometryForcingVelocity(snaxel,snakposition,forceparam,coeffstructure,volumefraction);
    [snaxel]=AssignVelocityToSnaxel(snaxel,snaxeltensvel);
    
    snaxelmodvel=snaxel;
end

function [snaxel]=AssignVelocityToSnaxel(snaxel,snaxeltensvel)
    % Assigns velocities to the snaxel structure
    
    
    snaxTensInd=[snaxeltensvel(:).index];
    
    for ii=1:length(snaxel)
        velTensSub=FindObjNum([],snaxel(ii).index,snaxTensInd);
        
        snaxel(ii).v=snaxeltensvel(velTensSub).forcevel;
        
        
    end
    
    
end

%% Geometry feature forcing
% In all code Dt is removed and assumed to be 1
function [snaxeltensvel,snakposition,velcalcinfostruct,sensSnax,forceparam]=...
        GeometryForcingVelocity(snaxel,snakposition,forceparam,coeffstructure,volumefraction)
    % Calculate a tensile velocity to influence the solution
    snakPosIndex=[snakposition(:).index];
    snaxIndex=[snaxel(:).index];
    
    if sum(snakPosIndex~=snaxIndex)>0
        error('Indices do not match')
    end
    % Legacy
    [snaxeltensvel,snakposition]=CalculateTensileVelocity(snaxel,snakposition,snakPosIndex);
    [snaxeltensvel,snakposition]=CalculateBendingVelocity(snaxel,snakposition,snakPosIndex,snaxeltensvel);
    bendingVelInfluence=forceparam.bendingVelInfluence;
    tensVelInfluence=forceparam.tensVelInfluence;
    maxForceVel=forceparam.maxForceVel;
    forceVelScaling=maxForceVel/(bendingVelInfluence+tensVelInfluence);
    tensCoeff=forceVelScaling*tensVelInfluence;
    bendCoeff=forceVelScaling*bendingVelInfluence;
    
    % current Linear
    [areaTargVec,areaConstrMat]=AreaConstraintMatrixAssignement(snaxel,coeffstructure,volumefraction);
%     [derivtenscalc]=CalculateTensileVelocity2(snaxel,snakposition,snakPosIndex);
%     [implicitMatTens,forceVec]=BuildImplicitMatrix(derivtenscalc);
%     [forcingVec,conditionMat]=BuildSolutionLaplacianMatrix(implicitMatTens,forceVec,areaTargVec,areaConstrMat);
    
    % Current SQP
    smearLengthEps=forceparam.lengthEpsilon;
    distEpsilon=forceparam.distEpsilon;
    dirEpsilon=forceparam.dirEpsilon;
    typeSmear=forceparam.typeSmear;
    lagMultiPast=forceparam.lagMulti;
    switch typeSmear
        case 'length'
            [derivtenscalc2]=ExtractDataForDerivatives_LengthSmear(snaxel,snakposition,snakPosIndex,smearLengthEps);
        case 'lengthD'
            [derivtenscalc2]=ExtractDataForDerivatives_LengthDSmear(snaxel,snakposition,snakPosIndex,smearLengthEps);
        case 'd'
            [derivtenscalc2]=ExtractDataForDerivatives_distanceSmear(snaxel,snakposition,snakPosIndex,smearLengthEps,distEpsilon);
        case 'dir'
            [derivtenscalc2]=ExtractDataForDerivatives_directionSmear(snaxel,snakposition,snakPosIndex,smearLengthEps,distEpsilon,dirEpsilon);
    end
    
    derivtenscalc2=MatchCellToderivtenscal(derivtenscalc2,coeffstructure,volumefraction);
    if isempty(lagMultiPast)
        lagMultiPast=zeros([numel(volumefraction) 1]);
    end
    [Df,Hf,HA]=BuildJacobianAndHessian(derivtenscalc2,volumefraction,lagMultiPast,numel(volumefraction));
    isFreeze=[snaxel(:).isfreeze];
    %   [Deltax]=SQPStep(Df,Hf,areaConstrMat',areaTargVec);
    warning('OFF','MATLAB:nearlySingularMatrix')
    if true
        % THis section needs a clean up
        [DeltaxisFreeze,~]=SQPStepFreeze(Df,Hf,areaConstrMat',areaTargVec,false(size(isFreeze)));
        [isFreezeRnd2]=VelocityThawing(isFreeze,DeltaxisFreeze);
        [Deltax,lagMulti,condMat]=SQPStepFreeze(Df,Hf,areaConstrMat',areaTargVec,isFreezeRnd2);
        finIsFreeze=isFreezeRnd2;
        if any(isnan(DeltaxisFreeze))
            DeltaxisFreeze(:)=0;
        end
        if any(isnan(Deltax))
            [Deltax,lagMulti,condMat]=SQPStepFreeze(Df,Hf,areaConstrMat',areaTargVec,(isFreeze==1));
            finIsFreeze=(isFreeze==1);
        end
        HL=Hf;
    else
        [DeltaxisFreeze,~]=SQPStepLagFreeze(Df,Hf,HA,areaConstrMat',areaTargVec,false(size(isFreeze)));
        [isFreezeRnd2]=VelocityThawing(isFreeze,DeltaxisFreeze);
        [Deltax,lagMulti,condMat]=SQPStepLagFreeze(Df,Hf,HA,areaConstrMat',areaTargVec,isFreezeRnd2);
        finIsFreeze=isFreezeRnd2;
        if any(isnan(DeltaxisFreeze))
            DeltaxisFreeze(:)=0;
        end
        if any(isnan(Deltax))
            [Deltax,lagMulti,condMat]=SQPStepLagFreeze(Df,Hf,HA,areaConstrMat',areaTargVec,(isFreeze==1));
            finIsFreeze=(isFreeze==1);
        end
        HL=Hf+HA;
    end
    
    warning('ON','MATLAB:nearlySingularMatrix')
    Deltax(finIsFreeze)=DeltaxisFreeze(finIsFreeze);
    [feasVal,optVal,posDefVal]=SQPOptim(Df,HL,areaConstrMat',areaTargVec,lagMultiPast);
    
    forceparam.lagMulti=lagMulti;
    
    %lagMulti(isFreeze)=lagMultiisFreeze(isFreeze);
    %     [DeltaxFin,lagMulti]=SQPStepFreeze(Df,Hf,areaConstrMat',areaTargVec,isFreeze);
    %[hessA]=BuildDAdd2(snaxel,coeffstructure,volumefraction,lagMulti,derivtenscalc2);
    %     [Deltax,lagMulti]=SQPStepFreeze(Df,(Hf+hessA),areaConstrMat',areaTargVec,isFreeze);
    
    maxForceMean=forceparam.maxForceMean;
    maxForceStd=forceparam.maxForceStd;
    scaleOscil=forceparam.scaleOscil;
    
    if (condMat<1e-12 && (mean(abs(Deltax))*maxForceStd<std(abs(Deltax)))) %|| posDefVal~=0
        fprintf(' ! DeltaX is Oscillatory std/mean: %.2e ! ',std(abs(Deltax))/mean(abs(Deltax)))
        Deltax=Deltax*scaleOscil/max(abs(Deltax));
    elseif mean(abs(Deltax))>maxForceMean
        fprintf(' ! DeltaX is Large %.2e ! ',mean(abs(Deltax)))
        Deltax=Deltax*maxForceMean/mean(abs(Deltax));
    end
    
    sensSnax=[];
    if forceparam.isLast
        [sensSnax]=CalculateSensitivity(Hf,HA,areaConstrMat,lagMulti);
        
    end
%     velcalcinfostruct.forcingVec=forcingVec;
%     velcalcinfostruct.conditionMat=conditionMat;
    velcalcinfostruct=[];
    %[tensVelVec]=GeometryForcingVelCalc(forcingVec,conditionMat,tensCoeff);
    
    for ii=1:length(snaxeltensvel)
        snaxeltensvel(ii).forcevel=Deltax(ii);
    end
    
end

function [isFreezeRnd2]=VelocityThawing(isFreeze,DxisFreeze)
    % function handles the thawing of vertices 
    
    isThaw=(isFreeze>1 & isFreeze<2 & DxisFreeze>=0) ...
        | (isFreeze>2 & isFreeze<3 & DxisFreeze<=0);
    
    isFreezeRnd2= logical(isFreeze) & ~isThaw;
end

function [derivtenscal]=MatchCellToderivtenscal(derivtenscal,coeffstructure,volumefraction)
    coeffSnaxInd=[coeffstructure(:).snaxelindex];
    
    allNewCells=[volumefraction(:).newCellInd];
    allOldCells=zeros(size(allNewCells));
    allOldCellsSub=allOldCells;
    kk=1;
    for ii=1:numel(volumefraction)
        nCurr=numel(volumefraction(ii).newCellInd);
        allOldCells(kk:kk+nCurr-1)=volumefraction(ii).oldCellInd;
        allOldCellsSub(kk:kk+nCurr-1)=ii;
        kk=kk+nCurr;
    end
    err=false;
    for ii=1:numel(derivtenscal)
        newCell=unique([coeffstructure(FindObjNum([],[derivtenscal(ii).index],coeffSnaxInd)).cellindex]);
        newCell2=unique([coeffstructure(FindObjNum([],[derivtenscal(ii).snaxprec],coeffSnaxInd)).cellindex]);
        subs=FindObjNum([],newCell,allNewCells);
        subs2=FindObjNum([],newCell2,allNewCells);
        oldCell=unique(allOldCells(subs));
        oldCell2=unique(allOldCells(subs2));
        [i2,i1]=find((ones([numel(oldCell2) 1])*oldCell)==(oldCell2'*ones([1 numel(oldCell)])));
        err=err || isempty(i1);
        derivtenscal(ii).cellprec=unique(allOldCells(subs(i1)));
        derivtenscal(ii).cellprecsub=unique(allOldCellsSub(subs(i1)));
        
    end
    
    if err
        error('snakes:connectivity:nonSharedCell',...
            'Snaxels do not share a cell despite connection \n connectivity information damaged')
        
    end
    
end

function [forcingVec,conditionMat]=BuildSolutionLaplacianMatrix(implicitMatTens,forceVec,areaTargVec,areaConstrMat)
    
    forcingVec=[forceVec;areaTargVec];
    
    [mT,nT]=size(implicitMatTens);
    [mA,nA]=size(areaConstrMat);
    
    if mT~=nA
        error('Number of velocities does not match')
    end
    
    zeroPad=zeros(mA);
    
    conditionMat=[[implicitMatTens;areaConstrMat],[areaConstrMat';zeroPad]];
    %conditionMat=[[implicitMatTens;areaConstrMat]];
    testZero=find(sum(abs([forcingVec,conditionMat]),2)==0);
    forcingVec(testZero)=[];
    conditionMat(testZero,:)=[];
    conditionMat(:,testZero)=[];
end

function [areaTargVec,areaConstrMat]=AreaConstraintMatrixAssignement(snaxel,coeffstructure,volumefraction)
    
    snaxInd=[snaxel(:).index];
    cellindexCoeff=[coeffstructure(:).cellindex];
    snaxIndCoeff=[coeffstructure(:).snaxelindex];
    coeffCellInd=[coeffstructure(:).cellindex];
    cellIndex=RemoveIdenticalEntries(cellindexCoeff);
    
    snaxIndSubCoeff=FindObjNum([],snaxIndCoeff,snaxInd);
    coeffIndSubCoeff=zeros(size(cellIndex));
    
    
    
    areaTargVec=-([volumefraction(:).targetfill]-[volumefraction(:).volumefraction]).*[volumefraction(:).totalvolume];
    areaTargVec=areaTargVec';
    areaConstrMat=zeros([length(volumefraction),length(snaxInd)]);
    
    for ii=1:length(volumefraction)
        coeffSubs=FindObjNum([],volumefraction(ii).newCellInd,coeffCellInd);
        coeffSubs(coeffSubs==0)=[];
        coeffIndSubCoeff(coeffSubs)=ii;
    end
    
    for ii=1:length(coeffstructure)
        areaConstrMat(coeffIndSubCoeff(ii),snaxIndSubCoeff(ii))=...
            areaConstrMat(coeffIndSubCoeff(ii),snaxIndSubCoeff(ii))...
            +coeffstructure(ii).value;
        
    end
    totVol=[volumefraction(:).totalvolume]';
    areaConstrMat=areaConstrMat./repmat(totVol,[1,size(areaConstrMat,2)]);
    areaTargVec=areaTargVec./totVol;
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

% % WARNING BENDING IS TURNED OFF
function [implicitMatTens,forceVec]=BuildImplicitMatrix(derivtenscalc)
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
    forceVec=[derivtenscalc(:).forcecoeff_i]';
    %% WARNING BENDING IS TURNED OFF
    %     for ii=1:length(derivtenscalc)
    %         % Base
    %         implicitMatBend(ii,ii)=derivtenscalc(ii).velcoeff_i;
    %         implicitMatBend(ii,snaxPrecSub(ii))=derivtenscalc(ii).velcoeff_m;
    %         implicitMatBend(ii,snaxNextSub(ii))=derivtenscalc(ii).velcoeff_p;
    %         %p1
    %         implicitMatBendp1(ii,snaxNextSub(ii))=derivtenscalc(snaxNextSub(ii)).velcoeff_i;
    %         implicitMatBendp1(ii,snaxPrecSub(snaxNextSub(ii)))=derivtenscalc(snaxNextSub(ii)).velcoeff_m;
    %         implicitMatBendp1(ii,snaxNextSub(snaxNextSub(ii)))=derivtenscalc(snaxNextSub(ii)).velcoeff_p;
    %         %m1
    %         implicitMatBendm1(ii,snaxPrecSub(ii))=derivtenscalc(snaxPrecSub(ii)).velcoeff_i;
    %         implicitMatBendm1(ii,snaxPrecSub(snaxPrecSub(ii)))=derivtenscalc(snaxPrecSub(ii)).velcoeff_m;
    %         implicitMatBendm1(ii,snaxNextSub(snaxPrecSub(ii)))=derivtenscalc(snaxPrecSub(ii)).velcoeff_p;
    %     end
    implicitMatBend=-2*implicitMatBend+implicitMatBendm1+implicitMatBendp1;
end

function [derivtenscalc]=CalculateTensileVelocity2(snaxel,snakposition,snakPosIndex)
    
    
    %CheckResults(1,unstructglobal,unstructglobal,snakposition,snaxel,0);
    
    
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
        
        [forcecoeff,velcoeff_i,velcoeff_p,velcoeff_m]=TensileVelDerivCoeff(derivtenscalc(ii));
        
        derivtenscalc(ii).forcecoeff_i=forcecoeff;
        derivtenscalc(ii).velcoeff_i=velcoeff_i;
        derivtenscalc(ii).velcoeff_p=velcoeff_p;
        derivtenscalc(ii).velcoeff_m=velcoeff_m;
        
        %quiver(derivtenscalc(ii).pos_i(1),derivtenscalc(ii).pos_i(2),snakposition(ii).tensVector(1),snakposition(ii).tensVector(2));
        %         quiver(derivtenscalc(ii).pos_i(1),derivtenscalc(ii).pos_i(2),bendingVec(1),bendingVec(2));
    end
    
end

function [forcecoeff,velcoeff_i,velcoeff_p,velcoeff_m]=TensileVelDerivCoeff(derivtenscalc)
    
    %     fields=fieldnames(derivtenscalc);
    %
    %     for ii=1:length(fields)
    %         eval([fields{ii},'=derivtenscalc.',fields{ii},';']);
    %     end
    
    dir_i=derivtenscalc.dir_i;
    dir_p=derivtenscalc.dir_p;
    dir_m=derivtenscalc.dir_m;
    pos_i=derivtenscalc.pos_i;
    pos_p=derivtenscalc.pos_p;
    pos_m=derivtenscalc.pos_m;
    
    forcecoeff=2*dot(dir_i,(2*pos_i-(pos_p+pos_m)));
    
    velcoeff_i=4*dot(dir_i,dir_i);
    
    velcoeff_p=-2*dot(dir_i,dir_p);
    
    velcoeff_m=-2*dot(dir_i,dir_m);
    
    
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

function [snaxeltensvel,snakposition]=CalculateBendingVelocity(snaxel,snakposition,snakPosIndex,snaxeltensvel)
    
    
    
    velRatio=8; % limits the velocity to the same range as the
    for ii=1:length(snakposition)
        neighbourSub=FindObjNum([],[snaxel(ii).snaxprec,snaxel(ii).snaxnext],snakPosIndex);
        tensVecNeighbour=vertcat(snakposition(neighbourSub).tensVector);
        
        snakposition(ii).bendVector=-sum(tensVecNeighbour)+2*snakposition(ii).tensVector;
        snaxeltensvel(ii).bendvel=dot(snakposition(ii).bendVector,snakposition(ii).vector)/(velRatio);
    end
    
end

function [tensVelVec]=GeometryForcingVelCalc(forcingVec,conditionMat,forcecoeff)
    %
    %     T=[derivtenscalc(:).tension]';
    %     B=[derivtenscalc(:).bending]';
    %
    %     F=T*tensCoeff+B*bendCoeff;
    
    
    
    forcingVec=forcingVec*forcecoeff;
    condMat=rcond(conditionMat);
    if condMat>1e-12 && condMat<1e12
        tensVelVec=(conditionMat)\(-forcingVec);
    else
        disp('Warning pseudo inverse used when calculating velocities')
        tensVelVec=pinv(conditionMat)*(-forcingVec);
        %tensVelVec=(conditionMat'*conditionMat)\(conditionMat')*(-forcingVec);
    end
    %tensVelVec=tensVelVec/mean(abs(tensVelVec));
    
    
    
end

%% Actual Profile length minimisation using SQP

function [derivtenscalc]=ExtractDataForDerivatives_LengthSmear(snaxel,snakposition,snakPosIndex,smearLengthEps)
    
    derivtenscalcTemplate=struct('index',[],...
        'snaxprec',[],...
        'precsub',[],...,...
        'cellprec',[],...
        'cellprecsub',[],...
        'Dg_i',[],...
        'Dg_m',[],...
        'g1_i',[],...
        'g1_m',[],...
        'd_i',[],...
        'd_m',[],...
        'p_i',[],...
        'p_m',[],...
        'normFi',[],...
        'a_i',[],...
        'a_m',[],...
        'a_im',[],...
        'b_i',[],...
        'b_m',[],...
        'c',[],...
        'dfiddi',[],...
        'dfiddm',[],...
        'd2fiddi2',[],...
        'd2fiddm2',[],...
        'd2fiddim',[]...
        );
    derivtenscalc=repmat(derivtenscalcTemplate,[1,length(snakposition)]);
    
    for ii=length(snakposition):-1:1
        neighSub=FindObjNum([],[snaxel(ii).snaxprec],snakPosIndex);
        
        derivtenscalc(ii).index=snakposition(ii).index;
        derivtenscalc(ii).snaxprec=snaxel(ii).snaxprec;
        derivtenscalc(ii).precsub=neighSub;
        % extracting data from preexisting arrays
        derivtenscalc(ii).Dg_i=snakposition(ii).vectornotnorm;
        derivtenscalc(ii).Dg_m=snakposition(neighSub).vectornotnorm;
        derivtenscalc(ii).g1_i=snakposition(ii).vertInit;
        derivtenscalc(ii).g1_m=snakposition(neighSub).vertInit;
        derivtenscalc(ii).p_i=snakposition(ii).coord;
        derivtenscalc(ii).p_m=snakposition(neighSub).coord;
        derivtenscalc(ii).d_i=snaxel(ii).d;
        derivtenscalc(ii).d_m=snaxel(neighSub).d;
        % calculating data
        
        derivtenscalc(ii).normFi=sqrt(smearLengthEps^2+sum(...
            (derivtenscalc(ii).p_i- derivtenscalc(ii).p_m).^2));
        
    end
    
    
    for ii=length(snakposition):-1:1
        [derivtenscalc(ii).a_i,...
            derivtenscalc(ii).a_m,...
            derivtenscalc(ii).a_im,...
            derivtenscalc(ii).b_i,...
            derivtenscalc(ii).b_m,...
            derivtenscalc(ii).c]=...
            Calc_LengthDerivCoeff(...
            derivtenscalc(ii).Dg_i,derivtenscalc(ii).Dg_m,...
            derivtenscalc(ii).g1_i,derivtenscalc(ii).g1_m);
        
        [derivtenscalc(ii)]=CalculateDerivatives(derivtenscalc(ii));
        
        %         if (abs(dot(derivtenscalc(ii).Dg_i,derivtenscalc(ii).Dg_m))<=1e-10) ...
        %                 && (sum(derivtenscalc(ii).g1_i==derivtenscalc(ii).g1_m)==numel(derivtenscalc(ii).g1_i))
        %             printf('We''d have a problem')
        %         end
        
    end
    testnan=find(isnan([derivtenscalc(:).d2fiddi2]));
    if ~isempty(testnan)
        testnan
    end
end

function [derivtenscalc]=ExtractDataForDerivatives_LengthDSmear(snaxel,snakposition,snakPosIndex,smearLengthEps)
    
    normVec=@(vec) sqrt(sum(vec.^2,2));
    
    derivtenscalcTemplate=struct('index',[],...
        'snaxprec',[],...
        'precsub',[],...,...
        'cellprec',[],...
        'cellprecsub',[],...
        'Dg_i',[],...
        'Dg_m',[],...
        'g1_i',[],...
        'g1_m',[],...
        'd_i',[],...
        'd_m',[],...
        'p_i',[],...
        'p_m',[],...
        'normFi',[],...
        'a_i',[],...
        'a_m',[],...
        'a_im',[],...
        'b_i',[],...
        'b_m',[],...
        'c',[],...
        'dfiddi',[],...
        'dfiddm',[],...
        'd2fiddi2',[],...
        'd2fiddm2',[],...
        'd2fiddim',[]...
        );
    derivtenscalc=repmat(derivtenscalcTemplate,[1,length(snakposition)]);
    
    for ii=length(snakposition):-1:1
        neighSub=FindObjNum([],[snaxel(ii).snaxprec],snakPosIndex);
        
        derivtenscalc(ii).index=snakposition(ii).index;
        derivtenscalc(ii).snaxprec=snaxel(ii).snaxprec;
        derivtenscalc(ii).precsub=neighSub;
        % extracting data from preexisting arrays
        derivtenscalc(ii).Dg_i=snakposition(ii).vectornotnorm;
        derivtenscalc(ii).Dg_m=snakposition(neighSub).vectornotnorm;
        derivtenscalc(ii).g1_i=snakposition(ii).vertInit;
        derivtenscalc(ii).g1_m=snakposition(neighSub).vertInit;
        derivtenscalc(ii).p_i=snakposition(ii).coord;
        derivtenscalc(ii).p_m=snakposition(neighSub).coord;
        derivtenscalc(ii).d_i=snaxel(ii).d;
        derivtenscalc(ii).d_m=snaxel(neighSub).d;
        % calculating data
        
        derivtenscalc(ii).normFi=sqrt(smearLengthEps^2*min(normVec(derivtenscalc(ii).Dg_i),normVec(derivtenscalc(ii).Dg_m))^2+sum(...
            (derivtenscalc(ii).p_i- derivtenscalc(ii).p_m).^2));
        
    end
    
    
    for ii=length(snakposition):-1:1
        [derivtenscalc(ii).a_i,...
            derivtenscalc(ii).a_m,...
            derivtenscalc(ii).a_im,...
            derivtenscalc(ii).b_i,...
            derivtenscalc(ii).b_m,...
            derivtenscalc(ii).c]=...
            Calc_LengthDerivCoeff(...
            derivtenscalc(ii).Dg_i,derivtenscalc(ii).Dg_m,...
            derivtenscalc(ii).g1_i,derivtenscalc(ii).g1_m);
        
        [derivtenscalc(ii)]=CalculateDerivatives(derivtenscalc(ii));
        
        %         if (abs(dot(derivtenscalc(ii).Dg_i,derivtenscalc(ii).Dg_m))<=1e-10) ...
        %                 && (sum(derivtenscalc(ii).g1_i==derivtenscalc(ii).g1_m)==numel(derivtenscalc(ii).g1_i))
        %             printf('We''d have a problem')
        %         end
        
    end
    testnan=find(isnan([derivtenscalc(:).d2fiddi2]));
    if ~isempty(testnan)
        testnan
    end
end

function [derivtenscalc]=ExtractDataForDerivatives_distanceSmear(snaxel,snakposition,snakPosIndex,smearLengthEps,distEpsilon)
    
    derivtenscalcTemplate=struct('index',[],...
        'snaxprec',[],...
        'precsub',[],...,...
        'cellprec',[],...
        'cellprecsub',[],...
        'Dg_i',[],...
        'Dg_m',[],...
        'g1_i',[],...
        'g1_m',[],...
        'd_i',[],...
        'd_m',[],...
        'p_i',[],...
        'p_m',[],...
        'normFi',[],...
        'a_i',[],...
        'a_m',[],...
        'a_im',[],...
        'b_i',[],...
        'b_m',[],...
        'c',[],...
        'dfiddi',[],...
        'dfiddm',[],...
        'd2fiddi2',[],...
        'd2fiddm2',[],...
        'd2fiddim',[]...
        );
    derivtenscalc=repmat(derivtenscalcTemplate,[1,length(snakposition)]);
    
    for ii=length(snakposition):-1:1
        neighSub=FindObjNum([],[snaxel(ii).snaxprec],snakPosIndex);
        
        derivtenscalc(ii).index=snakposition(ii).index;
        derivtenscalc(ii).snaxprec=snaxel(ii).snaxprec;
        derivtenscalc(ii).precsub=neighSub;
        % extracting data from preexisting arrays
        derivtenscalc(ii).Dg_i=snakposition(ii).vectornotnorm;
        derivtenscalc(ii).Dg_m=snakposition(neighSub).vectornotnorm;
        derivtenscalc(ii).g1_i=snakposition(ii).vertInit;
        derivtenscalc(ii).g1_m=snakposition(neighSub).vertInit;
        
        derivtenscalc(ii).d_i=(1-2*distEpsilon)*snaxel(ii).d+distEpsilon;
        derivtenscalc(ii).d_m=(1-2*distEpsilon)*snaxel(neighSub).d+distEpsilon;
        % calculating data
        derivtenscalc(ii).p_i=(derivtenscalc(ii).g1_i+...
            derivtenscalc(ii).Dg_i*derivtenscalc(ii).d_i);
        derivtenscalc(ii).p_m=(derivtenscalc(ii).g1_m+...
            derivtenscalc(ii).Dg_m*derivtenscalc(ii).d_m);
        derivtenscalc(ii).normFi=sqrt(smearLengthEps^2+sum(...
            (derivtenscalc(ii).p_i- derivtenscalc(ii).p_m).^2));
        
        %         derivtenscalc(ii).normFi=sqrt(smearLengthEps^2+...
        %             sum(((derivtenscalc(ii).g1_i+derivtenscalc(ii).Dg_i*...
        %             derivtenscalc(ii).d_i)-(derivtenscalc(ii).g1_m+...
        %             derivtenscalc(ii).Dg_m*derivtenscalc(ii).d_m)).^2));
    end
    
    for ii=length(snakposition):-1:1
        [derivtenscalc(ii).a_i,...
            derivtenscalc(ii).a_m,...
            derivtenscalc(ii).a_im,...
            derivtenscalc(ii).b_i,...
            derivtenscalc(ii).b_m,...
            derivtenscalc(ii).c]=...
            Calc_LengthDerivCoeff(...
            derivtenscalc(ii).Dg_i,derivtenscalc(ii).Dg_m,...
            derivtenscalc(ii).g1_i,derivtenscalc(ii).g1_m);
        
        [derivtenscalc(ii)]=CalculateDerivatives_d(derivtenscalc(ii),smearLengthEps);
        
    end
    testnan=find(isnan([derivtenscalc(:).d2fiddi2]));
    if ~isempty(testnan)
        testnan
    end
end

function [derivtenscalc]=ExtractDataForDerivatives_directionSmear...
        (snaxel,snakposition,snakPosIndex,smearLengthEps,distEpsilon,smearLengthDir)
    
    derivtenscalcTemplate=struct('index',[],...
        'snaxprec',[],...
        'precsub',[],...,...
        'cellprec',[],...
        'cellprecsub',[],...
        'Dg_i',[],...
        'Dg_m',[],...
        'g1_i',[],...
        'g1_m',[],...
        'd_i',[],...
        'd_m',[],...
        'p_i',[],...
        'p_m',[],...
        'normFi',[],...
        'a_i',[],...
        'a_m',[],...
        'a_im',[],...
        'b_i',[],...
        'b_m',[],...
        'c',[],...
        'dfiddi',[],...
        'dfiddm',[],...
        'd2fiddi2',[],...
        'd2fiddm2',[],...
        'd2fiddim',[]...
        );
    derivtenscalc=repmat(derivtenscalcTemplate,[1,length(snakposition)]);
    [snakposition]=ModifySnakposition(snaxel,snakposition,snakPosIndex,smearLengthDir);
    for ii=length(snakposition):-1:1
        neighSub=FindObjNum([],[snaxel(ii).snaxprec],snakPosIndex);
        
        derivtenscalc(ii).index=snakposition(ii).index;
        derivtenscalc(ii).snaxprec=snaxel(ii).snaxprec;
        derivtenscalc(ii).precsub=neighSub;
        % extracting data from preexisting arrays
        derivtenscalc(ii).Dg_i=snakposition(ii).vectornotnorm;
        derivtenscalc(ii).Dg_m=snakposition(neighSub).vectornotnorm;
        derivtenscalc(ii).g1_i=snakposition(ii).vertInit;
        derivtenscalc(ii).g1_m=snakposition(neighSub).vertInit;
        
        derivtenscalc(ii).d_i=(1-2*distEpsilon)*snaxel(ii).d+distEpsilon;
        derivtenscalc(ii).d_m=(1-2*distEpsilon)*snaxel(neighSub).d+distEpsilon;
        % calculating data
        derivtenscalc(ii).p_i=(derivtenscalc(ii).g1_i+...
            derivtenscalc(ii).Dg_i*derivtenscalc(ii).d_i);
        derivtenscalc(ii).p_m=(derivtenscalc(ii).g1_m+...
            derivtenscalc(ii).Dg_m*derivtenscalc(ii).d_m);
        derivtenscalc(ii).normFi=sqrt(smearLengthEps^2+sum(...
            (derivtenscalc(ii).p_i- derivtenscalc(ii).p_m).^2));
        
    end
    
    for ii=length(snakposition):-1:1
        [derivtenscalc(ii).a_i,...
            derivtenscalc(ii).a_m,...
            derivtenscalc(ii).a_im,...
            derivtenscalc(ii).b_i,...
            derivtenscalc(ii).b_m,...
            derivtenscalc(ii).c]=...
            Calc_LengthDerivCoeff(...
            derivtenscalc(ii).Dg_i,derivtenscalc(ii).Dg_m,...
            derivtenscalc(ii).g1_i,derivtenscalc(ii).g1_m);
        
        [derivtenscalc(ii)]=CalculateDerivatives_d(derivtenscalc(ii),smearLengthEps);
        
    end
    testnan=find(isnan([derivtenscalc(:).d2fiddi2]));
    if ~isempty(testnan)
        testnan
    end
end

function [snakposition]=ModifySnakposition(snaxel,snakposition,snakPosIndex,eps)
    
    for ii=1:length(snakposition)
        
        neighSub=FindObjNum([],[snaxel(ii).snaxprec,snaxel(ii).snaxnext],snakPosIndex)';
        
        dists=[snaxel([ii,neighSub]).d]';
        vector=vertcat(snakposition([ii,neighSub]).vector);
        
        
        edgeDist=sqrt(sum((vertcat(snakposition(neighSub).coord)-...
            ones([2,1])*snakposition(ii).coord).^2,2));
        edgeNorm=sqrt(sum((vertcat(snakposition(neighSub).vectornotnorm)+...
            ones([2,1])*snakposition(ii).vectornotnorm).^2,2));
        edgeDist=edgeDist./edgeNorm;
        % find intersection vector
        vertInit=vertcat(snakposition([ii,neighSub]).vertInit);
        vertEnd=vertcat(snakposition([ii,neighSub]).vertInit)...
            +vertcat(snakposition([ii,neighSub]).vectornotnorm);
        isVertInit=all([vertInit;vertEnd]==(ones([6,1])*vertInit(1,:)),2);
        isVertEnd=all([vertInit;vertEnd]==(ones([6,1])*vertEnd(1,:)),2);
        
        initMultipliers=isVertInit(1:3)-isVertInit(4:6);
        endMultipliers=isVertEnd(4:6)-isVertEnd(1:3);
        distCoeff=(eps-edgeDist)/eps;distCoeff(distCoeff<0)=0;distCoeff=[1;distCoeff];
        
        vecCoeff=distCoeff.*initMultipliers+distCoeff.*endMultipliers;
        vecCoeff(1)=1;
        
        newVec=(vector'*vecCoeff)'/sum(vecCoeff);
        newVecNoNorm=newVec*sqrt(sum(snakposition(ii).vectornotnorm.^2,2));
        snakposition(ii).vector=newVec;
        snakposition(ii).vectornotnorm=newVecNoNorm;
        %edgeDist
        %newVec-snakposition(ii).vector
    end
    
    % pause(0.5)
end

function [a_i,a_m,a_im,b_i,b_m,c]=Calc_LengthDerivCoeff(Dgi,Dgm,g1i,g1m)
    
    a_i=sum(Dgi.^2);
    a_m=sum(Dgm.^2);
    a_im=-2*dot(Dgi,Dgm);
    b_i=2*dot(Dgi,(g1i-g1m));
    b_m=-2*dot(Dgm,(g1i-g1m));
    c=sum((g1i-g1m).^2);
    
    
    
end

function [Df,Hf,HA]=BuildJacobianAndHessian(derivtenscalc,volumefraction,lagMulti,nCell)
    rot90=@(vec) ([0,1;-1,0]*(vec'))';
    n=length(derivtenscalc);
    Jacobian=zeros([n n]);
    Hessian=zeros([n n n]);
    HessianConst=zeros([n n nCell]);
    for ii=1:n
        neighSub=derivtenscalc(ii).precsub;
        Jacobian(ii,ii)=derivtenscalc(ii).dfiddi;
        Jacobian(neighSub,ii)=derivtenscalc(ii).dfiddm;
        
        Hessian(ii,ii,ii)=derivtenscalc(ii).d2fiddi2;
        Hessian(ii,neighSub,ii)=derivtenscalc(ii).d2fiddim;
        Hessian(neighSub,ii,ii)=derivtenscalc(ii).d2fiddim;
        Hessian(neighSub,neighSub,ii)=derivtenscalc(ii).d2fiddm2;
        
        HessConstr=0.5*(dot(derivtenscalc(ii).Dg_m,rot90(derivtenscalc(ii).Dg_i))...
            -dot(derivtenscalc(ii).Dg_i,rot90(derivtenscalc(ii).Dg_m)))...
            /volumefraction(derivtenscalc(ii).cellprecsub).totalvolume;
        HessianConst(ii,neighSub,derivtenscalc(ii).cellprecsub)=HessConstr;
        HessianConst(neighSub,ii,derivtenscalc(ii).cellprecsub)=HessConstr;
    end
    HA=sum(HessConstr.*repmat(reshape(lagMulti,[1,1,numel(lagMulti)]),[n,n]),3);
    Df=sum(Jacobian,2);
    Hf=sum(Hessian,3);
    
    %     if any(any(HA~=0))
    %         fprintf(' Hessian Not 0!! - ')
    %     end
    
end

function [Deltax]=SQPStep(Df,Hf,Dh,h_vec)
    
    rmvCol=find(sum(Dh~=0)==0);
    Dh(:,rmvCol)=[];
    h_vec(rmvCol)=[];
    
    Bkinv=(Hf)^(-1);
    matToInv=(Dh'*Bkinv*Dh);
    if rcond(matToInv)>1e-10
        u_kp1=matToInv\(h_vec-Dh'*Bkinv*Df);
    else
        u_kp1=pinv(matToInv)*(h_vec-Dh'*Bkinv*Df);
    end
    Deltax=-Bkinv*(Df+Dh*u_kp1);
    
end

function [optVal,feasVal,posDefErr]=SQPOptim(Df,HL,Dh,h_vec,lagMulti)
    rmvCol=find(sum(Dh~=0)==0);
    %     Dh(:,rmvCol)=[];
    %     h_vec(rmvCol)=[];
    %     lagMulti(rmvCol)=[];
    RMS=@(x) sqrt(mean(x.^2));
    
    feasVal=RMS(h_vec);
    optVal=RMS(Df+Dh*lagMulti);
    Z=null(Dh');
    %HL=HA+Hf;
    optTest=real(eig(Z'*HL*Z));
    posDefErr=RMS([optTest(optTest<0);zeros(size(optTest(optTest>=0)))]);
    fprintf(' feas:%.2e  opt:%.2e  posdef:%.2e - ',feasVal,optVal,posDefErr)
end

function [DeltaxFin,lagMulti,condMat]=SQPStepFreeze(Df,Hf,Dh,h_vec,isFreeze)
    
    rmvCol=find(isFreeze);
    Dh(rmvCol,:)=[];
    Hf(rmvCol,:)=[];
    Hf(:,rmvCol)=[];
    Df(rmvCol)=[];
    %h_vec(rmvCol)=[];
    % output vector for lagrange multiplier
    lagMulti=zeros([length(Dh(1,:)),1]);
    actCol=(sum(Dh~=0)~=0);
    
    rmvCol=find(sum(Dh~=0)==0);
    Dh(:,rmvCol)=[];
    h_vec(rmvCol)=[];
    
    
    
    Bkinv=(Hf)^(-1);
    matToInv=(Dh'*Bkinv*Dh);
    condMat=rcond(matToInv);
    fprintf(' cond: %.3e -', condMat)
    if true %rcond(matToInv)>1e-20
        u_kp1=matToInv\(h_vec-Dh'*Bkinv*Df);
    else
        u_kp1=pinv(matToInv)*(h_vec-Dh'*Bkinv*Df);
    end
    Deltax=-Bkinv*(Df+Dh*u_kp1);
    
    DeltaxFin=zeros(size(isFreeze));
    DeltaxFin(~isFreeze)=Deltax;
    lagMulti(actCol)=u_kp1;
end

function [DeltaxFin,lagMulti,condMat]=SQPStepLagFreeze(Df,Hf,HA,Dh,h_vec,isFreeze)
    
    rmvCol=find(isFreeze);
    
    Dh(rmvCol,:)=[];
    Hf(rmvCol,:)=[];
    Hf(:,rmvCol)=[];
    HA(rmvCol,:)=[];
    HA(:,rmvCol)=[];
    Df(rmvCol)=[];
    %h_vec(rmvCol)=[];
    % output vector for lagrange multiplier
    lagMulti=zeros([length(Dh(1,:)),1]);
    
    HL=Hf+HA;
    %DL=Df+Dh*pastLagMulti;
    
    actCol=(sum(Dh~=0)~=0);
    rmvCol=find(sum(Dh~=0)==0);
    Dh(:,rmvCol)=[];
    h_vec(rmvCol)=[];
    
    
    
    Bkinv=(HL)^(-1);
    matToInv=(Dh'*Bkinv*Dh);
    
    condMat=rcond(matToInv);
    fprintf(' cond: %.3e -', condMat)
    if true %rcond(matToInv)>1e-20
        u_kp1=matToInv\(h_vec-Dh'*Bkinv*Df);
    else
        u_kp1=pinv(matToInv)*(h_vec-Dh'*Bkinv*Df);
    end
    Deltax=-Bkinv*(Df+Dh*u_kp1);
    
    DeltaxFin=zeros(size(isFreeze));
    DeltaxFin(~isFreeze)=Deltax;
    lagMulti(actCol)=u_kp1;
end

function [DeltaxFin,lagMulti]=SQPStepFreeze_quadprog(Df,Hf,Dh,h_vec,isFreeze)
    
    rmvCol=find(isFreeze);
    Dh(rmvCol,:)=[];
    Hf(rmvCol,:)=[];
    Hf(:,rmvCol)=[];
    Df(rmvCol)=[];
    %h_vec(rmvCol)=[];
    % output vector for lagrange multiplier
    lagMulti=zeros([length(Dh(1,:)),1]);
    actCol=(sum(Dh~=0)~=0);
    
    rmvCol=find(sum(Dh~=0)==0);
    Dh(:,rmvCol)=[];
    h_vec(rmvCol)=[];
    
    [Deltax,~,~,~,lambda] = eval('quadprog(Hf,Df,[],[],Dh'',-h_vec);');
    u_kp1=lambda.eqlin;
    %     Bkinv=(Hf)^(-1);
    %     matToInv=(Dh'*Bkinv*Dh);
    %     if rcond(matToInv)>1e-10
    %         u_kp1=matToInv\(h_vec-Dh'*Bkinv*Df);
    %     else
    %         u_kp1=pinv(matToInv)*(h_vec-Dh'*Bkinv*Df);
    %     end
    %     Deltax=-Bkinv*(Df+Dh*u_kp1);
    
    DeltaxFin=zeros(size(isFreeze));
    DeltaxFin(~isFreeze)=Deltax;
    lagMulti(actCol)=u_kp1;
end

%% Derivative calculations - Length Smearing

function [derivtenscalcII]=CalculateDerivatives(derivtenscalcII)
    
    %     varExtract={'a_i','a_m','a_im','b_i','b_m','c','normFi','d_i','d_m'};
    %
    %     for ii=1:length(varExtract)
    %         eval([varExtract{ii},'=derivtenscalcII.(varExtract{ii});'])
    %     end
    
    a_i=derivtenscalcII.a_i;
    a_m=derivtenscalcII.a_m;
    a_im=derivtenscalcII.a_im;
    b_i=derivtenscalcII.b_i;
    b_m=derivtenscalcII.b_m;
    c=derivtenscalcII.c;
    normFi=derivtenscalcII.normFi;
    d_i=derivtenscalcII.d_i;
    d_m=derivtenscalcII.d_m;
    
    [derivtenscalcII.dfiddi]=Calc_DFiDdi(a_i,a_m,a_im,b_i,b_m,c,normFi,d_i,d_m);
    [derivtenscalcII.dfiddm]=Calc_DFiDdm(a_i,a_m,a_im,b_i,b_m,c,normFi,d_i,d_m);
    [derivtenscalcII.d2fiddi2]=Calc_D2FiDdi2(a_i,a_m,a_im,b_i,b_m,c,normFi,d_i,d_m);
    [derivtenscalcII.d2fiddm2]=Calc_DFiDdm2(a_i,a_m,a_im,b_i,b_m,c,normFi,d_i,d_m);
    [derivtenscalcII.d2fiddim]=Calc_D2FiDdim(a_i,a_m,a_im,b_i,b_m,c,normFi,d_i,d_m);
    
end

function [dfiddi]=Calc_DFiDdi(a_i,a_m,a_im,b_i,b_m,c,normFi,di,dm)
    
    
    dfiddi=(2*a_i*di+a_im*dm+b_i)/(2*(normFi));
    
    
end
function [dfiddm]=Calc_DFiDdm(a_i,a_m,a_im,b_i,b_m,c,normFi,di,dm)
    
    dfiddm=(2*a_m*dm+a_im*di+b_m)/(2*(normFi));
end
function [d2fiddi2]=Calc_D2FiDdi2(a_i,a_m,a_im,b_i,b_m,c,normFi,di,dm)
    
    
    d2fiddi2=((2*a_i*2*(normFi)^2)...
        -(2*a_i*di+a_im*dm+b_i)^2) ...
        /(4*((normFi)^3));
    
end
function [d2fiddm2]=Calc_DFiDdm2(a_i,a_m,a_im,b_i,b_m,c,normFi,di,dm)
    
    d2fiddm2=((2*a_m*2*(normFi)^2)...
        -(2*a_m*dm+a_im*di+b_m)^2) ...
        /(4*(normFi)^3);
end
function [d2fiddim]=Calc_D2FiDdim(a_i,a_m,a_im,b_i,b_m,c,normFi,di,dm)
    
    
    d2fiddim=((a_im*2*(normFi)^2)...
        -((2*a_m*dm+a_im*di+b_m)*(2*a_i*di+a_im*dm+b_i))) ...
        /(4*(normFi)^3);
    
    
end

%% Derivative calculations - Distance Smearing

function [derivtenscalcII]=CalculateDerivatives_d(derivtenscalcII,lSmear)
    
    %     varExtract={'a_i','a_m','a_im','b_i','b_m','c','normFi','d_i','d_m'};
    %
    %     for ii=1:length(varExtract)
    %         eval([varExtract{ii},'=derivtenscalcII.(varExtract{ii});'])
    %     end
    
    a_i=derivtenscalcII.a_i;
    a_m=derivtenscalcII.a_m;
    a_im=derivtenscalcII.a_im;
    b_i=derivtenscalcII.b_i;
    b_m=derivtenscalcII.b_m;
    c=derivtenscalcII.c;
    normFi=derivtenscalcII.normFi;
    d_i=derivtenscalcII.d_i;
    d_m=derivtenscalcII.d_m;
    
    [derivtenscalcII.dfiddi]=Calc_DFiDdi_d(a_i,a_m,a_im,b_i,b_m,c,normFi,d_i,d_m,lSmear);
    [derivtenscalcII.dfiddm]=Calc_DFiDdm_d(a_i,a_m,a_im,b_i,b_m,c,normFi,d_i,d_m,lSmear);
    [derivtenscalcII.d2fiddi2]=Calc_D2FiDdi2_d(a_i,a_m,a_im,b_i,b_m,c,normFi,d_i,d_m,lSmear);
    [derivtenscalcII.d2fiddm2]=Calc_DFiDdm2_d(a_i,a_m,a_im,b_i,b_m,c,normFi,d_i,d_m,lSmear);
    [derivtenscalcII.d2fiddim]=Calc_D2FiDdim_d(a_i,a_m,a_im,b_i,b_m,c,normFi,d_i,d_m,lSmear);
    
end
function [dfiddi]=Calc_DFiDdi_d(a_i,a_m,a_im,b_i,b_m,c,normFi,di,dm,lSmear)
    
    
    dfiddi=(2*a_i*di*(1-2*lSmear)+a_im*dm*(1-2*lSmear)+b_i*(1-2*lSmear))/(2*(normFi));
    
    
end
function [dfiddm]=Calc_DFiDdm_d(a_i,a_m,a_im,b_i,b_m,c,normFi,di,dm,lSmear)
    
    dfiddm=(2*a_m*dm*(1-2*lSmear)+a_im*di*(1-2*lSmear)+b_m*(1-2*lSmear))/(2*(normFi));
end
function [d2fiddi2]=Calc_D2FiDdi2_d(a_i,a_m,a_im,b_i,b_m,c,normFi,di,dm,lSmear)
    
    
    d2fiddi2=((2*a_i*(1-2*lSmear)*(1-2*lSmear)*2*(normFi)^2)...
        -(2*a_i*di*(1-2*lSmear)+a_im*dm*(1-2*lSmear)+b_i*(1-2*lSmear))^2) ...
        /(4*((normFi)^3));
    
end
function [d2fiddm2]=Calc_DFiDdm2_d(a_i,a_m,a_im,b_i,b_m,c,normFi,di,dm,lSmear)
    
    d2fiddm2=((2*a_m*(1-2*lSmear)*(1-2*lSmear)*2*(normFi)^2)...
        -(2*a_m*dm*(1-2*lSmear)+a_im*di*(1-2*lSmear)+b_m*(1-2*lSmear))^2) ...
        /(4*(normFi)^3);
end
function [d2fiddim]=Calc_D2FiDdim_d(a_i,a_m,a_im,b_i,b_m,c,normFi,di,dm,lSmear)
    
    
    d2fiddim=((a_im*(1-2*lSmear)*(1-2*lSmear)*2*(normFi)^2)...
        -((2*a_m*dm*(1-2*lSmear)+a_im*di*(1-2*lSmear)+b_m*(1-2*lSmear))*(2*a_i*di*(1-2*lSmear)+a_im*dm*(1-2*lSmear)+b_i*(1-2*lSmear)))) ...
        /(4*(normFi)^3);
    
    
end


%% Calculate position derivatives versus constraints

function [sensSnax]=SnaxelSensitivity(snaxel,coeffstructure,volumefraction,lagMultiplier,derivtenscalc,...
        Hf,areaConstrMat)
    %lagMultiplier=abs(lagMultiplier);
    [hessA]=BuildDAdd2(snaxel,coeffstructure,volumefraction,lagMultiplier,derivtenscalc);
    [sensSnax]=CalculateSensitivity(Hf,hessA,areaConstrMat,lagMultiplier);
    
end

function [hessA]=BuildDAdd2(snaxel,coeffstructure,volumefraction,lagMultiplier,derivtenscalc)
    
    [snaxtocell]=MatchSnaxtoCell(snaxel,coeffstructure,volumefraction);
    
    oldIndCell=[volumefraction(:).oldCellInd];
    n=length(snaxel);
    hessA=zeros(n);
    
    rotMat=[0 1;-1 0];
    
    for ii=1:n
        
        precSub=derivtenscalc(ii).precsub;
        cellPrec=ones(size(snaxtocell(ii).cellindex))'*snaxtocell(precSub).cellindex;
        cellCurr=snaxtocell(ii).cellindex'*ones(size(snaxtocell(precSub).cellindex));
        snaxCell=cellPrec(cellPrec==cellCurr);
        
        snaxCellSub=FindObjNum([],snaxCell,oldIndCell);
        lagMultiSnax=sum(lagMultiplier(snaxCellSub));
        
        Dg_i=derivtenscalc(ii).Dg_i;
        Dg_m=derivtenscalc(ii).Dg_m;
        
        ddAddSnax=lagMultiSnax*0.5*(dot(Dg_m,Dg_i)-dot(Dg_i,rotMat*Dg_m'));
        
        hessA(ii,precSub)=ddAddSnax;
        hessA(precSub,ii)=ddAddSnax;
        
    end
    
end

function [snaxtocell]=MatchSnaxtoCell(snaxel,coeffstructure,volumefraction)
    
    snaxInd=[snaxel(:).index];
    snaxIndCoeff=[coeffstructure(:).snaxelindex];
    newCellIndCoeff=[coeffstructure(:).cellindex];
    newCellIndVF=[volumefraction(:).newCellInd];
    
    oldCellIndVF=[volumefraction(:).oldCellInd];
    oldCellIndVFDistrib=zeros(size(newCellIndVF));
    kk=1;
    for ii=1:length(volumefraction)
        nNewCell=length(volumefraction(ii).newCellInd);
        oldCellIndVFDistrib(kk:kk+nNewCell-1)=oldCellIndVF(ii);
        kk=kk+nNewCell;
    end
    newCellSubCoeff=FindObjNum([],newCellIndCoeff,newCellIndVF);
    snaxSubCoeff=FindObjNum([],snaxIndCoeff,snaxInd);
    
    oldCellIndCoeff=oldCellIndVFDistrib(newCellSubCoeff);
    
    snaxtocell=struct('snaxindex',0,'cellindex',zeros([1,4]));
    snaxtocell=repmat(snaxtocell,[1, length(snaxel)]);
    
    for ii=1:length(snaxSubCoeff)
        
        snaxtocell(snaxSubCoeff(ii)).snaxindex=snaxIndCoeff(ii);
        snaxtocell(snaxSubCoeff(ii)).cellindex...
            (find(snaxtocell(snaxSubCoeff(ii)).cellindex==0,1))...
            =oldCellIndCoeff(ii);
        
    end
    for ii=1:length(snaxtocell)
        snaxtocell(ii).cellindex=RemoveIdenticalEntries(snaxtocell(ii).cellindex);
    end
    
    
end

function [sensSnax,sensLagMulti]=CalculateSensitivity(Hf,Ha,Ja_x,lagMulti)
    
    nSnax=length(Hf(1,:));
    nCond=length(Ja_x(:,1));
    HL=Hf+Ha;
    
    actCol=find(sum(abs(Ja_x),2)~=0);
    
    [Ja_xSmall,Ja_p,nCondAct]=TrimInactiveConditions(Ja_x,lagMulti,actCol);
    
    matToInv=[HL,Ja_xSmall';Ja_xSmall,zeros(nCondAct)];
    if rcond(matToInv)<1e-15
        actCol=find(abs(lagMulti)>1e-15);
        [Ja_xSmall,Ja_p,nCondAct]=TrimInactiveConditions(Ja_x,lagMulti,actCol);
        matToInv=[HL,Ja_xSmall';Ja_xSmall,zeros(nCondAct)];
    end
    matMultiplier=[zeros([nSnax,nCondAct]);Ja_p];
    
    resSens=-matToInv\matMultiplier;
    
    sensSnax=zeros([nSnax,nCond]);
    sensLagMulti=zeros([nCond,nCond]);
    
    sensSnax(:,actCol)=resSens(1:nSnax,:);
    signLagMulti=sign(lagMulti);
    signLagMulti(signLagMulti==0)=1;
    for ii=find(lagMulti)'
        
        sensSnax(:,ii)=sensSnax(:,ii)*signLagMulti(ii);
    end
    sensLagMulti(actCol,actCol)=resSens(nSnax+1:end,:);
    
    
end

function [Ja_x,Ja_p,nCondAct]=TrimInactiveConditions(Ja_x,lagMulti,actCol)
    Ja_x=Ja_x(actCol,:);
    lagMultiSmall=lagMulti(actCol);
    
    nCondAct=length(lagMultiSmall);
    
    Ja_p=eye(nCondAct);
    
    for ii=1:nCondAct
        Ja_p(ii,ii)=-(lagMultiSmall(ii));
    end
end




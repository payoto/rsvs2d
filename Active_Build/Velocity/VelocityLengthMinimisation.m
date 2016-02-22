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

function [snaxel,snakposition,snaxelmodvel,velcalcinfostruct]=VelocityLengthMinimisation(snaxel,snakposition,volumefraction,coeffstructure,forceparam)
    
    [snaxeltensvel,snakposition,velcalcinfostruct]=GeometryForcingVelocity(snaxel,snakposition,forceparam,coeffstructure,volumefraction);
    [snaxel]=AssignVelocityToSnaxel(snaxel,snaxeltensvel);
    
    snaxelmodvel=snaxel;
end

function [snaxel]=AssignVelocityToSnaxel(snaxel,snaxeltensvel)
    % Assigns velocities to the snaxel structure
    

    snaxTensInd=[snaxeltensvel(:).index];
    
    for ii=1:length(snaxel)
        velTensSub=FindObjNum([],snaxel(ii).index,snaxTensInd);
%         snaxelVel=[velstruct(velSub).averagevel]+[velstruct(velSub).deviationvel];
%         if numel(snaxelVel)>1
%             if abs(snaxelVel(1)-snaxelVel(2))>10^-12
%                 disp(['Difference in Vel is ',num2str(snaxelVel(1)-snaxelVel(2))])
%                 warning('Averaging difference in velocities')
%             end
%         end
%         snaxelvel(ii).index=snaxel(ii).index;
%         snaxelvel(ii).average=velstruct(velSub(1)).averagevel;
%         snaxelvel(ii).deviation=velstruct(velSub(1)).deviationvel;
%         snaxelvel(ii).force=snaxeltensvel(velTensSub).forcevel;
        %snaxel(ii).v=mean(snaxelVel)+snaxeltensvel(velTensSub).forcevel;
        snaxel(ii).v=snaxeltensvel(velTensSub).forcevel;
        
        
    end
    
    
end

%% Geometry feature forcing
% In all code Dt is removed and assumed to be 1
function [snaxeltensvel,snakposition,velcalcinfostruct]=GeometryForcingVelocity(snaxel,snakposition,forceparam,coeffstructure,volumefraction)
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
    
    % current
    [areaTargVec,areaConstrMat]=AreaConstraintMatrixAssignement(snaxel,coeffstructure,volumefraction);
    [derivtenscalc]=CalculateTensileVelocity2(snaxel,snakposition,snakPosIndex);
    [implicitMatTens,forceVec]=BuildImplicitMatrix(derivtenscalc);
    [forcingVec,conditionMat]=BuildSolutionLaplacianMatrix(implicitMatTens,forceVec,areaTargVec,areaConstrMat);
    
    
    velcalcinfostruct.forcingVec=forcingVec;
    velcalcinfostruct.conditionMat=conditionMat;
    
    [tensVelVec]=GeometryForcingVelCalc(forcingVec,conditionMat,tensCoeff);
    
    for ii=1:length(snaxeltensvel)
        snaxeltensvel(ii).forcevel=tensVelVec(ii);
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
        hold on
        %quiver(derivtenscalc(ii).pos_i(1),derivtenscalc(ii).pos_i(2),snakposition(ii).tensVector(1),snakposition(ii).tensVector(2));
%         quiver(derivtenscalc(ii).pos_i(1),derivtenscalc(ii).pos_i(2),bendingVec(1),bendingVec(2));
    end
    
end

function [forcecoeff,velcoeff_i,velcoeff_p,velcoeff_m]=TensileVelDerivCoeff(derivtenscalc)
    
    fields=fieldnames(derivtenscalc);
    
    for ii=1:length(fields)
        eval([fields{ii},'=derivtenscalc.',fields{ii},';']);
    end
    
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

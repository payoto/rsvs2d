%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%              University of Bristol
%      Department of Aerospace Engineering
%                     2015
%
%          Subdivision for Shape refinement
%             Alexandre Payot
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [newPoints,projPoints]=SubDivision(startPoints,nSteps,refineMethod)
    include_Utilities
    startPoints=RemoveIdenticalConsecutivePoints(startPoints);
    switch refineMethod
        case 'chaikin'
            [newPoints,projPoints]=SubSurfChainkin(startPoints,nSteps);
        case 'bspline'
            
            [newPoints,projPoints]=SubSurfBSpline(startPoints,nSteps);
        case 'test'
            figure
            hold on
            plot(startPoints(:,1),startPoints(:,2),'b--')
            [newChaik]=SubSurfChainkin(startPoints,nSteps);
            plot(newChaik(:,1),newChaik(:,2));
            [newSpline]=SubSurfBSpline(startPoints,nSteps);
            plot(newSpline(:,1),newSpline(:,2))
            [newInterp1]=SubSurfinterp1(startPoints,nSteps);
            plot(newInterp1(:,1),newInterp1(:,2))
            [newInterp2]=SubSurfinterp2(startPoints,nSteps);
            plot(newInterp2(:,1),newInterp2(:,2))
            [newcube]=SubSurfCubic(startPoints,nSteps);
            plot(newcube(:,1),newcube(:,2))
            newPoints.chaikin=newChaik;
            newPoints.bspline=newSpline;
            newPoints.newInterp1=newInterp1;
            newPoints.newInterp2=newInterp2;
            newPoints.newcube=newcube;
        case 'interp1'
            [newPoints,projPoints]=SubSurfinterp1(startPoints,nSteps);
        case 'interp2'
            [newPoints,projPoints]=SubSurfinterp2(startPoints,nSteps);
        case 'cubic'
            [newPoints,projPoints]=SubSurfCubic(startPoints,nSteps);
        case 'area'
            [newPoints,projPoints]=SubSurf_AreaConserv(startPoints,nSteps);
        otherwise
            
            error('Invalid method');
    end
    
    
    
end

%% Different subdivision processes

function [newPoints,projPoints]=SubSurfChainkin(startPoints,refineSteps)
    % Implements a Chainkin subdivision process
    TEisLeft=0;
    chainkinNoCorn=zeros([4,3]);
    chainkinNoCorn(1:4,2)=[0.25;0.75;0.75;0.25];
    
    chainkinCorn=zeros([5,3]);
    chainkinCorn(1:5,1)=[0;0.25;0;0;0];
    chainkinCorn(1:5,2)=[0.25;0.5;1;0.5;0.25];
    chainkinCorn(1:5,3)=[0;0;0;0.25;0];
    
    
    
    newPoints=startPoints;
    for nIter=1:refineSteps
        numPoints=length(startPoints(:,1));
        isCorner=DetectTrailingEdge(startPoints,TEisLeft);
        cumCorner=cumsum(isCorner);
        numNewPoints=(numPoints*2+cumCorner(end));
        subMask=zeros(numNewPoints,numPoints);
        
        for ii=0:numPoints-1
            iStart=ii-1;
            jStart=ii*2+cumCorner(ii+1)-isCorner(ii+1);
            
            if isCorner(ii+1)
                nJ=5;
                chainkinMask=chainkinCorn;
            else
                nJ=4;
                chainkinMask=chainkinNoCorn;
            end
            indX=zeros(1,3);
            indY=zeros(1,nJ);
            for iLoop=1:3
                indX(iLoop)=mod(iStart+(iLoop-1),numPoints)+1;
            end
            for jLoop=1:nJ
                indY(jLoop)=mod(jStart+(jLoop-1),numNewPoints)+1;
            end
            
            subMask(indY,indX)=chainkinMask+subMask(indY,indX);
        end
        newPoints=subMask*startPoints;
        startPointsCell{nIter}=startPoints;
        startPoints=newPoints;
        subMaskCell{nIter}=subMask;
    end
    
    limCurvMat=LimitCurve(subMask,4);
    [projPoints]=ProjectPoints(newPoints,limCurvMat,1);
end

function [newPoints,projPoints]=SubSurfBSpline(startPoints,refineSteps)
    % Implements a Chainkin subdivision process
    TEisLeft=0;
    bsplineNoCorn=zeros([7,7]);
    bsplineNoCorn(2:6,4)=[0.125;0.5;0.75;0.5;0.125];
    
    bsplineCorn=[[0.125;zeros(6,1)],[0.1875;0.125;zeros(5,1)],...
        [-0.3125;0;0;-0.125;zeros(3,1)],[0;0;0.5;1;0.5;0;0],[zeros(3,1);-0.125;0;0;-0.3125],...
        [zeros(5,1);0.125;0.1875],[zeros(6,1);0.125]];
    
    
    
    newPoints=startPoints;
    for nIter=1:refineSteps
        numPoints=length(startPoints(:,1));
        isCorner=DetectTrailingEdge(startPoints,TEisLeft)*0;
        cumCorner=cumsum(isCorner);
        numNewPoints=(numPoints*2);
        subMask=zeros(numNewPoints,numPoints);
        
        for ii=0:numPoints-1
            iStart=ii-3;
            jStart=ii*2;
            
            if isCorner(ii+1)
                nJ=7;
                bSplineMask=bsplineCorn;
            else
                nJ=7;
                bSplineMask=bsplineNoCorn;
            end
            indX=zeros(1,7);
            indY=zeros(1,nJ);
            for iLoop=1:7
                indX(iLoop)=mod(iStart+(iLoop-1),numPoints)+1;
            end
            for jLoop=1:nJ
                indY(jLoop)=mod(jStart+(jLoop-1),numNewPoints)+1;
            end
            
            subMask(indY,indX)=bSplineMask+subMask(indY,indX);
        end
        newPoints=subMask*startPoints;
        startPointsCell{nIter}=startPoints;
        startPoints=newPoints;
        
        subMaskCell{nIter}=subMask;
    end
    
    limCurvMat=LimitCurve(subMask,3);
    [projPoints]=ProjectPoints(newPoints,limCurvMat,1);
end

function [newPoints,projPoints]=SubSurfinterp1(startPoints,refineSteps)
    % Implements a Chainkin subdivision process
    TEisLeft=0;
    bsplineNoCorn=zeros([7,1]);
    testcoeff=1/16;
    bsplineNoCorn(1:7,1)=[-testcoeff;0;0.5+testcoeff;1;testcoeff+0.5;0;-testcoeff];
    
    bsplineCorn=bsplineNoCorn;
    
    
    
    newPoints=startPoints;
    for nIter=1:refineSteps
        numPoints=length(startPoints(:,1));
        isCorner=DetectTrailingEdge(startPoints,TEisLeft);
        cumCorner=cumsum(isCorner);
        numNewPoints=(numPoints*2);
        subMask=zeros(numNewPoints,numPoints);
        
        for ii=0:numPoints-1
            iStart=ii;
            jStart=ii*2;
            
            if isCorner(ii+1)
                nJ=7;
                bSplineMask=bsplineCorn;
            else
                nJ=7;
                bSplineMask=bsplineNoCorn;
            end
            indX=zeros(1,1);
            indY=zeros(1,nJ);
            for iLoop=1:1
                indX(iLoop)=mod(iStart+(iLoop-1),numPoints)+1;
            end
            for jLoop=1:nJ
                indY(jLoop)=mod(jStart+(jLoop-1),numNewPoints)+1;
            end
            
            subMask(indY,indX)=bSplineMask+subMask(indY,indX);
        end
        newPoints=subMask*startPoints;
        startPoints=newPoints;
        
    end
    limCurvMat=LimitCurve(subMask,5);
    [projPoints]=ProjectPoints(newPoints,limCurvMat,1);
end

function [newPoints,projPoints]=SubSurfCubic(startPoints,refineSteps)
    % Implements a Chainkin subdivision process
    TEisLeft=0;
    bsplineNoCorn=zeros([4,1]);
    
    bsplineNoCorn(1:4,1)=[-1/16;9/16;9/16;-1/16];
    
    bsplineCorn=bsplineNoCorn;
    
    
    
    newPoints=startPoints;
    for nIter=1:refineSteps
        numPoints=length(startPoints(:,1));
        isCorner=DetectTrailingEdge(startPoints,TEisLeft);
        cumCorner=cumsum(isCorner);
        numNewPoints=(numPoints*2);
        subMask=zeros(numNewPoints,numPoints);
        
        for ii=0:numPoints-1
            iStart=ii;
            jStart=ii*2;
            
            if isCorner(ii+1)
                nJ=4;
                bSplineMask=bsplineCorn;
            else
                nJ=4;
                bSplineMask=bsplineNoCorn;
            end
            indX=zeros(1,1);
            indY=zeros(1,nJ);
            for iLoop=1:1
                indX(iLoop)=mod(iStart+(iLoop-1),numPoints)+1;
            end
            for jLoop=1:nJ
                indY(jLoop)=mod(jStart+(jLoop-1),numNewPoints)+1;
            end
            
            subMask(indY,indX)=bSplineMask+subMask(indY,indX);
        end
        newPoints=subMask*startPoints;
        startPoints=newPoints;
        
    end
    limCurvMat=LimitCurve(subMask,5);
    [projPoints]=ProjectPoints(newPoints,limCurvMat,1);
end

function [newPoints,projPoints]=SubSurfinterp2(startPoints,refineSteps)
    % Implements a Chainkin subdivision process
    TEisLeft=0;
    bsplineNoCorn=zeros([7,1]);
    tcoeff=3/256;
    bsplineNoCorn(1:11,1)=...
        [tcoeff;0;-(1/16+3*tcoeff);0;9/16+2*tcoeff;1;9/16+2*tcoeff;0;...
        -(1/16+3*tcoeff);0;tcoeff];
    
    bsplineCorn=bsplineNoCorn;
    
    
    
    newPoints=startPoints;
    for nIter=1:refineSteps
        numPoints=length(startPoints(:,1));
        isCorner=DetectTrailingEdge(startPoints,TEisLeft);
        cumCorner=cumsum(isCorner);
        numNewPoints=(numPoints*2);
        subMask=zeros(numNewPoints,numPoints);
        
        for ii=0:numPoints-1
            iStart=ii;
            jStart=ii*2;
            
            if isCorner(ii+1)
                nJ=11;
                bSplineMask=bsplineCorn;
            else
                nJ=11;
                bSplineMask=bsplineNoCorn;
            end
            indX=zeros(1,1);
            indY=zeros(1,nJ);
            for iLoop=1:1
                indX(iLoop)=mod(iStart+(iLoop-1),numPoints)+1;
            end
            for jLoop=1:nJ
                indY(jLoop)=mod(jStart+(jLoop-1),numNewPoints)+1;
            end
            
            subMask(indY,indX)=bSplineMask+subMask(indY,indX);
        end
        newPoints=subMask*startPoints;
        startPoints=newPoints;
        
    end
    limCurvMat=LimitCurve(subMask,5);
    [projPoints]=ProjectPoints(newPoints,limCurvMat,1);
end

function [newPoints,projPoints]=SubSurf_AreaConserv(startPoints,refineSteps)
    
    newStencilInfo.varStencil=[0
        -1/40
        0
        (2^(1/2)*21^(1/2)*800^(1/2))/800 - 1/20
        21/40
        21/20 - (2^(1/2)*21^(1/2)*800^(1/2))/800
        21/20 - (2^(1/2)*21^(1/2)*800^(1/2))/800
        21/40
        (2^(1/2)*21^(1/2)*800^(1/2))/800 - 1/20
        0
        -1/40
        0];
    
    newStencilInfo.nNew=3;
    [newPoints,subMask]=SubSurfVarStencil_NoCorn(startPoints,refineSteps,newStencilInfo);
    limCurvMat=LimitCurve(subMask,7);
    [projPoints]=ProjectPoints(newPoints,limCurvMat,1/1.261436458122615);
end

%% General functions

function [newPoints,subMask]=SubSurfVarStencil_NoCorn(startPoints,refineSteps,newStencilInfo)
    % Implements a Chainkin subdivision process
    
    varStencil=newStencilInfo.varStencil;
    nNew=newStencilInfo.nNew;
    [nJ,nI]=size(varStencil);
    
    
    newPoints=startPoints;
    for nIter=1:refineSteps
        numPoints=length(startPoints(:,1));
        
        numNewPoints=(numPoints*nNew);
        subMask=zeros(numNewPoints,numPoints);
        
        for ii=0:numPoints-1
            iStart=ii;
            jStart=ii*nNew;
            bSplineMask=varStencil;
            
            indX=zeros(1,nI);
            indY=zeros(1,nJ);
            for iLoop=1:1
                indX(iLoop)=mod(iStart+(iLoop-1),numPoints)+1;
            end
            for jLoop=1:nJ
                indY(jLoop)=mod(jStart+(jLoop-1),numNewPoints)+1;
            end
            
            subMask(indY,indX)=bSplineMask+subMask(indY,indX);
        end
        newPoints=subMask*startPoints;
        startPointsCell{nIter}=startPoints;
        startPoints=newPoints;
        
        subMaskCell{nIter}=subMask;
    end
    
    limCurvMat=LimitCurve(subMask,7);
end

function isCorner=DetectTrailingEdge(coord,TEisLeft)
    TEisLeft=(TEisLeft-0.5)*2;
    testLocMin=((TEisLeft*(coord([2:end,1],1)-coord(:,1)))>0) ...
        & ((TEisLeft*(coord([end,1:end-1],1)-coord(:,1)))>0);
    isCorner=testLocMin;
end

function [limCurvMat,eigVal]=LimitCurve(subMask,nStencil)
    
    [nNew,nOld]=size(subMask);
    locStencil=cumsum(subMask~=0,2);
    iStart=sum(locStencil==0,2)+1;
    
    locStencil=cumsum(subMask(:,end:-1:1)~=0,2);
    iEnd=nOld-sum(locStencil==0,2);
    nStencLoc=iEnd-(iStart-1);
    [~,iCentre]=max(subMask,[],2);
    [~,iCentreRev]=max(subMask(:,end:-1:1),[],2);
    iCentreRev=nOld+1-iCentreRev;
    iCentreComp=abs(iCentre-iCentreRev);
    iCentreInverse=find(iCentreComp>2);
    for ii=1:length(iCentreInverse)
        iCentre(iCentreInverse(ii))=iCentreRev(iCentreInverse(ii));
    end
    nI=nStencil;
    nJ=nStencil;
    
    limCurvMat=zeros(nNew);
    eigVal=zeros(nNew,1);
    for ii=1:nNew
        iStart=-1+iCentre(ii)-floor(nStencil/2);
        jStart=-1+ii-floor(nStencil/2);
        
        indX=zeros(1,nI);
        indY=zeros(1,nJ);
        for iLoop=1:nI 
            indX(iLoop)=mod(iStart+(iLoop-1),nOld)+1;
        end
        for jLoop=1:nJ
            indY(jLoop)=mod(jStart+(jLoop-1),nNew)+1;
        end
        
        [w,d]=eig([subMask(indY,indX)]');
        [iEig,~]=find(d==1,1);
        if numel(iEig)==0
            [iEig,~]=find(1-d<1e-10);
        end
        limCurvMat(ii,indY)=w(:,iEig)';
         limCurvMat(ii,indY)= limCurvMat(ii,indY)/sum( limCurvMat(ii,indY));
        eigVal(ii)=d(iEig);
    end

end

function [projPoints]=ProjectPoints(points,eigMat,convFactor)
    
    projPoints=convFactor*eigMat*points;
    
end
%% OLD

%{
function [newPoints]=SubSurfChainkin2(startPoints,refineSteps)
    % Implements a Chainkin subdivision process
    
    
    chainkinMask=[0.25,0.75,0.75,0.25]';
    newPoints=startPoints;
    for nIter=1:refineSteps
        numPoints=length(startPoints(:,1));
        subMask=zeros((numPoints*2+2),numPoints);
        
        % create the right sized mask
        for ii=1:numPoints
            jj=(1+(ii-1)*2);
            subMask(jj:jj+3,ii)=chainkinMask;
        end
        
        newPoints=subMask*startPoints;
        newPoints(1:2,:)=[];
        newPoints(end-1:end,:)=[];
        startPoints=newPoints;
    end
end
%}

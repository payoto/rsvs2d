%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%              University of Bristol
%      Department of Aerospace Engineering
%                     2015
%
%          Subdivision for Shape refinement
%             Alexandre Payot
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function newPoints=SubDivision(startPoints,nSteps,refineMethod)
    
    
    switch refineMethod
        case 'chainkin'
            [newPoints]=SubSurfChainkin(startPoints,nSteps);
        case 'bspline'
            
            [newPoints]=SubSurfBSpline(startPoints,nSteps);
        case 'test'
            figure
            hold on
            plot(startPoints(:,1),startPoints(:,2),'b--')
            [newChaik]=SubSurfChainkin(startPoints,nSteps);
            plot(newChaik(:,1),newChaik(:,2),'mo-');
            [newSpline]=SubSurfBSpline(startPoints,nSteps);
            plot(newSpline(:,1),newSpline(:,2),'g+-')
            newPoints.chaikin=newChaik;
            newPoints.bspline=newSpline;
        case 'interp1'
            [newPoints]=SubSurfinterp1(startPoints,nSteps);
        case 'interp2'
            [newPoints]=SubSurfinterp2(startPoints,nSteps);
        otherwise
            
            error('Invalid method');
    end
    
    
    
end


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

function [newPoints]=SubSurfChainkin(startPoints,refineSteps)
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
        startPoints=newPoints;
        
    end
end

function [newPoints]=SubSurfBSpline(startPoints,refineSteps)
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
        isCorner=DetectTrailingEdge(startPoints,TEisLeft);
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
        startPoints=newPoints;
        
    end
end

function [newPoints]=SubSurfinterp1(startPoints,refineSteps)
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
end

function [newPoints]=SubSurfinterp2(startPoints,refineSteps)
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
end

function isCorner=DetectTrailingEdge(coord,TEisLeft)
    TEisLeft=(TEisLeft-0.5)*2;
    testLocMin=((TEisLeft*(coord([2:end,1],1)-coord(:,1)))>0) ...
        & ((TEisLeft*(coord([end,1:end-1],1)-coord(:,1)))>0);
    isCorner=testLocMin;
end

function [newPoints]=SubSurfVarStencil_NoCorn_STA(startPoints,refineSteps,newStencilInfo)
    % Implements a Chainkin subdivision process
    
    varStencil=newStencilInfo.varStencil;
    nNew=newStencilInfo.nNew;
    [nJ,nI]=size(varStencil);
    
    
    newPoints=startPoints;
    for nIter=1:refineSteps
        numPoints=length(startPoints(:,1));
        
        numNewPoints=(numPoints*nNew);
        numPoints
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
        startPoints=newPoints;
        
    end
end
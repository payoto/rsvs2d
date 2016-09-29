%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%              University of Bristol
%      Department of Aerospace Engineering
%                     2015
%
%          snake parametrisation
%      for Aerodynamic shape parametrisation
%           - Grid and Volume fraction matching-
%             Alexandre Payot
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


function []=GridMatching(origGrid,newGrid,refineOrig,refineNew)
    % This grid matching works on the implicit assumption that the grids
    % come from the same parent grid This builds connectivity information
    % between the two grids at a volume fraction level.
    % There is no edge matching
    
    
    
    BuildCellMatchTemplate(refineNew,refineOrig)
    
    
end



function []=BuildCellMatchTemplate(newSize,oldSize)
    % Accepts the split levels for each
    
    % perform checks
    if any(newSize>oldSize), 
        warning('New Mesh is finer than old refined mesh, Data will not pass well')
    end
    if any(mod(newSize,oldSize)),
        warning('Meshes do not match exactly data transfer will be approximate')
    end
    
    for ii=1:2
        [linearCoeffs{ii}]=FindRatios(newSize(ii),oldSize(ii));
    end
    indN=1:prod(newSize);
    indO=1:prod(oldSize);
    
    allCoeffs{1}=repmat(linearCoeffs{1},[newSize(2),oldSize(2)]);
    allCoeffs{2}=repmat(linearCoeffs{2},[newSize(1),oldSize(1)]);
    ordDim2N=mod(0:newSize(2):prod(newSize)*newSize(2)-1,prod(newSize)-1)+1;
    ordDim2N(end)=prod(newSize);
    ordDim2O=mod(0:oldSize(2):prod(oldSize)*oldSize(2)-1,prod(oldSize)-1)+1;
    ordDim2O(end)=prod(oldSize);
    allCoeffs{2}=allCoeffs{2}(ordDim2N,ordDim2O);
    
    coeffs=allCoeffs{1}.*allCoeffs{2};
end

function [coeffs]=FindRatios(sN,sO)
    
    lN=1/sN;
    lO=1/sO;
    
    coordN=0:lN:1;
    coordO=0:lO:1;
    endN=repmat(coordN(2:end),[sO,1])';
    endO=repmat(coordO(2:end),[sN,1]);
    startN=repmat(coordN(1:end-1),[sO,1])';
    startO=repmat(coordO(1:end-1),[sN,1]);
    
    starts=cat(3,startN,startO);
    ends=cat(3,endN,endO);
    
    coeffs=max(min(ends,[],3)-max(starts,[],3),0)/lO;
    
    
end




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


function [gridmatch,grids]=GridMatching(origGrid,newGrid,refineOrig,refineNew)
    % This grid matching works on the implicit assumption that the grids
    % come from the same parent grid This builds connectivity information
    % between the two grids at a volume fraction level.
    % There is no edge matching
    
    gridmatch.coeffs=BuildCellMatchTemplate(refineNew,refineOrig);
    gridmatch.origbreak=refineOrig;
    gridmatch.newbreak=refineNew;
    gridmatch.matchstruct=repmat(struct('newGridInd',[],'newvolume',[],'oldGridInd',[],...
        'oldvolume',[],'coeff',[]),[1 numel(newGrid.refined.cell)]);
    
    
    gridmatch.matchstruct=PopulateMatchStruct(origGrid,newGrid,...
        gridmatch.matchstruct,gridmatch.coeffs);
    grids=struct('origin',origGrid,'new',newGrid,'match',gridmatch);
end

function [coeffs]=BuildCellMatchTemplate(newSize,oldSize)
    % Accepts the split levels for each level
    
    % perform checks
    if any(newSize>oldSize), 
        warning('New Mesh is finer than old refined mesh, Data will not pass well')
    end
    if any(mod(oldSize,newSize)),
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

function [matchstruct]=PopulateMatchStruct(origGrid,newGrid,matchstruct,coeffs)
    
    oldIndsNewOrd=cell2mat(cellfun(@(new,old)old*ones([1,numel(new)]),...
        {newGrid.connec.cell(:).new},...
        {newGrid.connec.cell(:).old},'UniformOutput',false));
    
    oldIndsN=[origGrid.connec.cell(:).old];
    oldSubsNewOrd=FindObjNum([],oldIndsNewOrd,oldIndsN);
    newIndsNewOrd=[newGrid.connec.cell(:).new];
    newSubGridOrd=FindObjNum([],[newGrid.cellrefined(:).index],newIndsNewOrd);
    oldIndRef=[origGrid.cellrefined(:).index];
    
    for ii=1:numel(newGrid.cellrefined)
        matchstruct(ii).newGridInd=newGrid.cellrefined(ii).index;
        matchstruct(ii).oldGridInd=origGrid.connec.cell(oldSubsNewOrd(newSubGridOrd(ii))).new;
        matchstruct(ii).newvolume=newGrid.cellrefined(ii).volume;
        oldCellSub=FindObjNum([],matchstruct(ii).oldGridInd,oldIndRef);
        oldRefVec=vertcat(origGrid.cellrefined(oldCellSub).refineVec);
        
        [~,ordOldCell]=sort(oldRefVec(:,end));
        matchstruct(ii).oldGridInd=matchstruct(ii).oldGridInd(ordOldCell);
        matchstruct(ii).coeff=coeffs(newGrid.cellrefined(ii).refineVec(end),:);
        
        isact=matchstruct(ii).coeff~=0;
        matchstruct(ii).oldGridInd=matchstruct(ii).oldGridInd(isact);
        matchstruct(ii).coeff=matchstruct(ii).coeff(isact);
        matchstruct(ii).oldvolume=[origGrid.cellrefined(oldCellSub(ordOldCell(isact))).volume];
    end
    
    
    
end



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
    gridmatch.matchstruct=repmat(struct('newGridInd',[],'newvolume',[],'oldGridInd',[],...
        'oldvolume',[],'coeff',[]),[1 numel(newGrid.refined.cell)]);
    isRefineCell=[origGrid.base.cell(:).isrefine];
    isSnakRefCell=[origGrid.base.cell(:).snakref]+1;
    
    [origGrid.base.cell(:).isrefine]=deal(false);
    gridmatch.matchstruct=PopulateMatchStruct(origGrid,newGrid,...
        gridmatch.matchstruct,[],true);
    for kk = 1:size(refineOrig,1)
        for ii=1:size(refineNew,1)
            [gridmatch.dat(ii,kk).coeffs,gridmatch.dat(ii,kk).warn]=...
                BuildCellMatchTemplate(refineNew(ii,:),refineOrig(kk,:));
            gridmatch.dat(ii,kk).origbreak=refineOrig(kk,:);
            gridmatch.dat(ii,kk).newbreak=refineNew(ii,:);
            for jj=1:numel(origGrid.base.cell)
                origGrid.base.cell(jj).isrefine=isRefineCell(jj)==ii && isSnakRefCell(jj)==kk;
                newGrid.base.cell(jj).isrefine=isRefineCell(jj)==ii && isSnakRefCell(jj)==kk;
            end
            gridmatch.matchstruct=PopulateMatchStruct(origGrid,newGrid,...
                gridmatch.matchstruct,gridmatch.dat(ii,kk),false);
        end
    end
    
    grids=struct('origin',origGrid,'new',newGrid,'match',gridmatch);
end

function [coeffs,isWarn]=BuildCellMatchTemplate(newSize,oldSize)
    % Accepts the split levels for each level
    isWarn = 0;
    % perform checks
    if any(newSize>oldSize),
%         warning('New Mesh is finer than old refined mesh, Data will not pass well')
        isWarn = 1;
    end
    if any(mod(oldSize,newSize)),
%         warning('Meshes do not match exactly data transfer will be approximate')
        isWarn = 2;
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

function [matchstruct]=PopulateMatchStruct(origGrid,newGrid,matchstruct,matchcoeffs,isPreparation)
    
    
    oldIndsNewOrd=cell2mat(cellfun(@(new,old)old*ones([1,numel(new)]),...
        {newGrid.connec.cell(:).new},...
        {newGrid.connec.cell(:).old},'UniformOutput',false));
    
    oldIndsN=[origGrid.connec.cell(:).old];
    oldCoarseInds=[origGrid.base.cell(:).index];
    oldSubsNewOrd=FindObjNum([],oldIndsNewOrd,oldIndsN);
    newIndsNewOrd=[newGrid.connec.cell(:).new];
    newSubGridOrd=FindObjNum([],[newGrid.cellrefined(:).index],newIndsNewOrd);
    oldIndRef=[origGrid.cellrefined(:).index];
    warn = 0;
    for ii=1:numel(newGrid.cellrefined)
        
        oldCoarseSub=FindObjNum([],origGrid.connec.cell(oldSubsNewOrd(newSubGridOrd(ii))).old,oldCoarseInds);
        
        
        
        if origGrid.base.cell(oldCoarseSub).isrefine
            matchstruct(ii).newGridInd=newGrid.cellrefined(ii).index;
            matchstruct(ii).oldGridInd=origGrid.connec.cell(oldSubsNewOrd(newSubGridOrd(ii))).new;
            matchstruct(ii).newvolume=newGrid.cellrefined(ii).volume;
            oldCellSub=FindObjNum([],matchstruct(ii).oldGridInd,oldIndRef);
            oldRefVec=vertcat(origGrid.cellrefined(oldCellSub).refineVec);
            
            [~,ordOldCell]=sort(oldRefVec(:,end));
            
            matchstruct(ii).oldGridInd=matchstruct(ii).oldGridInd(ordOldCell);
            matchstruct(ii).coeff=matchcoeffs.coeffs(newGrid.cellrefined(ii).refineVec(end),:);
            
            warn = matchcoeffs.warn;
            isact=matchstruct(ii).coeff~=0;
            try
                matchstruct(ii).oldGridInd=matchstruct(ii).oldGridInd(isact);
            catch ME
                ii
                matchstruct(ii).oldGridInd
                isact
                matchstruct(ii).oldGridInd=matchstruct(ii).oldGridInd(isact);
            end
            matchstruct(ii).coeff=matchstruct(ii).coeff(isact);
            matchstruct(ii).oldvolume=[origGrid.cellrefined(oldCellSub(ordOldCell(isact))).volume];
        elseif isPreparation
            matchstruct(ii).newGridInd=newGrid.cellrefined(ii).index;
            matchstruct(ii).oldGridInd=origGrid.connec.cell(oldSubsNewOrd(newSubGridOrd(ii))).new;
            matchstruct(ii).newvolume=newGrid.cellrefined(ii).volume;
            oldCellSub=FindObjNum([],matchstruct(ii).oldGridInd,oldIndRef);
            oldRefVec=vertcat(origGrid.cellrefined(oldCellSub).refineVec);
            [~,ordOldCell]=sort(oldRefVec(:,end));
            matchstruct(ii).coeff=ones(size(matchstruct(ii).oldGridInd));
            matchstruct(ii).oldvolume=[origGrid.cellrefined(oldCellSub(ordOldCell)).volume];
        else
        end
        %         if ~origGrid.base.cell(oldCellSub(ordOldCell(isact))).isrefine
        %             matchstruct(ii).coeff=1;
        %         end
    end
    
    if warn==1
        warning('New Mesh is finer than old refined mesh, Data will not pass well')

    elseif warn==2
        warning('Meshes do not match exactly data transfer will be approximate')
    end
    
end



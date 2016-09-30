function [profloops,transformstruct]=RebuildIterOptim(iterstruct,newGrid,gridmatch,profloops)
    
    [transformstruct,~]=BuildMatrix(gridmatch);
    
    [profloops]=ConvertProfToFill(profloops,transformstruct);
end


function [transformstruct,coeffMat]=BuildMatrix(gridmatch)
    
    newGridIndsMulti=cell2mat(cellfun(@(new,old)old*ones([1,numel(new)]),...
        {gridmatch.matchstruct(:).oldGridInd},...
        {gridmatch.matchstruct(:).newGridInd},'UniformOutput',false));
    oldGridInds=[gridmatch.matchstruct(:).oldGridInd];
    oldGridCoeff=[gridmatch.matchstruct(:).coeff];
    newGridInds=[gridmatch.matchstruct(:).newGridInd];
    
    oldGridUniq=RemoveIdenticalEntries(sort(oldGridInds));
    oldGridSub=FindObjNum([],oldGridInds,oldGridUniq);
    newGridMultiSub=FindObjNum([],newGridIndsMulti,newGridInds);
    
    coeffMat=zeros([numel(newGridInds),numel(oldGridUniq)]);
    coeffMat(sub2ind(size(coeffMat),newGridMultiSub,oldGridSub))=oldGridCoeff;
    transformstruct.coeff=coeffMat;
    transformstruct.indNew=newGridInds;
    transformstruct.indOld=oldGridUniq;
    
end

function [profloops]=ConvertProfToFill(profloops,transformstruct)
    for ii=1:length(profloops)
        volSubs=FindObjNum([],profloops(ii).refinevolfrac.index,transformstruct.indOld);
        profloops(ii).newFracs=transformstruct.coeff*profloops(ii).refinevolfrac.fractionvol(volSubs)';
    end
end
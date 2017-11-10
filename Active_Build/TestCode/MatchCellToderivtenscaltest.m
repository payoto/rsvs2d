function [derivtenscal]=MatchCellToderivtenscaltest(derivtenscal,coeffstructure,volumefraction)
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
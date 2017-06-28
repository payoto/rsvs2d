function BP3_invdestoporestarts(refCrit,cornAct,aeroStr,numRat)
    MoveToDir('source',1)
    InitialiseSnakeFlow;
    numRat=floor((numRat-1)/5);
    disp(sprintf('invdestopo(''%s'',%i,''%s'',%i)',refCrit,cornAct,aeroStr,numRat))
    
    caseStr=sprintf('invdestopo(''%s'',%i,''%s'',%i)',refCrit,cornAct,aeroStr,numRat);
    caseStr=regexprep(regexprep(...
        regexprep(caseStr,'(\(|\)|,)','_'),'\.','_'),'''','');
    
    [restartPath]=IdentifyRestart('~/SnakVolParam/restarts',caseStr,'Restart_','.mat');
    
    ExecuteOptimisation(sprintf('invdestopo(''%s'',%i,''%s'',%i)',refCrit,cornAct,aeroStr,numRat),...
        {restartPath,{'conjgrad',1}});
    

end
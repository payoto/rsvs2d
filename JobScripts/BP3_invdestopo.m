function BP3_invdestopo(refCrit,cornAct,aeroStr,numRat)
    MoveToDir('source',1)
    InitialiseSnakeFlow;
    numRat=floor((numRat-1)/10);
    disp(sprintf('invdestopo(''%s'',%i,''%s'',%i)',refCrit,cornAct,aeroStr,numRat))
    ExecuteOptimisation(sprintf('invdestopo(''%s'',%i,''%s'',%i)',refCrit,cornAct,aeroStr,numRat));
    

end
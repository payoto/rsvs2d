function BP3_invdeslocaltest6(gridCase,airfoil,cornAct,numRat)
    MoveToDir('source',1)
    InitialiseSnakeFlow;
    numRat=floor((numRat-1)/12)+1;
    disp(sprintf('invdeslocal_test6(''%s'',''%s'',%i,%i)',gridCase,airfoil,cornAct,numRat))
    ExecuteOptimisation(sprintf('invdeslocal_test6(''%s'',''%s'',%i,%i)',gridCase,airfoil,cornAct,numRat));
    

end
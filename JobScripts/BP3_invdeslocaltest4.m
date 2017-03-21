function BP3_invdeslocaltest4(gridCase,airfoil,cornAct,numRat)
    MoveToDir('source',1)
    InitialiseSnakeFlow;
    numRat=floor((numRat-1)/12)+1;
    disp(sprintf('invdeslocal_test3(''%s'',''%s'',%i,%i)',gridCase,airfoil,cornAct,numRat))
    ExecuteOptimisation(sprintf('invdeslocal_test3(''%s'',''%s'',%i,%i)',gridCase,airfoil,cornAct,numRat));
    

end
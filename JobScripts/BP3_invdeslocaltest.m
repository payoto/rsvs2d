function BP3_invdeslocaltest(gridCase,airfoil,cornAct)
    MoveToDir('source',1)
    InitialiseSnakeFlow;

    disp(sprintf('invdeslocal_test(''%s'',''%s'',%i)',gridCase,airfoil,cornAct))
    ExecuteOptimisation(sprintf('invdeslocal_test(''%s'',''%s'',%i)',gridCase,airfoil,cornAct));
    

end
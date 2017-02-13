function BP3_invdeslocaltest(gridCase,airfoil)
    MoveToDir('source',1)
    InitialiseSnakeFlow;

    disp(sprintf('invdeslocal_test(''%s'',''%s'')',gridCase,airfoil))
    ExecuteOptimisation(sprintf('invdeslocal_test(''%s'',''%s'')',gridCase,airfoil));
    

end
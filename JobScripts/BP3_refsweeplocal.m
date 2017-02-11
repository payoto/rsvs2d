function BP3_refsweeplocal(gridCase,airfoil)
    MoveToDir('source',1)
    InitialiseSnakeFlow;

    disp(sprintf('refsweeplocal(''%s'',''%s'')',gridCase,airfoil))
    ExecuteOptimisation(sprintf('refsweeplocal(''%s'',''%s'')',gridCase,airfoil));
    

end
function BP3_refsweep(gridCase,airfoil,ii)
    MoveToDir('source',1)
    InitialiseSnakeFlow;

    disp(sprintf('refsweep(''%s'',''%s'',%i)',gridCase,airfoil,ii))
    ExecuteOptimisation(sprintf('refsweep(''%s'',''%s'',%i)',gridCase,airfoil,ii));
    

end
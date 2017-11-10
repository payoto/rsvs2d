function BP3_NACA0012Sweep(gridCase,ii,optim)
    MoveToDir('source',1)
    InitialiseSnakeFlow;

    disp(sprintf('NACA0012Sweep(''%s'',%i,''%s'')',gridCase,ii,optim))
    ExecuteOptimisation(sprintf('NACA0012Sweep(''%s'',%i,''%s'')',gridCase,ii,optim));
    

end
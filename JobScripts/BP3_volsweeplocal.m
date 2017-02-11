function BP3_volsweeplocal(e,gridCase)
    MoveToDir('source',1)
    InitialiseSnakeFlow;

    disp(sprintf('volsweeplocal(%e,''%s'')',e,gridCase))
    ExecuteOptimisation(sprintf('volsweeplocal(%e,''%s'')',e,gridCase));
    

end
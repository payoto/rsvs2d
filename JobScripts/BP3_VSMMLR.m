function BP3_VSMMLR(e)
    MoveToDir('source',1)
    InitialiseSnakeFlow;

    disp(sprintf('volsweepmmesh(''%s'')',e))
    ExecuteOptimisation(sprintf('volsweepmmesh(''%s'')',e));
    

end
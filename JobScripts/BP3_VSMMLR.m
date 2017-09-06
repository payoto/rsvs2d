function BP3_VSMMLR(e)
    MoveToDir('source',1)
    InitialiseSnakeFlow;

    disp(sprintf('volsweepmmesh(%.1e)',e))
    ExecuteOptimisation(sprintf('volsweepmmesh(%.1e)',e));
    

end
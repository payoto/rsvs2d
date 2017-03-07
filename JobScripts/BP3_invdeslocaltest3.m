function BP3_invdeslocaltest3(offset)
    MoveToDir('source',1)
    InitialiseSnakeFlow;

    disp(sprintf('RefLocOut(%.1e)',offset))
    ExecuteOptimisation(sprintf('RefLocOut(%.1e)',offset));
    

end
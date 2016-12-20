function BP3_volsweeprefine(e,gridCase,ii)
    MoveToDir('source',1)
    InitialiseSnakeFlow;

    disp(sprintf('volsweeprefine(%e,''%s'',%i)',e,gridCase,ii))
    ExecuteOptimisation(sprintf('volsweeprefine(%e,''%s'',%i)',e,gridCase,ii));
    

end
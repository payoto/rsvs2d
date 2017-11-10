function BP3_N12MeshMotionLocalRefine(gridCase,optim)
    MoveToDir('source',1)
    InitialiseSnakeFlow;

    disp(sprintf('N12_LRef_MMesh(''%s'',''%s'')',gridCase,optim))
    ExecuteOptimisation(sprintf('N12_LRef_MMesh(''%s'',''%s'')',gridCase,optim));
    

end
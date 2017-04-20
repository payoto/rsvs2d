function BP3_N12MeshMotionLocalRefine_re1
    MoveToDir('source',1)
    InitialiseSnakeFlow;
    gridCase='uo';
    optim='BFGS';
    disp(sprintf('N12_LRef_MMesh(''%s'',''%s'')',gridCase,optim))
    ExecuteOptimisation(sprintf('N12_LRef_MMesh(''%s'',''%s'')',gridCase,optim),{'Restart_uo_BFGS__1_2.mat',{'conjgrad',true}});
    

end
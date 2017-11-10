function BP3_N12MeshMotionLocalRefine_re2
    MoveToDir('source',1)
    InitialiseSnakeFlow;
    gridCase='uo';
    optim='conjgrad';
    disp(sprintf('N12_LRef_MMesh(''%s'',''%s'')',gridCase,optim))
    ExecuteOptimisation(sprintf('N12_LRef_MMesh(''%s'',''%s'')',gridCase,optim),{'Restart_uo_conjgrad__1_2.mat',{'conjgrad',true}});
    

end
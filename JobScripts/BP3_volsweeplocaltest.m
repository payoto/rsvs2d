function BP3_volsweeplocaltest(gridCase)
    MoveToDir('source',1)
    InitialiseSnakeFlow;
% 
%     disp(sprintf('volsweeplocal_test(%e,''%s'')',e,gridCase))
%     ExecuteOptimisation(sprintf('volsweeplocal_test(%e,''%s'')',e,gridCase));
    disp(sprintf('volsweeplocal_test(''%s'')',gridCase))
    ExecuteOptimisation(sprintf('volsweeplocal_test(''%s'')',gridCase));
    

end
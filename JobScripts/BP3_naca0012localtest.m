function BP3_naca0012localtest(gridCase,refmet)
    MoveToDir('source',1)
    InitialiseSnakeFlow;
% 
%     disp(sprintf('volsweeplocal_test(%e,''%s'')',e,gridCase))
%     ExecuteOptimisation(sprintf('volsweeplocal_test(%e,''%s'')',e,gridCase));
    disp(sprintf('NACA0012local(''%s'',''%s'')',gridCase,refmet))
    ExecuteOptimisation(sprintf('NACA0012local(''%s'',''%s'')',gridCase,refmet));
    

end
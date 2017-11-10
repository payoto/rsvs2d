function BP3_LocalConstr(ii)
    MoveToDir('source',1)
    InitialiseSnakeFlow;

    disp(['LocConsTopo_prof(',num2str(ii,'%e'),')'])
    ExecuteOptimisation(['LocConsTopo_prof(',num2str(ii,'%e'),')']);
    



end


